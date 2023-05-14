/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include "geom3d.hpp"
#include "graph_gn.hpp"
#include <fstream>
#include <iostream>

using namespace opt4imu;

struct Con6D
{
    using Mat6D = Eigen::Matrix<double, 6, 6>;
    using Vec6D = Eigen::Matrix<double, 6, 1>;
    int id1, id2;
    Vec6D t;
    Mat6D info;
    void setCovarianceMatrix( const Mat6D &mat )
    {
        info = mat.inverse();
    }
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

using Con6DVec = std::vector<Con6D>;
using RVVec = std::vector<RotVec>;
using Node = Eigen::Matrix<double, 6, 1>;
using Jac  = Eigen::Matrix<double, 6, 6>;
using NodeList = Eigen::Matrix<double, 6, Eigen::Dynamic>;

class GraphGNSE3 : public GraphGaussNewtonOptimizer<double, 6, Con6D>
{
    bool debug_numerical_grad = false;
    static Mat3D Rx() {
        Mat3D mat;
        mat << 0, 0, 0, 0, 0,-1, 0, 1, 0;
        return mat;
    }
    static Mat3D Ry() {
        Mat3D mat;
        mat << 0, 0, 1, 0, 0, 0, -1, 0, 0;
        return mat;
    }
    static Mat3D Rz() {
        Mat3D mat;
        mat << 0,-1, 0, 1, 0, 0, 0, 0, 0;
        return mat;
    }

    static Eigen::VectorXd errorFunc(const Pose3D &pa, const Pose3D &pb, const Vec6D &con)
    {
        Pose3D ba = pb.ominus(pa);
        Pose3D err = ba.ominus(Pose3D(con));
        return err.vec();
    }

    Eigen::VectorXd calcError(const Node &ina, const Node &inb, const Con6D &edge, Jac &Ja, Jac &Jb)
    {
        Pose3D pa(ina);
        Pose3D pb(inb);
        Pose3D con(edge.t);

        Ja = Jb = Mat6D::Identity();
        Eigen::VectorXd e0 = errorFunc(pa, pb, edge.t);
        Mat3D RaInv = pa.rv.toRotationMatrix().transpose();
        Mat3D Rb = pb.rv.toRotationMatrix();
        Mat3D RcInv = con.rv.inverted().toRotationMatrix();
        RotVec rve(e0.tail<3>());
        Mat3D Re = rve.toRotationMatrix();

        Ja.block<3,3>(0, 0) = -RcInv;
        Jb.block<3,3>(0, 0) = RotVec(e0.tail<3>()).toRotationMatrix();

        // Rc' * {Rx' Ry' Rz'} * R0' * d
        Vec3D d;
        d << pb.x - pa.x, pb.y - pa.y, pb.z - pa.z;
        Ja.block<3,1>(0,3) = -RcInv * Rx() * RaInv * d;
        Ja.block<3,1>(0,4) = -RcInv * Ry() * RaInv * d;
        Ja.block<3,1>(0,5) = -RcInv * Rz() * RaInv * d;

        // rotation part: Re = Rc' * Ra' * Rb;
        double buf[9 * 3];
        Eigen::Map<Eigen::Matrix<double, 9, 3>> dRedRa(buf);
        Mat3D dx = -RcInv * Rx() * RaInv * Rb;
        Mat3D dy = -RcInv * Ry() * RaInv * Rb;
        Mat3D dz = -RcInv * Rz() * RaInv * Rb;
        memcpy(buf   , dx.data(), sizeof(double)*9);
        memcpy(buf+9 , dy.data(), sizeof(double)*9);
        memcpy(buf+18, dz.data(), sizeof(double)*9);
        Ja.block<3,3>(3,3) = jacobianR(Re) * dRedRa;
        Jb.block<3,3>(3,3) = rightJacobianInvRV(rve);
        if (debug_numerical_grad) {
            double eps = 1e-5;
            Mat6D JaN = Mat6D::Identity();
            Mat6D JbN = Mat6D::Identity();
            for(size_t i=0; i<6; ++i) {
                double d[6] = {0, 0, 0, 0, 0, 0};
                d[i] = eps;
                RotVec dr(d[3], d[4], d[5]);
                Pose3D pa1 = pa.oplus(Pose3D(d[0], d[1], d[2], dr));
                Pose3D pb1 = pb.oplus(Pose3D(d[0], d[1], d[2], dr));
                JaN.block<6,1>(0,i) = (errorFunc(pa1, pb, edge.t) - e0)/eps;
                JbN.block<6,1>(0,i) = (errorFunc(pa, pb1, edge.t) - e0)/eps;
            }
            std::cout << "----- Ja, JaN" << std::endl;
            std::cout << Ja << std::endl;
            std::cout << JaN << std::endl;

            std::cout << "----- Jb, JbN" << std::endl;
            std::cout << Jb << std::endl;
            std::cout << JbN << std::endl;
            Ja = JaN;
            Jb = JbN;
        }
        return e0;
    }

    static NodeList updatefunc(const NodeList &n, const Eigen::Matrix<double, Eigen::Dynamic, 1> &d)
    {
        NodeList out(6, n.cols());
        for(unsigned int i=0;i<n.cols();i++) {
            Pose3D p(n.col(i));
            Pose3D dpose(d.segment(i*6,6));
            out.col(i) = p.oplus(dpose).vec();
        }
        return out;
    }

public:
    GraphGNSE3() {
        using namespace std::placeholders;
        error_func = std::bind(&GraphGNSE3::calcError, this, _1, _2, _3, _4, _5);
        update_func = &GraphGNSE3::updatefunc;
        //oplus_func = &Pose3D::concat;
    }
};
