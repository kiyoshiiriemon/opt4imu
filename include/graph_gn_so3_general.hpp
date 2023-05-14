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

namespace opt4imu {

struct ConSO3WithErrorFunc
{
    using Node = Eigen::Matrix<double, 3, 1>;
    using Jac  = Eigen::Matrix<double, Eigen::Dynamic, 3>;
    int type = 0;
    int id1;
    int id2;
    Vec3D t; // SO3 or Euler angles (Z-Y-X)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> info;
    std::function<Eigen::VectorXd (const Node& , const Node &, const ConSO3WithErrorFunc &con, Jac &, Jac &)> error_func;
    static Vec3D errorFunc(const RotVec &pa, const RotVec &pb, const RotVec &con)
    {
        RotVec err = con.inverted().concat(pa.inverted()).concat(pb);
        return err.vec();
    }
    static Eigen::VectorXd calcErrorRelSO3(const Node &ina, const Node &inb, const ConSO3WithErrorFunc &con, Jac &Ja, Jac &Jb)
    {
        Mat3D Rx, Ry, Rz;
        Rx << 0, 0, 0,
                0, 0,-1,
                0, 1, 0;
        Ry << 0, 0, 1,
                0, 0, 0,
                -1, 0, 0;
        Rz << 0,-1, 0,
                1, 0, 0,
                0, 0, 0;
        RotVec  pa(ina[0], ina[1], ina[2]);
        RotVec  pb(inb[0], inb[1], inb[2]);
        RotVec rel = con.t;
        //std::cout << pa.vec().transpose() << " " << pb.vec().transpose() << " " << con.t.transpose() << std::endl;
        Eigen::VectorXd e0 = errorFunc(pa, pb, rel);
        Mat3D RaInv = pa.toRotationMatrix().transpose();
        Mat3D Rb = pb.toRotationMatrix();
        Mat3D RcInv = rel.toRotationMatrix().transpose();
        RotVec rve(e0.tail<3>());
        Mat3D Re = rve.toRotationMatrix();

        // rotation part: Re = Rc' * Ra' * Rb;
        double buf[9 * 3];
        Eigen::Map<Eigen::Matrix<double, 9, 3>> dRedRa(buf);
        Mat3D dx = -RcInv * Rx * RaInv * Rb;
        Mat3D dy = -RcInv * Ry * RaInv * Rb;
        Mat3D dz = -RcInv * Rz * RaInv * Rb;
        memcpy(buf   , dx.data(), sizeof(double)*9);
        memcpy(buf+9 , dy.data(), sizeof(double)*9);
        memcpy(buf+18, dz.data(), sizeof(double)*9);
        //std::cout << "--- dRedRa, rightJacobianInv ---" << std::endl;
        //std::cout << dRedRa << std::endl;
        //std::cout << rightJacobianInvRV(rve) << std::endl;
        Ja = jacobianR(Re) * dRedRa;

        Jb = rightJacobianInvRV(rve);
        return e0;
    }
    static Eigen::VectorXd calcErrorGravity(const Node &ina, const Node &inb, const ConSO3WithErrorFunc &con, Jac &Ja, Jac &Jb)
    {
        double rx0, ry0, rz0;
        double ix, iy, iz;
        RotVec rv0(inb);
        rv0.getEulerZYX(&rx0, &ry0, &rz0);
        RotVec acc_rv(con.t);
        acc_rv.getEulerZYX(&ix, &iy, &iz);
        //std::cout << rx0 << " " << ry0 << " " << rz0 << std::endl;
        //std::cout << ix << " " << iy << " " << iz << std::endl;
        double eps = 1e-5;
        Ja = Eigen::MatrixXd(2, 3);
        //Ja.resize(2,3);
        Ja.setZero();
        Jb = Eigen::MatrixXd(2, 3);
        //Jb.resize(2,3);
        for(int i=0; i<3; ++i) {
            double d[3] = {0, 0, 0};
            d[i] = eps;
            double rx, ry, rz;
            RotVec rv = rv0.concat(RotVec(d[0], d[1], d[2]));
            rv.getEulerZYX(&rz, &ry, &rx);
            Jb(0, i) = (rx-rx0)/eps;
            Jb(1, i) = (ry-ry0)/eps;
        }
        return Eigen::Vector2d(rx0 - ix, ry0 - iy);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class GraphGNSO3General : public GraphGaussNewtonOptimizer<double, 3, ConSO3WithErrorFunc>
{
public:
    using Con3DVec = std::vector<ConSO3WithErrorFunc>;
    using RVVec = std::vector<RotVec>;
    using Node = Eigen::Matrix<double, 3, 1>;
    using Jac  = Eigen::Matrix<double, Eigen::Dynamic, 3>;
    using NodeList = Eigen::Matrix<double, 3, Eigen::Dynamic>;

    static Eigen::VectorXd calcError2(const Node &ina, const Node &inb, const ConSO3WithErrorFunc &con, Jac &Ja, Jac &Jb)
    {
        return con.error_func(ina, inb, con, Ja, Jb);
    }

    static NodeList updatefunc(const NodeList &n, const Eigen::Matrix<double, Eigen::Dynamic, 1> &d)
    {
        NodeList out(3, n.cols());
        for(unsigned int i=0;i<n.cols();i++) {
            RotVec p(n.col(i));
            out.col(i) = p.concat(RotVec(d.segment(i*3, 3))).vec();
        }
        return out;
    }

public:
    GraphGNSO3General(){
        error_func = &GraphGNSO3General::calcError2;
        update_func = &GraphGNSO3General::updatefunc;
        //oplus_func = &RotVec::concat_rv;
    }
};

} // opt4imu
