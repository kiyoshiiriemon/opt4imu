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

using ConSO3 = ConstraintND<double, 3>;

class GraphGNSO3 : public GraphGaussNewtonOptimizer<double, 3, ConSO3>
{
public:
    using Con3DVec = std::vector<ConSO3>;
    using RVVec = std::vector<opt4imu::RotVec>;
    using Node = Eigen::Matrix<double, 3, 1>;
    //using Jac  = Eigen::Matrix<double, 3, 3>;
    using Jac  = Eigen::Matrix<double, -1, 3>;
    using NodeList = Eigen::Matrix<double, 3, Eigen::Dynamic>;

    static opt4imu::Vec3D errorFunc(const opt4imu::RotVec &pa, const opt4imu::RotVec &pb, const opt4imu::RotVec &con)
    {
        opt4imu::RotVec err = con.inverted().concat(pa.inverted()).concat(pb);
        return err.vec();
    }

    static Eigen::VectorXd calcError2(const Node &ina, const Node &inb, const ConSO3 &con, Jac &Ja, Jac &Jb)
    {
        using Mat3D = Eigen::Matrix<double, 3, 3>;
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
        opt4imu::RotVec  pa(ina[0], ina[1], ina[2]);
        opt4imu::RotVec  pb(inb[0], inb[1], inb[2]);
        opt4imu::RotVec rel(con.t[0], con.t[1], con.t[2]);

        Eigen::VectorXd e0 = errorFunc(pa, pb, rel);
        Mat3D RaInv = pa.toRotationMatrix().transpose();
        Mat3D Rb = pb.toRotationMatrix();
        Mat3D RcInv = rel.toRotationMatrix().transpose();
        opt4imu::RotVec rve(e0.tail<3>());
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
        Ja = opt4imu::jacobianR(Re) * dRedRa;

        Jb = rightJacobianInvRV(rve);
        return e0;
    }

    static NodeList updatefunc(const NodeList &n, const Eigen::Matrix<double, Eigen::Dynamic, 1> &d)
    {
        NodeList out(3, n.cols());
        for(unsigned int i=0;i<n.cols();i++) {
            opt4imu::RotVec p(n.col(i));
            out.col(i) = p.concat(opt4imu::RotVec(d.segment(i*3, 3))).vec();
        }
        return out;
    }

public:
    GraphGNSO3(){
        error_func = &GraphGNSO3::calcError2;
        update_func = &GraphGNSO3::updatefunc;
        //oplus_func = &RotVec::concat_rv;
    }
};
