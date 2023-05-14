/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include <imutrj_optimizer.hpp>
#include "graph_gn_so3.hpp"
#include "graph_gn_so3_general.hpp"
#include "graph_gn_lin3.hpp"

using namespace opt4imu;

ImuTrjOptimizer::ImuTrjOptimizer() {

}

void ImuTrjOptimizer::correctBySagawa() {
    size_t n = imulog.size();
    Vec3D diff = RotVec(-att(0)).concat(RotVec(att(n-1))).vec();
    // correct att. error
    for(int i=1; i<n; ++i) {
        att(i) = RotVec(att(i)).concat(RotVec(-(diff/n) * i)).vec();
    }
    // recalculate vel
    for(int i=0; i<n-1; ++i) {
        double dt = imulog.duration(i, i+1);
        Vec3D accrot = RotVec(att(i)).toRotationMatrix() * imulog.acc(i);
        accrot(2) = accrot(2) - gravity;
        vel(i+1) = vel(i) + accrot * dt;
    }
    // correct vel. error
    Vec3D dv = vel(n-1) - vel(0);
    for(int i=1; i<n; ++i) {
        vel(i) = vel(i) - (dv/n) * i;
    }
    // recalculate pos
    for(int i=0; i<n-1; ++i) {
        double dt = imulog.duration(i, i+1);
        pos(i + 1) = pos(i) + vel(i) * dt;
    }
    // correct att. error
    Vec3D dp = pos(n-1) - pos(0);
    for(int i=1; i<n; ++i) {
        pos(i) = pos(i) - (dp/n) * i;
    }
}

void ImuTrjOptimizer::estimateByIntegration() {
    for(int i=0; i<imulog.size()-1; ++i) {
        double dt = imulog.duration(i, i+1);
        att(i+1) = RotVec::concat_rv(att(i), imulog.gyro(i) * dt);
        Vec3D accrot = RotVec(att(i)).toRotationMatrix() * imulog.acc(i);
        accrot(2) = accrot(2) - gravity;
        vel(i+1) = vel(i) + accrot * dt;
        pos(i+1) = pos(i) + vel(i) * dt;
        //std::cout << att(i).transpose() << "\t" << vel(i).transpose() << "\t" << pos(i).transpose() << std::endl;
    }
    std::cout << "naive integration final pos: " << pos(trj.cols()-1).transpose() << std::endl;
}

void ImuTrjOptimizer::initializeByStill(double dur) {
    trj = Trj9D::Zero(trj.rows(), imulog.size());
    if (dur > 0) {
        Vec3D gyro_bias = imulog.getInitialMeanGyro(dur);
#if 0
        Vec3D acc_bias(0.2281, -0.1486, -0.6369);
        Vec3D gyro_gain(1-0.1999/100., 1-4.1312/100., 1+0.6296/100.);
        Vec3D acc_gain(1+0.7928/100., 1+7.708/100., 1-4.8846/100.);
        imulog.applyCalibrationParams(acc_bias, gyro_bias, acc_gain, gyro_gain);
#endif
        imulog.applyCalibrationParams(Vec3D::Zero(), gyro_bias, Vec3D::Ones(), Vec3D::Ones());

        Vec3D acc0 = imulog.getInitialMeanAcc(dur);
        Vec3D att0 = calcTilt(0, acc0).vec();
        std::cout << "gyro    init: " << gyro_bias.transpose() << std::endl;
        std::cout << "acc     init: " << acc0.transpose() << std::endl;
        std::cout << "initial tilt: " << att0.transpose() << std::endl;
        att(0) = att0;
    }
}

void ImuTrjOptimizer::writeTrj(const char *fname) {
    std::ofstream ofs(fname);
    if (!ofs) return;
    ofs << trj.transpose() << std::endl;
}

void ImuTrjOptimizer::correctByPosloop(const Con3DList &conlist) {
    int n = imulog.size();
    GraphGNSO3::Con3DVec vec;
    for(int i=0; i<n-1; ++i) {
       ConSO3 con;
       con.id1 = i;
       con.id2 = i+1;
       con.t = imulog.gyro(i) * imulog.duration(i, i+1);
       //con.t = RotVec::concat_rv(-att(i), att(i+1));
       con.info = Mat3D::Identity();
       vec.push_back(con);
    }
    vec.push_back(ConSO3{0, n-1, Vec3D::Zero(), Mat3D::Identity()});

    GraphGNSO3 so3opt;
    so3opt.verbose = true;
    GraphGNSO3::NodeList newatt = so3opt.optimize(att(), vec);
    att() = newatt;
    reintegrateVel();

    GraphGNLin3::Con3DVec vcvec;
    for(int i=0; i<1000; ++i) {
        ConLin3D con;
        con.id1 = 0;
        con.id2 = i;
        con.t = vel(i);
        con.info = Mat3D::Identity() * 1/(1000 * 0.002 * (i+1))/(1000 * 0.002 * (i+1));
        vcvec.push_back(con);
    }
    for(int i=0; i<n-1; ++i) {
        ConLin3D con;
        con.id1 = i;
        con.id2 = i+1;
        con.t = vel(i+1)-vel(i);
        con.info = Mat3D::Identity();
        vcvec.push_back(con);
    }
    vcvec.push_back(ConLin3D{0, n-1, Vec3D::Zero(), Mat3D::Identity()});
    GraphGNLin3 posopt;
    posopt.verbose = true;
    GraphGNLin3::NodeList newvel = posopt.optimize(vel(), vcvec);
    vel() = newvel;
    reintegratePos();

    GraphGNLin3::Con3DVec pcvec;
    for(const auto &c : conlist) {
        pcvec.push_back(ConLin3D{c.id1, c.id2, c.t, c.info});
    }
    for(int i=0; i<n-1; ++i) {
        ConLin3D con;
        con.id1 = i;
        con.id2 = i+1;
        con.t = pos(i+1)-pos(i);
        con.info = Mat3D::Identity();
        pcvec.push_back(con);
    }
    for(int i=0; i<1000; ++i) {
        ConLin3D con;
        con.id1 = 0;
        con.id2 = i;
        con.t = pos(i);
        con.info = Mat3D::Identity() * 1/(1000 * 0.002 * (i+1))/(1000 * 0.002 * (i+1));
        pcvec.push_back(con);
    }
    pcvec.push_back(ConLin3D{0, n-1, Vec3D::Zero(), Mat3D::Identity()});
    GraphGNLin3::NodeList newpos = posopt.optimize(pos(), pcvec);
    pos() = newpos;
}

void ImuTrjOptimizer::correctVel(const Con3DList &conlist) {
    int n = imulog.size();
    GraphGNLin3::Con3DVec vcvec;
#if 0
    for(int i=0; i<1000; ++i) {
        ConLin3D con;
        con.id1 = 0;
        con.id2 = i;
        con.t = vel(i);
        con.info = Mat3D::Identity() * 1/(1000 * 0.002 * (i+1))/(1000 * 0.002 * (i+1));
        vcvec.push_back(con);
    }
#endif
    for(int i=0; i<n-1; ++i) {
        ConLin3D con;
        con.id1 = i;
        con.id2 = i+1;
        con.t = vel(i+1)-vel(i);
        con.info = Mat3D::Identity();
        vcvec.push_back(con);
    }
    for(const auto &c : conlist) {
        ConLin3D con;
        if (c.id1 < n && c.id2 < n) {
            con.id1 = c.id1;
            con.id2 = c.id2;
            con.t = c.t;
            con.info = c.info;
            vcvec.push_back(con);
        }
    }
    GraphGNLin3 posopt;
    posopt.verbose = true;
    GraphGNLin3::NodeList newvel = posopt.optimize(vel(), vcvec);
    vel() = newvel;
    reintegratePos();
}

void ImuTrjOptimizer::correctVelCombined(const Con3DList &conlist, ImuTrjOptimizer &opt2, const Con3DList &conlist2)
{
    int n = imulog.size();
    GraphGNLin3::Con3DVec vcvec;
    for(int i=0; i<n-1; ++i) {
        ConLin3D con;
        con.id1 = i;
        con.id2 = i+1;
        con.t = vel(i+1)-vel(i);
        con.info = Mat3D::Identity();
        vcvec.push_back(con);
    }
    for(const auto &c : conlist) {
        ConLin3D con;
        if (c.id1 < n && c.id2 < n) {
            con.id1 = c.id1;
            con.id2 = c.id2;
            con.t = c.t;
            con.info = c.info;
            vcvec.push_back(con);
        }
    }
    for(int i=0; i<n-1; ++i) {
        ConLin3D con;
        con.id1 = i + n;
        con.id2 = i+1+n;
        con.t = opt2.vel(i+1)-opt2.vel(i);
        con.info = Mat3D::Identity();
        vcvec.push_back(con);
    }
    for(const auto &c : conlist2) {
        ConLin3D con;
        if (c.id1 < n && c.id2 < n) {
            con.id1 = c.id1 + n;
            con.id2 = c.id2 + n;
            con.t = c.t;
            con.info = c.info;
            vcvec.push_back(con);
        }
    }
    for(int i=0; i < n; ++i) {
        if (imulog.gyro(i).norm() > 0.05) continue;
        ConLin3D con;
        con.id1 = i;
        con.id2 = i + n;
        con.t = Vec3D::Zero();
        con.info = Mat3D::Identity() * 0.00001;
        vcvec.push_back(con);
    }
    std::cout << "no. constraints: " << vcvec.size() << std::endl;
    GraphGNLin3::NodeList vel_combined(3, imulog.size()*2);
    vel_combined << vel(), opt2.vel();
    GraphGNLin3 posopt;
    posopt.verbose = true;
    GraphGNLin3::NodeList newvel_combined = posopt.optimize(vel_combined, vcvec);

    std::cout << "newvel: " << newvel_combined.rows() << " " << newvel_combined.cols() << std::endl;
    vel() = newvel_combined.block(0, 3, n, 3);
    reintegratePos();
}

void ImuTrjOptimizer::correctAttByGravity(const std::vector<int> &stationary)
{
    using namespace std::placeholders;
    int n = imulog.size();
    GraphGNSO3General attopt;
    attopt.verbose = true;
    GraphGNSO3General::Con3DVec cvec;
    for(int i=0; i<n-1; ++i) {
        ConSO3WithErrorFunc con;
        con.id1 = i;
        con.id2 = i+1;
        //con.t = (imulog.gyro(i) + imulog.gyro(i+1))* imulog.duration(i, i+1) / 2; // may be better
        con.t = imulog.gyro(i) * imulog.duration(i, i+1);
        //auto rel = RotVec::concat_rv(-att(i), att(i+1));
        //std::cout << "rel :" << con.t.transpose() << " | " << rel.transpose() << std::endl;
        //std::cout << "err :" << RotVec::concat_rv(-con.t, rel) << std::endl;
        con.info = Eigen::Matrix<double, -1, -1>(3, 3);
        con.info << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        con.error_func = ConSO3WithErrorFunc::calcErrorRelSO3;
        cvec.emplace_back(con);
    }
    for(int i : stationary ) {
        ConSO3WithErrorFunc con;
        con.id1 = 0;
        con.id2 = i;
        if (i >= n) continue;
        con.t = calcTilt(0, imulog.acc(i)).vec(); // TODO: low pass
        //std::cout << "calctilt " << imulog.acc(i).transpose() << " " << con.t.transpose() << std::endl;
        con.info = Eigen::Matrix<double, -1, -1>(2, 2);
        //con.info << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        con.info << 1, 0, 0, 1;
        con.info *= 1e-4;
        con.error_func = ConSO3WithErrorFunc::calcErrorGravity;
        cvec.emplace_back(con);
    }
    attopt.max_steps = 20;
    GraphGNSO3::NodeList newatt = attopt.optimize(att(), cvec);
    att() = newatt.block(0, 0, n, 3);
    reintegrateVel();
    reintegratePos();
}

void ImuTrjOptimizer::reintegrateVel() {
    for(int i=0; i<imulog.size()-1; ++i) {
        double dt = imulog.duration(i, i+1);
        Vec3D accrot = RotVec(att(i)).toRotationMatrix() * imulog.acc(i);
        accrot(2) = accrot(2) - gravity;
        vel(i+1) = vel(i) + accrot * dt;
    }
}

void ImuTrjOptimizer::reintegratePos() {
    for(int i=0; i<imulog.size()-1; ++i) {
        double dt = imulog.duration(i, i+1);
        //std::cout << "dt: " << dt << std::endl;
        pos(i+1) = pos(i) + vel(i) * dt;
    }
}

void ImuTrjOptimizer::setInitAtt(const RotVec &rv0) {
    att(0) = rv0.vec();
}
