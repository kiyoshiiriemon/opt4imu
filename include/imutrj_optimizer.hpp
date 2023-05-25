/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include "imulog.hpp"
#include <Eigen/Dense>

namespace opt4imu {

struct Con3D {
   int id1; // from
   int id2; // to
   Vec3D t;
   Mat3D info;
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
using Con3DList = std::vector<Con3D>;

class ImuTrjOptimizer
{
public:
    using Trj9D = Eigen::Matrix<double, 9, Eigen::Dynamic>;
    ImuTrjOptimizer();
    void setImuLog(const IMULog &in_log, double init_dur = -1) {
        imulog = in_log;
        initializeByStill(init_dur);
    }
    void setInitAtt(const RotVec &rv0);
    void initializeByStill(double dur);
    void correctByPosloop(const Con3DList &conlist);
    void correctVel(const Con3DList &conlist);
    void correctAttByGravity(const std::vector<int> &stationary);
    void correctVelCombined(const Con3DList &conlist, ImuTrjOptimizer &opt2, const Con3DList &conlist2);
    void estimateByIntegration();
    void correctBySagawa();
    void reintegrateVel();
    void reintegratePos();
    void writeTrj(const char *fname);

    static RotVec calcTilt(double yaw, const Eigen::Vector3d &acc)
    {
        double pitch = asin(acc(0) / acc.norm());
        double roll = atan2(-acc(1), acc(2));
        return RotVec::fromXYZ(roll, pitch, yaw).inverted();
    }

    Trj9D trj;
    IMULog imulog;
    double gravity = 9.8021;
    auto att() { return trj.block<3,-1>(0,0, 3, trj.cols()); }
    auto vel() { return trj.block<3,-1>(3,0, 3, trj.cols()); }
    auto pos() { return trj.block<3,-1>(6,0, 3, trj.cols()); }
    auto att(int i) { return trj.block<3,1>(0,i); }
    auto vel(int i) { return trj.block<3,1>(3,i); }
    auto pos(int i) { return trj.block<3,1>(6,i); }
};

} //opt4imu

