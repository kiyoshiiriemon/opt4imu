/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include <graph_gn_so3.hpp>
#include <graph_gn_lin3.hpp>
#include <imulog.hpp>
#include <imutrj_optimizer.hpp>

using namespace opt4imu;

static RotVec calcTilt(double yaw, const Eigen::Vector3d &acc)
{
    double pitch = asin(acc(0) / acc.norm());
    double roll = atan2(-acc(1), acc(2));
    return RotVec::fromXYZ(roll, pitch, yaw).inverted();
}

void process_log(const char *fname, const Con3DList &conlist)
{
    IMULog imulog;
    imulog.loadImuFile(fname, -1, -1);
    if (imulog.data[0].header.stamp.sec == 0)
        imulog.setUniformTimeStamp(2*1000);
    std::cout << "load " << imulog.data.size() << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    ImuTrjOptimizer opt;
    opt.setImuLog(imulog);
    opt.estimateByIntegration();
    //opt.writeTrj("naive.txt");
    opt.correctByPosloop(conlist);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast< std::chrono::microseconds> (t1-t0);
    std::cout << " time: " << elapsed.count()*1e-6 << "s" << std::endl;
    opt.writeTrj("posloop_single.txt");
}

static bool loadLoops(const char *fname, Con3DList &list)
{
    std::ifstream is(fname);
    if (!is) return false;

    while(is){
        Con3D con;
        char buf[1024];
        is.getline(buf,1024);
        if (!is) return true;
        std::istringstream sstrm(buf);
        std::string tag;
        sstrm >> con.id1 >> con.id2;
        con.t = Vec3D::Zero();
        con.info = Mat3D::Identity();
        list.push_back(con);
    }
    return true;
}

int main(int argc, char *argv[])
{
    if (argc > 1) {
        Con3DList conlist;
        if (argc > 2) {
           loadLoops(argv[2], conlist);
        }
        process_log(argv[1], conlist);
    }
    return 0;
}
