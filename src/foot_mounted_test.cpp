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

void process_log(const char *fname, const std::vector<int> &stationary, const Con3DList &con_vel)
{
    IMULog imulog;
    imulog.loadImuFile(fname, -1, -1);
    if (imulog.data[0].header.stamp.sec == 0)
        imulog.setUniformTimeStamp(2*1000);
    std::cout << "load " << imulog.data.size() << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    ImuTrjOptimizer opt;
    opt.setImuLog(imulog);
    opt.setInitAtt(imulog.data[0].orientation);
    opt.estimateByIntegration();
    opt.writeTrj("naive.txt");
    opt.correctAttByGravity(stationary);
    opt.correctVel(con_vel);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast< std::chrono::microseconds> (t1-t0);
    std::cout << " time: " << elapsed.count()*1e-6 << "s" << std::endl;
    opt.writeTrj("est_foot_mount.txt");
}

static std::vector<int> load_stationary_flags(const char *fname, Con3DList &con_vel)
{
    std::vector<int> ret;
    std::ifstream is(fname);

    int id=0;
    while(is){
        int is_stationary;
        char buf[1024];
        is.getline(buf,1024);
        if (!is) break;
        std::istringstream sstrm(buf);
        sstrm >> is_stationary;
        if (is_stationary) {
            Con3D con;
            ret.push_back(is_stationary);
            if (id > 0) {
                con.id1 = 0;
                con.id2 = id;
                con.t = Vec3D::Zero();
                con.info = Mat3D::Identity();
                con_vel.push_back(con);
            }
        }
        id++;
    }
    return ret;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << argv[0] << "<imufile> <stationary>" << std::endl;
        return 0;
    }
    Con3DList con_vel;
    std::vector<int> stationary_ind = load_stationary_flags(argv[2], con_vel);
    process_log(argv[1], stationary_ind, con_vel);
    return 0;
}
