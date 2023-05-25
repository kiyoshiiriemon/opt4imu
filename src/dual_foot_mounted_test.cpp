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

std::vector<Imu> movingAverageAccGyro(const std::vector<Imu> &data, int window_size) {
    std::vector<Imu> ret;
    ret.reserve(data.size());
    for (int i = 0; i < data.size(); ++i) {
        Imu imu = data[i];
        if (i >= window_size) {
            imu.linear_acceleration = Eigen::Vector3d::Zero();
            imu.angular_velocity = Eigen::Vector3d::Zero();
            for (int j = 0; j < window_size; ++j) {
                imu.linear_acceleration += data[i - j].linear_acceleration;
                imu.angular_velocity += data[i - j].angular_velocity;
            }
            imu.linear_acceleration /= window_size;
            imu.angular_velocity /= window_size;
        }
        ret.push_back(imu);
    }
    return ret;
}

// detect stationary from acc and gyro
static void detect_stationary(const IMULog &imulog, std::vector<int> &out_stationary, Con3DList &out_con_vel)
{
    out_con_vel.clear();
    out_con_vel.clear();
    auto imu_movavg = movingAverageAccGyro(imulog.data, 50);

    std::ofstream ofs("stationary.txt");
    for (int i=0; i<imu_movavg.size(); ++i) {
        const auto &imu = imu_movavg[i];
        bool is_stationary = false;
        std::cout << imu.linear_acceleration << std::endl;
        if (imu.angular_velocity.norm() < 0.5 && fabs(imu.linear_acceleration.norm()-9.8) < 0.5) {
            is_stationary = true;
        } else {
            is_stationary = false;
        }
        out_stationary.push_back(is_stationary);
        if (is_stationary) {
            Con3D con;
            if (i > 0) {
                con.id1 = 0;
                con.id2 = i;
                con.t = Vec3D::Zero();
                con.info = Mat3D::Identity();
                out_con_vel.push_back(con);
            }
        }
        ofs << is_stationary << std::endl;
    }
}

ImuTrjOptimizer::Trj9D process_single_feet(const char *fname)
{
    IMULog imulog;
    imulog.loadImuFile(fname, -1, -1);
    if (imulog.data[0].header.stamp.sec == 0)
        imulog.setUniformTimeStamp(2*1000);
    std::cout << "load " << imulog.data.size() << std::endl;
    std::vector<int> stationary;
    Con3DList con_vel;
    detect_stationary(imulog, stationary, con_vel);

    auto t0 = std::chrono::high_resolution_clock::now();
    ImuTrjOptimizer opt;
    opt.setImuLog(imulog);
    opt.setInitAtt(imulog.data[0].orientation);
    opt.estimateByIntegration();
    opt.correctAttByGravity(stationary);
    opt.correctVel(con_vel);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast< std::chrono::microseconds> (t1-t0);
    std::cout << " time: " << elapsed.count()*1e-6 << "s" << std::endl;
    opt.writeTrj("est.txt");
    return opt.trj;
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
        std::cout << argv[0] << "<imufile_right> <imufile_left>" << std::endl;
        return 0;
    }
    Con3DList con_vel;
    ImuTrjOptimizer::Trj9D trj_r = process_single_feet(argv[1]);
    ImuTrjOptimizer::Trj9D trj_l = process_single_feet(argv[2]);
    return 0;
}
