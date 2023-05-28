/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include <graph_gn_so3.hpp>
#include <graph_gn_se3.hpp>
#include <imulog.hpp>
#include <imutrj_optimizer.hpp>

using namespace opt4imu;

static Eigen::VectorXd extract_upper_right(const Mat6D &info)
{
    Eigen::VectorXd ret(21);
    int idx = 0;
    for (int i=0; i<6; ++i) {
        for (int j=i; j<6; ++j) {
            ret[idx++] = info(i, j);
        }
    }
    return ret;
}

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
        //std::cout << imu.linear_acceleration << std::endl;
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

ImuTrjOptimizer process_single_feet(const char *fname)
{
    IMULog imulog;
    if (!imulog.loadImuFile(fname, -1, -1)) {
        std::cerr << "failed to load " << fname << std::endl;
        exit(1);
    }
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
    return opt;
}

NodeList process_dual_foot_mounted(ImuTrjOptimizer &opt_r, ImuTrjOptimizer &opt_l)
{
    Pose3DVec p3dnodes;
    Con6DVec constraints;
    size_t len_r = opt_r.trj.cols();
    size_t len_l = opt_l.trj.cols();
    std::cout << "len: " << len_r << "  " << len_l << std::endl;
    size_t len_min = std::min<size_t>(len_r, len_l);

    p3dnodes.reserve(len_r+len_l);
    p3dnodes.emplace_back(0, 0, 0, 0, 0, 0); // node id=0
    std::ofstream ofs("dual_foot.g2o");
    for(int i=0; i<len_min; ++i) {
        Pose3D node(opt_r.pos(i), RotVec(opt_r.att(i)));
        p3dnodes.push_back(node);
        Eigen::Quaterniond q = node.rv.toQuaternion();
        ofs << "VERTEX_SE3:QUAT " << i << " " << node.x << " " << node.y << " " << node.z << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << std::endl;
    }
    for(int i=0; i<len_min; ++i) {
        Pose3D node(opt_l.pos(i), RotVec(opt_l.att(i)));
        p3dnodes.push_back(node);
        Eigen::Quaterniond q = node.rv.toQuaternion();
        ofs << "VERTEX_SE3:QUAT " << i+len_min << " " << node.x << " " << node.y << " " << node.z << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << std::endl;
    }

    for(int i=1; i<len_min; ++i) {
        Con6D con;
        con.id1 = i;
        con.id2 = i + 1;
        con.t = p3dnodes[con.id2].ominus(p3dnodes[con.id1]).vec();
        con.info = Mat6D::Identity();
        constraints.push_back(con);

        con.id1 = i + len_min;
        con.id2 = i + 1 + len_min;
        con.t = p3dnodes[con.id2].ominus(p3dnodes[con.id1]).vec();
        con.info = Mat6D::Identity();
        constraints.push_back(con);
    }
    // left-right relative pose constraint
    auto rel_r2l = Pose3D(0, -0.15, 0, p3dnodes[len_min].rv).ominus(Pose3D(0, 0, 0, p3dnodes[1].rv)).vec();
    for(int i=1; i<len_min+1; ++i) {
        Con6D con;
        con.id1 = i;
        con.id2 = i+len_min;
        con.t = rel_r2l;
        con.info = Mat6D::Identity()*1e-9;
        constraints.push_back(con);
    }
    {
        // initial pose constraint
        Con6D con;
        con.id1 = 0;
        con.id2 = 1;
        con.t = Pose3D(0, 0, 0, p3dnodes[1].rv).vec();
        con.info = Mat6D::Identity()*1e-9;
    }
    for(const auto &c : constraints) {
        auto q = RotVec(c.t.tail(3)).toQuaternion();
        ofs << "EDGE_SE3:QUAT " << c.id1 << " " << c.id2 << " ";
        ofs << c.t.head(3).transpose() << " ";
        ofs << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << " ";
        ofs << extract_upper_right(c.info).transpose() << std::endl;
    }
    GraphGNSE3 opt;
    opt.verbose = true;
    NodeList nodes(6, p3dnodes.size());
    for(int i=0; i<nodes.cols(); ++i) {
        nodes.col(i) = p3dnodes[i].vec();
    }
    opt.max_steps = 10;
    NodeList result = opt.optimize(nodes, constraints);
    return result;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << argv[0] << "<imufile_right> <imufile_left>" << std::endl;
        return 0;
    }
    Con3DList con_vel;
    ImuTrjOptimizer opt_r = process_single_feet(argv[1]);
    ImuTrjOptimizer opt_l = process_single_feet(argv[2]);
    auto nodelist = process_dual_foot_mounted(opt_r, opt_l);
    {
        std::ofstream ofs_r("est_dual_foot_r.txt");
        std::ofstream ofs_l("est_dual_foot_l.txt");
        size_t len = (nodelist.cols() - 1) / 2;
        ofs_r << nodelist(Eigen::all, Eigen::seq(1, len)).transpose() << std::endl;
        ofs_l << nodelist(Eigen::all, Eigen::seq(len + 1, len * 2)).transpose() << std::endl;
    }

    return 0;
}
