/****************************************************************************
 * Copyright (C) 2023-2024 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include <random>
#include <fstream>
#include <Eigen/Dense>
#include "geom3d.hpp"

class GaussNewtonOptimizer
{
    // calculate least-squares
    // E = err(T)'*err(T);
public:
    using State = Eigen::VectorXd;
    using Jac = Eigen::MatrixXd;
    using ErrFunc = std::function<Eigen::VectorXd(size_t, const State&, Jac&)>;
    using ErrFuncNoJac = std::function<Eigen::VectorXd(size_t, const State&)>;
    using OplusFunc = std::function<State(const State&, const State&)>;

    int num = 0;
    int max_steps = 1000;
    double stop_thre = 1e-9;
    bool verbose = false;
    ErrFunc error_func;
    ErrFuncNoJac error_func_nojac;
    OplusFunc oplus_func;

    double optimizeOneStep( State &s_out, const State &s_in)
    {
        double ret = 0;
        size_t dim = s_in.size();
        Eigen::MatrixXd mat(dim, dim);
        Eigen::VectorXd bf(dim);
        mat.setZero();
        bf.setZero();
        Jac jac;
        for (size_t i=0; i<num; ++i) {
            Eigen::VectorXd err;
            if(error_func) {
                err = error_func(i, s_in, jac);
                mat.noalias() += jac.transpose() * jac;
                bf.noalias()  += jac.transpose() * err;
#if 0 // debug jecobian matrix
                std::cout << "=== jac, numjac ===" << std::endl;
                std::cout << jac << std::endl;
                double eps = 1e-4;
                err = error_func_nojac(i, s_in);
                jac.resize(err.rows(), dim);
                for(size_t j=0; j<dim; ++j) {
                    State d(dim); d.setZero();
                    d(j) += eps;
                    auto s = oplus_func(s_in, d);
                    opt4imu::Vec3D v = (error_func_nojac(i, s) - err)/eps;
                    jac.block(0,j,err.rows(),1) = v;
                }
                std::cout << jac << std::endl;
                mat.noalias() += jac.transpose() * jac;
                bf.noalias()  += jac.transpose() * err;
#endif
            } else {
                double eps = 1e-4;
                err = error_func_nojac(i, s_in);
                //std::cout << err << std::endl;
                jac.resize(err.rows(), dim);
                for(size_t j=0; j<dim; ++j) {
                    State d(dim); d.setZero();
                    d(j) += eps;
                    auto s = oplus_func(s_in, d);
                    opt4imu::Vec3D v = (error_func_nojac(i, s) - err)/eps;
                    //std::cout << v << std::endl;
                    jac.block(0,j,err.rows(),1) = v;
                }
                mat.noalias() += jac.transpose() * jac;
                bf.noalias()  += jac.transpose() * err;
            }
            ret += err.transpose()*err;
        }
        Eigen::LDLT<Eigen::MatrixXd> solver;
        solver.compute(mat);
        State x = solver.solve(-bf);

        if (oplus_func) {
            s_out = oplus_func(s_in, x);
        } else {
            s_out = x + s_in;
        }
        return ret;
    }

    State optimize(const State &init, size_t in_num)
    {
        num = in_num;
        State xret, x = init;
        double prevres = std::numeric_limits<double>::max();
        for(int i=1; i<=max_steps; i++) {
            auto t0 = std::chrono::high_resolution_clock::now();
            double res = optimizeOneStep(xret, x);
            auto t1 = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast< std::chrono::microseconds> (t1-t0);
            double convergence = (xret - x).norm() / x.norm();
            if (verbose)
                std::cout << "step " << i << ": " << x.transpose() << ", convergence=" << res << " time: " << elapsed.count()*1e-6 << "s" << std::endl;
            x = xret;

            if (i >= 2 && convergence < stop_thre) {
                if (verbose)
                    std::cout << "converged: " << prevres - res << " < " << stop_thre << std::endl;
                break;
            }
            prevres = res;
        }
        return x;
    }
};

auto generate_points(size_t npoints)
{
    Eigen::Matrix3Xd ret(3, npoints);
    std::default_random_engine engine(0);
    std::normal_distribution<> dist(0.0, 1.0);
    for (size_t i = 0; i < npoints; ++i) {
        ret(0, i) = dist(engine);
        ret(1, i) = dist(engine);
        ret(2, i) = dist(engine);
    }
    return ret;
}

Eigen::Matrix3Xd transform_points(const Eigen::Matrix3Xd &points, const opt4imu::Pose3D &trans)
{
    Eigen::Matrix4Xd ext_points(4, points.cols());
    ext_points.topRows(3) = points;
    ext_points.row(3).setOnes();
    Eigen::Matrix4d T = trans.toMatrix();
    Eigen::Matrix4Xd trans_points = T * ext_points;
    return trans_points.topRows(3);
}

void register_points(const Eigen::Matrix3Xd &org, const Eigen::Matrix3Xd &dest)
{
    GaussNewtonOptimizer optim;
    opt4imu::Mat3D Rx, Ry, Rz;
    Rx << 0, 0, 0, 0, 0,-1, 0, 1, 0;
    Ry << 0, 0, 1, 0, 0, 0,-1, 0, 0;
    Rz << 0,-1, 0, 1, 0, 0, 0, 0, 0;
    optim.error_func_nojac = [=](size_t i, const Eigen::VectorXd &state)->opt4imu::Vec3D {
        Eigen::Vector3d t(state.head(3));
        opt4imu::RotVec rv(state.tail(3));
        auto p = org.col(i);
        opt4imu::Vec3D ret = rv.toRotationMatrix()*p + t - dest.col(i);
        return ret;
    };
#if 1
    optim.error_func = [=](size_t i, const Eigen::VectorXd &state, Eigen::MatrixXd &jac)->opt4imu::Vec3D {
        Eigen::Vector3d t(state.head(3));
        opt4imu::RotVec rv(state.tail(3));
        auto R0 = rv.toRotationMatrix();
        auto p = org.col(i);
        jac.resize(3,6);
        jac << R0, R0*Rx*p, R0*Ry*p, R0*Rz*p;
        opt4imu::Vec3D ret = R0 * p + t - dest.col(i);
        return ret;
    };
#endif
    optim.oplus_func = [](const Eigen::VectorXd &l, const Eigen::VectorXd &r) {
        auto p = opt4imu::Pose3D(l).oplus(opt4imu::Pose3D(r)).vec();
        return p;
    };
    optim.verbose = true;
    auto ret = optim.optimize(opt4imu::Vec6D(0,0,0, 0, 0, 0), org.cols());
    std::cout << "result: " << ret.transpose() << std::endl;
}

void test_registration()
{
    Eigen::Matrix3Xd org = generate_points(10000);
    opt4imu::Pose3D p3d(1.7, -2.2, 3.57, -0.24, 0.05, 0.18);
    Eigen::Matrix3Xd trans = transform_points(org, p3d) + 0.1*Eigen::Matrix3Xd::Random(3, org.cols());
    //std::cout << org << std::endl;
    //std::cout << trans << std::endl;
    register_points(org, trans);
}

int main() {
    test_registration();
    return 0;
}
