/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace opt4imu
{
// === 2D geometry utility classes and functions ===

using Vec3D = Eigen::Matrix<double,3,1>;
using Mat3D = Eigen::Matrix<double,3,3>;

static inline double normalize_rad_pi_mpi(double rad)
{
    #ifndef M_PI
    static const double M_PI = 3.141592653589793238;
    #endif
    double val = fmod(rad, 2*M_PI);
    if (val > M_PI) {
        val -= 2*M_PI;
    } else if ( val < -M_PI) {
        val += 2*M_PI;
    }
    return val;
}

// === 3D utility classes and functions ===

using Vec6D = Eigen::Matrix<double,6,1>;
using Mat6D = Eigen::Matrix<double,6,6>;

struct RotVec
{
    double ax, ay, az;
    RotVec() : ax(0), ay(0), az(0) {}
    RotVec(const Vec3D &vec) : ax(vec(0)), ay(vec(1)), az(vec(2)) {}
    RotVec(double x, double y, double z) : ax(x), ay(y), az(z) {}
    RotVec(const Eigen::Quaterniond &in_q) {
        Eigen::Quaterniond q = in_q.normalized();
        if (q.w() < 0) q = Eigen::Quaterniond(-q.w(), -q.x(), -q.y(), -q.z());
        double x = q.x();
        double y = q.y();
        double z = q.z();
        double norm_im = sqrt(x*x+y*y+z*z);
        if (norm_im < 1e-7) {
            ax = 2*x;
            ay = 2*y;
            az = 2*z;
        } else {
            double th = 2 * atan2(norm_im, q.w());
            ax = x / norm_im * th;
            ay = y / norm_im * th;
            az = z / norm_im * th;
        }
    }
    static void rotmatFromXYZ(Mat3D &rot, double rotx, double roty, double rotz)
    {
        rot(0,0) = cos(roty)*cos(rotz);
        rot(0,1) = -cos(roty)*sin(rotz);
        rot(0,2) = sin(roty);
        rot(1,0) = sin(rotx)*sin(roty)*cos(rotz)+cos(rotx)*sin(rotz);
        rot(1,1) = -sin(rotx)*sin(roty)*sin(rotz)+cos(rotx)*cos(rotz);
        rot(1,2) = -sin(rotx)*cos(roty);
        rot(2,0) = -cos(rotx)*sin(roty)*cos(rotz)+sin(rotx)*sin(rotz);
        rot(2,1) = cos(rotx)*sin(roty)*sin(rotz)+sin(rotx)*cos(rotz);
        rot(2,2) = cos(rotx)*cos(roty);
    }
    static RotVec fromXYZ(double rotx, double roty, double rotz) {
        Mat3D rot;
        rotmatFromXYZ(rot, rotx, roty, rotz);
        return RotVec(Eigen::Quaterniond(rot));
    }
    double norm() const { return sqrt(ax*ax + ay*ay + az*az); }
    RotVec inverted() const { return RotVec(-ax, -ay, -az); }
    Eigen::Quaterniond toQuaternion() const {
        double v = sqrt(ax*ax + ay*ay + az*az);
        if (v < 1e-6) {
            return Eigen::Quaterniond(1, 0, 0, 0);
        }
        return Eigen::Quaterniond(cos(v/2), sin(v/2)*ax/v, sin(v/2)*ay/v, sin(v/2)*az/v);
    }
    Mat3D toRotationMatrix() const {
        return toQuaternion().toRotationMatrix();
    }
    RotVec concat(const RotVec &rv) const {
        return RotVec(toQuaternion() * rv.toQuaternion());
    }
    void getEulerZYX(double *rotx, double *roty, double *rotz) const {
        Mat3D rot = toRotationMatrix();
        double sin_pitch = -rot(2,0);
        if (sin_pitch <= -1.0f) {
            *roty = -M_PI / 2;
        } else if (sin_pitch >= 1.0f) {
            *roty = M_PI / 2;
        } else {
            *roty = asin(sin_pitch);
        }
        if (sin_pitch > 0.99999f) {
            *rotz = 0.0f;
            *rotx = atan2(rot(0, 1), rot(0, 2));//need to verify
        } else {
            *rotz  = atan2(rot(1, 0), rot(0, 0));
            *rotx   = atan2(rot(2, 1), rot(2, 2));
        }
    }
    Vec3D vec() const {
        Vec3D v;
        v << ax, ay, az;
        return v;
    }
    static Vec3D concat_rv(const Vec3D &a, const Vec3D &b)
    {
        return RotVec(a).concat(RotVec(b)).vec();
    }
};

inline std::ostream& operator<<(std::ostream& os, const RotVec &p)
{
    os << p.ax << ' ' << p.ay << ' ' << p.az << ' ';
    return os;
}

static inline Eigen::Matrix<double, 4, 3> dQuat_dRV(const RotVec &rv)
{
    Eigen::Matrix<double, 4, 3> dqu;
    double u1 = rv.ax, u2 = rv.ay, u3 = rv.az;
    if (rv.norm() < 1e-6) {
        dqu << -rv.ax, -rv.ay, -rv.az,
                    2,      0,      0,
                    0,      2,      0,
                    0,      0,      2;

        return dqu * 0.25;
    }
    double v = sqrt(u1*u1+u2*u2+u3*u3);
    double vd = v*2;
    double v2 = v*v;
    double v3 = v*v*v;

    double S = sin(v/2); double C = cos(v/2);
    dqu << -u1 * S/vd, -u2*S/vd, -u3*S/vd,
    S/v + u1*u1*C/(2*v2) - u1*u1*S/v3, u1*u2*(C/(2*v2)-S/v3), u1*u3*(C/(2*v2)-S/v3),
    u1*u2*(C/(2*v2)-S/v3), S/v+u2*u2*C/(2*v2)-u2*u2*S/v3, u2*u3*(C/(2*v2)-S/v3),
    u1*u3*(C/(2*v2)-S/v3), u2*u3*(C/(2*v2)-S/v3), S/v+u3*u3*C/(2*v2)-u3*u3*S/v3;

    return dqu;
}

static inline void dR_dRV(const RotVec &rv, Mat3D &dux, Mat3D &duy, Mat3D &duz)
{
    Eigen::Quaterniond q = rv.toQuaternion();
    double qw = q.w(), qx = q.x(), qy = q.y(), qz = q.z();
    Mat3D dRdqw, dRdqx, dRdqy, dRdqz;
    dRdqw << qw,-qz, qy,
             qz, qw,-qx,
            -qy, qx, qw; dRdqw *= 2;
    dRdqx << qx, qy, qz,
             qy,-qx,-qw,
             qz, qw,-qx; dRdqx *= 2;
    dRdqy <<-qy, qx, qw,
             qx, qy, qz,
            -qw, qz,-qy; dRdqy *= 2;
    dRdqz <<-qz,-qw, qx,
             qw,-qz, qy,
             qx, qy, qz; dRdqz *= 2;
    Eigen::Matrix<double, 4, 3> dqdu = dQuat_dRV(rv);
    dux = dRdqw * dqdu(0,0) + dRdqx * dqdu(1, 0) + dRdqy * dqdu(2, 0) + dRdqz * dqdu(3, 0);
    duy = dRdqw * dqdu(0,1) + dRdqx * dqdu(1, 1) + dRdqy * dqdu(2, 1) + dRdqz * dqdu(3, 1);
    duz = dRdqw * dqdu(0,2) + dRdqx * dqdu(1, 2) + dRdqy * dqdu(2, 2) + dRdqz * dqdu(3, 2);
}

static inline Eigen::Matrix<double, 3, 4> dRV_dQuat(Eigen::Quaterniond q)
{
    if (q.w() < 0) q = Eigen::Quaterniond(-q.w(), -q.x(), -q.y(), -q.z());
    if (1-q.w()*q.w() < 1e-7) {
        Eigen::Matrix<double, 3, 4> ret;
        ret << 0, 2, 0, 0,
               0, 0, 2, 0,
               0, 0, 0, 2;
        return ret;
    }
    double c = 1/(1-q.w()*q.w());
    double d = acos((double)q.w())/(sqrt(1-q.w()*q.w()));
    //std::cout << "dRV_dQuat: " << c << "\t" << d << std::endl;
    Eigen::Matrix<double, 3, 4> ret;
    ret << 2*c*q.x()*(d*q.w()-1),2*d,  0,  0,
           2*c*q.y()*(d*q.w()-1),  0,2*d,  0,
           2*c*q.z()*(d*q.w()-1),  0,  0,2*d;
    return ret;
}

static inline Eigen::Matrix<double,4,4> QMat(Eigen::Quaterniond q)
{
    if (q.w() < 0) q = Eigen::Quaterniond(-q.w(), -q.x(), -q.y(), -q.z());
    double qw = q.w(), qx = q.x(), qy = q.y(), qz = q.z();
    Eigen::Matrix<double,4,4> Q;
    Q << qw,-qx,-qy,-qz,
         qx, qw,-qz, qy,
         qy, qz, qw,-qx,
         qz,-qy, qx, qw;
    return Q;
}

static inline Eigen::Matrix<double,4,4> QMatBar(Eigen::Quaterniond q)
{
    if (q.w() < 0) q = Eigen::Quaterniond(-q.w(), -q.x(), -q.y(), -q.z());
    double qw = q.w(), qx = q.x(), qy = q.y(), qz = q.z();
    Eigen::Matrix<double,4,4> Q;
    Q << qw,-qx,-qy,-qz,
         qx, qw, qz,-qy,
         qy,-qz, qw, qx,
         qz, qy,-qx, qw;
    return Q;
}

static inline Mat3D rightJacobianInvRV(const RotVec &rv)
{
    Mat3D S;
    S <<      0, -rv.az,   rv.ay,
          rv.az,      0,  -rv.ax,
         -rv.ay,  rv.ax,       0;
    double n0 = rv.norm();
    if (n0 < 1e-7) return Mat3D::Identity() + S/2;
    return Mat3D::Identity() + S/2 + (1/(n0*n0)-(1+cos(n0))/(2*n0*sin(n0)))*S*S;
}

static inline Eigen::Matrix<double,3,9> jacobianR(const Mat3D &R)
{
    double cost = ((R(0,0) + R(1,1) + R(2,2))-1)/2.;
    double t = acos(cost);
    double sint = sqrt(1-cost*cost);
    Eigen::Matrix<double,3,9> J;
    J <<
     0,    0,    0,    0,    0,    1,    0,   -1,    0,
     0,    0,   -1,    0,    0,    0,    1,    0,    0,
     0,    1,    0,   -1,    0,    0,    0,    0,    0;
    if (fabs(1-cost) < 1e-8) {
        J = J * 0.5;
    } else {
        Vec3D a;
        a << R(2,1)-R(1,2), R(0,2)-R(2,0), R(1,0)-R(0,1);
        a = a * (t*cost-sint)/(4*sint*sint*sint);
        double b = t/(2*sint);

        J = J * b;
        J.block<3,1>(0,0) = a;
        J.block<3,1>(0,4) = a;
        J.block<3,1>(0,8) = a;
    }
    return J;
}

// === 3D pose-graph optimizer ===

struct Pose3D
{
    double x, y, z;
    RotVec rv;
    Pose3D() : x(0), y(0), z(0) {}
    Pose3D(Vec3D in_t, const RotVec &in_rv) : x(in_t[0]), y(in_t[1]), z(in_t[2]), rv(in_rv) {}
    Pose3D(double in_x, double in_y, double in_z, const RotVec &in_rv) : x(in_x), y(in_y), z(in_z), rv(in_rv) {}
    Pose3D(double in_x,  double in_y,  double in_z,
           double ax, double ay, double az) : x(in_x), y(in_y), z(in_z), rv(ax, ay, az) {}
    Pose3D(const Eigen::Matrix<double, 6, 1> &v) : x(v[0]), y(v[1]), z(v[2]), rv(v[3], v[4], v[5]) {}
    Vec6D vec() const {
        Vec6D v;
        v << x, y, z, rv.ax, rv.ay, rv.az;
        return v;
    }
    Vec3D pos() const {
        Vec3D v;
        v << x, y, z;
        return v;
    }
    Pose3D oplus(const Pose3D &rel) const {
        Vec3D t = rv.toRotationMatrix()*rel.pos() + this->pos();
        RotVec rv2(rv.toQuaternion() * rel.rv.toQuaternion());
        return Pose3D(t[0], t[1], t[2], rv2);
    }
    Pose3D ominus(const Pose3D &base) const {
        Vec3D t = base.rv.toRotationMatrix().transpose()*(this->pos() - base.pos());
        RotVec rv2(base.rv.toQuaternion().conjugate() * rv.toQuaternion());
        return Pose3D(t[0], t[1], t[2], rv2);
    }
    static Vec6D concat(const Vec6D &a, const Vec6D &b)
    {
        return Pose3D(a).oplus(Pose3D(b)).vec();
    }
    Eigen::Matrix<double,4,4> toMatrix() {
        Eigen::Matrix<double,4,4> mat;
        mat.setIdentity();
        mat.block<3,3>(0,0) = rv.toRotationMatrix();
        mat(0,3) = x;
        mat(1,3) = y;
        mat(2,3) = z;
        return mat;
    }
};

using Pose3DVec = std::vector<Pose3D>;

} //namespace opt4imu
