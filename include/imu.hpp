/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include "msg_header.hpp"
#include "geom3d.hpp"

namespace opt4imu {

struct Imu
{
    opt4imu::Header header;
    Eigen::Quaterniond orientation;
    opt4imu::Vec3D angular_velocity;
    opt4imu::Vec3D linear_acceleration;
    std::string to_str()
    {
        std::ostringstream ss;
        ss << "[Imu] seq=" << header.seq << " " << header.stamp << std::endl;
        ss << " orientation(x, y, z, w): " << orientation.x() << " " <<
                                              orientation.y() << " " <<
                                              orientation.z() << " " <<
                                              orientation.w() << std::endl;
        ss << " acc  (x, y, z): " << linear_acceleration(0) << " "
           << linear_acceleration(1) << " "
           << linear_acceleration(2) << std::endl;
        ss << " gyro (x, y, z): " << angular_velocity(0) << " "
            << angular_velocity(1) << " "
            << angular_velocity(2) << std::endl;
        return ss.str();
    }
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // opt4imu
