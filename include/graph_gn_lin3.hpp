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

namespace opt4imu {

using ConLin3D = ConstraintND<double, 3>;

class GraphGNLin3 : public GraphGaussNewtonOptimizer<double, 3, ConLin3D>
{
public:
    using Con3DVec = std::vector<ConLin3D>;
    using RVVec = std::vector<RotVec>;
    using Node = Eigen::Matrix<double, 3, 1>;
    //using Jac  = Eigen::Matrix<double, 3, 3>;
    using Jac  = Eigen::Matrix<double, -1, 3>;
    using NodeList = Eigen::Matrix<double, 3, Eigen::Dynamic>;

    static Eigen::VectorXd calcError(const Node &pa, const Node &pb, const ConLin3D &con, Jac &Ja, Jac &Jb)
    {
        Jb = Ja = Eigen::Matrix<double, -1, 3>(3, 3);
        Jb.setIdentity();
        Ja = -Jb;
        return pb-pa-con.t;
    }

public:
    GraphGNLin3() {
        error_func = &GraphGNLin3::calcError;
    }
};

} //opt4imu

