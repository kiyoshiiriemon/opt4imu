/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include "geom3d.hpp"
#include "graph_gn.hpp"
#include "graph_gn_so3.hpp"
#include <fstream>
#include <iostream>

using namespace opt4imu;

using RVVec = std::vector<RotVec>;
using ConSO3Vec = std::vector<ConSO3>;

static bool load_g2o_file(const char *fname, RVVec &nodes, ConSO3Vec &constraints)
{
    std::ifstream is(fname);
    if (!is) return false;

    nodes.clear();
    constraints.clear();
    int id;
    double rx, ry, rz;
    while(is){
        char buf[1024];
        is.getline(buf,1024);
        std::istringstream sstrm(buf);
        std::string tag;
        sstrm >> tag;
        if (tag=="VERTEX_SO3:RV"){
            sstrm >> id >> rx >> ry >> rz; 
            nodes.push_back(RotVec(rx, ry, rz));
        } else if (tag=="EDGE_SO3:RV"){
            ConSO3 con;
            sstrm >> con.id1 >> con.id2 >> rx >> ry >> rz;
            con.t = Vec3D(rx, ry, rz);
            for(int i=0; i<3; i++) {
                for(int j=i; j<3; j++) {
                    double val;
                    sstrm >> val;
                    con.info(i,j) = val;
                    con.info(j,i) = val;
                }
            }
            constraints.push_back(con);
        }
    }
    return true;
}

int main()
{
    using NodeList = GraphGNSO3::NodeList;
    RVVec rvnodes;
    ConSO3Vec constraints;
    load_g2o_file("hoge.g2o", rvnodes, constraints);

    NodeList nodes(rvnodes.size(),3);
    for(int i=0; i<rvnodes.size(); ++i) {
        nodes.row(i) = rvnodes[i].vec();
    }

    GraphGNSO3 opt;
    opt.verbose = true;
    NodeList result = opt.optimize(nodes, constraints);
    for (int i=0; i<result.rows(); ++i) {
        std::cout << result.block<1,3>(i,0) << std::endl;
    }

    return 0;
}

