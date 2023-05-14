/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include "graph_gn_se3.hpp"
#include <iostream>
#include <iterator>
#include <fstream>

bool load_g2o_file( const char *g2ofile, Pose3DVec &nodes, Con6DVec &constraints )
{
    std::ifstream is(g2ofile);
    if (!is) return false;

    nodes.clear();
    constraints.clear();
    int id;
    double x, y, z, qx, qy, qz, qw;
    while(is){
        char buf[1024];
        is.getline(buf,1024);
        std::istringstream sstrm(buf);
        std::string tag;
        std::vector<std::string> tok = {std::istream_iterator<std::string>{sstrm}, std::istream_iterator<std::string>{}};
        std::cout << tok.size() << std::endl;
        if (tok.size() < 9) break;
        int i = 0;
        if (tok[i]=="VERTEX_SE3:QUAT") {
            ++i;
            id = std::stoi(tok[i++]);
            x  = std::stod(tok[i++]);
            y  = std::stod(tok[i++]);
            z  = std::stod(tok[i++]);
            qx = std::stod(tok[i++]);
            qy = std::stod(tok[i++]);
            qz = std::stod(tok[i++]);
            qw = std::stod(tok[i++]);
            nodes.push_back(Pose3D(x, y, z, Eigen::Quaterniond(qw, qx, qy, qz)));
        } else if (tok[i]=="EDGE_SE3:QUAT"){
            ++i;
            Con6D con;
            con.info.setIdentity();
            con.id1 = std::stoi(tok[i++]);
            con.id2 = std::stoi(tok[i++]);
            x  = std::stod(tok[i++]);
            y  = std::stod(tok[i++]);
            z  = std::stod(tok[i++]);
            qx = std::stod(tok[i++]);
            qy = std::stod(tok[i++]);
            qz = std::stod(tok[i++]);
            qw = std::stod(tok[i++]);
            con.t = Pose3D(x, y, z, Eigen::Quaterniond(qw, qx, qy, qz)).vec();
            if (tok.size() == 9 + 21) {
                for (int ii = 0; ii < 6; ii++) {
                    for (int jj = ii; jj < 6; jj++) {
                        double val = std::stod(tok[i++]);
                        con.info(ii, jj) = val;
                        con.info(jj, ii) = val;
                    }
                }
                std::cout << "read information matrix" << std::endl;
            }
            if (con.id1 > con.id2) {
                con.t = Pose3D().ominus(con.t).vec();
                int id0 = con.id1;
                con.id1 = con.id2;
                con.id2 = id0;
            }
            constraints.push_back(con);
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    std::vector<Pose3D> p3dnodes;
    Con6DVec constraints;
    if (argc < 2) return 0;
    load_g2o_file(argv[1], p3dnodes, constraints);

    NodeList nodes(6, p3dnodes.size());
    for(int i=0; i<nodes.cols(); ++i) {
        nodes.col(i) = p3dnodes[i].vec();
    }
    //std::cout << nodes << std::endl;
    GraphGNSE3 opt;
    opt.verbose = true;
    NodeList result = opt.optimize(nodes, constraints);
    //std::cout << result.rows() << " " << result.cols() << std::endl;
    std::cout << result.transpose() << std::endl;

    return 0;
}

