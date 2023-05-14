/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#include <fstream>
#include <iostream>
#include "geom3d.hpp"
#include "graph_gn_lin3.hpp"

using namespace opt4imu;

using Con3D = ConLin3D;

using Vec3DList = std::vector<Vec3D>;
using Con3DList = std::vector<Con3D>;

static bool load_g2o_file(const char *fname, Vec3DList &nodes, Con3DList &constraints)
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
        if (tag=="VERTEX_R3"){
            sstrm >> id >> rx >> ry >> rz; 
            Vec3D v;
            v << rx, ry, rz;
            nodes.push_back(v);
        } else if (tag=="EDGE_R3"){
            Con3D con;
            sstrm >> con.id1 >> con.id2 >> rx >> ry >> rz;
            con.t << rx, ry, rz;
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

Vec3D calcError(const Vec3D &pa, const Vec3D &pb, const Con3D &con, Mat3D &Ja, Mat3D &Jb, Mat3D &Haa, Mat3D &Hab, Mat3D &Hbb)
{
    //Ja = con.info(pa + pb - con.t);
    //Jb = con.info(pb - pb - con.t);
    Jb.setIdentity();
    Ja = -Jb;
    Haa = con.info;
    Hab =-con.info;
    Hbb = con.info;
    return pb-pa-con.t;
}

double optimizeOneStep( Vec3DList &out_nodes, const Vec3DList &graphnodes, const Con3DList &constraints )
{
    double total_cost = 0;
    static const int dim = 3;
    size_t indlist[] = {0, 1, 2};
    size_t numnodes = graphnodes.size();
    Eigen::Matrix<double, Eigen::Dynamic, 1> bf = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(numnodes * dim);

    typedef Eigen::Triplet<double> opt4imu_triplet;
    std::vector<opt4imu_triplet> tripletList;
    tripletList.clear();
    tripletList.reserve(constraints.size() * dim * dim * 4);

    for(const Con3D &con: constraints) {
        int ida = con.id1;
        int idb = con.id2;
        assert(0 <= ida && ida < graphnodes.size());
        assert(0 <= idb && idb < graphnodes.size());
        const Vec3D &pa = graphnodes[ida];
        const Vec3D &pb = graphnodes[idb];

        Mat3D Ja, Jb;
        Mat3D Haa, Hab, Hbb;
        Vec3D r = calcError(pa, pb, con, Ja, Jb, Haa, Hab, Hbb);
        Mat3D info = con.info;
        total_cost += r.transpose() * info * r;

        Mat3D trJaInfo = Ja.transpose() * info;
        Mat3D trJaInfoJa = trJaInfo * Ja;
        Mat3D trJbInfo = Jb.transpose() * info;
        Mat3D trJbInfoJb = trJbInfo * Jb;
        Mat3D trJaInfoJb = trJaInfo * Jb;

        for(size_t k: indlist) {
            for(size_t m: indlist) {
                tripletList.push_back(opt4imu_triplet(ida*dim+k, ida*dim+m, Haa(k,m)));
                tripletList.push_back(opt4imu_triplet(idb*dim+k, idb*dim+m, Hbb(k,m)));
                tripletList.push_back(opt4imu_triplet(ida*dim+k, idb*dim+m, Hab(k,m)));
                tripletList.push_back(opt4imu_triplet(idb*dim+k, ida*dim+m, Hab(m,k)));
            }
        }
        bf.segment(ida*dim, dim) += trJaInfo * r;
        bf.segment(idb*dim, dim) += trJbInfo * r;
    }
    for(size_t k: indlist) {
        tripletList.push_back(opt4imu_triplet(k, k, 1e10));
    }
    double lambda = 1e-9;
    for(size_t i=0; i<dim*numnodes; ++i) {
        tripletList.push_back(opt4imu_triplet(i, i, lambda));
    }

    Eigen::SparseMatrix<double> mat(numnodes*dim, numnodes*dim);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(mat);
    Eigen::Matrix<double, Eigen::Dynamic, 1> x = solver.solve(-bf);

    out_nodes.clear();
    for(unsigned int i=0;i<graphnodes.size();i++) {
        int u_i = i * dim;
        const Vec3D &p = graphnodes[i];
        Vec3D d;
        d << x[u_i], x[u_i+1], x[u_i+2];
        out_nodes.push_back(p + d);
    }
    return total_cost;
}

Vec3DList optimize(const Vec3DList &in_nodes, const Con3DList &constraints)
{
    int max_steps = 10;
    double stop_thre = 1e-9;
    Vec3DList graphnodes = in_nodes;
    Vec3DList ret;
    double prevres = std::numeric_limits<double>::max();
    for(int i=1; i<=max_steps; i++) {
        auto t0 = std::chrono::high_resolution_clock::now();
        double res = optimizeOneStep(ret, graphnodes, constraints);
        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast< std::chrono::microseconds> (t1-t0);
        std::cout << "step " << i << ": " << res << " time: " << elapsed.count()*1e-6 << "s" << std::endl;
        graphnodes = ret;

        if (i >= 2 && prevres - res < stop_thre) {
            std::cout << "converged: " << prevres - res << " < " << stop_thre << std::endl;
            break;
        }
        prevres = res;
    }
    return graphnodes;
}

void test()
{
    const Eigen::Matrix<double, Eigen::Dynamic, 1> d(9);
    Eigen::Matrix<double, Eigen::Dynamic, 1> mat = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::RowMajor> (d.data(), 3);
}

int main()
{
    test();
    Vec3DList nodes_v;
    Con3DList constraints;
    load_g2o_file("pos.g2o", nodes_v, constraints);
    Vec3DList result2 = optimize(nodes_v, constraints);

    GraphGNLin3::NodeList nodes(nodes_v.size(),3);
    std::vector<Con3D> edges;
    for(int i=0; i<nodes_v.size(); ++i) {
        nodes.block<1,3>(i,0) = nodes_v[i];
    }

    GraphGNLin3 opt;
    opt.verbose = true;
    GraphGNLin3::NodeList result = opt.optimize(nodes, constraints);
    for (int i=0; i<result.rows(); ++i) {
        std::cout << result.block<1,3>(i,0) << std::endl;
    }

    return 0;
}

