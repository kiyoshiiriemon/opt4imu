/****************************************************************************
 * Copyright (C) 2023 Kiyoshi Irie
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this file,
 * You can obtain one at https://mozilla.org/MPL/2.0/.
 ****************************************************************************/

#pragma once
#include <vector>
#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Sparse>

template<typename myfloat_t, size_t D, class Edge>
class GraphGaussNewtonOptimizer
{
public:
    using Node = Eigen::Matrix<myfloat_t, D, 1>;
    using NodeList = Eigen::Matrix<myfloat_t,D,Eigen::Dynamic>;
    using EdgeList = std::vector<Edge>;
    using Info = Eigen::Matrix<myfloat_t, -1, -1>;
    //using Jac  = Eigen::Matrix<myfloat_t, D, D>;
    using Jac  = Eigen::Matrix<myfloat_t, -1, D>;
    using Vector = Eigen::Matrix<myfloat_t, Eigen::Dynamic, 1>;
    using ErrFunc = std::function<Vector(const Node&, const Node&, const Edge&, Jac&, Jac&)>;
    using UpdateFunc = std::function<NodeList(const NodeList &, const Eigen::Matrix<myfloat_t, Eigen::Dynamic, 1> &)>;
    using OplusFunc = std::function<Node(const Node&, const Node&)>;

    int max_steps = 1000;
    double stop_thre = 1e-9;
    bool verbose = false;
    ErrFunc error_func;
    UpdateFunc update_func = &GraphGaussNewtonOptimizer::update_linear;
    OplusFunc oplus_func;
    std::vector<Eigen::Triplet<myfloat_t>> tripletList;

    static NodeList update_linear(const NodeList &n, const Eigen::Matrix<myfloat_t, Eigen::Dynamic, 1> &d)
    {
        NodeList dmat = Eigen::Map<const Eigen::Matrix<myfloat_t,D,Eigen::Dynamic>, Eigen::RowMajor>(d.data(), D, d.rows()/D);
        return n + dmat;
    }

    double optimizeOneStep( NodeList &out_n, const NodeList &nodes, const std::vector<Edge> &edges)
    {
        double total_cost = 0;
        constexpr static int dim = D;
        std::array<int, dim> indlist;
        std::iota(indlist.begin(), indlist.end(), 0);
        size_t numnodes = nodes.cols();
        Eigen::Matrix<myfloat_t, Eigen::Dynamic, 1> bf = Eigen::Matrix<myfloat_t, Eigen::Dynamic, 1>::Zero(numnodes * dim);

        tripletList.clear();
        tripletList.reserve(edges.size() * dim * dim * 4);

        for(const Edge &con: edges) {
            int ida = con.id1;
            int idb = con.id2;
            assert(0 <= ida && ida < numnodes);
            assert(0 <= idb && idb < numnodes);
            const Node &pa = nodes.col(ida);
            const Node &pb = nodes.col(idb);

            //Info Ja, Jb;
            Eigen::Matrix<double, -1, D> Ja, Jb;
            Vector r = error_func(pa, pb, con, Ja, Jb);
            Info info = con.info;
            //std::cout << "info : " << info.rows() << " " << info.cols() << " " << r.rows() << std::endl;
            total_cost += r.transpose() * info * r;

            Info trJaInfo = Ja.transpose() * info;
            Info trJaInfoJa = trJaInfo * Ja;
            Info trJbInfo = Jb.transpose() * info;
            Info trJbInfoJb = trJbInfo * Jb;
            Info trJaInfoJb = trJaInfo * Jb;

            for(size_t k: indlist) {
                for(size_t m: indlist) {
                    tripletList.push_back(Eigen::Triplet<myfloat_t>(ida*dim+k, ida*dim+m, trJaInfoJa(k,m)));
                    tripletList.push_back(Eigen::Triplet<myfloat_t>(idb*dim+k, idb*dim+m, trJbInfoJb(k,m)));
                    tripletList.push_back(Eigen::Triplet<myfloat_t>(ida*dim+k, idb*dim+m, trJaInfoJb(k,m)));
                    tripletList.push_back(Eigen::Triplet<myfloat_t>(idb*dim+k, ida*dim+m, trJaInfoJb(m,k)));
                }
            }
            bf.segment(ida*dim, dim) += trJaInfo * r;
            bf.segment(idb*dim, dim) += trJbInfo * r;
        }
        for(size_t k: indlist) {
            tripletList.push_back(Eigen::Triplet<myfloat_t>(k, k, 1e10));
        }
        double lambda = 1e-9;
        for(size_t i=0; i<dim*numnodes; ++i) {
            tripletList.push_back(Eigen::Triplet<myfloat_t>(i, i, lambda));
        }

        Eigen::SparseMatrix<myfloat_t> mat(numnodes*dim, numnodes*dim);
        mat.setFromTriplets(tripletList.begin(), tripletList.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<myfloat_t> > solver;
        solver.compute(mat);
        Eigen::Matrix<myfloat_t, Eigen::Dynamic, 1> x = solver.solve(-bf);

        if (oplus_func) {
            out_n.resize(nodes.rows(),D);
            for(unsigned int i=0;i<nodes.cols();i++) {
                out_n.col(i) = oplus_func(nodes.col(i), x.segment(i*D, D));
            }
        } else {
            out_n = update_func(nodes, x);
        }
        return total_cost;
    }

    NodeList optimize(const NodeList &in_nodes, const EdgeList &edges)
    {
        NodeList graphnodes = in_nodes;
        NodeList ret(in_nodes.rows(), in_nodes.cols());
        myfloat_t prevres = std::numeric_limits<myfloat_t>::max();
        for(int i=1; i<=max_steps; i++) {
            auto t0 = std::chrono::high_resolution_clock::now();
            myfloat_t res = optimizeOneStep(ret, graphnodes, edges);
            auto t1 = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast< std::chrono::microseconds> (t1-t0);
            myfloat_t convergence = (ret-graphnodes).norm()/graphnodes.norm(); 
            if (verbose)
                std::cout << "step " << i << ": " << convergence << " time: " << elapsed.count()*1e-6 << "s" << std::endl;
            graphnodes = ret;

            if (i >= 2 && convergence < stop_thre) {
                if (verbose)
                    std::cout << "converged: " << prevres - res << " < " << stop_thre << std::endl;
                break;
            }
            prevres = res;
        }
        return graphnodes;
    }
};

// same dimensions for node and error
template<typename myfloat_t, size_t D>
struct ConstraintND
{
    int id1, id2;
    Eigen::Matrix<myfloat_t, D, 1> t;
    Eigen::Matrix<myfloat_t, D, D> info;
    void setCovarianceMatrix( const Eigen::Matrix<myfloat_t, D, D> &mat )
    {
        info = mat.inverse();
    }
EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
