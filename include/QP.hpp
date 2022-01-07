#pragma once
#include "Eigen/Dense"

#include "utilities.hpp"
#include <iostream>
class QP
{
private:

    //Problem Data
    //Eigen::SparseMatrix<double> Q;
    Eigen::MatrixXd Q;
    Eigen::VectorXd q;
    //Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    //Eigen::SparseMatrix<double> G;
    Eigen::MatrixXd G;
    Eigen::VectorXd h;

    //Optimization variables
    Eigen::VectorXd x;
    Eigen::VectorXd s;
    Eigen::VectorXd z;
    Eigen::VectorXd y;

    //KKT System
    //Eigen::SparseMatrix<double> KKT;
    Eigen::MatrixXd KKT;
    Eigen::VectorXd rhs_a;
    Eigen::VectorXd rhs_c;
    Eigen::VectorXd p_a;
    Eigen::VectorXd p_c;

    IDX idx;

    DELTA delta;

public:
    void solve();
    void solve(bool verb);

    QP(Eigen::MatrixXd Qi, Eigen::VectorXd qi, Eigen::MatrixXd Ai, Eigen::VectorXd bi, Eigen::MatrixXd Gi, Eigen::VectorXd hi) : Q(Qi), q(qi), A(Ai), b(bi), G(Gi), h(hi)
    {

        //Variable lengths
        unsigned int nx = q.size();
        unsigned int ns = h.size();
        unsigned int nz = h.size();
        unsigned int ny = b.size();
        unsigned int N = nx + ns + nz + ny;

        IDX idx(nx,ns,nz,ny,N);

        //KKT System
        KKT.setZero(N,N);
        rhs_a.setZero(N,1);
        rhs_c.setZero(N,1);
        p_a.setZero(N,1);
        p_c.setZero(N,1);

        //Optimization variables
        x.setZero(nx, 1);
        s.setZero(ns, 1);
        z.setZero(nz, 1);
        y.setZero(ny, 1);

        DELTA delta(nx,ns,nz,ny);

        std::cout<<Q<<std::endl;
    }
};
