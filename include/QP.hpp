#pragma once
#include "Eigen/Dense"

#include "utilities.hpp"
#include <iostream>
class QP
{
private:
    //Variable lengths
    unsigned int nx;
    unsigned int ns;
    unsigned int nz;
    unsigned int ny;

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

    void rhs_kkt_a();
    void index_sol_a();
    void rhs_kkt_c(double sig, double mu);
    void index_sol_c();
    double linesearch(Eigen::VectorXd x,Eigen::VectorXd dx);
    void centering_params(double &sig, double &mu);
    void combine_deltas();
    void update_vars(double a);
    void initialize_kkt();
    void update_kkt();
    void logging();

public:
    unsigned int N;
    QP(Eigen::MatrixXd Qi, Eigen::VectorXd qi, Eigen::MatrixXd Ai, Eigen::VectorXd bi, Eigen::MatrixXd Gi, Eigen::VectorXd hi);
    DELTA delta;
    void solve();
    void solve(bool verbose);
};
