#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "utilities.hpp"
class QP
{
private:
    Eigen::SparseMatrix<double> Q;
    Eigen::VectorXd q;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Eigen::SparseMatrix<double> G;
    Eigen::VectorXd h;

    //DELTA delta;
    //IDX idx;

public:
    QP(Eigen::SparseMatrix<double> Q, Eigen::VectorXd q, Eigen::SparseMatrix<double> A, Eigen::VectorXd b, Eigen::SparseMatrix<double> G, Eigen::VectorXd h)
    {
        
    }
};
