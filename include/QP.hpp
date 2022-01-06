#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "utilities.hpp"
class QP
{
private:
    Eigen::SparseMatrix<double> Q;
    Eigen::VectorXd q;
    Eigen::SparseMatrix A;
    Eigen::VectorXd b;
    Eigen::SparseMatrix G;
    Eigen::VectorXd h;

    DELTA delta;
    IDX idx;

public:
    QP::QP(Eigen::SparseMatrix<double> Q, Eigen::VectorXd q, Eigen::VectorXd b, Eigen::SparseMatrix G, Eigen::VectorXd h)
    {
        
    }
};
