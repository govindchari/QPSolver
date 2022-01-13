#pragma once
#include "Eigen/Dense"

struct DELTA
{
    Eigen::VectorXd x_a;
    Eigen::VectorXd s_a;
    Eigen::VectorXd z_a;
    Eigen::VectorXd y_a;

    Eigen::VectorXd x_c;
    Eigen::VectorXd s_c;
    Eigen::VectorXd z_c;
    Eigen::VectorXd y_c;

    Eigen::VectorXd x;
    Eigen::VectorXd s;
    Eigen::VectorXd z;
    Eigen::VectorXd y;

    DELTA() {}

    DELTA(unsigned int nx, unsigned int ns, unsigned int nz, unsigned int ny)
    {
        x_a.setZero(nx, 1);
        s_a.setZero(ns, 1);
        z_a.setZero(nz, 1);
        y_a.setZero(ny, 1);

        x_c.setZero(nx, 1);
        s_c.setZero(ns, 1);
        z_c.setZero(nz, 1);
        y_c.setZero(ny, 1);

        x.setZero(nx, 1);
        s.setZero(ns, 1);
        z.setZero(nz, 1);
        y.setZero(ny, 1);
    }
};

struct SOLUTION
{
    Eigen::VectorXd x;
    Eigen::VectorXd s;
    Eigen::VectorXd z;
    Eigen::VectorXd y;
    double obj;
    SOLUTION(){}
};

struct TOLERANCE
{
    double constraint;
    double gap;
    double cost;
    TOLERANCE(){
        constraint = 1e-6;
        gap = 1e-6;
        cost = 1e-4;
    }
};