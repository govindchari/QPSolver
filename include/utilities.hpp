#pragma once
#include "Eigen/Dense"

struct DELTA
{
    Eigen::VectorXd x_a;
    Eigen::VectorXd s_a;
    Eigen::VectorXd z_a;
    Eigen::VectorXd y_a;

    Eigen::VectorXd x_cc;
    Eigen::VectorXd s_cc;
    Eigen::VectorXd z_cc;
    Eigen::VectorXd y_cc;

    Eigen::VectorXd x;
    Eigen::VectorXd s;
    Eigen::VectorXd z;
    Eigen::VectorXd y;

    DELTA(){}

    DELTA(unsigned int nx, unsigned int ns, unsigned int nz, unsigned int ny)
    {
        x_a.setZero(nx, 1);
        s_a.setZero(ns, 1);
        z_a.setZero(nz, 1);
        y_a.setZero(ny, 1);

        x_cc.setZero(nx, 1);
        s_cc.setZero(ns, 1);
        z_cc.setZero(nz, 1);
        y_cc.setZero(ny, 1);

        x.setZero(nx, 1);
        s.setZero(ns, 1);
        z.setZero(nz, 1);
        y.setZero(ny, 1);
    }
};
struct IDX
{
    unsigned int nx;
    unsigned int ns;
    unsigned int nz;
    unsigned int ny;
    unsigned int N;

    IDX(){}

    IDX(int numx, int nums, int numz, int numy, int numN) : nx(numx), ns(nums), nz(numz), ny(numy), N(numN) {}
};

struct SOLUTION {
    Eigen::VectorXd x;
    Eigen::VectorXd s;
    Eigen::VectorXd z;
    Eigen::VectorXd y;
    double obj;
};
