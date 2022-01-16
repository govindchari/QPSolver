#pragma once
#include "catch.hpp"
#include "QPSolver.hpp"
#include "epigraph.hpp"
#include <iostream>
#include <chrono>

using namespace cvx;
using namespace std::chrono;

TEST_CASE("Randomized Markowitz Portfolio QP")
{
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "Randomized Markowitz Portfolio" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    // Setup problem data
    int n = 100; //Number of assets
    Eigen::MatrixXd A(1, n);
    Eigen::MatrixXd Q(n, n); //Return Covariance
    Eigen::MatrixXd G(n, n);
    Eigen::VectorXd q(n); //Expected Return
    Eigen::VectorXd b(1);
    Eigen::VectorXd h(n);

    
    //Generate random return covariance matrix (PSD matrix)
    Q.setRandom();
    Q = Q.transpose() * Q;
    Q += n * Eigen::MatrixXd::Identity(n, n);

    //Generate random expected return
    q.setRandom();

    //No shorting (non-negativity constraint)
    G = -Eigen::MatrixXd::Identity(n, n);
    h.setZero(n);

    A.setConstant(1);
    b << 1;

    //Set up QPSolver solver
    QP qp(Q, q, A, b, G, h);
    std::cout << "=============================QPSolver===============================" << std::endl;
    auto start1 = high_resolution_clock::now();
    qp.solve(true);
    auto stop1 = high_resolution_clock::now();
    auto time1 = duration_cast<microseconds>(stop1 - start1);

    //Set up OSQP solver
    OptimizationProblem qp_epi;
    VectorX x = qp_epi.addVariable("x", n);

    qp_epi.addConstraint(equalTo(par(A) * (x), par(b)));
    for (int i = 0; i < G.rows(); i++)
    {
        auto Gi = G.row(i);
        qp_epi.addConstraint(lessThan(par(Gi).dot(x), par(h(i))));
    }
    qp_epi.addCostTerm(par(0.5) * x.transpose() * par(Q) * x + par(q).dot(x));
    osqp::OSQPSolver solver(qp_epi);
    std::cout << "===============================OSQP=================================" << std::endl;
    auto start2 = high_resolution_clock::now();
    solver.solve(true);
    auto stop2 = high_resolution_clock::now();
    auto time2 = duration_cast<microseconds>(stop2 - start2);

    //Speed Report
    std::cout << "===========================Speed Results=============================" << std::endl;
    std::cout << "QPSolver Run Time (us):" << std::endl;
    std::cout << time1.count() << std::endl;
    std::cout << "OSQP Run Time (us):" << std::endl;
    std::cout << time2.count() << std::endl;

    //Compare OSQP and QPSolver solutions
    Eigen::VectorXd QPSolver_sol = qp.solution.x;
    Eigen::VectorXd OSQP_sol = eval(x);

    double tol = 1e-3;
    REQUIRE((QPSolver_sol - OSQP_sol).lpNorm<Eigen::Infinity>() <= tol);
}