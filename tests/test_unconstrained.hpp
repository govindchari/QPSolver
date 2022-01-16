#pragma once
#include "catch.hpp"
#include "QPSolver.hpp"
#include "epigraph.hpp"
#include <iostream>
#include <chrono>

using namespace cvx;
using namespace std::chrono;

TEST_CASE("Unconstrained Optimization")
{
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "Unconstrained Optimization" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    // Setup problem data
    int n = 100; 
    Eigen::MatrixXd A;
    Eigen::MatrixXd Q(n, n);    
    Eigen::MatrixXd G;
    Eigen::VectorXd q(n);       
    Eigen::VectorXd b;
    Eigen::VectorXd h;

    Q.setRandom();
    Q = 0.5 * (Q.transpose() + Q);
    Q += n * Eigen::MatrixXd::Identity(n, n);
    q.setRandom();

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