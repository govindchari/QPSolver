#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "QP.hpp"
#include "Eigen/Dense"
#include "utilities.hpp"

TEST_CASE("Test Case") {
    int n = 10;
    Eigen::MatrixXd A;
    Eigen::MatrixXd Q(n,n);
    Eigen::MatrixXd G(2*n,n);
    Eigen::VectorXd q(n);
    Eigen::VectorXd b;
    Eigen::VectorXd h(2*n);

    Eigen::VectorXd temp(n);
    temp.setConstant(1);

    Q = Eigen::MatrixXd::Identity(n, n);

    q.setConstant(1);
    h.setZero(2*n,1);
    h(Eigen::seq(0, n-1)) = temp;

    G(Eigen::seq(0, n - 1),Eigen::seq(0, n - 1)) = Eigen::MatrixXd::Identity(n, n);
    G(Eigen::seq(n, 2*n - 1),Eigen::seq(0, n- 1)) = -Eigen::MatrixXd::Identity(n, n);

    QP qp(Q,q,A,b,G,h);
    qp.solve();
    REQUIRE(2==1);
}