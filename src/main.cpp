#include <iostream>
#include "QP.hpp"
#include "utilities.hpp"
#include "Eigen/Dense"

int main(){
    Eigen::MatrixXd A(2,2);
    Eigen::MatrixXd Q(2,2);
    Eigen::MatrixXd G(2,2);
    Eigen::VectorXd q(2);
    Eigen::VectorXd b(2);
    Eigen::VectorXd h(2);

    A << 1,2,3,4;
    Q << 1,2,3,4;
    G << 1,2,3,4;
    q << 1,2;
    b << 1,2;
    h << 1,2;

    QP qp(Q,q,A,b,G,h);
    qp.solve();
    return 0;
}