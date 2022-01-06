#include <iostream>
#include <QP.hpp>
#include "Eigen/Dense"

int main(){
    Eigen::Vector3f vec;
    vec = {1.0,1.0,1.0};
    std::cout << vec << std::endl;
    return 0;
}