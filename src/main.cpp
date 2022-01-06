#include <iostream>
#include "QP.hpp"
#include "utilities.hpp"
#include "Eigen/Dense"

int main(){
    DELTA delta(3,2,3,4);
    std::cout<<delta.x_a<<std::endl;
    return 0;
}