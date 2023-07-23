#include <iostream>

#include <eigen3/Eigen/Dense>

int main() {
    double val[10];
    for (int i = 0 ; i < 10 ; ++i) 
        val[i] = i;

    Eigen::Map<Eigen::MatrixXd> testmat(val,2,5);
    std::cout << "testmat size = " <<  testmat.rows() << ", " << testmat.cols() << std::endl;
    std::cout << testmat << std::endl;

    return 0;
}