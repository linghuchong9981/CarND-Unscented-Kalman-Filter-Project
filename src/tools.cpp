#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector <VectorXd> &estimations,
                              const vector <VectorXd> &ground_truth) {
    /**
    TODO:
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        std::cout << "CalculateRMSE Estimation size is zero!!!" << std::endl;
        return rmse;
    }

    for (int i = 0; i < estimations.size(); i++) {
//        std::cout << "CalculateRMSE index: " << i << ",truth:" << ground_truth[i].transpose() << ",estimation:" << estimations[i].transpose() << std::endl;
        VectorXd diff = ground_truth[i] - estimations[i];
        diff = diff.array() * diff.array();
        rmse += diff;
    }



    rmse /= estimations.size();
//    std::cout << "CalculateRMSE Estimation is " << rmse.transpose() << ",size:" << estimations.size() << std::endl;
    rmse = rmse.array().sqrt();
    return rmse;
}

double Tools::NormalizationAngle(double angle) {
    angle = fmod(angle, 2. * M_PI);
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle < -M_PI) angle += 2. * M_PI;
    return angle;
}
