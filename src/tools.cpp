#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    if (estimations.size() != ground_truth.size() || estimations.size()==0) {
        std::cout << "Estimation size is zero!!!" << std::endl;
        return rmse;
    }
    
    for (int i=0; i<estimations.size(); i++) {
        VectorXd diff = ground_truth[i] - estimations[i];
        diff = diff.array() * diff.array();
        rmse += diff;
    }
    
    rmse /= estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

double Tools::NormalizationAngle(double angle){
    angle =  fmod(angle,2.*M_PI);
    while (angle > M_PI) angle -= 2.*M_PI;
    while (angle<-M_PI) angle += 2.*M_PI;
    return angle;
}
