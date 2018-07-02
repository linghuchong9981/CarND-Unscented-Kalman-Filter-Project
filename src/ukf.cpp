#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 6;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_aug_;
    Xsig_pred_ = MatrixXd(n_x_,2 * n_aug_ + 1);
    
    weights_ = VectorXd(2*n_aug_+1);
    double weight_0 = lambda_/(lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
    

  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    cout << "ProcessMeasurement begin,sensor type " << meas_package.sensor_type_ <<  endl;
    if (!use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
        return;
    }

    if (!use_radar_ && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
        return;
    }
    
    //init
    if (!is_initialized_) {
        InitState(meas_package);
        return;
    }
    
    
    double delta_t = meas_package.timestamp_ - time_us_;
    Prediction(delta_t);
    time_us_ = meas_package.timestamp_;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
        return;
    }
    
    UpdateLidar(meas_package);
}

//初始化CTRV状态及协方差矩阵
void UKF::InitState(MeasurementPackage meas_package) {
    double px;
    double py;
    x_.fill(0.0);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        double ro = meas_package.raw_measurements_(0);
        double theta = meas_package.raw_measurements_(1);
        double ro_dot = meas_package.raw_measurements_(2);
        px = ro * cos(theta);
        py = ro * sin(theta);
    } else {
        px =  meas_package.raw_measurements_(0);
        py =  meas_package.raw_measurements_(1);
    }
    
    if (px < 0.001) {
        px = 0.001;
    }
    
    if (py < 0.001) {
        py = 0.001;
    }
    
    x_(0) = px;
    x_(1) = py;
    
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    
    cout << "Init statu " << x_.transpose() << endl;
    cout << "Init statu P " << P_ << endl;
    cout << "Init statu weights " << weights_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    //for test data
//    delta_t = 0.1;
//    x_ <<   5.7441,
//    1.3800,
//    2.2049,
//    0.5015,
//    0.3528;
//
//    P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
//    -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
//    0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
//    -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
//    -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
//
//    std_a_ = 0.2;
//    std_yawdd_ = 0.2;
    
    
    cout << "Prediction start " << x_.transpose() << endl;
    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    x_aug.fill(0.0);
    x_aug.head(5) = x_;
    cout << "Prediction x_aug " << x_aug.transpose() << endl;
    
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_ * std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    cout << "Prediction P_aug " << P_aug << endl;
    
    MatrixXd L = P_aug.llt().matrixL();
    
    cout << "Prediction L " << L << endl;
    
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    
    cout << "Prediction Xsig_aug " << Xsig_aug << endl;
    
//    cout << "Prediction weights " << Xsig_aug << endl;

//    cout << "Prediction Xsig aug finish " << x_.transpose() << endl;
    
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
    
    cout << "Prediction Xsig_pred_ " << Xsig_pred_ << endl;
    
    //for test data
//    Xsig_pred_ <<
//    5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
//    1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
//    2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
//    0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
//    0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
    
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_+ weights_(i) * Xsig_pred_.col(i);
    }
    
    cout << "Predict x " << x_.transpose() << endl;
    
    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        x_diff(3) = tools_.NormalizationAngle(x_diff(3));
        cout << " predict x diff " << x_diff.transpose() << endl;
        MatrixXd plus = weights_(i) * x_diff * x_diff.transpose();
        cout << " predict plus " << plus << endl;
        P_ = P_ + plus ;
    }
    
    
    cout << "Predict P_" << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
     TODO:
     You'll also need to calculate the lidar NIS.
     */
    
    cout << "Update By Lidar start " << meas_package.raw_measurements_.transpose() << endl;
//    cout << "Update By Lidar x " << x_.transpose() << endl;
    VectorXd z = meas_package.raw_measurements_;
    Eigen::MatrixXd H = MatrixXd(2,5);
    H <<   1,0,0,0,0,
    0,1,0,0,0;
    VectorXd z_pred = H * x_;
    cout << "Update By Lidar z_pred " << z_pred.transpose() << endl;
    VectorXd diff = z - z_pred;
    MatrixXd Ht = H.transpose();
    
    cout << "Update By Lidar after diff " << diff << "Ht" << Ht << endl;
    
    MatrixXd R = MatrixXd(2,2);
    R <<    std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;
    
    cout << "Update By Lidar P_ " << P_ << "R" << R << endl;
    
    MatrixXd S = H * P_ * Ht + R;
    cout << "Update By Lidar after Keman S " << S << endl;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    cout << "Update By Lidar after Keman gain " << K << endl;
    
    //new estimate
    x_ = x_ + (K * diff);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;
    cout << "Update By Lidar " << x_.transpose() << endl;
    cout << "Update By Lidar P_" << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
     TODO:
     You'll also need to calculate the radar NIS.
     */

    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    
    std::cout << "Zsig: "  << Zsig << std::endl;
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        z_diff(1) = tools_.NormalizationAngle(z_diff(1));
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(3,3);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;
    
    std::cout << "S: "  << S << std::endl;
    
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        z_diff(1) = tools_.NormalizationAngle(z_diff(1));
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        x_diff(3) = tools_.NormalizationAngle(x_diff(3));
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    std::cout << "Tc: "  << Tc << std::endl;
    
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    std::cout << "K: "  << K << std::endl;
    
    std::cout << "K: " << std::endl << K << std::endl;
    
    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
    
    //angle normalization
    z_diff(1) = tools_.NormalizationAngle(z_diff(1));
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    cout << "Update By Radar " << x_.transpose() << endl;
}

