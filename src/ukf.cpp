#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

 
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
  
  //set state dimension (from lecture 7-18)
  int n_x_ = 5; 
  //set augmented dimension (from lecture 7-18)
  int n_aug_ = 7;
  //define spreading parameter (from lecture 7-18)
  double lambda_ = 3 - n_aug_;
  int n_seg_ = 2 * n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_, n_seg_);
  
  //set vector for weights (from lecture 7-30)
  VectorXd weights = VectorXd(n_seg_);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }
  //set NIS Lidar
  double NIS_LIDAR_ = 0;
  //set NIS Radar
  double NIS_RADAR_ = 0;  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  //initialization like for EKF:
  if (!is_initialized_) {
    //initialize timestamp
    time_us_ = meas_package.timestamp_;
    //initialize state(x_) and covariance matrix(P_) (from lecture 7-17 + 7-29)
    x_ << 0.0,   //ekf x   ukf px 
          0.0,   //ekf y   ukf py
          0.0,   //ekf dx  ukf v
          0.0,   //ekf dy  ukf yaw
          0.0;   //ekf na  ukf yawd
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0, //expected error for derivatives a little larger
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 10.0, 0.0, 0.0, 
          0.0, 0.0, 0.0, 10.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.1;
    
    //set values for RADAR and LIDAR:
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rhodot = meas_package.raw_measurements_(2);
      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = rhodot;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    //set initialization:
    is_initialized_ = true;
    }
  else {
    // set timestep
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
    //if RADAR
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
    }
    //if LIDAR
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
    }
    //if NONE
    //else {
    //  return;
    //}
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /* AUGMENTATION 
  ***************/
  //create augmented mean vector (from lecture 7-17)
  VectorXd x_aug = VectorXd(7);
  //create augmented state covariance (from lecture 7-17)
  MatrixXd P_aug = MatrixXd(7, 7);
  //create sigma point matrix (from lecture 7-18)
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //create augmented mean state (from lecutre 7-17)
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  //create augmented covariance matrix (from lecture 7-17)
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_ = 0.2; //from Lecture

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_ = 0.2; //from lecture
  
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  //create square root matrix (from lecture 7-17)
  MatrixXd L = P_aug.llt().matrixL();
  //create augmented sigma points (from lecture 7-17)
  
  //try something different
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  /*
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd B = sqrt(lambda_ + n_aug_) * A;
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.leftCols(n_aug_ + 1).rightCols(n_aug_) = x_aug.replicate(1, B.cols()) + B;
  Xsig_aug.rightCols(n_aug_) = x_aug.replicate(1, B.cols()) - B;
  */
  
 
  /* SIGMA POINT PREIDCTION
  *************************/
  
  // all from lecture 7-21
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
    //avoid multiple calculations (from EKF)
    double d_t2 = delta_t * delta_t;
    double yawddt = yawd * delta_t;
    //avoid division by zero and make calculation faster
    //write predicted sigma point directly into right column
    if (fabs(yawd) > 0.001) {
        Xsig_pred_(0,i) = p_x + v/yawd * ( sin (yaw + yawddt) - sin(yaw)) + (0.5 * nu_a * d_t2 * cos(yaw));
        Xsig_pred_(1,i) = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawddt) ) + (0.5 * nu_a * d_t2 * sin(yaw));
    }
    else {
        Xsig_pred_(0,i) = p_x + v*delta_t*cos(yaw) + (0.5 * nu_a * d_t2 * cos(yawd));
        Xsig_pred_(1,i) = p_y + v*delta_t*sin(yaw) + (0.5 * nu_a * d_t2 * sin(yawd));
    }
    Xsig_pred_(2,i) = v + nu_a * delta_t;
    Xsig_pred_(3,i) = yaw + yawddt + (0.5 * nu_yawdd * d_t2);
    Xsig_pred_(4,i) = yawd + nu_yawdd * delta_t;
  }

    /* PREDICTED MEAN AND COVARIANCE
  ********************************/
  // from lecture 7-23

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization: call from separate function located in tools.cpp (from EKF)
    x_diff(3,i) = NormalizeAngles(x_diff(3,i));
    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    /* PREDICT LIDAR
  ******************/
    // copied and adapted from RADAR PREDICTION (see below)
  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;  
  S = S + weights(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;
    
    /* UPDATE LIDAR
  *****************/
  //create Vector for incoming lidar measurement
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  NIS_LIDAR_ = z_diff.transpose() * S.inverse() * z_diff;
  
  // prevents P from diverging (check source)
  if (NIS_LIDAR_ > 100.0) {
  is_initialized_ = false;
  }
  
  else {
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    /* PREDICT RADER
  ******************/
  // from lecture 7-26
  int n_z = 3;
  //create matrix for sigma points in measurement space
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
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1,i) = NormalizeAngles(z_diff(1,i));
    S = S + weights(i) * z_diff * z_diff.transpose();
  }
   
  //add measurement noise covariance matrix
  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;

    /* UPDATE RADAR
  *****************/
  //create Vector for incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1,i) = NormalizeAngles(z_diff(1,i));
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3,i) = NormalizeAngles(x_diff(3,i));
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = NormalizeAngles(z_diff(1));
  NIS_RADAR_ = z_diff.transpose() * S.inverse() * z_diff;
  
  // prevents P from diverging (check source)
  if (NIS_RADAR_ > 100.0) {
  is_initialized_ = false;
  }
  else {
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
  }
  
}


double UKF::NormalizeAngles(double angle) {
  auto times = round( fabs( angle / (2.0 * M_PI) ) );
  
  if (angle > M_PI) {
  angle -= times * 2.0 * M_PI;
  }
  
  if (angle < -M_PI) {
  angle += times * 2.0 * M_PI;
  }
  
  return angle;
}
