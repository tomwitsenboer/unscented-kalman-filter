#include "ukf.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = n_x_ + 2; // new dimensions are std_a_ and std_yawdd

  // sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;//true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  // Flag indicating whether the first measurement processed
  is_initialized_ = false;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // weights of sigma-points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // number of a measurement being processed
  counter_ = 0;
}

UKF::~UKF() = default;

void UKF::NormalizeAngle(double *angle) {
  auto times = (unsigned long long) fabs(*angle / (2.0 * M_PI));

  while (*angle > M_PI) {
    *angle -= times * 2.0 * M_PI;
    times = 1;
  }

  while (*angle < -M_PI) {
    *angle += times * 2.0 * M_PI;
    times = 1;
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    return;
  }

  if (!is_initialized_) {
    // the first measurement is handled here

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // convert polar from polar to cartesian coordinate system
      double meas_rho     = meas_package.raw_measurements_[0];
      double meas_phi     = meas_package.raw_measurements_[1];
      double meas_px      = meas_rho     * cos(meas_phi);
      double meas_py      = meas_rho     * sin(meas_phi);
      double meas_v       = 0; // not enough measurements to determine
      double meas_yaw     = M_PI_2 - meas_phi; // phi is measured relative to Y axis while yaw is relative to X axis
      double meas_yawd    = 0; // not enough measurements to determine

      // initial state in case the first measurement comes from radar sensor
      x_ << meas_px,
            meas_py,
            meas_v,
            meas_yaw,
            meas_yawd;

      P_.fill(0.0);
      VectorXd diag(n_x_);
      diag << 1, 1, 10, 1, 10; // we are uncertain about the velocity magnitude and orientation vector change rate
      P_.diagonal() = diag;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // initial state in case the first measurement comes from lidar sensor
      //TODO
    } else {
      std::cerr << "unknown sensor type of measurement; skipping initialization" << std::endl;
      return;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
  } else {
    // difference in seconds between the current measurement and the previous one
    double delta_t_s = (meas_package.timestamp_ - time_us_) / 1000000.0; // us -> s

    // the first step in the UKF processing chain (prediction)
    Prediction(delta_t_s);

    // the second step in the UKF processing chain (update)
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    } if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER){
      // TODO
    }
  }

  // initialize the time at which the measurement was taken
  time_us_ = meas_package.timestamp_;
  ++counter_; // increase counter of processed measurements

  // print statistics and intermediate results
  std::cout << "====================\n"
            << "Measurement number = " << counter_ << "\n"
            << "Measurement type   = " << meas_package.sensor_type_ << "\n"
            << "Time (us)          = " << time_us_ << "\n"
            << "x                  =\n" << x_ << "\n"
            << "P                  =\n" << P_ << "\n"
            << "====================\n" << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // create sigma point matrix which rows are the n_x_ vector
  // with two extra dimensions std_a_ and str_yawdd_ (to be initialized)
  MatrixXd Xsig_aug(n_aug_, 2 * n_aug_ + 1);

  // initialize the matrix (columns are the selected sigma-points in the augmented space)
  GenerateAugmentedSigmaPoints(Xsig_aug);

  // for each sigma-point apply non-linear transformation and store the results in the Xsig_pred_ matrix
  PredictSigmaPoints(Xsig_aug, delta_t);

  // based on the predicted sigma-points, calculate the mean state x_ and covariance matrix P_
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // actual measurements
  VectorXd z = meas_package.raw_measurements_.head(n_z);

  // predicted measurement mean (to be initialized)
  VectorXd z_pred(n_z);

  // measurement covariance matrix S
  MatrixXd S(n_z,n_z);

  // matrix for sigma points in measurement space (to be initialized)
  MatrixXd Zsig(n_z, 2 * n_aug_ + 1);

  // predict the radar measurement and initialize z_pred, S, and Zsig
  PredictRadarMeasurement(z_pred, S, Zsig, n_z);

  UpdateState(z, z_pred, S, Zsig, n_z);
}

void UKF::GenerateAugmentedSigmaPoints(Ref<MatrixXd> Xsig_aug) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(n_aug_-2, n_aug_-2) = std_a_ * std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
}

void UKF::PredictSigmaPoints(const Ref<const MatrixXd> Xsig_aug, double delta_t) {
  double delta_t_2 = delta_t * delta_t; // delta_t_2 will be used in several places later

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
      px_p = p_x + v / yawd * (sin (yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t_2 * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t_2 * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t_2;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    NormalizeAngle(&x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::PredictRadarMeasurement(Ref<VectorXd> z_pred, Ref<MatrixXd> S, Ref<MatrixXd> Zsig, int n_z) {
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

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

  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    NormalizeAngle(&z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0,std_radrd_*std_radrd_;
  S = S + R;
}

void UKF::UpdateState(const Ref<const VectorXd> z, const Ref<const VectorXd> z_pred,
                      const Ref<const MatrixXd> S, const Ref<const MatrixXd> Zsig, int n_z) {
  // cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    NormalizeAngle(&z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    NormalizeAngle(&x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  NormalizeAngle(&z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
