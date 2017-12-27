#include "ukf.h"
#include <iostream>
#include <iomanip>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = n_x_ + 2; // new dimensions are std_a_ and std_yawdd

  // measurement dimension, radar can measure r, phi, and r_dot
  n_r_ = 3;

  // measurement dimension, lidar can measure px and py
  n_l_ = 2;

  // sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //
  // Based on the "Parkin, J. and Rotheram, J. (2010) Design speeds and acceleration characteristics of bicycle
  // traffic for use in planning, design and appraisal. Transport Policy, 17 (5). pp. 335-341. ISSN 0967-070X
  // Available from: http://eprints.uwe.ac.uk/20767" Table 1 "Summary of gradient, speed and acceleration data",
  // the maximum acceleration is 0.71. According to the rule of thumb, the process noise standard deviation
  // longitudinal acceleration should be a half of the maximum acceleration value, so we set it to 0.355.
  std_a_ = 0.355;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
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

  // Measurement noise covariance matrix for radar
  R_ = MatrixXd(n_r_, n_r_);
  R_ <<    std_radr_ * std_radr_,                         0,                       0,
                               0, std_radphi_ * std_radphi_,                       0,
                               0,                         0, std_radrd_ * std_radrd_;

  // Measurement noise covariance matrix for lidar
  L_ = MatrixXd(n_l_, n_l_);
  L_ << std_laspx_ * std_laspx_,                       0,
                              0, std_laspy_ * std_laspy_;
  
  // Flag indicating whether the first measurement processed
  is_initialized_ = false;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // weights of sigma-points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // number of a measurement being processed
  counter_ = 0;

  // instance of Tools class with is a collection of helper methods
  tools_ = Tools();

  // normalized innovation squared for lidar (set to -1 to indicate that NIS has not been measured yet)
  NIS_l_ = -1;

  // normalized innovation squared for radar (set to -1 to indicate that NIS has not been measured yet)
  NIS_r_ = -1;

  // file to write NIS values into for latter analysis
  NIS_data_file_.open( "NIS_data.csv", ios::out );
  if (!NIS_data_file_.is_open()) {
    std::cerr << "failed to open NIS_data.csv file" << std::endl;
    exit(1);
  }

  // write column headers
  NIS_data_file_ << "SENSOR_TYPE NIS" << std::endl;
}

/**
 * Closes opened file.
 */
UKF::~UKF() {
  NIS_data_file_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    // the first measurement is handled here
    double px, py, v, yaw, yawd; // these variables will be initialized in radar or lidar conditional branch

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // convert from polar to cartesian coordinate system
      double meas_rho     = meas_package.raw_measurements_[0];
      double meas_phi     = meas_package.raw_measurements_[1];

      px   = meas_rho     * cos(meas_phi);
      py   = meas_rho     * sin(meas_phi);
      v    = 0;                 // not enough measurements to determine
      yaw  = M_PI_2 - meas_phi; // phi is measured relative to Y axis while yaw is relative to X axis
      yawd = 0;                 // not enough measurements to determine
    }
    else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // initial state in case the first measurement comes from lidar sensor
      px   = meas_package.raw_measurements_[0];
      py   = meas_package.raw_measurements_[1];
      v    = 0; // not enough measurements to determine
      yaw  = 0; // not enough measurements to determine
      yawd = 0; // not enough measurements to determine
    } else {
      std::cerr << "unknown sensor type of measurement or all sensors are turned off"
                << std::endl;
      exit(2);
    }

    // initial state vector
    x_ << px, py, v, yaw, yawd;

    // initial state covariance matrix
    P_.fill(0.0);
    P_.diagonal().setOnes();

    // done initializing, no need to predict or update
    is_initialized_ = true;
  } else {
    // difference in seconds between the current measurement and the previous one
    double delta_t_s = (meas_package.timestamp_ - time_us_) / 1000000.0; // us -> s

    // the first step in the UKF processing chain (prediction)
    // the second step in the UKF processing chain (update)
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      Prediction(delta_t_s);
      UpdateRadar(meas_package);
    } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      Prediction(delta_t_s);
      UpdateLidar(meas_package);
    } else {
      std::cerr << "unknown sensor type of measurement or all sensors are turned off"
                << std::endl;
      exit(3);
    }

    // write normalized innovation squared value to file
    NIS_data_file_ << meas_package.sensor_type_ << ' '
                   << (meas_package.sensor_type_ == MeasurementPackage::LASER ? NIS_l_ : NIS_r_)
                   << std::endl;
  }

  // we do not want to do any updates in cases <use_laser=false, sensor=LASER> and <use_radar=false, sensor=RADAR>
  if ( (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) ||
       (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) ) {
    // initialize the time at which the measurement was taken
    time_us_ = meas_package.timestamp_;
    ++counter_; // increase counter of processed measurements

    // print statistics and intermediate results
    std::cout << "====================\n"
              << "Measurement number = " << counter_ << "\n"
              << "Measurement type   = " <<
                  (meas_package.sensor_type_ == MeasurementPackage::LASER ? "LASER": "RADAR") << "\n"
              << "Time (us)          = " << time_us_ << "\n"
              << "NIS " << (meas_package.sensor_type_ == MeasurementPackage::LASER ?
                            "lidar          = " : "radar          = ")
                        << (meas_package.sensor_type_ == MeasurementPackage::LASER ? NIS_l_ : NIS_r_) << "\n"
              << "x                  =\n" << x_ << "\n"
              << "P                  =\n" << P_ << "\n"
              << "====================\n" << std::endl;
  }
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
  // actual measurements
  VectorXd z = meas_package.raw_measurements_.head(n_l_);

  // predicted measurement mean (to be initialized)
  VectorXd z_pred(n_l_);

  // measurement covariance matrix S
  MatrixXd S(n_l_,n_l_);

  // matrix for sigma points in measurement space (to be initialized)
  MatrixXd Zsig(n_l_, 2 * n_aug_ + 1);

  // predict the lidar measurement and initialize z_pred, S, and Zsig
  PredictMeasurement(z_pred, S, Zsig, L_, n_l_, false);

  // finish the update step with the update of x_ and P_
  UpdateState(z, z_pred, S, Zsig, NIS_l_, n_l_, false);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // actual measurements
  VectorXd z = meas_package.raw_measurements_.head(n_r_);

  // predicted measurement mean (to be initialized)
  VectorXd z_pred(n_r_);

  // measurement covariance matrix S
  MatrixXd S(n_r_,n_r_);

  // matrix for sigma points in measurement space (to be initialized)
  MatrixXd Zsig(n_r_, 2 * n_aug_ + 1);

  // predict the radar measurement and initialize z_pred, S, and Zsig
  PredictMeasurement(z_pred, S, Zsig, R_, n_r_, true);

  // finish the update step with the update of x_ and P_
  UpdateState(z, z_pred, S, Zsig, NIS_r_, n_r_, true);
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
    tools_.NormalizeAngle(&x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::PredictMeasurement(Ref<VectorXd> z_pred, Ref<MatrixXd> S, Ref<MatrixXd> Zsig,
                             const Ref<const MatrixXd> M_meas_noise, int n_z, bool is_radar) {
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    if (is_radar) {
      // radar
      double v  = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);

      double v1 = cos(yaw)*v;
      double v2 = sin(yaw)*v;

      Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                            // rho
      Zsig(1, i) = atan2(p_y, p_x);                                        // phi
      Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);    // rho_dot
    } else {
      // lidar
      Zsig(0,i) = p_x;
      Zsig(1,i) = p_y;
    }
  }

  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (is_radar) {
      //angle normalization
      tools_.NormalizeAngle(&z_diff(1));
    }

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + M_meas_noise;
}

void UKF::UpdateState(const Ref<const VectorXd> z, const Ref<const VectorXd> z_pred,
                      const Ref<const MatrixXd> S, const Ref<const MatrixXd> Zsig,
                      double &nis, int n_z, bool is_radar) {
  // cross correlation matrix
  MatrixXd Tc(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (is_radar) {
      //angle normalization
      tools_.NormalizeAngle(&z_diff(1));
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    if (is_radar) {
      //angle normalization
      tools_.NormalizeAngle(&x_diff(3));
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  if (is_radar) {
    //angle normalization
    tools_.NormalizeAngle(&z_diff(1));
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // update normalized innovation squared value
  nis = z_diff.transpose() * S.inverse() * z_diff;
}
