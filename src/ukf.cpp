#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // Set initialization to false initially
  is_initialized_  = false;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3.0 - n_aug_;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // initialize the sigma points matrix with zeroes
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  ///* time when the state is true, in us
  previous_timestamp_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.1;

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

  ///* Weights of sigma points
  weights_ = VectorXd::Zero(2*n_aug_+1);

  ///* the current NIS for radar
  NIS_radar_ = 0.0;

  ///* the current NIS for laser
  NIS_laser_ = 0.0;


  /**
  See ukf.h for other member properties.
  one or more values initialized above might be wildly off...
  */

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  if (!is_initialized_) {

    //  Initialize the state x_,state covariance matrix P_
    SetInitialValues(meas_package);
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;

  }

  float delta_t =  (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  Prediction(delta_t);


  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);

  } else {
    // Laser updates
    UpdateLidar(meas_package);


  }



}


/**
 * @param {MeasurementPackage} meas_package. Initializes the state and
 * the covariance matrix P
 */
void UKF::SetInitialValues(MeasurementPackage meas_package) {

  // initialize state covariance matrix P 
  P_ = MatrixXd(5, 5);
  P_ <<   10, 0,  0,  0,  0,
          0,  10, 0,  0,  0,
          0,  0,  10, 0,  0,
          0,  0,  0,  10, 0,
          0,  0,  0,  0,  10;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    
    // Convert radar from polar to cartesian coordinates and initialize state.
    float rho = meas_package.raw_measurements_(0);
    float phi = meas_package.raw_measurements_(1);
    float rhodot = meas_package.raw_measurements_(2);

    float px = rho * cos(phi);
    float py = rho * sin(phi);
    float vx = rhodot * cos(phi);
    float vy = rhodot * sin(phi);

    x_ << px,py,sqrt(vx*vx +vy*vy),0,0;

  } 

  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

    // Initialize state.
    x_ << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0,0,0;

  }

  return;

}

 /**
 * @param MatrixXd &Xsig_out. Computes augmented sigma points
 */
void UKF::AugmentedSigmaPoints(MatrixXd &Xsig_out){

    //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = pow(std_a_,2);
  P_aug(6,6) = pow(std_yawdd_,2);

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  Xsig_out = Xsig_aug;

  return;  
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = MatrixXd(2*n_aug_ + 1, n_aug_);
  AugmentedSigmaPoints(Xsig_aug); // compute augmented sigma points
  


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
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
