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

  //* radar measurement dimension
  n_zrad_ = 3;

  //* radar measurement dimension
  n_zlas_ = 2;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

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
  std_yawdd_ = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.35;

  ///* Weights of sigma points
  weights_ = VectorXd::Zero(2*n_aug_+1);

  ///* the current NIS for radar
  NIS_radar_ = 0.0;

  ///* the current NIS for laser
  NIS_laser_ = 0.0;

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  return;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage meas_package) {
  
  if (!is_initialized_) {

    //  Initialize the state x_,state covariance matrix P_
    SetInitialValues(meas_package);
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }


  double delta_t =  (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  while (delta_t > 0.1)
  {
//  cout<< "delta_t: " << delta_t << endl;
//  cout << "Press ENTER to continue...";
//  cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
  
  const double dt = 0.05;
  Prediction(dt);
  delta_t -= dt;
  }

  Prediction(delta_t);


  if ( meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ) {
    // Radar updates 
    // rho, phi, rho_dot

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_zrad_);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_zrad_,n_zrad_);
    // cross-correlation matrix Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_zrad_);
    // get predictions for x,S and Tc in RADAR space
    PredictRadarMeasurement(z_pred, S, Tc);  
    // update the state using the RADAR measurement
    UpdateRadar(meas_package, z_pred, Tc, S);
    // update the time
    previous_timestamp_ = meas_package.timestamp_;

  } 

  else if ( meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ ) {
    // Laser updates
    // px, py

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_zlas_);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_zlas_,n_zlas_);
    // cross-correlation matrix Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_zlas_);
    // get predictions for x,S and Tc in Lidar space
    PredictLidarMeasurement(z_pred, S, Tc);
    // update the state using the LIDAR measurement
    UpdateLidar(meas_package, z_pred, Tc, S);
    // update the time
    previous_timestamp_ = meas_package.timestamp_;

  }

  return;

}

/**
 * @param {MeasurementPackage} meas_package. Initializes the state and
 * the covariance matrix P
 */
void UKF::SetInitialValues(const MeasurementPackage meas_package) {

  // initialize state covariance matrix P 
  P_ = MatrixXd(5, 5);
  P_ <<   1, 0,  0,  0,  0,
          0,  1, 0,  0,  0,
          0,  0,  1, 0,  0,
          0,  0,  0,  1, 0,
          0,  0,  0,  0,  1;

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
 * @param MatrixXd &Xsig_out. Computes augmented sigma points in state space
 */
void UKF::AugmentedSigmaPoints(MatrixXd &Xsig_out){

  //create augmented mean vector
  static VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  static MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  static MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

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
 * @param MatrixXd &Xsig_out. predicts augmented sigma points
 */
void UKF::SigmaPointPrediction(const MatrixXd &Xsig_aug, const double delta_t, MatrixXd  &Xsig_out){
  
    //create matrix with predicted sigma points as columns
  static MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

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
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  Xsig_out = Xsig_pred;

  return;

}


void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd &x_out, MatrixXd &P_out) {
 
  //create vector for predicted state
  static VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  static MatrixXd P = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;

  }

  //print result
  //std::cout << "Predicted state px,py,v,psi,psi_dot" << std::endl;
  //std::cout << x << std::endl;
  //std::cout << "Predicted covariance matrix P " << std::endl;
  //std::cout << P << std::endl;

  //write result
  x_out = x;
  P_out = P;

  return;
}

void UKF::PredictRadarMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Tc_out) {

  //radar can measure r, phi, and r_dot

  //create matrix for sigma points in measurement space
  static MatrixXd Zsig = MatrixXd(n_zrad_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    // Avoid division by zero
    if(fabs(p_x) <= 0.0001){
            p_x = 0.0001;
    }
    if(fabs(p_y) <= 0.0001){
            p_y = 0.0001;
    }


    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  static VectorXd z_pred = VectorXd(n_zrad_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  static MatrixXd S = MatrixXd(n_zrad_,n_zrad_);
  S.fill(0.0);
  //create matrix for cross correlation Tc
  static MatrixXd Tc = MatrixXd(n_x_, n_zrad_);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }

  //add measurement noise covariance matrix
  static MatrixXd R = MatrixXd(n_zrad_,n_zrad_);
  R <<    pow(std_radr_,2), 0, 0,
          0, pow(std_radphi_,2), 0,
          0, 0, pow(std_radrd_,2);
  S = S + R;

  //write result
  z_out = z_pred;
  S_out = S;
  Tc_out = Tc;

  return;

}


void UKF::PredictLidarMeasurement(VectorXd &z_out, MatrixXd &S_out, MatrixXd &Tc_out) {

  //create matrix for sigma points in measurement space
  static MatrixXd Zsig = MatrixXd(n_zlas_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // measurement model
    Zsig(0,i) = Xsig_pred_(0,i);          //px
    Zsig(1,i) = Xsig_pred_(1,i);          //py
    
  }

  //mean predicted measurement
  static VectorXd z_pred = VectorXd(n_zlas_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  static MatrixXd S = MatrixXd(n_zlas_,n_zlas_);
  S.fill(0.0);
  //create matrix for cross correlation Tc
  static MatrixXd Tc = MatrixXd(n_x_, n_zlas_);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }

  //add measurement noise covariance matrix
  static MatrixXd R = MatrixXd(n_zlas_,n_zlas_);
  R <<    pow(std_laspx_,2), 0,
          0, pow(std_laspy_,2);
  S = S + R;

  //write result
  z_out = z_pred;
  S_out = S;
  Tc_out = Tc;

  return;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimates the object's location. Predicts the state
  vector, x_, the sigma points, then the augmented sigma points, and the state covariance matrix.
  All results are in the state vector coordinate system
  */
  static MatrixXd Xsig_aug = MatrixXd(2*n_aug_ + 1, n_aug_);
 
  //create example matrix with predicted sigma points
  static MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // compute augmented sigma points
  AugmentedSigmaPoints(Xsig_aug); 
  // predict augmented sigma points
  SigmaPointPrediction(Xsig_aug, delta_t, Xsig_pred);

  static VectorXd x_pred = VectorXd(n_x_);
  static MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  PredictMeanAndCovariance(Xsig_pred,x_pred, P_pred);

  x_ = x_pred;
  P_ = P_pred;
  Xsig_pred_ = Xsig_pred;

  return;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &Tc, MatrixXd &S) {
  /**
  Uses lidar data to update the belief about the object's
  position. Updates the state vector, x_, and covariance, P_, and lidar NIS.
  */

  //mean predicted measurement
  static VectorXd z = VectorXd::Zero(n_zlas_);
  z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);

  //Kalman gain K;
  static MatrixXd K = MatrixXd::Zero(n_x_,n_zlas_);
  K = Tc * S.inverse();

  //residual
  static VectorXd z_diff = VectorXd::Zero(n_zlas_);
  z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  return;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &Tc, MatrixXd &S) {
  /**
  Uses radar data to update the belief about the object's
  position. Updates the state vector, x_, and covariance, P_ and radar NIS.
  */

  //mean predicted measurement  rho, phi, rho_dot
  static VectorXd z = VectorXd::Zero(n_zrad_);
  z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),meas_package.raw_measurements_(2);

  //Kalman gain K;
  static MatrixXd K = MatrixXd::Zero(n_x_,n_zrad_);
  K = Tc * S.inverse();

  //residual
  static VectorXd z_diff = VectorXd::Zero(n_zrad_);
  z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  return;

}
