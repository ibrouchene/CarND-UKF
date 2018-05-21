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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5; // 30 m/s^2 is way too high for an acceleration noise value

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.01; // 30 rad/s^2 is way too high for a yaw rate rate value
  
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
  
  n_x_ = 5; 

  n_aug_ = 7;

  n_augmented_mat_dim = (2 * n_aug_) + 1;

  lambda_ = 3 - n_x_; // Suggested value from the lesson

  Xsig_pred_ = MatrixXd(n_x_, n_augmented_mat_dim);

  weights_ = VectorXd(n_augmented_mat_dim);
  is_initialized_ = false;

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
  double radial_distance, azimuth, dT;
  /* Is the tracker initialized? */
  if (!is_initialized_) {
    /* Initialize the covariance matrix as the identity matrix */
    P_ = MatrixXd::Identity(5, 5);
    /* Depending on the first received input data, set the initial state vector [px, py, v, yaw, yawrate]*/
    x_ << 1, 1, 1, 1, 0;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      radial_distance = meas_package.raw_measurements_(0);
      azimuth = meas_package.raw_measurements_(1);
      x_(0) = radial_distance * cos(azimuth);
      x_(1) = radial_distance * sin(azimuth);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    /* Initialization is over */
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
  }
  else {

    /* Get delta time */
    dT = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0F;

    /* ***** Predict ***** */
    Prediction(dT);

    /* ***** Update ***** */
    /* Step 1- Predict Measurement */
    /* Step 2- Update State */


  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /* Step 1- Generate Sigma Points */
  MatrixXd X_sigma_points = MatrixXd(n_aug_, n_augmented_mat_dim);
  GenerateSigmaPoints(&X_sigma_points);
  /* Step 2- Predict Sigma Points */
  PredictSigmaPoints(delta_t, X_sigma_points);
  /* Step 3- Predict Mean and Covariance */
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
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}


/**
 * Generates the sigma points, with state augmentation.
 * @param {X_sigma_points} pointer to the matrix in which the sigma points are stored
 */
void UKF::GenerateSigmaPoints(MatrixXd* X_sigma_points) {
  /* Create augmented vectors and matrices */  
  VectorXd x_augmented = VectorXd(7);
  MatrixXd P_augmented = MatrixXd(7, 7);

  /* Create augmented state vector */
  x_augmented.head(5) = x_;
  x_augmented(5) = 0;
  x_augmented(6) = 0;

  /* Create augmented covariance matrix */
  P_augmented.fill(0.0);
  P_augmented.topLeftCorner(5, 5) = P_;
  P_augmented(5, 5) = std_a_*std_a_;
  P_augmented(6, 6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_augmented.llt().matrixL();

  //create augmented sigma points
  X_sigma_points->col(0) = x_augmented;
  for (int i = 0; i< n_aug_; i++)
  {
    X_sigma_points->col(i + 1) = x_augmented + sqrt(lambda_ + n_aug_) * L.col(i);
    X_sigma_points->col(i + 1 + n_aug_) = x_augmented - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

/**
* Generates the sigma points, with state augmentation.
* @param {delta_time} elapsed time between the current cycle and the previous one
* @param {X_augmented_sigma_points} the generated sigma points for the augmented state
*/
void UKF::PredictSigmaPoints(double delta_time, MatrixXd X_augmented_sigma_points){

  for (int i = 0; i< n_augmented_mat_dim; i++)
  {
    //extract values for better readability
    double p_x = X_augmented_sigma_points(0, i);
    double p_y = X_augmented_sigma_points(1, i);
    double v = X_augmented_sigma_points(2, i);
    double yaw = X_augmented_sigma_points(3, i);
    double yawd = X_augmented_sigma_points(4, i);
    double nu_a = X_augmented_sigma_points(5, i);
    double nu_yawdd = X_augmented_sigma_points(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.0001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd*delta_time) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_time));
    }
    else {
      px_p = p_x + v*delta_time*cos(yaw);
      py_p = p_y + v*delta_time*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_time;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_time*delta_time * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_time*delta_time * sin(yaw);
    v_p = v_p + nu_a*delta_time;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_time*delta_time;
    yawd_p = yawd_p + nu_yawdd*delta_time;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

}
/**
 * Predicts the mean and covariance using the predicted sigma points.
 */
void UKF::PredictMeanAndCovariance(void) {
  
  /* Set weights */
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i<n_augmented_mat_dim; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  /* Predict the state mean using the sigma points */
  x_.fill(0.0);
  for (int i = 0; i < n_augmented_mat_dim; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  /* Predict the state covariance matrix using the sigma points */
  P_.fill(0.0);
  for (int i = 0; i < n_augmented_mat_dim; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (!((x_diff(3) < M_PI) && (x_diff(3) > -M_PI))) { // Put angle in range
      if (x_diff(3) > M_PI) {
        x_diff(3) -= (2.0F*M_PI);
      }
      else if (x_diff(3) < -M_PI) {
        x_diff(3) += (2.0F*M_PI);
      }
    }

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}