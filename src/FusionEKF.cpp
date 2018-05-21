#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
#if TWO_STEP_INIT
  is_initialized_radar_ = false;
  is_initialized_lidar_ = false;
  initial_state_ = VectorXd(4);
  initial_state_ << 0, 0, 0, 0;
#endif

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.H_ = MatrixXd(2, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  double radial_distance;
  double azimuth;
  double dT, dT2, dT3, dT4;
  double noise_ax = 9;
  double noise_ay = 9;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.P_ = MatrixXd(4, 4);
#if TWO_STEP_INIT
    ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 100, 0, // Initial covariance lower than with the one step init, since we do set an estimated value for the speeds.
      0, 0, 0, 100;
#else
    ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;
#endif


#if TWO_STEP_INIT
    if ((is_initialized_lidar_ == false) && (is_initialized_radar_ == false)) { // this is the first data we have, store it in the temporary initial state vector
      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        radial_distance = measurement_pack.raw_measurements_(0);
        azimuth = measurement_pack.raw_measurements_(1);
        initial_state_(0) = radial_distance * cos(azimuth);
        initial_state_(1) = radial_distance * sin(azimuth);
        is_initialized_radar_ = true;
      }
      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */
        initial_state_(0) = measurement_pack.raw_measurements_(0);
        initial_state_(1) = measurement_pack.raw_measurements_(1);
        is_initialized_lidar_ = true;
      }
    }
    else { // We already got one set of measurements from a sensor, with the next one we can initialize all states in the state vector, including the velocities
      dT = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0F;
      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        radial_distance = measurement_pack.raw_measurements_(0);
        azimuth = measurement_pack.raw_measurements_(1);
        ekf_.x_(0) = radial_distance * cos(azimuth);
        ekf_.x_(1) = radial_distance * sin(azimuth);
        ekf_.x_(2) = (ekf_.x_(0) - initial_state_(0)) / dT;
        ekf_.x_(3) = (ekf_.x_(1) - initial_state_(1)) / dT;
      }
      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */
        ekf_.x_(0) = measurement_pack.raw_measurements_(0);
        ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        ekf_.x_(2) = (ekf_.x_(0) - initial_state_(0)) / dT;
        ekf_.x_(3) = (ekf_.x_(1) - initial_state_(1)) / dT;
      }
      // done initializing, no need to predict or update
      is_initialized_ = true;
    }

#else
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      radial_distance = measurement_pack.raw_measurements_(0);
      azimuth = measurement_pack.raw_measurements_(1);
      ekf_.x_(0) = radial_distance * cos(azimuth);
      ekf_.x_(1) = radial_distance * sin(azimuth);
  }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
#endif

    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /* Find delta time in seconds, and for the sake of simplicity compute it's powers for the process noise matrix */
  dT = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0F;
  dT2 = dT * dT;
  dT3 = 0.5F * (dT2 * dT);
  dT4 = 0.25F * (dT2 * dT2);
  /* Set the state transition matrix for the current measurement */
  ekf_.F_ <<
    1, 0, dT, 0,
    0, 1, 0, dT,
    0, 0, 1, 0,
    0, 0, 0, 1;
  /* Set the process noise matrix for the current measurement */
  ekf_.Q_ << (dT4 * noise_ax), 0, (dT3 * noise_ax), 0,
    0, (dT4 * noise_ay), 0, (dT3 * noise_ay),
    (dT3 * noise_ax), 0, (dT2 * noise_ax), 0,
    0, (dT3 * noise_ay), 0, (dT2 * noise_ay);


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = MatrixXd(3, 3);
    ekf_.R_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
      0, 1, 0, 0;
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ << 0.0225, 0,
               0, 0.0225;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
