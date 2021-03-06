#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  long long previous_timestamp_;

  int n_augmented_mat_dim_;

  std::ofstream output_file_;

  double nis_radar_;
  double nis_lidar_;

  Tools tools_;
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
  * Generates the sigma points, with state augmentation.
  * @param {X_sigma_points} pointer to the matrix in which the sigma points are stored
  */
  void GenerateSigmaPoints(MatrixXd* X_sigma_points);

  /**
  * Generates the sigma points, with state augmentation.
  * @param {delta_time} elapsed time between the current cycle and the previous one
  * @param {X_augmented_sigma_points} the generated sigma points for the augmented state
  */
  void PredictSigmaPoints(double delta_time, MatrixXd X_augmented_sigma_points);

  /**
  * Predicts the mean and covariance using the predicted sigma points.
  */
  void PredictMeanAndCovariance(void);

  /**
  * Predicts radar measurements
  * @param {z_out} Pointer to the vector storing the mean predicted measurement
  * @param {S_out} Pointer to the matrix storing the measurement covariance matrix S
  * @param {Z_sigma_points} Pointer to the matrix storing the sigma points transformed in the measurement space
  */
  void PredictRadarMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Z_sigma_points);

  /**
  * Predicts lidar measurements
  * @param {z_out} Pointer to the vector storing the mean predicted measurement
  * @param {S_out} Pointer to the matrix storing the measurement covariance matrix S
  * @param {Z_sigma_points} Pointer to the matrix storing the sigma points transformed in the measurement space
  */
  void PredictLidarMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Z_sigma_points);
};

#endif /* UKF_H */
