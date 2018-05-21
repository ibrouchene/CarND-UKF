#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  Eigen::MatrixXd I = MatrixXd::Identity(4, 4); // 4 states: x, y, vx, vy
}

void KalmanFilter::Predict() {

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  MatrixXd I = MatrixXd::Identity(4, 4); // 4 states: x, y, vx, vy

  //new state
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  Tools tools;
  /* convert predicted state x_ from cartesian to polar coordinate system */
  MatrixXd h = MatrixXd(3, 1);
  h(0) = (x_(0)*x_(0) + x_(1)*x_(1));
  h(0) = sqrt(h(0));
  if (h(0) < 0.0001) { // division by zero check
    return;
  }
  h(1) = atan2(x_(1), x_(0));
  h(2) = (x_(0)*x_(2) + x_(1)*x_(3));
  h(2) = h(2) / h(0);
  
  VectorXd y = z - h;
  while (!((y(1) < M_PI) && (y(1) > -M_PI))) { // Put azimuth in range
    if (y(1) > M_PI) {
      y(1) -= (2.0F*M_PI);
    }
    else if (y(1) < -M_PI) {
      y(1) += (2.0F*M_PI);
    }
  }


  MatrixXd Hj = tools.CalculateJacobian(x_);

  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  MatrixXd I = MatrixXd::Identity(4, 4); // 4 states: x, y, vx, vy

                                         //new state
  x_ = x_ + (K * y);
  P_ = (I - K * Hj) * P_;
}

