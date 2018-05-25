#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd ret(4);
  ret << 0, 0, 0, 0;
  for (unsigned int i = 0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    ret += residual;
  }

  //calculate the mean
  ret = ret / estimations.size();

  //calculate the squared root
  ret = ret.array().sqrt();

	return ret;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3, 4);
	//recover state parameters
	float px = (float)x_state(0);
	float py = (float)x_state(1);
	float vx = (float)x_state(2);
	float vy = (float)x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

	return Hj;
}

VectorXd Tools::ComputeNISStats(const vector<double> &nis_vector) {
  double sum = 0;
  double mean = 0;
  double residual = 0;
  double std = 0;
  double stat = 0;
  double count = 0;
  VectorXd ret(3);
  ret << 0, 0, 0;
  if (nis_vector.size() > 10) {
    for (unsigned int i = 0; i < nis_vector.size(); ++i) {
      sum += nis_vector[i];
      if (nis_vector[i] > 5.99) { // Change value to 7.8 for radar!
        count += 1;
      }
    }
    mean = sum / nis_vector.size();
    ret(1) = mean;
    for (unsigned int i = 0; i < nis_vector.size(); ++i) {
      residual += ((nis_vector[i] - mean)*(nis_vector[i] - mean));
    }
    std = sqrt(residual / mean);
    ret(2) = std;
    stat = count / (double)nis_vector.size();
    ret(0) = stat;
  }
  return ret;
}

double Tools::NormalizeAngle(double angle)
{
  while (angle> M_PI) angle -= 2.*M_PI;
  while (angle<-M_PI) angle += 2.*M_PI;
  return angle;
}
