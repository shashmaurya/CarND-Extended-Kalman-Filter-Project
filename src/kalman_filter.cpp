#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

#define PI 3.14159265

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
}


void KalmanFilter::Predict() {
    /**
     * Function to predict the state
     */
    x_ = F_ * x_;  
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}


void KalmanFilter::Update(const VectorXd &z) {
    /**
     * LASER: Update the state by using Kalman Filter equations
     */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;  
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    // New estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}


void KalmanFilter::UpdateEKF(const VectorXd &z){ //, const VectorXd &Hj_, const VectorXd &R_radar_) {
    /**
     * RADAR: Update the state by using Extended Kalman Filter equations
     */

    // Set R matrix
    MatrixXd R_radar_(3,3);
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

  // Get position and velocity for cartesian coordinates
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    VectorXd z_pred = VectorXd(3);  
    MatrixXd Hj(3,4);
  
  	// Use a tools instance and calculate Jacobian Matrix
    Tools tools;
    Hj = tools.CalculateJacobian(x_);

    // Calculate predicted state from cartesian coordinates  
    z_pred << pow((pow(px, 2) + pow(py, 2)), 0.5),
              atan2(py, px),
              (px*vx + py* vy)/(pow((pow(px, 2) + pow(py, 2)), 0.5)); 

    VectorXd y = z - z_pred;
    
    // Ensure theta is between -PI and PI
    float theta = y(1);
    if (theta<-PI){
        while(theta < -PI){
            theta = theta + (2 * PI);        
        }
    }
    else if (theta>PI){
        while(theta > PI){
            theta = theta - (2 * PI);        
        }
    }
    y(1) =theta;

    // Matrix operations
    MatrixXd Ht = Hj.transpose();
    MatrixXd S = Hj * P_ * Ht + R_radar_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;  

    // New Estimates
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj) * P_;  

}

