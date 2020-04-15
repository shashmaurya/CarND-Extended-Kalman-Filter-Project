#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "math.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

#define PI 3.14159265

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    // Set flag to intialize on the first call
    is_initialized_ = false;

    // Initialize previous timestamp
    previous_timestamp_ = 0;

    // Initializing matrices
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

    H_laser_ << 1.0, 0, 0, 0,
                0, 1.0, 0, 0;

    Hj_ << 0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0;

}




/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}




void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /**
    * Initialization
    */

    // Set the acceleration noise
    float noise_ax = 9.0;
    float noise_ay = 9.0;

    // Variable to hold con converted coordinates
    VectorXd cart_coord = VectorXd(4);

    if (!is_initialized_) {
      /**
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Convert radar from polar to cartesian coordinates.
         */

        // first measurement
        cout << "EKF: Execution Started" << endl;
        ekf_.x_ = VectorXd(4);

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // Convert radar from polar to cartesian coordinates 
            // and initialize state.
            // Set the state with the initial location and zero velocity
            float rho = measurement_pack.raw_measurements_[0];
            float theta =  measurement_pack.raw_measurements_[1];


            // Normalize theta  
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
            
            // Convert polar coordinates to cartesian for initialization
            cart_coord << rho* cos(theta),
                          rho* sin(theta),
                          0,
                          0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            // Store value to cartesian coordinate variable
            cart_coord << measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1], 
            0,
            0;
        }

        // Initialize state with first measurement 
        ekf_.x_ = cart_coord;

        // Intialize F Matrix
        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ << 1.0, 0, 0, 0,
                   0, 1.0, 0, 0,
                   0, 0, 1.0, 0,
                   0, 0, 0, 1.0;

        ekf_.H_ = MatrixXd(2, 4);
        ekf_.H_ << 1.0, 0, 0, 0,
                   0, 1.0, 0, 0;

        ekf_.R_ = MatrixXd(2, 2);
        ekf_.R_ = R_laser_;

        // Initialize the timestamp
        previous_timestamp_ = measurement_pack.timestamp_;

        // Create the process covariance matrix
        ekf_.Q_ = MatrixXd(4, 4);
        ekf_.Q_ << 1.0, 0, 0, 0,
                   0, 1.0, 0, 0,
                   0, 0, 1.0, 0,
                   0, 0, 0, 1.0;

        // Create the process covariance matrix
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ << 1.0, 0, 0, 0,
                   0, 1.0, 0, 0,
                   0, 0, 1000.0, 0,
                   0, 0, 0, 1000.0;

        // Done initializing, no need to predict or update this time
        is_initialized_ = true;
        return;
    }


    /**
     * Prediction
     */
  
    // Time elapsed between the current and previous measurements
    // dt - expressed in seconds
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  
    // Update the timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // Modify the F matrix so that the time is integrated
    ekf_.F_<< 1.0, 0, dt, 0,
              0, 1.0, 0, dt,
              0, 0, 1.0, 0,
              0, 0, 0, 1.0;

    // Update process covariance matrix
    ekf_.Q_ << (pow(dt, 4)*(noise_ax)/4.0), 0, (pow(dt, 3)*(noise_ax)/2.0), 0,
               0, (pow(dt, 4)*(noise_ay)/4), 0, (pow(dt, 3)*(noise_ay)/2.0),
               (pow(dt, 3)*(noise_ax)/2.0), 0, (pow(dt, 2)*(noise_ax)), 0,
               0, (pow(dt, 3)*(noise_ay)/2.0), 0, (pow(dt, 2)*(noise_ay));
    
    
    // Call function to predict the state
    ekf_.Predict();

  
    /**
     * Update
     */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } 
    else {
        // Laser updates
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}
