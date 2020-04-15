#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
     * Calculate the RMSE values
     */
  
    // Initialize rmse matrix
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // Include error checking
    if (estimations.size() == 0){
        cout<<"vector size is zero"<<endl;
    }
    else if (estimations.size() != ground_truth.size()){
        cout<<"vector sizes do not match"<< endl;
    } 
    else{
      
        // Initialize sum matrix; holds summation of residual squares
        VectorXd sum = rmse;

        // Accumulate squared residuals
        for (int i=0; i < estimations.size(); ++i) {
            VectorXd res = estimations[i] - ground_truth[i];
            res =  res.array()*res.array();
            sum = sum.array() + res.array();
        }
        // Calculate the mean
        VectorXd mean = sum.array()/ estimations.size();

        // Calculate the squared root
        rmse =  mean.array().sqrt();
    }
    // return the result
    return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /**
     * Calculate Jacobian
     */
    MatrixXd Hj(3,4);
    // Recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // Intial Hj matrix 
    Hj << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;  

    // check division by zero
    if ((px == 0) && (py == 0)){
        // Do nothing
        cout << "Divide by 0" << endl;
    } else{
        Hj(0,0) = px/pow((pow(px, 2) + pow(py, 2)), 0.5);
        Hj(0,1) = py/pow((pow(px, 2) + pow(py, 2)), 0.5);
        Hj(1,0) = -py/(pow(px, 2) + pow(py, 2));
        Hj(1,1) = px/(pow(px, 2) + pow(py, 2));
        Hj(2,0) = (py*((vx* py) - (vy*px)))/pow((pow(px, 2) + pow(py, 2)), 1.5);
        Hj(2,1) = (px*(- (vx* py) + (vy*px)))/pow((pow(px, 2) + pow(py, 2)), 1.5);
        Hj(2,2) = px/pow((pow(px, 2) + pow(py, 2)), 0.5);
        Hj(2,3) = py/pow((pow(px, 2) + pow(py, 2)), 0.5);

    }
    // Return the Jacobian matrix
    return Hj; 
}
