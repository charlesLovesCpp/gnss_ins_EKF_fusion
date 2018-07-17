#pragma once

#include "global.h"
#include <ros/ros.h>

namespace Fusion {
   
class State {
public:
    
    // State
    Eigen::Matrix<double, 3, 1> p_;         	///< position (IMU centered)          (0-2 / 0-2)
    Eigen::Matrix<double, 3, 1> v_;        	///< velocity                         (3- 5 / 3- 5)
    Eigen::Quaterniond q_;        		///< attitude                         (6- 9 / 6- 8)
    Eigen::Matrix<double, 3, 1> b_a_;      	///< acceleration biases              (10-12 / 9-11)
    Eigen::Matrix<double, 3, 1> b_w_;       	///< gyro biases                      (13-15 / 12-14)
    
    // IMU measurement
    Eigen::Matrix<double, 3, 1> f_b_ib;		///< acceleration measurement
    Eigen::Matrix<double, 3, 1> w_b_ib;		///< angular rate measurement
    
    // Process Covariance
    Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE> P_;
    
    // Time
    double time_;
    ros::Time ros_time;
    
    State():
	time_(-1)
    {}
    
    void reset() {
	p_.setZero();
	v_.setZero();
	q_.setIdentity();
	b_w_.setZero();
	b_a_.setZero();
	
	f_b_ib.setZero();
	w_b_ib.setZero();
	
	P_.setZero();
	time_ = -1;	
    }
};
    
}