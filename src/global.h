#pragma once

#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>  
#include <unistd.h>

#include<algorithm> 
#include <stdio.h>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <queue>
#include <assert.h>
#include <cmath>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <sophus/se3.h>
#include <sophus/so3.h>



namespace global {
    
const static int NUM_FULL_STATE_ = 15;
const static int QUALITY_THRES_ = 1e3; 
const static int NUM_OF_STATE = 15;
const static int NUM_OF_GNSS_MEAS = 6;
const static int NUM_STATE_BUFFER = 512;

const static Eigen::Matrix<double, 3, 1> G(0, 0, 9.81);	///< 重力
const static Eigen::Matrix<double, 3, 1> L_b_ba(0, 0, 0); 	///< 杆臂

// R
const static double POS_N = 10;
const static double VEL_N = 0.5;

// Q
const static double ACC_N = 0.2;
const static double ACC_W = 0.0002 ;
const static double GYR_N = 0.02;
const static double GYR_W = 2.0e-5;







// acc_n: 0.2          # accelerometer measurement noise standard deviation. #0.2
// gyr_n: 0.02         # gyroscope measurement noise standard deviation.     #0.05
// acc_w: 0.0002         # accelerometer bias random work noise standard deviation.  #0.02
// gyr_w: 2.0e-5       # gyroscope bias random work noise standard deviation.     #4.0e-5
// g_norm: 9.81007     # gravity magnitude

// const static double ACC_N = 200 * micro_g_to_meters_per_second_squared;
// const static double ACC_W = 1.0e-7;
// const static double GYR_N = 0.02 * DEG2RAD / 60;
// const static double GYR_W = 2.0e-12;

}

