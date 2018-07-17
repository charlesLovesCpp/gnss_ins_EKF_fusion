#pragma once

#include <ros/ros.h>
#include <stdio.h>
#include <iostream>
#include <cstring>

#include <sensor_msgs/Imu.h>
#include <std_msgs/String.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/PointCloud.h>
#include <geometry_msgs/PointStamped.h>
#include <integrated_nav/GNSS_ECEF.h>
#include <cv_bridge/cv_bridge.h>

#include "global.h"
#include "state.h"	
#include "utility.h"

namespace Fusion {
    
class FusionHandler
{
public:
    typedef Eigen::Matrix<double, global::NUM_OF_STATE, 1> ErrorState;
    typedef Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE> ErrorStateCov;
    typedef const Eigen::Matrix<double, 4, 4> ConstMatrix4;
    typedef const Eigen::Matrix<double, 3, 3> ConstMatrix3;
    typedef const Eigen::Matrix<double, 3, 1> ConstVector3;
    typedef Eigen::Matrix<double, 4, 4> Matrix4;
    
    FusionHandler(ros::NodeHandle& nh, ros::NodeHandle& pnh);
    void run();
private:
    
    ros::NodeHandle it;
    ros::Subscriber sub_state_;
    ros::Subscriber sub_imu_;
    ros::Subscriber sub_image_;
    ros::Subscriber sub_gnss_;
    ros::Publisher pub_state_;     
    
    Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE> Fd_;
    Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE> Gd_;
    Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE> Qd_;
    Eigen::Matrix<double, global::NUM_OF_GNSS_MEAS, global::NUM_OF_GNSS_MEAS> CONST_R_;
    Eigen::Matrix<double, 3, 1> l_b_ba; 
    Eigen::Matrix<double, 3, 1> last_am;
    State state_buffer_[global::NUM_STATE_BUFFER]; 
    std::list<State> state_buffer;
    size_t idx_state_;
    size_t idx_P_;
    
    bool isInitialized_;
    bool use_fixed_covariance_;
    

    
    void initialize(
	const double& time,
	const ros::Time& ros_time,
	const Eigen::Vector3d& p0,
	const Eigen::Vector3d& v0);
    void initByGNSS(
	ros::Time time, 
	const Eigen::Matrix<double, 3, 1>& p,
	const Eigen::Matrix<double, 3, 1>& v,
	const double& pos_var,
	const double& vel_var);
    void propagateState(const double& dt);
    void predictProcessCovariance(const double& dt);
    bool correctState(const size_t& idx, const ErrorState& delta);
    size_t getClosestState(State& state, double tstamp, double delay);
    void propagatePtoIdx(const size_t& idx);
    template<class H_type, class Res_type, class R_type>
    bool applyMeasurement(
	const size_t& idx,
	const Eigen::MatrixBase<H_type>& H, 
	const Eigen::MatrixBase<Res_type> & res, 
	const Eigen::MatrixBase<R_type>& R);    
    
    void imuCallback(const sensor_msgs::ImuConstPtr& msg);
    void imageCallback(const sensor_msgs::ImageConstPtr& msg);
    void gnssCallback(const integrated_nav::GNSS_ECEF::ConstPtr& msg);
};    
    
    
}




 
//     State cur_state;
//     Eigen::Matrix<double, NUM_OF_STATE, NUM_OF_STATE> P_0;
//     Eigen::Matrix<double, NUM_OF_STATE, NUM_OF_STATE> Q_0;
//     Eigen::Matrix<double, NUM_OF_STATE, NUM_OF_STATE> R_0;
//     double gyro_noise_PSD;
//     double accel_noise_PSD;
//     double accel_bias_PSD;
//     double gyro_bias_PSD;
    			 ///< gravity vector