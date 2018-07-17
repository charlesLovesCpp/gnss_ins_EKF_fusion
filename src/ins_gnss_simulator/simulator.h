#include <ros/ros.h>

#include <stdio.h>
#include <iostream>
#include <cstring>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>

#include <sensor_msgs/Imu.h>
#include <std_msgs/String.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/PointCloud.h>
#include <geometry_msgs/PointStamped.h>
#include <integrated_nav/GNSS_ECEF.h>
#include <cv_bridge/cv_bridge.h>

#include "../global.h"
#include "../utility.h"

namespace Simulator { 

///< vector< time, pos+vel>
typedef std::vector< std::pair<int, std::vector<double>> > GnssInput;
    
struct GnssConfig {
    double epoch_interval;
    double r_os;
    double inclination;
    double const_delta_lambda;
    double const_delta_t;
    double mask_angle;
    double SIS_err_SD;
    double zenith_iono_err_SD;
    double zenith_trop_err_SD;
    double code_track_err_SD;
    double rate_track_err_SD;
    double rx_clock_offset;
    double rx_clock_drift;
    int num_sat;
    Eigen::Vector3d init_est_r_ea_e;
    
    GnssConfig():
	epoch_interval(0.5),
	r_os(2.656175e7),
	inclination(55),
	const_delta_lambda(0),
	const_delta_t(0),
	mask_angle(10),
	SIS_err_SD(1),
	zenith_iono_err_SD(2),
	zenith_trop_err_SD(0.2),
	code_track_err_SD(1),
	rate_track_err_SD(0.02),
	rx_clock_offset(10000),
	rx_clock_drift(100),
	num_sat(30)
    {
       init_est_r_ea_e.setZero();
    }   
    
    void setDefault() {
	
    }
};

struct ImuParameter {
    Eigen::Vector3d b_a;
    Eigen::Vector3d b_g;
    Eigen::Matrix3d M_a;
    Eigen::Matrix3d M_g;
    Eigen::Matrix3d G_g;
    double accel_noise_root_PSD;
    double gyro_noise_root_PSD;
    double accel_quant_level;
    double gyro_quant_level;
    
    ImuParameter():
	b_a(Utility::micro_g_to_meters_per_second_squared * Eigen::Vector3d(900, -1300, 800)),
	b_g(Utility::DEG2RAD / 3600 * Eigen::Vector3d(-9, 13, -8)),
	accel_noise_root_PSD(100 * Utility::micro_g_to_meters_per_second_squared),
	gyro_noise_root_PSD(0.01 * Utility::DEG2RAD / 60),
	accel_quant_level(1e-2),
	gyro_quant_level(2e-4)
    {
	M_a << 500, -300, 200,
	      -150, -600, 250,
	      -250,  100, 450;
	M_a *= 1e-6;
	      
	M_g << 400, -300,  250,
		0,  -300, -150,
		0,    0,  -350;
	M_g *= 1e-6;
	
	G_g << 0.9, -1.1, -0.6,
	      -0.5,  1.9, -1.6,
	       0.3,  1.1, -1.3;      
	G_g *= Utility::DEG2RAD / (3600 * 9.80665);  
    }
    
    void setDefault() {
	
    }
};

struct GnssMeasurement {
    double range;
    double range_rate;
    Eigen::Vector3d sat_p;
    Eigen::Vector3d sat_v;
};

struct GnssResult {
    double time;
    Eigen::Matrix3d R;
    Eigen::Vector3d p;
    Eigen::Vector3d v;
    
    void clear() {
	R.setZero();
	p.setZero();
	v.setZero();
    }    
};

struct ImuMsg {
    bool status;
    Eigen::Vector3d f_ib_b;
    Eigen::Vector3d omega_ib_b;   
    
    void clear() {
	status = 0;
	f_ib_b.setZero();
	omega_ib_b.setZero();
    }
    
    void set(const Eigen::Vector3d& _f_ib_b, const Eigen::Vector3d& _omega_ib_b) 
    {
	status = 1;
	f_ib_b = _f_ib_b;
	omega_ib_b = _omega_ib_b;
    }
};

struct GnssMsg {
    bool status;
    GnssResult ecef;
    double position_acc; 
    double speed_acc;
    int num_SV;
    int p_DOP;
    int gps_fix_status;   
    
    void clear() {
	status = 0;
	ecef.clear();
	position_acc = 0;
	speed_acc = 0;
	num_SV = 0;
	p_DOP = 0;
	gps_fix_status = 0;
    }   
    
    void set(
	const GnssResult& _ecef, 
	const double& _position_acc, 
	const double& _speed_acc, 
	const int& _num_SV,
	const int& _p_DOP = 0,
	const int& _gps_fix_status = 3) 
    {
	status = 1;
	ecef = _ecef;
	position_acc = _position_acc;
	speed_acc = _speed_acc;
	num_SV = _num_SV;
	p_DOP = _p_DOP;
	gps_fix_status = _gps_fix_status;
    }
};

class GnssSimulator
{
public:
    
    unsigned seed;
    std::default_random_engine generator_;	       
    
    ImuParameter imu_parameter_;
    GnssConfig gnss_config_;
    
    std::vector<GnssResult> ned_;
    std::vector<GnssResult> ecef_;
    std::vector<GnssResult> est_ecef_;
    std::vector<Eigen::Vector2d> est_clock_;  
    
    ImuMsg imu_msg_;
    GnssMsg gnss_msg_;
    bool hasImu;
    bool hasGnss;
    
    GnssSimulator(ros::NodeHandle& nh, ros::NodeHandle& pnh);
    
    void run();
    
    ///< Read .csv for ground true [time, position, velocity]
    bool readFile(const std::string& file_name);
    
    void ned2ecef(const GnssResult& ned, GnssResult& ecef);
    
    void neds2ecefs(const std::vector<GnssResult>& ned, std::vector<GnssResult>& ecef);
    
    void pubImu(const double& time, const Eigen::Vector3d& _f_ib_b, const Eigen::Vector3d& _omega_ib_b);
    
    void pubGnss(
	const GnssResult& ecef, 
	const double& position_acc, 
	const double& speed_acc, 
	const int& num_SV,
	const int& p_DOP = 0,
	const int& gps_fix_status = 3);
    
    void imuCallback(const sensor_msgs::ImuConstPtr& msg);
    
    void reset();
    
    ///< Compute the Satellite positions at specific time
    bool solveSvPV(
	const double& time, 
	const GnssConfig& config,
	std::vector<Eigen::Vector3d>& sats_p, 
	std::vector<Eigen::Vector3d>& sats_v);
    
    
    void generateGnssMeasurements(    
	const double& time, 
	const GnssConfig& config,   
	const std::vector<Eigen::Vector3d>& sats_p, 
	const std::vector<Eigen::Vector3d>& sats_v,
	const Eigen::Vector3d& r_ea_e,
	const Eigen::Vector3d& v_ea_e,
	const double& L_a,
	const double& lambda_a,
	const std::vector<double>& gnss_biases,
	std::vector<GnssMeasurement>& gnss_measurements,
	int& num_gnss_meas);

    void initialGnssBias(	
	const GnssConfig& config,   
	const std::vector<Eigen::Vector3d>& sats_p, 
	const Eigen::Vector3d& r_ea_e,
	const double& L_a,
	const double& lambda_a, 
	std::vector<double>& gnss_biases);
    
    
    
    bool solveReceiverPV(
	const std::vector<GnssMeasurement>& gnss_measurements,
	const int& num_gnss_meas,
	const Eigen::Vector3d& predicted_r_ea_e,
	const Eigen::Vector3d& predicted_v_ea_e,
	Eigen::Vector3d& est_r_ea_e,
	Eigen::Vector3d& est_v_ea_e,
	Eigen::Vector2d& est_clock);
    
    ///< Derive true IMU measurement by ground true using IMU kinematic model
    void deriveImuMeasurement(
	const double& tor_i, 
	const Eigen::Matrix<double, 3, 3>& C_b_e,
	const Eigen::Matrix<double, 3, 3>& old_C_b_e,
	const Eigen::Vector3d& v_eb_e,
	const Eigen::Vector3d& old_v_eb_e,
	const Eigen::Vector3d& r_eb_e,
	Eigen::Vector3d& f_ib_b,
	Eigen::Vector3d& omega_ib_b
    );
    
    void perturbeImuMeasurement(
	const double& tor_i, 
	const ImuParameter& imu_parameter,
	const Eigen::Vector3d& true_f_ib_b,
	const Eigen::Vector3d& true_omega_ib_b,
	const Eigen::Matrix<double, 6, 1>& old_quant_residuals,
	Eigen::Vector3d& meas_f_ib_b,
	Eigen::Vector3d& meas_omega_ib_b,
	Eigen::Matrix<double, 6, 1>& quant_residuals);
    
    // 用求来的卫星位置,初始化GNSS的误差
    // 用这一误差,生成伪距测量值
    // 用LS法计算当前的ECEF位置
    // 将结果转换成NED参考系表示,作为结果输出
    // 将file里的真实NED坐标,转换成ECEF坐标,然后用它生成真实的IMU测量值
    // 将真实IMU值加上噪声干扰,得到测量值
    // 加上ros::Time, 发送出去
    // 当时间间隔超过一定门限时,判断当前卫星的位置,然后根据此生成伪距测量值
    // 用LS法计算当前ECEF位置,加上ros::Time,发送出去
    
    
private:
    ros::NodeHandle it;
    ros::Publisher pub_imu_;
    ros::Publisher pub_gnss_;
    ros::Subscriber sub_imu_;
    
};   

}
