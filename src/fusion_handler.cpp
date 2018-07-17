#include "fusion_handler.h"

namespace Fusion {

using namespace std;
using namespace Eigen;
    
FusionHandler::FusionHandler(ros::NodeHandle& nh, ros::NodeHandle& pnh):
    it(nh),
    Fd_(Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE>::Identity()),
    Qd_(Eigen::Matrix<double, global::NUM_OF_STATE, global::NUM_OF_STATE>::Zero()),
    idx_state_(0),
    idx_P_(0),
    isInitialized_(false),
    use_fixed_covariance_(false)
{
    sub_imu_ = it.subscribe("/imu0", 2000, &FusionHandler::imuCallback, this);
    sub_image_ = it.subscribe("/image", 1000, &FusionHandler::imageCallback, this);	
    sub_gnss_ = it.subscribe("/gnss/ecef", 1000, &FusionHandler::gnssCallback, this);
    pub_state_ = it.advertise<std_msgs::String> ("state_out", 1);  
}

void FusionHandler::run()
{

}

static int imu_count = 0;
static int gnss_count = 0;  
void FusionHandler::imuCallback(const sensor_msgs::ImuConstPtr& msg)
{
    if(!isInitialized_) {
	ROS_WARN_STREAM("imuCallback: FusionHandler has not been initialized yet. Wait for GNSS measurement");
	return;
    }
    imu_count ++;
    if (idx_state_ == global::NUM_STATE_BUFFER || idx_P_ == global::NUM_STATE_BUFFER ) {
	ROS_WARN_STREAM("imuCallback: State buffer is full");
	return;
    }
    state_buffer_[idx_state_].time_ = msg->angular_velocity_covariance[0];
    state_buffer_[idx_state_].ros_time = msg->header.stamp;
    state_buffer_[idx_state_].f_b_ib << msg->linear_acceleration.x,
					msg->linear_acceleration.y,
					msg->linear_acceleration.z;
    state_buffer_[idx_state_].w_b_ib << msg->angular_velocity.x,
					msg->angular_velocity.y,
					msg->angular_velocity.z;  
    
					
//     cout << "imu_count: " << imu_count << ", f_b_ib: " << state_buffer_[idx_state_].f_b_ib.transpose() << endl;

    if (state_buffer_[idx_state_].f_b_ib.norm() > 50) {
	state_buffer_[idx_state_].f_b_ib = last_am;
    } else {
	last_am = state_buffer_[idx_state_].f_b_ib;
    }
    
    double dt = state_buffer_[idx_state_].time_ - state_buffer_[idx_state_ - 1].time_;
//     cout << "current idx_state: " << idx_state_ << ", dt: " << dt << endl;
    propagateState(dt);	
    predictProcessCovariance(dt);

    cout << "later IMU time: " <<  state_buffer_[idx_state_- 1].time_ << ", ECEF: " << state_buffer_[idx_state_- 1].p_.transpose() << endl;
    
//     cout << " state_buffer_.q_: \n" <<  state_buffer_[idx_state_- 1].q_.coeffs().transpose() << endl;
//     cout << " state_buffer_.p_: \n" <<  state_buffer_[idx_state_- 1].p_.transpose() << endl;
//     cout << " state_buffer_.v_: \n" <<  state_buffer_[idx_state_- 1].v_.transpose() << endl;
    
    if(!Utility::checkForNumeric((double*)(&state_buffer_[idx_state_ - 1].p_[0]), 3, "prediction p")) {
	ROS_WARN_STREAM("checkForNumeric: Numeric problem");
    }
    
}

void FusionHandler::imageCallback(const sensor_msgs::ImageConstPtr& msg)
{
    
}

void FusionHandler::gnssCallback(const integrated_nav::GNSS_ECEF::ConstPtr& msg)
{
    gnss_count++;
    // Measured time
    ros::Time ros_time = msg->header.stamp;
    double time = msg->time;
    // ECEF position measurement
    Eigen::Matrix<double, 3, 1> pos_ECEF(msg->ecef_x, msg->ecef_y, msg->ecef_z);
    // ECEF velocity measurement
    Eigen::Matrix<double, 3, 1> vel_ECEF(msg->ecefV_x, msg->ecefV_y, msg->ecefV_z);  
    // Position Accuracy Estimate (3D)
    double p_var = pow(msg->position_acc, 2); 
    // Velocity Accuracy Estimate (3D)
    double v_var = pow(msg->speed_acc, 2); 
    // Position DOP
    double p_DOP = msg->p_DOP; 
    // Number of used satellites
    int num_SV = msg->num_SV;
    // Positioning status
    int gps_fix_status = msg->gps_fix_status;
    

    cout << "gnss_count: " << gnss_count << endl;
    cout << "pos_ECEF: " << pos_ECEF.transpose() <<endl; 
    cout << "vel_ECEF: " << vel_ECEF.transpose() <<endl; 
    cout << endl;
    
    if(!isInitialized_) {
	initialize(time, ros_time, pos_ECEF, vel_ECEF);
	return;
    }    
  
    Eigen::Matrix<double, global::NUM_OF_GNSS_MEAS, global::NUM_OF_STATE> H;
    Eigen::Matrix<double, global::NUM_OF_GNSS_MEAS, 1> res;
    Eigen::Matrix<double, global::NUM_OF_GNSS_MEAS, global::NUM_OF_GNSS_MEAS> R;
    H.setZero();
    res.setZero();
    R.setZero();
    
    if (!use_fixed_covariance_) {
	R.block<3, 3>(0, 0) = p_var * Eigen::Matrix<double, 3, 3>::Identity();
	R.block<3, 3>(3, 3) = v_var * Eigen::Matrix<double, 3, 3>::Identity();	
    } else {
	R = CONST_R_;
    }

//     cout << " *** current GNSS time: " <<  time << endl;
    
    // Return the state with closest time
    State state_measured;
    size_t idx = getClosestState(state_measured, time, 0);
    if (idx == 0) {
	ROS_WARN_STREAM("gnssCallback: No matched state");
	return;
    }
    
//     cout << "idx_state_: " << idx_state_ << ", gnss idx: " << idx << endl;
	
    Eigen::Matrix<double, 3, 1> p = state_measured.p_;
    Eigen::Matrix<double, 3, 1> v = state_measured.v_;
    
    H.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity();
    H.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
    
    res.block<3, 1>(0, 0) = pos_ECEF - p;
    res.block<3, 1>(3, 0) = vel_ECEF - v;
    
    if (!applyMeasurement(idx, H, res, R)) {
	ROS_WARN_STREAM("Wrong measurement!");
    }
    
}

void FusionHandler::initialize(
    const double& time,
    const ros::Time& ros_time,
    const Eigen::Vector3d& p0,
    const Eigen::Vector3d& v0
)
{
    isInitialized_ = false;
    for (size_t i = 0; i < global::NUM_STATE_BUFFER; i++) {
	state_buffer_[i].reset();
	
// 	cout << " state_buffer_.q_: \n" <<  state_buffer_[i].q_.coeffs().transpose() << endl;
// 	cout << " state_buffer_.p_: \n" <<  state_buffer_[i].p_.transpose() << endl;
// 	cout << " state_buffer_.v_: \n" <<  state_buffer_[i].v_.transpose() << endl;
    }    
    
    idx_state_ = 0;
    idx_P_ = 0;  

    ConstVector3 ypr(0, 0, 0);
    State& state = state_buffer_[idx_state_];   
    
    state.time_ = time;
    state.ros_time = ros_time;
    
    state.p_ = p0;
    state.v_ = v0;
    Eigen::Matrix3d R_test = Utility::ypr2R(ypr);
    cout << "R_test: \n" << R_test << endl;
    state.q_ = Utility::ypr2R(ypr);
    state.b_a_.setZero();
    state.b_w_.setZero();   

    const double init_pos_unc = 10;
    const double init_vel_unc = 0.1;
    const double init_att_unc = 1 * Utility::DEG2RAD;
    const double init_b_a_unc = 1000 * Utility::micro_g_to_meters_per_second_squared;
    const double init_b_g_unc = 10 * Utility::DEG2RAD / 3600;

    state.P_.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Identity() * init_pos_unc * init_pos_unc;
    state.P_.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity() * init_vel_unc * init_vel_unc;
    state.P_.block<3, 3>(6, 6) = Eigen::Matrix<double, 3, 3>::Identity() * init_att_unc * init_att_unc;
    state.P_.block<3, 3>(9, 9) = Eigen::Matrix<double, 3, 3>::Identity() * init_b_a_unc * init_b_a_unc;
    state.P_.block<3, 3>(12, 12) = Eigen::Matrix<double, 3, 3>::Identity() * init_b_g_unc * init_b_g_unc;    

    const double accel_noise_PSD = pow(200 * Utility::micro_g_to_meters_per_second_squared, 2);
    const double gyro_noise_PSD = pow(0.02 * Utility::DEG2RAD / 60, 2);
    const double accel_bias_PSD = 1.0e-7;	
    const double gyro_bias_PSD = 2.0e-12;
    
    Qd_.setZero();
    Qd_.block<3, 3>(3, 3) = accel_noise_PSD * Eigen::Matrix<double, 3, 3>::Identity();
    Qd_.block<3, 3>(6, 6) = gyro_noise_PSD * Eigen::Matrix<double, 3, 3>::Identity();
    Qd_.block<3, 3>(9, 9) = accel_bias_PSD * Eigen::Matrix<double, 3, 3>::Identity();
    Qd_.block<3, 3>(12, 12) = gyro_bias_PSD * Eigen::Matrix<double, 3, 3>::Identity();     
    
    const double pos_meas_SD = 2.5;
    const double vel_meas_SD = 0.1;
    
    CONST_R_.setZero();
    CONST_R_.block<3, 3>(0, 0) = pos_meas_SD * pos_meas_SD * Eigen::Matrix<double, 3, 3>::Identity();
    CONST_R_.block<3, 3>(3, 3) = vel_meas_SD * vel_meas_SD * Eigen::Matrix<double, 3, 3>::Identity();
    
    l_b_ba = ConstVector3(0, 0, 0);
    last_am = Eigen::Matrix<double, 3, 1>(0, 0, 0);
    
    idx_state_++;
    idx_P_++;
    
    use_fixed_covariance_ = true;
    isInitialized_ = true;
}


void FusionHandler::propagateState(const double& dt)
{
    ConstMatrix3 eye3 = Eigen::Matrix<double, 3, 3> ::Identity();
    ConstMatrix3 w_e_ie = Utility::skewSymmetric(Eigen::Matrix<double, 3, 1>(0, 0, Utility::W_IE));
    double alpha_ie = Utility::W_IE * dt;
    Eigen::Matrix<double, 3, 3> R_e_i;
    R_e_i << cos(alpha_ie), sin(alpha_ie), 0,
	    -sin(alpha_ie), cos(alpha_ie), 0,
                       0,             0,  1;   
		       
    State & cur_state = state_buffer_[idx_state_];
    State & prev_state = state_buffer_[idx_state_-1];    
    Eigen::Matrix<double, 3, 1> d_v;
    Eigen::Matrix<double, 3, 1> g = Utility::Earth::getGravityECEF(prev_state.p_);
    
//     cout << " propagateState: cur_state.q_: " <<  cur_state.q_.coeffs().transpose() << endl;
//     cout << " propagateState: prev_state.q_: " <<  prev_state.q_.coeffs().transpose() << endl;
    
    cur_state.b_a_ = prev_state.b_a_; 
    cur_state.b_w_ = prev_state.b_w_;

    ConstVector3 est_a0 = prev_state.f_b_ib - prev_state.b_a_;  
    ConstVector3 est_w0 = prev_state.w_b_ib - prev_state.b_w_;
    ConstVector3 est_a1 = cur_state.f_b_ib - cur_state.b_a_;
    ConstVector3 est_w1 = cur_state.w_b_ib - cur_state.b_w_;
    
//     // Principle of GNSS version
//     // (5.9) 姿态增量
//     Eigen::Matrix<double, 3, 1> alpha_b_ib;
//     Eigen::Matrix<double, 3, 3> Alpha_b_ib;
//     Eigen::Matrix<double, 3, 3> C_b_b;
//     Eigen::Matrix<double, 3, 1> f_e_ib;
//     alpha_b_ib = (est_w0 + est_w1) / 2 * dt;
//     Alpha_b_ib = Utility::skewSymmetric(alpha_b_ib);
// 
//     if (1) {
// 	C_b_b = Utility::Earth::getRodreigueMatrix(alpha_b_ib);
//     } else {
// 	C_b_b = eye3 + Alpha_b_ib;
//     }
//     
//     // (5.75)
//     // cur_state.q_ = prev_state.q_.toRotationMatrix() * C_b_b - w_e_ie * prev_state.q_.toRotationMatrix() * dt; // approximiate version
//     cur_state.q_ = R_e_i * prev_state.q_.toRotationMatrix() * C_b_b; // precise version
// 
//     Eigen::Matrix<double, 3, 3> R_e_b_bar;
//     if (1) {
// 	double mag_alpha = alpha_b_ib.norm();
// 	double mag_alpha2 = alpha_b_ib.squaredNorm();
// 	Eigen::Matrix<double, 3, 3> R_tmp = eye3 + (1 - cos(mag_alpha)) / mag_alpha2 * Alpha_b_ib + (1 - sin(mag_alpha) / mag_alpha) / mag_alpha2 * Alpha_b_ib * Alpha_b_ib;
// 	R_e_b_bar = prev_state.q_.toRotationMatrix() * R_tmp - 0.5 * w_e_ie * dt * prev_state.q_.toRotationMatrix();
//     } else {
// 	R_e_b_bar = prev_state.q_.toRotationMatrix() - 0.5 * w_e_ie * dt * prev_state.q_.toRotationMatrix();
//     }
// 
//     f_e_ib = R_e_b_bar * est_a1;
//     //TODO (5.36), 用(2.113), (2.150)变换在n系下算出的重力 
//     d_v = (f_e_ib + g) - 2 * w_e_ie * prev_state.v_;
//     cur_state.v_ = prev_state.v_ + d_v * dt;
//     cur_state.p_ = prev_state.p_ + (cur_state.v_ + prev_state.v_) / 2 * dt;
    
    
    // ETH version
    ConstMatrix4 Omega = Utility::omegaMatJPL(est_w1);
    ConstMatrix4 OmegaOld = Utility::omegaMatJPL(est_w0);
    Matrix4 OmegaMean = Utility::omegaMatJPL((est_w1 + est_w0) / 2);
    
    int div = 1;
    Matrix4 MatExp;
    MatExp.setIdentity();
    OmegaMean *= 0.5 * dt;
    for (int i = 1; i < 5; i++) {
	div *= i;
	MatExp = MatExp + OmegaMean / div;
	OmegaMean *= OmegaMean;
    }
    
    // first oder quat integration matrix
    ConstMatrix4 quat_int = MatExp + 1.0 / 48.0 * (Omega * OmegaOld - OmegaOld * Omega) * dt * dt;//(kinematics 225)
    cur_state.q_.coeffs() = quat_int * prev_state.q_.coeffs();
    cur_state.q_.normalize();
    
    d_v = (cur_state.q_.toRotationMatrix() * est_a1 + prev_state.q_.toRotationMatrix() * est_a0) / 2;
    
    //TODO G要去掉
    cur_state.v_ = prev_state.v_ + (d_v - g) * dt;
    cur_state.p_ = prev_state.p_ + ((cur_state.v_ + prev_state.v_) / 2 * dt);
    
    idx_state_++;
}

void FusionHandler::predictProcessCovariance(const double& dt)
{
    ConstMatrix3 w_e_ie = Utility::skewSymmetric(Eigen::Matrix<double, 3, 1>(0, 0, Utility::W_IE));
    ConstMatrix3 eye3 = Eigen::Matrix<double, 3, 3>::Identity();
    ConstMatrix3 est_R_e_eb = state_buffer_[idx_P_].q_.toRotationMatrix();
    ConstVector3 est_p_e_eb = state_buffer_[idx_P_].p_;
    const double est_L_b = 20;
    ConstVector3 est_f_b_ib = state_buffer_[idx_P_].f_b_ib - state_buffer_[idx_P_].b_a_;
    ConstVector3 v_e_ib = Utility::Earth::getGravityECEF(est_p_e_eb); 
    const double r_e_eS = Utility::Earth::getEarthRadius(est_L_b); // geocentric_radius

    ConstMatrix3 F_21 =  - 2 * v_e_ib / r_e_eS * state_buffer_[idx_P_].p_.transpose() / state_buffer_[idx_P_].p_.norm();
    ConstMatrix3 F_22 = -2 * w_e_ie;
    ConstMatrix3 F_23 = - est_R_e_eb * Utility::skewSymmetric(est_f_b_ib);
    ConstMatrix3 F_24 = - est_R_e_eb;
    ConstMatrix3 F_33 = -w_e_ie;
    ConstMatrix3 F_35 = - est_R_e_eb;
    
    Fd_.setIdentity();
    Fd_.block<3, 3> (0, 3) += eye3 * dt;
    Fd_.block<3, 3> (3, 0) += F_21 * dt;
    Fd_.block<3, 3> (3, 3) += F_22 * dt;
    Fd_.block<3, 3> (3, 6) += F_23 * dt;
    Fd_.block<3, 3> (3, 9) += F_24 * dt;
    Fd_.block<3, 3> (6, 6) += F_33 * dt;
    Fd_.block<3, 3> (6, 12) += F_35 * dt;    
    
    
    Gd_.setZero();
    Gd_.block<3, 3> (3, 3) += -est_R_e_eb;
    Gd_.block<3, 3> (6, 6) += -est_R_e_eb;
    Gd_.block<3, 3> (9, 9) += eye3;
    Gd_.block<3, 3> (12, 12) += eye3;
    
    // Principle of GNSS version (3.46W)
    state_buffer_[idx_P_].P_ = Fd_ * (state_buffer_[idx_P_-1].P_ + Qd_ / 2) * Fd_.transpose() + Qd_ / 2;

    // ETH version
    //  state_buffer_[idx_P_].P_ = Fd_ * state_buffer_[idx_P_-1].P_ * Fd_.transpose() + Gd_ * Qd_ * Gd_.transpose();
    
    idx_P_++;
}

template<class H_type, class Res_type, class R_type>
bool FusionHandler::applyMeasurement(
    const size_t& idx, 
    const Eigen::MatrixBase< H_type >& H, 
    const Eigen::MatrixBase< Res_type >& res, 
    const Eigen::MatrixBase< R_type >& R)
{
    
    // For idx_P_ < idx < idx_state_, propagate idx_P_ to idx
    propagatePtoIdx(idx);
    
    R_type S;
    Eigen::Matrix<double, global::NUM_OF_STATE, R_type::RowsAtCompileTime> K;
    Eigen::Matrix<double, global::NUM_OF_STATE, 1> delta;
    ErrorStateCov& P = state_buffer_[idx].P_;
    S = H * P * H.transpose() + R;
    //TODO Use QR and so on to speed up matrix inverse
    K = P * H.transpose() * S.inverse();
    delta = K * res;
    const ErrorStateCov I_KH = (ErrorStateCov::Identity() - K * H);
    P = I_KH * P * I_KH.transpose() + K * R * K.transpose();
    P = 0.5 * (P + P.transpose());
    
    return correctState(idx, delta);
    
}

void FusionHandler::propagatePtoIdx(const size_t& idx)
{
    // propagate cov matrix until idx. The smallest idx is 1.
    // idx_P_ 指的是下一个P,所以对应的state_buffer一定是空的
    if (idx < idx_state_ && idx_P_<= idx)	
	while (idx != idx_P_ - 1)
	    predictProcessCovariance(state_buffer_[idx_P_].time_-state_buffer_[idx_P_-1].time_);
}


bool FusionHandler::correctState(const size_t& idx, const FusionHandler::ErrorState& delta)
{
    State& state_measured = state_buffer_[idx];
    
    // Correct state
    state_measured.p_ += delta.block<3, 1>(0, 0);
    state_measured.v_ += delta.block<3, 1>(3, 0);
    
//     // Kinematics version
//     state_measured.q_ = Utility::deltaQ(delta.block<3, 1>(6, 0)) * state_measured.q_; // Global perturbation
//     state_measured.q_.normalize();
    
    // ETH version
    Eigen::Quaternion<double> qbuff_q = Utility::quaternionFromSmallAngle(delta.block<3, 1> (6, 0));
    state_measured.q_ = state_measured.q_ * qbuff_q;
    state_measured.q_.normalize();    
    
    state_measured.b_a_ += delta.block<3, 1>(9, 0);
    state_measured.b_w_ += delta.block<3, 1>(12, 0);
    
    size_t idx_latest_state = idx_state_;
    idx_state_ = idx + 1;
    idx_P_ = idx + 1;
    
    while (idx_state_ != idx_latest_state) {
	propagateState(state_buffer_[idx_state_].time_ - state_buffer_[(unsigned char)(idx_state_ - 1)].time_);
	predictProcessCovariance(state_buffer_[idx_P_].time_ - state_buffer_[(unsigned char)(idx_P_ - 1)].time_);
    
    }
    
    if(!Utility::checkForNumeric(&delta[0], global::NUM_OF_STATE, "update")) {
	
    }
    
    
    return 1;
}

// 要看 1. 这个state是否有东西 2. 尤其在收尾两个极端时,注意是不是Buffer里根本没有相应state 3. 搜索时,如果测量的时间太后面,远远超过前面了,该返回失败

size_t FusionHandler::getClosestState(State& state, double tstamp, double delay)
{
    if (idx_state_ == 1) {
	ROS_WARN_STREAM("getClosestState: No enough IMU measurement");
	return 0;
    }  
    
    size_t idx = idx_state_ - 1;
    double closest_time = 1.0;
    for (; idx > 0; idx--) {
	double time_dist = fabs(tstamp - state_buffer_[idx].time_);
	if (time_dist > 1.0) {
	    ROS_WARN_STREAM("getClosestState: Measurement delays too much");
	    return 0;
	}
	if (time_dist < closest_time) 
	    closest_time = time_dist;
	else {	    
	    break;
	}
	    
    }
    idx++;
    
//     double timedist = 1e100;
//     double time_measured = tstamp.toSec() - delay;    
//     
//     while (fabs(time_measured - state_buffer_[idx].time_) < timedist) {
// 	timedist = fabs(time_measured - state_buffer_[idx].time_);
// 	if (idx == 0) {
// 	    ROS_WARN_STREAM("timedist is too small! ");
// 	    break;
// 	}  
// 	idx--;
//     }  
//     idx++;
    
    if (state_buffer_[idx].time_ == -1) {
	ROS_INFO_STREAM("getClosestState: state needed to be initialized");
	return 0;
    }
	
    
    if (idx == idx_state_ -1) {
	ROS_INFO_STREAM("getClosestState: Correspond to latest state");
    }
    
//     propagatePtoIdx(idx);
    
    state = state_buffer_[idx];
    
    return idx;
}


void FusionHandler::initByGNSS(
    ros::Time time, 
    const Matrix< double, int(3), int(1) >& p, 
    const Matrix< double, int(3), int(1) >& v, 
    const double& pos_var, 
    const double& vel_var)
{
//     isInitialized_ = false;
//     for (size_t i = 0; i < global::NUM_STATE_BUFFER; i++) {
// 	state_buffer_[i].reset();
//     }    
//     
//     idx_state_ = 0;
//     idx_P_ = 0;  
//     
//     // 初始化状态
//     State& state = state_buffer_[idx_state_];    
//     state.p_ = p;
//     state.v_ = v;
//     //TODO 方位角用速度向量初始化
//     state.q_; 
//     state.b_a_.setZero();
//     state.b_w_.setZero();
// 
//     // 初始化矩阵P_0, 代表的是初始状态数值的不确定性.跟你怎么赋值有关系
//     state.P_.block<3, 3>(0, 0) = (Eigen::Matrix<double, 3, 1>() << pos_var, pos_var, pos_var*2).finished().asDiagonal();
//     state.P_.block<3, 3>(3, 3) = (Eigen::Matrix<double, 3, 1>() << vel_var, vel_var, vel_var*2).finished().asDiagonal();
//     //TODO 方位角怎么初始化?一会儿要写初始化方位角的程序,用弧度表示roll和pitch的不确定性,公式在(14.93)
//     state.P_.block<3, 3>(9, 9) = Eigen::Matrix<double, 3, 3>::Identity() * global::ACC_W * global::ACC_W;
//     state.P_.block<3, 3>(12, 12) = Eigen::Matrix<double, 3, 3>::Identity() * global::GYR_W * global::GYR_W;
//     
//     Qd_.setZero();
//     Qd_.block<3, 3>(3, 3) = (global::ACC_N * global::ACC_N) * Eigen::Matrix<double, 3, 3>::Identity();
//     Qd_.block<3, 3>(6, 6) = (global::GYR_N * global::GYR_N) * Eigen::Matrix<double, 3, 3>::Identity();
//     Qd_.block<3, 3>(9, 9) = (global::ACC_W * global::ACC_W) * Eigen::Matrix<double, 3, 3>::Identity();
//     Qd_.block<3, 3>(12, 12) = (global::GYR_W * global::GYR_W) * Eigen::Matrix<double, 3, 3>::Identity();  
// 
//     l_b_ba = global::L_b_ba;
//     
//     idx_state_++;
//     idx_P_++;
//     
//     isInitialized_ = true;  
}  
    
}

