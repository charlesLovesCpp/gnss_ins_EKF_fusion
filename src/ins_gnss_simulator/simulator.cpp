#include "simulator.h"

namespace Simulator {

using namespace std;
using namespace Eigen;

GnssSimulator::GnssSimulator(ros::NodeHandle& nh, ros::NodeHandle& pnh):
    it(nh),
    hasImu(false),
    hasGnss(false)
{
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator_ = std::default_random_engine(seed); 
    
    pub_imu_ = nh.advertise<sensor_msgs::Imu> ("/imu0", 1);  
    pub_gnss_ = nh.advertise<integrated_nav::GNSS_ECEF> ("/gnss/ecef", 1); 
}	



void GnssSimulator::pubGnss(	
	const GnssResult& ecef, 
	const double& position_acc, 
	const double& speed_acc, 
	const int& num_SV,
	const int& p_DOP,
	const int& gps_fix_status)
{
    cout << "pub GNSS!" << endl;
    integrated_nav::GNSS_ECEF gnss_msg;
    gnss_msg.header.stamp = ros::Time::now();
    gnss_msg.header.frame_id = "gnss_link";
    gnss_msg.time = ecef.time;
    gnss_msg.ecef_x = ecef.p.x();
    gnss_msg.ecef_y = ecef.p.y();
    gnss_msg.ecef_z = ecef.p.z();
    gnss_msg.ecefV_x = ecef.v.x();
    gnss_msg.ecefV_y = ecef.v.y();
    gnss_msg.ecefV_z = ecef.v.z();
    gnss_msg.position_acc = position_acc;
    gnss_msg.speed_acc = speed_acc;
    gnss_msg.p_DOP = p_DOP;
    gnss_msg.num_SV = num_SV;
    gnss_msg.gps_fix_status = gps_fix_status;
    pub_gnss_.publish(gnss_msg);
}
    
void GnssSimulator::pubImu(const double& time, const Eigen::Vector3d& f_ib_b, const Eigen::Vector3d& omega_ib_b)
{
    sensor_msgs::Imu imu_msg;
    imu_msg.header.stamp = ros::Time::now();
    imu_msg.header.frame_id = "base_link";    
    imu_msg.linear_acceleration.x = f_ib_b.x();
    imu_msg.linear_acceleration.y = f_ib_b.y();
    imu_msg.linear_acceleration.z = f_ib_b.z();
    imu_msg.angular_velocity.x = omega_ib_b.x();
    imu_msg.angular_velocity.y = omega_ib_b.y();
    imu_msg.angular_velocity.z = omega_ib_b.z();  
    
    imu_msg.angular_velocity_covariance[0] = time;
    
    pub_imu_.publish(imu_msg);
}

void GnssSimulator::reset()
{
    ned_.clear();
    ecef_.clear();    
    est_ecef_.clear();
    est_clock_.clear();
}

void GnssSimulator::run()
{
    Eigen::Vector3d idtt;
    idtt.setZero();
	    
    reset();

    readFile("/home/charles/nav_ws/src/integrated_nav/data/Ground_true_1.csv");
    
    imu_parameter_.setDefault();
    gnss_config_.setDefault();    

    int num_gnss_epoch = 0;
    
    // First GNSS measurement
    GnssResult& res_ecef0 = ecef_.at(0);
    GnssResult& res_ned0 = ned_.at(0);

    double last_gnss_time = res_ecef0.time;
    
    // Generate constellation information
    std::vector<Eigen::Vector3d> sats_p0;
    std::vector<Eigen::Vector3d> sats_v0;
    solveSvPV(last_gnss_time, gnss_config_, sats_p0, sats_v0);     
    
    // Generate pseudo range (rate) bias
    std::vector<double> gnss_biases;
    initialGnssBias(	
	gnss_config_,   
	sats_p0, 
	res_ecef0.p,
	res_ned0.p.x(),
	res_ned0.p.y(), 
	gnss_biases);

    // Set original IMU quantity residuals
    Eigen::Matrix<double, 6, 1> quant_residuals, old_quant_residuals;
    quant_residuals.setZero(); 
    old_quant_residuals.setZero();
    
    std::vector<GnssMeasurement> gnss_measurements0;
    int num_gnss_meas0 = 0;
    generateGnssMeasurements(    
	last_gnss_time, 
	gnss_config_,   
	sats_p0, 
	sats_v0,
	res_ecef0.p,
	res_ecef0.v,
	res_ned0.p.x(),
	res_ned0.p.y(),
	gnss_biases,
	gnss_measurements0,
	num_gnss_meas0);    
    
    GnssResult est_ecef0;
    Eigen::Vector2d est_clock0;
    est_ecef0.time = res_ecef0.time;
    
    bool result = false;
    result = solveReceiverPV(
	gnss_measurements0,
	num_gnss_meas0,
	idtt,
	idtt,
	est_ecef0.p,
	est_ecef0.v,
	est_clock0);

    est_ecef_.push_back(est_ecef0);
    est_clock_.push_back(est_clock0);
   
    double time_delay = 0.1;
    // Publish first position&velocity
    pubGnss(est_ecef0, 2.5, 0.1, num_gnss_meas0);
    ros::Duration(0.5).sleep();
    pubGnss(est_ecef0, 2.5, 0.1, num_gnss_meas0);
    ros::Duration(0.5).sleep();
    num_gnss_epoch++;
/*
    for (int i = 0; i < sats_p0.size(); ++i) {
	cout << "sats_p " << i << " : "<< sats_p0[i].transpose() << endl;
    }
    
    
    for (int i = 0; i < gnss_measurements0.size(); ++i) {
	cout << "range " << i << " : "<< gnss_measurements0[i].range << endl;
    }
    */
    
	    
    cout <<"ECEF true P: "<< res_ecef0.p.transpose() << "     V: " << res_ecef0.v.transpose() << endl;  
    cout <<"ECEF estimate P: "<< est_ecef0.p.transpose() << "     V: " << est_ecef0.v.transpose() << "\n" << endl;  
    cout <<"estimate clock offset: "<< est_clock0.x() << "     estimate clock drift: " << est_clock0.y() << endl; 
    
    
    for (size_t i = 1; i < 512; ++i) {
	GnssResult& cur_ecef = ecef_.at(i);
	GnssResult& cur_ned = ned_.at(i);
	GnssResult& prev_ecef = ecef_.at(i - 1);
	
	double tor_i = cur_ecef.time - prev_ecef.time;
	//TODO Generate IMU measurement
	Eigen::Vector3d f_ib_b;
	Eigen::Vector3d omega_ib_b;	
	deriveImuMeasurement(
	    tor_i,
	    cur_ecef.R, prev_ecef.R,
	    cur_ecef.v, prev_ecef.v,
	    cur_ecef.p,
	    f_ib_b, omega_ib_b);

// 	cout << f_ib_b.transpose() << ", " << omega_ib_b.transpose() << endl;
	
	Eigen::Vector3d meas_f_ib_b;
	Eigen::Vector3d meas_omega_ib_b;
	perturbeImuMeasurement(
	    tor_i, 
	    imu_parameter_,
	    f_ib_b,
	    omega_ib_b,
	    old_quant_residuals,
	    meas_f_ib_b,
	    meas_omega_ib_b,
	    quant_residuals);
	
// 	cout << meas_f_ib_b.transpose() << ", " << meas_omega_ib_b.transpose() << endl;
	
	old_quant_residuals = quant_residuals;
	
	// Publish IMU measurements
// 	pubImu(cur_ecef.time, meas_f_ib_b, meas_omega_ib_b);
// 	ros::Duration(time_delay).sleep();
	
	// When GNSS output is accessible
	if ((cur_ecef.time - last_gnss_time) > gnss_config_.epoch_interval) {
	    //TODO Generate GNSS measurement
	    last_gnss_time = cur_ecef.time;
	    
	    std::vector<Eigen::Vector3d> sats_p;
	    std::vector<Eigen::Vector3d> sats_v;
	    solveSvPV(last_gnss_time, gnss_config_, sats_p, sats_v); 
	    
// 	    for (int i = 0; i < sats_p.size(); ++i) {
// 		cout << "sats_p " << i << " : "<< sats_p[i].transpose() << endl;
// 	    }
	    
	    std::vector<GnssMeasurement> gnss_measurements;
	    int num_gnss_meas = 0;
	    generateGnssMeasurements(    
		last_gnss_time, 
		gnss_config_,   
		sats_p, 
		sats_v,
		cur_ecef.p,
		cur_ecef.v,
		cur_ned.p.x(),
		cur_ned.p.y(),
		gnss_biases,
		gnss_measurements,
		num_gnss_meas);
	    
// 	    for (int i = 0; i < gnss_measurements.size(); ++i) {
// 		cout << "range " << i << " : "<< gnss_measurements[i].range << endl;
// 	    }
	    
	    GnssResult est_ecef;
	    Eigen::Vector2d est_clock;
	    est_ecef.time = cur_ecef.time;
	    
	    bool result = false;
	    result = solveReceiverPV(
		gnss_measurements,
		num_gnss_meas,
		prev_ecef.p,
		prev_ecef.v,
		est_ecef.p,
		est_ecef.v,
		est_clock);

	    est_ecef_.push_back(est_ecef);
	    est_clock_.push_back(est_clock);
cout <<"ECEF true P:     "<< ecef_[i].p.transpose() << endl;
cout <<"ECEF estimate P: "<< est_ecef.p.transpose() << endl;
cout <<"ECEF true V:     "<< ecef_[i].v.transpose() << endl;
cout <<"ECEF estimate V: " << est_ecef.v.transpose() << endl;  
cout <<"estimate clock offset: "<< est_clock.x() << "     estimate clock drift: " << est_clock.y() << endl; 
	    
	    Eigen::Vector3d p_err = est_ecef.p - ecef_[i].p;
	    Eigen::Vector3d v_err = est_ecef.v - ecef_[i].v;
	    cout <<"position error: "<< p_err.transpose() << ", norm: " << p_err.norm() << endl;
	    cout <<"velocity error: "<< v_err.transpose() << ", norm: " << v_err.norm() << "\n" << endl;

// 	    pubGnss(est_ecef, 2.5, 0.1, num_gnss_meas);    
// 	    ros::Duration(time_delay).sleep();
	    num_gnss_epoch++;
	}
    }
    
    cout << num_gnss_epoch << " GNSS msgs published successfully!" << endl;
}


bool GnssSimulator::solveSvPV(
    const double& time, 
    const GnssConfig& config,
    std::vector<Eigen::Vector3d>& sats_p, 
    std::vector<Eigen::Vector3d>& sats_v)
{
    double mu = Utility::MU;
    double omega_ie = Utility::W_IE;
    double pi = M_PI;
    sats_p.clear();
    sats_v.clear();
    
    // Convert inclination angle to degrees
    double inclination = Utility::DEG2RAD * config.inclination; 
    
    // Determine orbital angular rate using (8.8)
    double omega_is = sqrt(mu / pow(config.r_os, 3));    
    
    // Determine constellation time
    double const_time = time + config.const_delta_t;    
    
    for (int i = 0; i < config.num_sat; ++i) {
	
	Eigen::Vector3d sat_r_es_e;
	Eigen::Vector3d sat_v_es_e;
	Eigen::Vector3d r_os_o, v_os_o;
	
	double u_os_o = 2 * pi * (i-1)/config.num_sat + omega_is*const_time;
	double Omega = (pi * (i % 6) / 3 + Utility::DEG2RAD * config.const_delta_lambda) - omega_ie * const_time; // Longitude of the ascending node from (8.16)
	
	// Satellite position in the orbital frame from (8.14)
	r_os_o << config.r_os*cos(u_os_o), config.r_os*sin(u_os_o), 0;
	
	// ECEF Satellite Position from (8.19)
	sat_r_es_e << r_os_o(0)*cos(Omega) - r_os_o(1)*cos(inclination)*sin(Omega), 
		    r_os_o(0)*sin(Omega) + r_os_o(1)*cos(inclination)*cos(Omega),
		    r_os_o(1)*sin(inclination);

	// Satellite velocity in the orbital frame from (8.25), noting that with
	// a circular orbit r_os_o is constant and the time derivative of u_os_o
	// is omega_is.
	v_os_o = config.r_os * omega_is * Eigen::Vector3d(-sin(u_os_o), cos(u_os_o), 0);   
	
	// ECEF Satellite velocity from (8.26)
	sat_v_es_e << v_os_o(0)*cos(Omega) - v_os_o(1)*cos(inclination)*sin(Omega) + omega_ie*sat_r_es_e(1),
		    v_os_o(0)*sin(Omega) + v_os_o(1)*cos(inclination)*cos(Omega) - omega_ie*sat_r_es_e(0),
		    v_os_o(1)*sin(inclination);
	
	sats_p.push_back(sat_r_es_e);
	sats_v.push_back(sat_v_es_e);
    }      
    return 1;
}

void GnssSimulator::generateGnssMeasurements(
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
    int& num_gnss_meas)
{
    std::normal_distribution<double> distribution(0, 1); 
	
    double c = Utility::C;
    double omega_ie = Utility::W_IE;
    Eigen::Matrix<double, 3, 3> Omega_ie = Utility::skewSymmetric(Eigen::Matrix<double, 3, 1>(0, 0, Utility::W_IE));
	
    gnss_measurements.clear();    
    
    double cos_lat = cos(L_a);
    double sin_lat = sin(L_a);
    double cos_long = cos(lambda_a);
    double sin_long = sin(lambda_a);
    
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
	     -sin_long,            cos_long,        0,
	     -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    num_gnss_meas = 0;
    for (int i = 0; i < config.num_sat; ++i) {
	// Determine ECEF line-of-sight vector using (8.41)
	Eigen::Vector3d delta_r = sats_p.at(i) - r_ea_e;
	double approx_range = delta_r.norm();
	Eigen::Vector3d u_as_e = delta_r / approx_range;
	
	// Convert line-of-sight vector to NED using (8.39) and determine
	// elevation using (8.57)
	double elevation = -asin(C_e_n.block<1, 3>(2, 0) * u_as_e);

	if (elevation >= Utility::DEG2RAD * config.mask_angle) {   
	 
	    Eigen::Matrix3d C_e_I;
	    
	    // Calculate frame rotation during signal transit time using (8.36)
	    C_e_I << 1, omega_ie * approx_range / c, 0,
		    -omega_ie * approx_range / c, 1, 0,
		     0, 0, 1;
	    // Calculate range using (8.35)
	    delta_r = C_e_I * sats_p.at(i) - r_ea_e;
	    double range = delta_r.norm();
	    
	    // Calculate range rate using (8.44)
	    double range_rate = u_as_e.transpose() * ( C_e_I * (sats_v.at(i) + Omega_ie * sats_p.at(i)) - 
				(v_ea_e + Omega_ie * r_ea_e) );
	    
	    // Calculate pseudo-range measurement
	    double range_measurement = range + gnss_biases.at(i) + config.rx_clock_offset + config.rx_clock_drift * time + config.code_track_err_SD * distribution(generator_);
	    
// 	    double range_measurement = range +  config.rx_clock_offset + config.rx_clock_drift * time +  config.code_track_err_SD * distribution(generator_);
	    
	    // Calculate pseudo-range rate measurement
	    double range_rate_measurement = range_rate + config.rx_clock_drift + config.rate_track_err_SD * distribution(generator_);
	    
	    GnssMeasurement gnss_meas;
	    gnss_meas.range = range_measurement;
	    gnss_meas.range_rate = range_rate_measurement;
	    gnss_meas.sat_p = sats_p.at(i);
	    gnss_meas.sat_v = sats_v.at(i);
	    
	    gnss_measurements.push_back(gnss_meas);
	    num_gnss_meas++;	 
	}
    }
	     
    
}

void GnssSimulator::initialGnssBias(
    const GnssConfig& config, 
    const vector< Vector3d >& sats_p, 
    const Vector3d& r_ea_e, 
    const double& L_a, 
    const double& lambda_a, 
    std::vector< double >& gnss_biases)
{
    std::normal_distribution<double> distribution(0, 1); 
	
    double cos_lat = cos(L_a);
    double sin_lat = sin(L_a);
    double cos_long = cos(lambda_a);
    double sin_long = sin(lambda_a);
    
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
	     -sin_long,            cos_long,        0,
	     -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
	     
    for (int i = 0; i < config.num_sat; ++i) {
	// Determine ECEF line-of-sight vector using (8.41)
	Eigen::Vector3d delta_r = sats_p.at(i) - r_ea_e;
	double approx_range = delta_r.norm();
	Eigen::Vector3d u_as_e = delta_r / approx_range;

	// Convert line-of-sight vector to NED using (8.39) and determine
	// elevation using (8.57)
	double elevation = -asin(C_e_n.block<1, 3>(2, 0) * u_as_e);
	
	// Limit the minimum elevation angle to the masking angle
	elevation = max(elevation, Utility::DEG2RAD*config.mask_angle);

	// Calculate ionosphere and troposphere error SDs using (9.79) and (9.80)	
	double iono_SD = config.zenith_iono_err_SD / sqrt(1 - 0.899 * cos(elevation) * cos(elevation));
	double trop_SD = config.zenith_trop_err_SD / sqrt(1 - 0.998 * cos(elevation) * cos(elevation));
	// Determine range bias
	gnss_biases.push_back(config.SIS_err_SD*distribution(generator_) 
			    + iono_SD * distribution(generator_) 
			    + trop_SD * distribution(generator_));
    }
}


bool GnssSimulator::solveReceiverPV(
    const vector<GnssMeasurement>& gnss_measurements, 
    const int& num_gnss_meas, 
    const Vector3d& predicted_r_ea_e, 
    const Vector3d& predicted_v_ea_e, 
    Vector3d& est_r_ea_e, 
    Vector3d& est_v_ea_e, 
    Vector2d& est_clock)
{
    double c = Utility::C;
    double omega_ie = Utility::W_IE;
    Eigen::Matrix<double, 3, 3> Omega_ie = Utility::skewSymmetric(Eigen::Matrix<double, 3, 1>(0, 0, Utility::W_IE));
    
    int max_iter = 1000;
    Eigen::Matrix<double, 4, 1> x_pred, x_est;
    double test_convergence;
    
    if (gnss_measurements.size() < 4) {
	ROS_WARN_STREAM("solveReceiverPV: No enough visible satellite");
	return 0;
    }
    
    Eigen::VectorXd res_final(num_gnss_meas, 1);

    ///*********Solve P*********//
    
    x_pred.head<3>() = predicted_r_ea_e;
    x_pred.tail<1>()(0) = 0;

    test_convergence = 1;
    while (test_convergence > 0.0001 && max_iter != 0) {

	Eigen::VectorXd res(num_gnss_meas, 1);
	Eigen::MatrixXd H(num_gnss_meas, 4);
	res.setZero();
	H.setZero();
	
	for (int i = 0; i < num_gnss_meas; ++i) {
	    // Predict approx range 
	    
	    const GnssMeasurement& gnss_meas = gnss_measurements.at(i);
	    Eigen::Vector3d  delta_r = gnss_meas.sat_p - x_pred.head<3>() ;
	    double approx_range = delta_r.norm();

	    // Calculate frame rotation during signal transit time using (8.36)
	    Eigen::Matrix3d C_e_I;
	    C_e_I << 1, omega_ie * approx_range / c, 0,
		    -omega_ie * approx_range / c, 1, 0,
		     0, 0, 1;
	    // Predict pseudo-range using (9.143)
	    delta_r = C_e_I * gnss_meas.sat_p - x_pred.head<3>() ;
	    double range = delta_r.norm();
	    Eigen::Vector3d u_aS_e = delta_r / range;
	    
	    // Resisual
	    res(i) = gnss_meas.range - (range + x_pred(3));
	    
	    // Predict line of sight and deploy in measurement matrix, (9.144)
	    H(i, 0) = -u_aS_e.x();
	    H(i, 1) = -u_aS_e.y();
	    H(i, 2) = -u_aS_e.z();
	    H(i, 3) = 1;
	}
	
	// Unweighted least-squares solution, (9.35)/(9.141)
	Eigen::Matrix<double, 4, 4> S = H.transpose() * H;
	x_est = x_pred +  S.inverse() * H.transpose() * res;
	
	// Test convergence   
	test_convergence = (x_est - x_pred).norm();
	
	// Set predictions to estimates for next iteration
	x_pred = x_est;	
	
	max_iter--;
	
	res_final = res;
    } 
    
    if (max_iter == 0) {
	ROS_WARN_STREAM("solveReceiverPV: P Does not converge");
	return 0;
    }
    
    est_r_ea_e = x_est.head<3>();
    est_clock.x() = x_est.tail<1>()(0);
    
/*    cout << "p test_convergence: " << test_convergence << endl;
    cout << "p res_final: " << res_final.transpose() << endl; */   
    
    ///*********Solve V*********//
    
    x_pred.head<3>() = predicted_v_ea_e;
    x_pred.tail<1>()(0) = 0;
    
    max_iter = 1000;
    test_convergence = 1;
    while (test_convergence > 0.0001 && max_iter != 0) {

	Eigen::VectorXd res(num_gnss_meas, 1);
	Eigen::MatrixXd H(num_gnss_meas, 4);
	res.setZero();
	H.setZero();

	for (int i = 0; i < num_gnss_meas; ++i) {
	    const GnssMeasurement& gnss_meas = gnss_measurements.at(i);
	    
	    // Predict approx range 
	    Eigen::Vector3d  delta_r = gnss_meas.sat_p - est_r_ea_e;
	    double approx_range = delta_r.norm();

	    // Calculate frame rotation during signal transit time using (8.36)
	    Eigen::Matrix3d C_e_I;
	    C_e_I << 1, omega_ie * approx_range / c, 0,
		    -omega_ie * approx_range / c, 1, 0,
		     0, 0, 1;
	    // Predict pseudo-range using (9.143)
	    delta_r = C_e_I * gnss_meas.sat_p - est_r_ea_e;
	    double range = delta_r.norm();
	    Eigen::Vector3d u_aS_e = delta_r / range;
	    
	    double range_rate = u_aS_e.transpose() * (C_e_I * (gnss_meas.sat_v + Omega_ie * gnss_meas.sat_p)
				- ( x_pred.head<3>() + Omega_ie * est_r_ea_e));
	    
	    // Resisual
	    res(i) = gnss_meas.range_rate - (range_rate + x_pred(3));
	    
	    // Predict line of sight and deploy in measurement matrix, (9.144)
	    H(i, 0) = -u_aS_e.x();
	    H(i, 1) = -u_aS_e.y();
	    H(i, 2) = -u_aS_e.z();
	    H(i, 3) = 1;
	}
	
	// Unweighted least-squares solution, (9.35)/(9.141)
	Eigen::Matrix<double, 4, 4> S = H.transpose() * H;
	x_est = x_pred +  S.inverse() * H.transpose() * res;
	
	// Test convergence   
	test_convergence = (x_est - x_pred).norm();
	
	// Set predictions to estimates for next iteration
	x_pred = x_est;	
	
	max_iter--;
	
	res_final = res;
    } 
//     cout << "v test_convergence: " << test_convergence << endl;
//     cout << "v res_final: " << res_final.transpose() << endl;
    
    
    if (max_iter == 0) {
	ROS_WARN_STREAM("solveReceiverPV: V Does not converge");
	return 0;
    }    
    
    
    
    est_v_ea_e = x_est.head<3>();
    est_clock.y() = x_est.tail<1>()(0);
    
    return 1;
}

void GnssSimulator::deriveImuMeasurement(
    const double& tor_i, 
    const Matrix< double, int(3), int(3) >& C_b_e, 
    const Matrix< double, int(3), int(3) >& old_C_b_e, 
    const Vector3d& v_eb_e, 
    const Vector3d& old_v_eb_e, 
    const Vector3d& r_eb_e, 
    Vector3d& f_ib_b, 
    Vector3d& omega_ib_b)
{
    const Eigen::Matrix<double, 3, 3> eye3 = Eigen::Matrix<double, 3, 3>::Identity();
    double omega_ie = Utility::W_IE;
    Eigen::Matrix<double, 3, 3> Omega_ie = Utility::skewSymmetric(Eigen::Matrix<double, 3, 1>(0, 0, Utility::W_IE));

    if (tor_i > 0) {
	// From (2.145) determine the Earth rotation over the update interval
	// C_Earth = C_e_i' * old_C_e_i
	double alpha_ie = omega_ie * tor_i;
	
	Eigen::Matrix<double, 3, 3> C_Earth;
	C_Earth << cos(alpha_ie), sin(alpha_ie),   0,
		    -sin(alpha_ie), cos(alpha_ie), 0,
                           0,             0,       1;
                           
	// Obtain coordinate transformation matrix from the old attitude (w.r.t.
	// an inertial frame) to the new
	Eigen::Matrix<double, 3, 3> C_old_new = C_b_e.transpose() * C_Earth * old_C_b_e;                           

	// Calculate the approximate angular rate w.r.t. an inertial frame
	Eigen::Vector3d alpha_ib_b;
	alpha_ib_b(0) = 0.5 * (C_old_new(1,2) - C_old_new(2,1));
	alpha_ib_b(1) = 0.5 * (C_old_new(2,0) - C_old_new(0,2));
	alpha_ib_b(2) = 0.5 * (C_old_new(0,1) - C_old_new(1,0));
    
	// Calculate and apply the scaling factor
	double temp = acos(0.5 * (C_old_new(0,0) + C_old_new(1,1) + C_old_new(2,2) - 1.0));	
	
	if (temp>2e-5)  // scaling is 1 if temp is less than this
	    alpha_ib_b = alpha_ib_b * temp/sin(temp);
	
	// Calculate the angular rate
	omega_ib_b = alpha_ib_b / tor_i;

	// Calculate the specific force resolved about ECEF-frame axes
	// From (5.36),
	
	Vector3d g = Utility::Earth::getGravityECEF(r_eb_e);
// 	cout << "gravity: " << g.transpose() << endl;
	Vector3d f_ib_e = ((v_eb_e - old_v_eb_e) / tor_i) - g + 2 * Omega_ie * old_v_eb_e;

	// Calculate the average body-to-ECEF-frame coordinate transformation
	// matrix over the update interval using (5.84) and (5.85)
	double mag_alpha = alpha_ib_b.norm();
	double mag_alpha2 = alpha_ib_b.squaredNorm();
	Eigen::Matrix3d Alpha_ib_b = Utility::skewSymmetric(alpha_ib_b); 

	Eigen::Matrix3d ave_C_b_e;
	if (mag_alpha > 1.e-8) {
	    ave_C_b_e = old_C_b_e * (eye3 + (1 - cos(mag_alpha))/mag_alpha2*Alpha_ib_b 
					  + (1 - sin(mag_alpha) / mag_alpha)/mag_alpha2*Alpha_ib_b*Alpha_ib_b)
					  - 0.5*Utility::skewSymmetric(Eigen::Vector3d(0, 0, alpha_ie)) * old_C_b_e;
	} else {
	    ave_C_b_e = old_C_b_e - 0.5 * Utility::skewSymmetric(Eigen::Vector3d(0, 0, alpha_ie)) * old_C_b_e;
	}
	f_ib_b = ave_C_b_e.inverse() * f_ib_e;
	
    } else {
	omega_ib_b.setZero();
	f_ib_b.setZero();
    }
}

void GnssSimulator::perturbeImuMeasurement(
    const double& tor_i, 
    const ImuParameter& imu_parameter_, 
    const Eigen::Vector3d& true_f_ib_b, 
    const Eigen::Vector3d& true_omega_ib_b, 
    const Eigen::Matrix< double, int(6), int(1) >& old_quant_residuals, 
    Eigen::Vector3d& meas_f_ib_b, 
    Eigen::Vector3d& meas_omega_ib_b,
    Eigen::Matrix<double, 6, 1>& quant_residuals)
{
    std::normal_distribution<double> distribution(0, 1); 
    const Eigen::Matrix<double, 3, 3> eye3 = Eigen::Matrix<double, 3, 3>::Identity();
    
    Eigen::Vector3d randn_vec1, randn_vec2;
    randn_vec1 << distribution(generator_), distribution(generator_), distribution(generator_);
    randn_vec2 << distribution(generator_), distribution(generator_), distribution(generator_);
    
    Eigen::Vector3d accel_noise, gyro_noise;
    if (tor_i > 0) {
	accel_noise = randn_vec1 * imu_parameter_.accel_noise_root_PSD / sqrt(tor_i);  	
	gyro_noise = randn_vec2 * imu_parameter_.gyro_noise_root_PSD / sqrt(tor_i);  
    } else {
	accel_noise.setZero();
	gyro_noise.setZero();
    }
    
    // Calculate accelerometer and gyro outputs using (4.16) and (4.17)
    Eigen::Vector3d uq_f_ib_b = imu_parameter_.b_a + (eye3 + imu_parameter_.M_a) * true_f_ib_b + accel_noise;
    Eigen::Vector3d uq_omega_ib_b = imu_parameter_.b_g + (eye3 + imu_parameter_.M_g) * true_omega_ib_b + imu_parameter_.G_g * true_f_ib_b + gyro_noise;

    if (imu_parameter_.accel_quant_level > 0) {
	Eigen::Vector3d tmp = (uq_f_ib_b + old_quant_residuals.block<3, 1>(0, 0)) / imu_parameter_.accel_quant_level;
	meas_f_ib_b = imu_parameter_.accel_quant_level * Utility::roundVector3(tmp);
	quant_residuals.block<3, 1>(0, 0) = uq_f_ib_b + old_quant_residuals.block<3, 1>(0, 0) - meas_f_ib_b;
    } else {
	meas_f_ib_b = uq_f_ib_b;
	quant_residuals.block<3, 1>(0, 0) << 0, 0, 0;
    }
    
    if (imu_parameter_.gyro_quant_level > 0) {
	Eigen::Vector3d tmp = (uq_omega_ib_b + old_quant_residuals.block<3, 1>(3, 0)) / imu_parameter_.gyro_quant_level;
	meas_omega_ib_b = imu_parameter_.gyro_quant_level * Utility::roundVector3(tmp);
	quant_residuals.block<3, 1>(3, 0) = uq_omega_ib_b + old_quant_residuals.block<3, 1>(3, 0) - meas_omega_ib_b;
    } else {
	meas_omega_ib_b = uq_omega_ib_b;
	quant_residuals.block<3, 1>(3, 0) << 0, 0, 0;
    }
     
}


void splitLine(const string& str,const char& separator, std::vector<double>& result) {
    istringstream sin(str);
    string str_data;
    double data;
    while (getline(sin, str_data, separator)) {
	stringstream ss;
	ss << str_data;
	ss >> data; 
	result.push_back(data);
    }
}

void GnssSimulator::neds2ecefs(const std::vector<GnssResult>& neds, std::vector<GnssResult>& ecefs)
{
    if (neds.size() == 0)
	return;
    ecefs.clear();
    
    for (size_t i = 0; i < neds.size(); ++i) {
	const GnssResult& res_ned = neds.at(i);
	GnssResult res_ecef;
	ned2ecef(res_ned, res_ecef);
	ecefs.push_back(res_ecef);
    }
}


void GnssSimulator::ned2ecef(const GnssResult& ned, GnssResult& ecef)
{
    ecef.time = ned.time;
    Utility::Earth::NED2ECEF(ned.R, ecef.R, ned.p, ecef.p, ned.v, ecef.v);
}


bool GnssSimulator::readFile(const string& file_name)
{
    string value;
    ifstream fin;
    fin.open(file_name, ios::in);
    if (!fin.is_open())
    {
        cout << "Read file error: " << file_name << endl;
    }
    
    ned_.clear();
    string line;
    while (getline(fin, line)) {
	vector<double> data;
	splitLine(line, ',', data);
	GnssResult res_ned;
	res_ned.time = data[0];
	res_ned.p << Utility::DEG2RAD*data[1], Utility::DEG2RAD*data[2], data[3];
	res_ned.v << data[4], data[5], data[6];
	res_ned.R = Utility::ypr2R(Eigen::Vector3d(data[7], data[8], data[9]));
	ned_.push_back(res_ned);
	
	GnssResult res_ecef;
	ned2ecef(res_ned, res_ecef);
	ecef_.push_back(res_ecef);
	
// 	cout <<"NED result: P: "<< res_ned.p.transpose() << "     V: " << res_ned.v.transpose() << endl; 
// 	cout <<"ECEF result: P: "<< res_ecef.p.transpose() << "    V: " << res_ecef.v.transpose() << endl; 
    } 
    
    
    
    return 1;
}
    /* 
    pubImu(Eigen::Vector3d(11, 11, 11), Eigen::Vector3d(4, 5, 6));
    ros::Duration(1).sleep();
    
    for (int i = 0; i < 10; ++i) {
	cout << "!!" << endl;
	pubImu(Eigen::Vector3d(i*1, i*2, i*3), Eigen::Vector3d(4, 5, 6));
	ros::Duration(1).sleep();
    }
    }*/

}