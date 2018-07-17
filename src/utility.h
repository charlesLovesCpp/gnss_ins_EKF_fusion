#pragma once

#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <eigen3/Eigen/Dense>



namespace Utility {

const static double micro_g_to_meters_per_second_squared = 9.80665e-6;
const static double DEG2RAD = 0.01745329252;
const static double RAD2DEG = 1/DEG2RAD;

const static double C = 299792458; 		///< Speed of light in m/s
const static double W_IE = 7.292115E-5;  	///< Earth rotation rate in rad/s
const static double R_0 = 6378137; 		///< WGS84 Equatorial radius in meters
const static double R_P = 6356752.31425; 	
const static double E = 0.0818191908425; 	///< WGS84 eccentricity
const static double F = 1 / 298.257223563;
const static double MU = 3.986004418e14; 	///< WGS84 Earth gravitational constant (m^3 s^-2)
const static double J_2 = 1.082627e-3; 	///< WGS84 Earth's second gravitational constant
const static Eigen::Matrix<double, 3, 1> G(0, 0, 9.81);	///< Gravity

template <typename Derived>
static Eigen::Quaternion<typename Derived::Scalar> deltaQ(const Eigen::MatrixBase<Derived> &theta)
{
    typedef typename Derived::Scalar Scalar_t;

    Eigen::Quaternion<Scalar_t> dq;
    Eigen::Matrix<Scalar_t, 3, 1> half_theta = theta;
    half_theta /= static_cast<Scalar_t>(2.0);
    dq.w() = static_cast<Scalar_t>(1.0);
    dq.x() = half_theta.x();
    dq.y() = half_theta.y();
    dq.z() = half_theta.z();
    return dq;
}

template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymmetric(const Eigen::MatrixBase<Derived> &q)
{
    Eigen::Matrix<typename Derived::Scalar, 3, 3> ans;
    ans << typename Derived::Scalar(0), -q(2), q(1),
	q(2), typename Derived::Scalar(0), -q(0),
	-q(1), q(0), typename Derived::Scalar(0);
    return ans;
}

template <typename Derived>
static Eigen::Quaternion<typename Derived::Scalar> positify(const Eigen::QuaternionBase<Derived> &q)
{
    //printf("a: %f %f %f %f", q.w(), q.x(), q.y(), q.z());
    //Eigen::Quaternion<typename Derived::Scalar> p(-q.w(), -q.x(), -q.y(), -q.z());
    //printf("b: %f %f %f %f", p.w(), p.x(), p.y(), p.z());
    //return q.template w() >= (typename Derived::Scalar)(0.0) ? q : Eigen::Quaternion<typename Derived::Scalar>(-q.w(), -q.x(), -q.y(), -q.z());
    return q;
}

template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 4, 4> Qleft(const Eigen::QuaternionBase<Derived> &q)
{
    Eigen::Quaternion<typename Derived::Scalar> qq = positify(q);
    Eigen::Matrix<typename Derived::Scalar, 4, 4> ans;
    ans(0, 0) = qq.w(), ans.template block<1, 3>(0, 1) = -qq.vec().transpose();
    ans.template block<3, 1>(1, 0) = qq.vec(), ans.template block<3, 3>(1, 1) = qq.w() * Eigen::Matrix<typename Derived::Scalar, 3, 3>::Identity() + skewSymmetric(qq.vec());
    return ans;
}

template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 4, 4> Qright(const Eigen::QuaternionBase<Derived> &p)
{
    Eigen::Quaternion<typename Derived::Scalar> pp = positify(p);
    Eigen::Matrix<typename Derived::Scalar, 4, 4> ans;
    ans(0, 0) = pp.w(), ans.template block<1, 3>(0, 1) = -pp.vec().transpose();
    ans.template block<3, 1>(1, 0) = pp.vec(), ans.template block<3, 3>(1, 1) = pp.w() * Eigen::Matrix<typename Derived::Scalar, 3, 3>::Identity() - skewSymmetric(pp.vec());
    return ans;
}

static Eigen::Vector3d R2ypr(const Eigen::Matrix3d &R)
{
    Eigen::Vector3d n = R.col(0);
    Eigen::Vector3d o = R.col(1);
    Eigen::Vector3d a = R.col(2);

    Eigen::Vector3d ypr(3);
    double y = atan2(n(1), n(0));
    double p = atan2(-n(2), n(0) * cos(y) + n(1) * sin(y));
    double r = atan2(a(0) * sin(y) - a(1) * cos(y), -o(0) * sin(y) + o(1) * cos(y));
    ypr(0) = y;
    ypr(1) = p;
    ypr(2) = r;

    return ypr / M_PI * 180.0;
}

template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 3, 3> ypr2R(const Eigen::MatrixBase<Derived> &ypr)
{
    typedef typename Derived::Scalar Scalar_t;

    Scalar_t y = ypr(0) / 180.0 * M_PI;
    Scalar_t p = ypr(1) / 180.0 * M_PI;
    Scalar_t r = ypr(2) / 180.0 * M_PI;

    Eigen::Matrix<Scalar_t, 3, 3> Rz;
    Rz << cos(y), -sin(y), 0,
	sin(y), cos(y), 0,
	0, 0, 1;

    Eigen::Matrix<Scalar_t, 3, 3> Ry;
    Ry << cos(p), 0., sin(p),
	0., 1., 0.,
	-sin(p), 0., cos(p);

    Eigen::Matrix<Scalar_t, 3, 3> Rx;
    Rx << 1., 0., 0.,
	0., cos(r), -sin(r),
	0., sin(r), cos(r);

    return Rz * Ry * Rx;
}

template<class Derived>
  inline Eigen::Matrix<typename Derived::Scalar, 4, 4> omegaMatJPL(const Eigen::MatrixBase<Derived> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
    return (
        Eigen::Matrix<typename Derived::Scalar, 4, 4>() <<
        0, vec[2], -vec[1], vec[0],
        -vec[2], 0, vec[0], vec[1],
        vec[1], -vec[0], 0, vec[2],
        -vec[0], -vec[1], -vec[2], 0
        ).finished();
  }

/// returns a matrix with angular velocities used for quaternion derivatives/integration with the Hamilton notation
/**
 The quaternion to be multiplied with this matrix has to be in the order x y z w !!!
 \param <vec> {3D vector with angular velocities}
 \return {4x4 matrix for multiplication with the quaternion}
 */
template<class Derived>
  inline Eigen::Matrix<typename Derived::Scalar, 4, 4> omegaMatHamilton(const Eigen::MatrixBase<Derived> & vec)
  {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
    return (
        Eigen::Matrix<typename Derived::Scalar, 4, 4>() <<
        0, -vec[2], vec[1], vec[0],
        vec[2], 0, -vec[0], vec[1],
        -vec[1], vec[0], 0, vec[2],
        -vec[0], -vec[1], -vec[2], 0
        ).finished();
  }
  
/// computes a quaternion from the 3-element small angle approximation theta
template<class Derived>
  Eigen::Quaternion<typename Derived::Scalar> quaternionFromSmallAngle(const Eigen::MatrixBase<Derived> & theta)
  {
  typedef typename Derived::Scalar Scalar;
  EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);
  const Scalar q_squared = theta.squaredNorm() / 4.0;

    if ( q_squared < 1)
    {
      return Eigen::Quaternion<Scalar>(sqrt(1 - q_squared), theta[0] * 0.5, theta[1] * 0.5, theta[2] * 0.5);
    }
    else
    {
      const Scalar w = 1.0 / sqrt(1 + q_squared);
      const Scalar f = w*0.5;
      return Eigen::Quaternion<Scalar>(w, theta[0] * f, theta[1] * f, theta[2] * f);
    }
  }
  
/// debug output to check misbehavior of Eigen
template<class T>
static bool checkForNumeric(const T & vec, int size, const std::string & info)
{
    for (int i = 0; i < size; i++)
    {
	if (std::isnan(vec[i]))
	{
	    std::cerr << "=== ERROR ===  " << info << ": NAN at index " << i << std::endl;
	    return false;
	}
	if (std::isinf(vec[i]))
	{
	    std::cerr << "=== ERROR ===  " << info << ": INF at index " << i << std::endl;
	    return false;
	}
    }
    return true;
}



static Eigen::Matrix<double, 3, 1> roundVector3(Eigen::Matrix<double, 3, 1> v) {
    Eigen::Matrix<double, 3, 1> r_vec(round(v.x()), round(v.y()), round(v.z()));
//     std::cout << v.transpose() << " | " << r_vec.transpose()w << std::endl;
    return r_vec;
}

template <typename Derived>
static int sign(const Derived& _val)
{
    if (_val > 0)
	return 1;
    else
	return -1;
}

namespace Earth {

static void ECEF2ECI(
    const Eigen::Matrix3d& _R_e_b, Eigen::Matrix3d& _R_i_b, 
    const Eigen::Vector3d& _p_e_eb, Eigen::Vector3d& _p_i_ib, 
    const Eigen::Vector3d& _v_e_eb, Eigen::Vector3d& _v_i_ib, 
    const int& _t)
{
    Eigen::Matrix3d R_i_e;
    R_i_e << cos(Utility::W_IE * _t), -sin(Utility::W_IE * _t), 0, 
	     sin(Utility::W_IE * _t), cos(Utility::W_IE * _t), 0,
	     0, 0, 1;
    _p_i_ib = R_i_e * _p_e_eb;
    _v_i_ib = R_i_e * (_v_e_eb + Utility::W_IE * Eigen::Vector3d(-_p_e_eb.y(), _p_e_eb.x(), 0));
    _R_i_b = R_i_e * _R_e_b;
}

static void ECEF2NED(
    const Eigen::Matrix3d& _R_e_b, Eigen::Matrix3d& _R_n_b, 
    const Eigen::Vector3d& _p_e_eb, Eigen::Vector3d& _p_b, 
    const Eigen::Vector3d& _v_e_eb, Eigen::Vector3d& _v_n_eb)
{
    double lambda_b = atan2(_p_e_eb.y(), _p_e_eb.x());	// latitude
    double k1 = sqrt(1 - pow(Utility::E, 2)) * abs(_p_e_eb.z());
    double k2 = pow(Utility::E, 2) * Utility::R_0;
    double beta = sqrt(pow(_p_e_eb.x(), 2) + pow(_p_e_eb.y(), 2));
    double E = (k1 - k2) / beta;
    double F = (k1 + k2) / beta;    
    double P = 4/3 * (E*F + 1);
    double Q = 2 * (pow(E, 2) - pow(F, 2));
    double D = pow(P, 3) + pow(Q, 2);
    double V = pow(sqrt(D) - Q, 1/3) - pow(sqrt(D) + Q, 1/3);
    double G = 0.5 * (sqrt(pow(E, 2) + V) + E);
    double T = sqrt(pow(G, 2) + (F - V * G) / (2 * G - E)) - G;
    double L_b = sign(_p_e_eb.z()) * atan((1 - pow(T, 2)) / (2 * T * sqrt (1 - pow(Utility::E, 2))));
    double h_b = (beta - Utility::R_0 * T) * cos(L_b) + (_p_e_eb.z() - sign(_p_e_eb.z()) * Utility::R_0 * sqrt(1 - pow(Utility::E, 2))) * sin (L_b);
    
    double cos_lat = cos(L_b);
    double sin_lat = sin(L_b);
    double cos_long = cos(lambda_b);
    double sin_long = sin(lambda_b);
    
    Eigen::Matrix3d R_n_e;
    R_n_e << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
	     -sin_long, cos_long, 0,
	     -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    _v_n_eb = R_n_e * _v_e_eb;
    _R_n_b = R_n_e * _R_e_b;
    _p_b.x() = L_b;
    _p_b.y() = lambda_b;
    _p_b.z() = h_b;
}

static void ECI2ECEF(const Eigen::Matrix3d& _R_i_b, Eigen::Matrix3d& _R_e_b, 
	      const Eigen::Vector3d& _p_i_ib, Eigen::Vector3d& _p_e_eb, 
	      const Eigen::Vector3d& _v_i_ib, Eigen::Vector3d& _v_e_eb, 
	      const int& _t)
{
    Eigen::Matrix3d R_i_e;
    R_i_e << cos(Utility::W_IE * _t), -sin(Utility::W_IE * _t), 0, 
	     sin(Utility::W_IE * _t), cos(Utility::W_IE * _t), 0,
	     0, 0, 1;
    _p_e_eb = R_i_e.transpose() * _p_i_ib;
    _v_e_eb = R_i_e.transpose() * _v_i_ib;
    _R_e_b = R_i_e.transpose() * _R_i_b;
}

static void NED2ECEF(
    const Eigen::Matrix3d& _R_n_b, Eigen::Matrix3d& _R_e_b, 
    const Eigen::Vector3d& _p_b, Eigen::Vector3d& _p_e_eb, 
    const Eigen::Vector3d& _v_n_eb, Eigen::Vector3d& _v_e_eb)
{
    double L_b = _p_b.x();
    double lambda_b = _p_b.y();
    double h_b = _p_b.z();
    double R_E = Utility::R_0 / sqrt(1 - pow(Utility::E*sin(L_b), 2));
    double cos_lat = cos(L_b);
    double sin_lat = sin(L_b);
    double cos_long = cos(lambda_b);
    double sin_long = sin(lambda_b);    
    _p_e_eb.x() = (R_E + h_b) * cos_lat * cos_long;
    _p_e_eb.y() = (R_E + h_b) * cos_lat * sin_long;
    _p_e_eb.z() = ((1 - pow(Utility::E, 2)) * R_E + h_b) * sin_lat;
    Eigen::Matrix3d R_n_e;
    R_n_e << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
	    -sin_long, cos_long, 0,
	    -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    _v_e_eb = R_n_e.transpose() * _v_n_eb;
    _R_e_b = R_n_e.transpose() * _R_n_b;
}

static void ECI2NED(const Eigen::Matrix3d& _R_i_b, Eigen::Matrix3d& _R_n_b, 
	     const Eigen::Vector3d& _p_i_ib, Eigen::Vector3d& _p_b, 
	     const Eigen::Vector3d& _v_i_ib, Eigen::Vector3d& _v_n_eb, 
	     const int& _t)
{
    Eigen::Matrix3d R_e_b;
    Eigen::Vector3d p_e_eb;
    Eigen::Vector3d v_e_eb;
    
    ECI2ECEF(_R_i_b, R_e_b, _p_i_ib, p_e_eb, _v_i_ib, v_e_eb, _t);
    ECEF2NED(R_e_b, _R_n_b, p_e_eb, _p_b, v_e_eb, _v_n_eb);
}

static void NED2ECI(const Eigen::Matrix3d& _R_n_b, Eigen::Matrix3d& _R_i_b, 
	     const Eigen::Vector3d& _p_b, Eigen::Vector3d& _p_i_ib, 
	     const Eigen::Vector3d& _v_n_eb, Eigen::Vector3d& _v_i_ib, 
	     const int& _t)
{
    Eigen::Matrix3d R_e_b;
    Eigen::Vector3d p_e_eb;
    Eigen::Vector3d v_e_eb;
    
    NED2ECEF(_R_n_b, R_e_b, _p_b, p_e_eb, _v_n_eb, v_e_eb);
    ECEF2ECI(R_e_b, _R_i_b, p_e_eb, _p_i_ib, v_e_eb, _v_i_ib, _t);
}    

static double getRe(const double& L_b) {
    // (2.106)
    return R_0 / sqrt(1 - (E * sin(L_b)) * (E * sin(L_b)));
}
 
static double getEarthRadius(const double& L_b) {
    // (2.137)
    double Re = getRe(L_b);
    return Re * sqrt(cos(L_b)*cos(L_b) + pow((1 - E*E),2) * sin(L_b) * sin(L_b));    
}

// (5.73) 右乘的姿态增量
static Eigen::Matrix<double, 3, 3> getRodreigueMatrix(const Eigen::Matrix<double, 3, 1>& a_b_ib) {
    double n = a_b_ib.norm();
    double nn = a_b_ib.squaredNorm();
    Eigen::Matrix<double, 3, 3> skew =  Utility::skewSymmetric(a_b_ib);
    return Eigen::Matrix<double, 3, 3>::Identity() + sin(n) / n * skew + (1 - cos(n)) / nn * skew * skew;
}
 

 
static Eigen::Matrix<double, 3, 1> getGravityECEF(const Eigen::Matrix<double, 3, 1>& p) {
    Eigen::Matrix<double, 3, 1> g;
    Eigen::Matrix<double, 3, 1> gamma;
    g.setZero();
    
    double z_scale;
    double mag_r = p.norm();
    if (mag_r == 0)
	return G;
    
    z_scale = 5 * pow(p.z() / mag_r, 2);
    gamma = -MU / pow(mag_r, 3) *(p + 1.5 * J_2 * pow(R_0 / mag_r, 2) * 
		    Eigen::Vector3d((1 - z_scale) * p.x(), (1 - z_scale) * p.y(), (3 - z_scale) * p.z()));
    g.head<2>() = gamma.head<2>() + pow(W_IE, 2) * p.head<2>();
    g.z() = gamma.z();
    
//     if ( (g - G).transpose() * (g - G) > 10) {
// 	ROS_WARN_STREAM("getGravityECEF: Wrong g");
// 	return G;
//     }
	 
    return g;
}

static Eigen::Matrix<double, 3, 1> getGravityNED(const Eigen::Vector3d& p_b) {
    Eigen::Vector3d g;
    double L_b = p_b.x(); 
    double h_b = p_b.z();
    double sinsqL = pow(sin(L_b), 2);
    double g0 = 9.7803253359 * (1 + 0.001931853 * sinsqL) / sqrt(1 - pow(Utility::E, 2) * sinsqL);
    g.x() = -8.08E-9 * h_b * sin(2 * L_b);
    g.y() = 0;
    g.z() = g0 * (1 - (2 / Utility::R_0) * (1 + Utility::F * (1 - 2 * sinsqL)
	    + (pow(Utility::W_IE, 2) * pow(Utility::R_0, 2) * Utility::R_P / Utility::MU)) * h_b 
	    + (3 * pow(h_b, 2) / pow(Utility::R_0, 2)));
    return g;
}

static Eigen::Matrix<double, 3, 1> getGravityECI(const Eigen::Vector3d& p_i_ib) {
    Eigen::Vector3d g;
    Eigen::Vector3d gamma;
    g.setZero();
    
    double z_scale;
    double mag_r = p_i_ib.norm();

    if (mag_r == 0) {
	return G;
    } 
    
    z_scale = 5 * pow(p_i_ib.z() / mag_r, 2);
    gamma = -Utility::MU / pow(mag_r, 3) *(p_i_ib + 1.5 * Utility::J_2 * pow(Utility::R_0 / mag_r, 2) * 
			Eigen::Vector3d((1 - z_scale) * p_i_ib.x(), (1 - z_scale) * p_i_ib.y(), (3 - z_scale) * p_i_ib.z()));
    g = gamma;
    
    return g;   
}

} // end namespace EARTH

 
} // end namespace Utility

#endif