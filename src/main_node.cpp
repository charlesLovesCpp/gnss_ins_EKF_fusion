#include <ros/ros.h>
#include "global.h"
#include "fusion_handler.h"

using namespace Fusion;

ros::Publisher pub_imu;

int main(int argc, char** argv) 
{
    ros::init(argc, argv, "main_node");
    ros::NodeHandle nh;
    ros::NodeHandle pnh("~");
    FusionHandler fh(nh, pnh);
    ROS_INFO("initialization complete...Looping");
    //----------
    ros::Rate loop_rate(1);
    while (ros::ok()) {
        ros::spinOnce();   
        loop_rate.sleep();
    }
    return(0);
}




//     pub_imu = nh.advertise<sensor_msgs::Imu> ("/imu1", 1);  

/*	Eigen::Vector3d f_ib_b(1, 2, 3);
	Eigen::Vector3d omega_ib_b(4, 5, 6);
	sensor_msgs::Imu imu_msg;
	imu_msg.header.stamp = ros::Time::now();
	imu_msg.header.frame_id = "base_link";    
	imu_msg.linear_acceleration.x = f_ib_b.x();
	imu_msg.linear_acceleration.y = f_ib_b.y();
	imu_msg.linear_acceleration.z = f_ib_b.z();
	imu_msg.angular_velocity.x = omega_ib_b.x();
	imu_msg.angular_velocity.y = omega_ib_b.y();
	imu_msg.angular_velocity.z = omega_ib_b.z();   
	pub_imu.publish(imu_msg);*/	