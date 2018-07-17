#include <ros/ros.h>
#include "../global.h"
#include "simulator.h"

using namespace Simulator;

int main(int argc, char** argv) 
{
    ros::init(argc, argv, "ins_gnss_simulator_node");
    ros::NodeHandle nh;
    ros::NodeHandle pnh("~");
    ROS_INFO("Runing...");
    GnssSimulator simulator(nh, pnh);
    simulator.run();
    ros::spin();  

    return(0);
}


      
//     simulator.imu_msg_.set(Eigen::Vector3d(1, 2, 3), Eigen::Vector3d(4, 5, 6));
//     simulator.pubImu();    
/*
    //----------
    ros::Rate loop_rate(1);
    while (ros::ok()) {
        ros::spinOnce();   
	
// 	simulator.pubImu();
	
	//TODO
// 	simulator.pubImu();
// 	for (int i = 0; i < 100; ++i) {
// 	    simulator.pubImu();
// 	}
// 	if (simulator.hasImu) {
// 	    simulator.pubImu();
// // 	    simulator.hasImu = false;
// 	}
	    


	    

        loop_rate.sleep();
    }    
    ROS_INFO("End!");*/