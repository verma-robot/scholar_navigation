#include <ros/ros.h>

#include "scholar_gmapping.h"

int main(int argc, char** argv)
{
  ros::init(argc, argv, "scholar_gmapping");

  ScholarGMapping gn;
  gn.startLiveSlam();
  ros::spin();

  return(0);
}

