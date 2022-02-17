#include "ros/ros.h"
#include "sensor_msgs/LaserScan.h"
#include "std_msgs/Float64.h"
#include "nav_msgs/GetMap.h"
#include "tf/transform_listener.h"
#include "tf/transform_broadcaster.h"
#include "message_filters/subscriber.h"
#include "tf/message_filter.h"

#include "gmapping/gridfastslam/scholar_gridslamprocessor.h"
#include "gmapping/sensor/sensor_base/sensor.h"
#include <gmapping/utils/point.h>

#include <boost/thread.hpp>

class ScholarGMapping
{
  public:
    ScholarGMapping();
    ScholarGMapping(ros::NodeHandle& nh, ros::NodeHandle& pnh);
    ScholarGMapping(unsigned long int seed, unsigned long int max_duration_buffer);
    ~ScholarGMapping();

    void init();
    void startLiveSlam();
    void publishTransform();
  
    void laserCallback(const sensor_msgs::LaserScan::ConstPtr& scan);
   
    void publishLoop(double transform_publish_period);

  private:
    ros::NodeHandle node_;
    ros::Publisher sst_;
    ros::Publisher sstm_;

    tf::TransformListener tf_;
    message_filters::Subscriber<sensor_msgs::LaserScan>* scan_filter_sub_;

    tf::MessageFilter<sensor_msgs::LaserScan>* scan_filter_;

    tf::TransformBroadcaster* tfB_;

    Scholar_GMapping::GridSlamProcessor* gsp_;
    Scholar_GMapping::RangeSensor* gsp_laser_;

    std::vector<double> laser_angles_;
    tf::Stamped<tf::Pose> centered_laser_pose_;
   
    bool do_reverse_range_;
    unsigned int gsp_laser_beam_count_;
    Scholar_GMapping::OdometrySensor* gsp_odom_;

    bool got_first_scan_;

    bool got_map_;
    nav_msgs::GetMap::Response map_;

    ros::Duration map_update_interval_;
    tf::Transform map_to_odom_;
    boost::mutex map_to_odom_mutex_;
    boost::mutex map_mutex_;

    int laser_count_;
    int throttle_scans_;

    boost::thread* transform_thread_;

    std::string base_frame_;
    std::string laser_frame_;
    std::string map_frame_;
    std::string odom_frame_;

    void updateMap(const sensor_msgs::LaserScan& scan);
    bool getOdomPose(Scholar_GMapping::OrientedPoint& gmap_pose, const ros::Time& t);
    bool initMapper(const sensor_msgs::LaserScan& scan);
    bool addScan(const sensor_msgs::LaserScan& scan, Scholar_GMapping::OrientedPoint& gmap_pose);

    double maxRange_;
    double maxUrange_;
    double maxrange_;
    double minimum_score_;
    double sigma_;
    int kernelSize_;
    double lstep_;
    double astep_;
    int iterations_;
    double lsigma_;
    double ogain_;
    int lskip_;

    double srr_;
    double srt_;
    double str_;
    double stt_;

    double linearUpdate_;
    double angularUpdate_;
    double temporalUpdate_;
    double resampleThreshold_;
    int particles_;
    double xmin_;
    double ymin_;
    double xmax_;
    double ymax_;
    double delta_;
    double totalmap_xmin;
    double totalmap_ymin;
    double totalmap_xmax;
    double totalmap_ymax;

    double occ_thresh_;
    double llsamplerange_;
    double llsamplestep_;
    double lasamplerange_;
    double lasamplestep_;
    
    ros::NodeHandle private_nh_;
    
    unsigned long int seed_;
    
    double transform_publish_period_;
    double tf_delay_;
};
