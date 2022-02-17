#include "scholar_gmapping.h"

#include <iostream>

#include <time.h>

#include "ros/ros.h"
#include "ros/console.h"
#include "nav_msgs/MapMetaData.h"

#include "gmapping/sensor/sensor_range/rangesensor.h"
#include "gmapping/sensor/sensor_odometry/odometrysensor.h"

#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define MAP_IDX(sx, i, j) ((sx) * (j) + (i))

ScholarGMapping::ScholarGMapping():
  map_to_odom_(tf::Transform(tf::createQuaternionFromRPY( 0, 0, 0 ), tf::Point(0, 0, 0 ))),
  laser_count_(0), private_nh_("~"), scan_filter_sub_(NULL), scan_filter_(NULL), transform_thread_(NULL)
{
  seed_ = time(NULL);
  init();
}

ScholarGMapping::ScholarGMapping(ros::NodeHandle& nh, ros::NodeHandle& pnh):
  map_to_odom_(tf::Transform(tf::createQuaternionFromRPY( 0, 0, 0 ), tf::Point(0, 0, 0 ))),
  laser_count_(0),node_(nh), private_nh_(pnh), scan_filter_sub_(NULL), scan_filter_(NULL), transform_thread_(NULL)
{
  seed_ = time(NULL);
  init();
}

ScholarGMapping::ScholarGMapping(long unsigned int seed, long unsigned int max_duration_buffer):
  map_to_odom_(tf::Transform(tf::createQuaternionFromRPY( 0, 0, 0 ), tf::Point(0, 0, 0 ))),
  laser_count_(0), private_nh_("~"), scan_filter_sub_(NULL), scan_filter_(NULL), transform_thread_(NULL),
  seed_(seed), tf_(ros::Duration(max_duration_buffer))
{
  init();
}


void ScholarGMapping::init()
{
  
  gsp_laser_ = NULL;
  gsp_odom_ = NULL;

  got_first_scan_ = false;
  got_map_ = false;
  
  if(!private_nh_.getParam("throttle_scans", throttle_scans_))throttle_scans_ = 1;
  if(!private_nh_.getParam("base_frame", base_frame_))base_frame_ = "base_link";
  if(!private_nh_.getParam("map_frame", map_frame_))map_frame_ = "map";
  if(!private_nh_.getParam("odom_frame", odom_frame_))odom_frame_ = "odom";

  private_nh_.param("transform_publish_period", transform_publish_period_, 0.05);

  double tmp;
  if(!private_nh_.getParam("map_update_interval", tmp))tmp = 5.0;
    map_update_interval_.fromSec(tmp);
  
  maxUrange_ = 0.0;  maxRange_ = 0.0; 
  if(!private_nh_.getParam("minimumScore", minimum_score_))minimum_score_ = 0;
  if(!private_nh_.getParam("sigma", sigma_))sigma_ = 0.05;
  if(!private_nh_.getParam("kernelSize", kernelSize_))kernelSize_ = 1;
  if(!private_nh_.getParam("lstep", lstep_))lstep_ = 0.05;
  if(!private_nh_.getParam("astep", astep_))astep_ = 0.05;
  if(!private_nh_.getParam("iterations", iterations_))iterations_ = 5;
  if(!private_nh_.getParam("lsigma", lsigma_))lsigma_ = 0.075;
  if(!private_nh_.getParam("ogain", ogain_))ogain_ = 3.0;
  if(!private_nh_.getParam("lskip", lskip_))lskip_ = 0;
  if(!private_nh_.getParam("srr", srr_))srr_ = 0.1;
  if(!private_nh_.getParam("srt", srt_))srt_ = 0.2;
  if(!private_nh_.getParam("str", str_))str_ = 0.1;
  if(!private_nh_.getParam("stt", stt_))stt_ = 0.2;
  if(!private_nh_.getParam("linearUpdate", linearUpdate_))linearUpdate_ = 1.0;
  if(!private_nh_.getParam("angularUpdate", angularUpdate_))angularUpdate_ = 0.5;
  if(!private_nh_.getParam("temporalUpdate", temporalUpdate_))temporalUpdate_ = -1.0;
  if(!private_nh_.getParam("resampleThreshold", resampleThreshold_))resampleThreshold_ = 0.5;
  if(!private_nh_.getParam("particles", particles_))particles_ = 30;

  if(!private_nh_.getParam("totalmap_xmin", totalmap_xmin))totalmap_xmin = -100.0;
  if(!private_nh_.getParam("totalmap_ymin", totalmap_ymin))totalmap_ymin = -100.0;
  if(!private_nh_.getParam("totalmap_xmax", totalmap_xmax))totalmap_xmax = 100.0;
  if(!private_nh_.getParam("totalmap_ymax", totalmap_ymax))totalmap_ymax = 100.0;


  if(!private_nh_.getParam("xmin", xmin_))xmin_ = -10.0;
  if(!private_nh_.getParam("ymin", ymin_))ymin_ = -10.0;
  if(!private_nh_.getParam("xmax", xmax_))xmax_ = 10.0;
  if(!private_nh_.getParam("ymax", ymax_))ymax_ = 10.0;
  if(!private_nh_.getParam("delta", delta_))delta_ = 0.05;
  if(!private_nh_.getParam("occ_thresh", occ_thresh_))occ_thresh_ = 0.25;
  if(!private_nh_.getParam("llsamplerange", llsamplerange_))llsamplerange_ = 0.01;
  if(!private_nh_.getParam("llsamplestep", llsamplestep_))llsamplestep_ = 0.01;
  if(!private_nh_.getParam("lasamplerange", lasamplerange_))lasamplerange_ = 0.005;
  if(!private_nh_.getParam("lasamplestep", lasamplestep_))lasamplestep_ = 0.005;    
  if(!private_nh_.getParam("tf_delay", tf_delay_)) tf_delay_ = transform_publish_period_;

  gsp_ = new Scholar_GMapping::GridSlamProcessor();
  ROS_ASSERT(gsp_);

  tfB_ = new tf::TransformBroadcaster();
  ROS_ASSERT(tfB_);

}

void ScholarGMapping::startLiveSlam()
{
  sst_ = node_.advertise<nav_msgs::OccupancyGrid>("map", 1, true);
  sstm_ = node_.advertise<nav_msgs::MapMetaData>("map_metadata", 1, true);

  scan_filter_sub_ = new message_filters::Subscriber<sensor_msgs::LaserScan>(node_, "scan", 5);
  scan_filter_ = new tf::MessageFilter<sensor_msgs::LaserScan>(*scan_filter_sub_, tf_, odom_frame_, 5);

  scan_filter_->registerCallback(boost::bind(&ScholarGMapping::laserCallback, this, _1));

  transform_thread_ = new boost::thread(boost::bind(&ScholarGMapping::publishLoop, this, transform_publish_period_));
}


void ScholarGMapping::publishLoop(double transform_publish_period)
{
  if(transform_publish_period == 0)return;

  ros::Rate r(1.0 / transform_publish_period);
  while(ros::ok())
  {
    publishTransform();
    r.sleep();
  }
}

ScholarGMapping::~ScholarGMapping()
{
  if(transform_thread_)
  {
    transform_thread_->join();
    delete transform_thread_;
  }

  delete gsp_;
  if(gsp_laser_)delete gsp_laser_;
  if(gsp_odom_)delete gsp_odom_;
  if (scan_filter_)delete scan_filter_;
  if (scan_filter_sub_)delete scan_filter_sub_;
}

bool ScholarGMapping::getOdomPose(Scholar_GMapping::OrientedPoint& gmap_pose, const ros::Time& t)
{
  centered_laser_pose_.stamp_ = t;
  tf::Stamped<tf::Transform> odom_pose;
  try
  {
    tf_.transformPose(odom_frame_, centered_laser_pose_, odom_pose);
  }
  catch(tf::TransformException e)
  {
    ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
    return false;
  }
  double yaw = tf::getYaw(odom_pose.getRotation());

  gmap_pose = Scholar_GMapping::OrientedPoint(odom_pose.getOrigin().x(),   odom_pose.getOrigin().y(),  yaw);

  return true;
}

bool ScholarGMapping::initMapper(const sensor_msgs::LaserScan& scan)
{
  laser_frame_ = scan.header.frame_id;
  tf::Stamped<tf::Pose> ident;
  tf::Stamped<tf::Transform> laser_pose;
  ident.setIdentity();
  ident.frame_id_ = laser_frame_;
  ident.stamp_ = scan.header.stamp;
  try
  {
    tf_.transformPose(base_frame_, ident, laser_pose);
  }
  catch(tf::TransformException e)
  {
    ROS_WARN("Failed to compute laser pose, aborting initialization (%s)", e.what());
    return false;
  }

  tf::Vector3 v;
  v.setValue(0, 0, 1 + laser_pose.getOrigin().z());
  tf::Stamped<tf::Vector3> up(v, scan.header.stamp,  base_frame_);
  try
  {
    tf_.transformPoint(laser_frame_, up, up);
    ROS_DEBUG("Z-Axis in sensor frame: %.3f", up.z());
  }
  catch(tf::TransformException& e)
  {
    ROS_WARN("Unable to determine orientation of laser: %s",e.what());
    return false;
  }
  
  if (fabs(fabs(up.z()) - 1) > 0.001)
  {
    ROS_WARN("Laser has to be mounted planar! Z-coordinate has to be 1 or -1, but gave: %.5f",up.z());
    return false;
  }

  gsp_laser_beam_count_ = scan.ranges.size();

  double angle_center = (scan.angle_min + scan.angle_max)/2;

  if (up.z() > 0)
  {
    do_reverse_range_ = scan.angle_min > scan.angle_max;
    centered_laser_pose_ = tf::Stamped<tf::Pose>(tf::Transform(tf::createQuaternionFromRPY(0,0,angle_center),   tf::Vector3(0,0,0)), ros::Time::now(), laser_frame_);
    ROS_INFO("Laser is mounted upwards.");
  }
  else
  {
    do_reverse_range_ = scan.angle_min < scan.angle_max;
    centered_laser_pose_ = tf::Stamped<tf::Pose>(tf::Transform(tf::createQuaternionFromRPY(M_PI,0,-angle_center),  tf::Vector3(0,0,0)), ros::Time::now(), laser_frame_);
    ROS_INFO("Laser is mounted upside down.");
  }

  laser_angles_.resize(scan.ranges.size());
  double theta = - std::fabs(scan.angle_min - scan.angle_max)/2;
  for(unsigned int i=0; i<scan.ranges.size(); ++i)
  {
    laser_angles_[i]=theta;
    theta += std::fabs(scan.angle_increment);
  }

  Scholar_GMapping::OrientedPoint gmap_pose(0, 0, 0);

  ros::NodeHandle private_nh_("~");
  if(!private_nh_.getParam("maxRange", maxRange_))maxRange_ = scan.range_max - 0.01;
  if(!private_nh_.getParam("maxUrange", maxUrange_))maxUrange_ = maxRange_;
  gsp_laser_ = new Scholar_GMapping::RangeSensor("FLASER", gsp_laser_beam_count_,  fabs(scan.angle_increment), gmap_pose, 0.0,   maxRange_);
  ROS_ASSERT(gsp_laser_);

  Scholar_GMapping::SensorMap smap;
  smap.insert(make_pair(gsp_laser_->getName(), gsp_laser_));
  gsp_->setSensorMap(smap);

  gsp_odom_ = new Scholar_GMapping::OdometrySensor(odom_frame_);
  ROS_ASSERT(gsp_odom_);

  Scholar_GMapping::OrientedPoint initialPose;
  if(!getOdomPose(initialPose, scan.header.stamp))
  {
    ROS_WARN("Unable to determine inital pose of laser! Starting point will be set to zero.");
    initialPose = Scholar_GMapping::OrientedPoint(0.0, 0.0, 0.0);
  }

  gsp_->setMatchingParameters(maxUrange_, maxRange_, sigma_, kernelSize_, lstep_, astep_, iterations_, lsigma_, ogain_, lskip_);
  gsp_->setMotionModelParameters(srr_, srt_, str_, stt_);
  gsp_->setUpdateDistances(linearUpdate_, angularUpdate_, resampleThreshold_);
  gsp_->setUpdatePeriod(temporalUpdate_);
  gsp_->setgenerateMap(false);
  gsp_->GridSlamProcessor::init(particles_, xmin_, ymin_, xmax_, ymax_, delta_, initialPose);

  gsp_->setllsamplerange(llsamplerange_);
  gsp_->setllsamplestep(llsamplestep_);
  gsp_->setlasamplerange(lasamplerange_);
  gsp_->setlasamplestep(lasamplestep_);
  gsp_->setminimumScore(minimum_score_);

  Scholar_GMapping::sampleGaussian(1,seed_);

  ROS_INFO("Initialization complete");

  return true;
}

bool ScholarGMapping::addScan(const sensor_msgs::LaserScan& scan, Scholar_GMapping::OrientedPoint& gmap_pose)
{
  if(!getOdomPose(gmap_pose, scan.header.stamp)) return false;
  
  if(scan.ranges.size() != gsp_laser_beam_count_) return false;

  double* ranges_double = new double[scan.ranges.size()];
  if (do_reverse_range_)
  {
    int num_ranges = scan.ranges.size();
    for(int i=0; i < num_ranges; i++)
    {
      if(scan.ranges[num_ranges - i - 1] < scan.range_min) ranges_double[i] = scan.range_max;
      else   ranges_double[i] = (double)scan.ranges[num_ranges - i - 1];
    }
  } 
  else 
  {
    for(unsigned int i=0; i < scan.ranges.size(); i++)
    {
      if(scan.ranges[i] < scan.range_min)ranges_double[i] = scan.range_max;
      else  ranges_double[i] = (double)scan.ranges[i];
    }
  }

  Scholar_GMapping::RangeReading reading(scan.ranges.size(), ranges_double,  gsp_laser_, scan.header.stamp.toSec());
  delete[] ranges_double;
  reading.setPose(gmap_pose);

  return gsp_->processScan(reading);
}

void ScholarGMapping::laserCallback(const sensor_msgs::LaserScan::ConstPtr& scan)
{

  sensor_msgs::LaserScan lase = *scan;
  Scholar_GMapping::OrientedPoint mpose;

  laser_count_++;//
  if ((laser_count_ % throttle_scans_) != 0) return;

  if(!got_first_scan_)
  {
    if(!initMapper(*scan))return;//    
    got_first_scan_ = true;
  }

  Scholar_GMapping::OrientedPoint odom_pose;
  float delt_time = 0;

  ros::Time now_time = ros::Time::now();

  bool add_scan_ = addScan(*scan, odom_pose );

  ros::Time last_time = ros::Time::now();

  int laser_i = 0;
  if(add_scan_)//处理成功
  {

    ROS_DEBUG("scan processed");

    mpose= gsp_->getParticles()[gsp_->getBestParticleIndex()].pose;
  
    tf::Transform laser_to_map = tf::Transform(tf::createQuaternionFromRPY(0, 0, mpose.theta), tf::Vector3(mpose.x, mpose.y, 0.0)).inverse();
    tf::Transform odom_to_laser = tf::Transform(tf::createQuaternionFromRPY(0, 0, odom_pose.theta), tf::Vector3(odom_pose.x, odom_pose.y, 0.0));

    map_to_odom_mutex_.lock();
    map_to_odom_ = (odom_to_laser * laser_to_map).inverse();
    map_to_odom_mutex_.unlock();

    updateMap(*scan);   //更新地图
    
  } 
  else 
  {
  }
  
}

void ScholarGMapping::updateMap(const sensor_msgs::LaserScan& scan)
{
  boost::mutex::scoped_lock map_lock (map_mutex_);
  Scholar_GMapping::ScanMatcher matcher;

  matcher.setLaserParameters(scan.ranges.size(), &(laser_angles_[0]),  gsp_laser_->getPose());

  matcher.setlaserMaxRange(maxRange_);
  matcher.setusableRange(maxUrange_);
  matcher.setgenerateMap(true);

  Scholar_GMapping::GridSlamProcessor::Particle best =  gsp_->getParticles()[gsp_->getBestParticleIndex()];  
  if(!got_map_) 
  {
    map_.map.info.resolution = delta_;
    map_.map.info.origin.position.x = 0.0;
    map_.map.info.origin.position.y = 0.0;
    map_.map.info.origin.position.z = 0.0;
    map_.map.info.origin.orientation.x = 0.0;
    map_.map.info.origin.orientation.y = 0.0;
    map_.map.info.origin.orientation.z = 0.0;
    map_.map.info.origin.orientation.w = 1.0;
  } 

  Scholar_GMapping::Point center;
  center.x=(totalmap_xmin + totalmap_xmax) / 2.0;
  center.y=(totalmap_ymin + totalmap_ymax) / 2.0;

  Scholar_GMapping::ScanMatcherMap smap(center, totalmap_xmin, totalmap_ymin, totalmap_xmax, totalmap_ymax,  delta_);

  ROS_DEBUG("Trajectory tree:");
  for(Scholar_GMapping::GridSlamProcessor::TNode* n = best.node;  n;  n = n->parent)
  {

    ROS_DEBUG("  %.3f %.3f %.3f",    n->pose.x,        n->pose.y,          n->pose.theta);
    if(!n->reading)
    {
      ROS_DEBUG("Reading is NULL");
      continue;
    }
    matcher.invalidateActiveArea();
    matcher.computeActiveArea(smap, n->pose, &((*n->reading)[0]));
    matcher.registerScan(smap, n->pose, &((*n->reading)[0]));
  }

  if(map_.map.info.width != (unsigned int) smap.getMapSizeX() || map_.map.info.height != (unsigned int) smap.getMapSizeY()) 
  {

    Scholar_GMapping::Point wmin = smap.map2world(Scholar_GMapping::IntPoint(0, 0));
    Scholar_GMapping::Point wmax = smap.map2world(Scholar_GMapping::IntPoint(smap.getMapSizeX(), smap.getMapSizeY()));
    
    ROS_DEBUG("map size is now %dx%d pixels (%f,%f)-(%f, %f)", smap.getMapSizeX(), smap.getMapSizeY(),  xmin_, ymin_, xmax_, ymax_);

    map_.map.info.width = smap.getMapSizeX();
    map_.map.info.height = smap.getMapSizeY();
    map_.map.info.origin.position.x = wmin.x;
    map_.map.info.origin.position.y = wmin.y;
    map_.map.data.resize(map_.map.info.width * map_.map.info.height);

    ROS_DEBUG("map origin: (%f, %f)", map_.map.info.origin.position.x, map_.map.info.origin.position.y);
  }

  for(int x=0; x < smap.getMapSizeX(); x++)
  {
    for(int y=0; y < smap.getMapSizeY(); y++)
    {
      Scholar_GMapping::IntPoint p(x, y);
      double occ=smap.cell(p);
      assert(occ <= 1.0);
      if(occ < 0) map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = -1;
      else if(occ > occ_thresh_)
      {
        map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = 100;
      }
      else   map_.map.data[MAP_IDX(map_.map.info.width, x, y)] = 0;
    }
  }
  got_map_ = true;

  map_.map.header.stamp = ros::Time::now();
  map_.map.header.frame_id = tf_.resolve( map_frame_ );

  sst_.publish(map_.map);
  sstm_.publish(map_.map.info);

  map_lock.unlock();
}

void ScholarGMapping::publishTransform()
{
  map_to_odom_mutex_.lock();
  ros::Time tf_expiration = ros::Time::now() + ros::Duration(tf_delay_);
  tfB_->sendTransform( tf::StampedTransform (map_to_odom_, tf_expiration, map_frame_, odom_frame_));
  map_to_odom_mutex_.unlock();
}
