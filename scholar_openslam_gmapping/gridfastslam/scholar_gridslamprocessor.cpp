
#include "gmapping/gridfastslam/scholar_gridslamprocessor.h"


namespace Scholar_GMapping {

 
using namespace std;


  GridSlamProcessor::GridSlamProcessor()
  {
    
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;

  };

  GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor& gsp) :last_update_time_(0.0), m_particles(gsp.m_particles)
  {
    
    period_ = 5.0;

    m_obsSigmaGain=gsp.m_obsSigmaGain;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_minimumScore=gsp.m_minimumScore;    

    m_beams=gsp.m_beams;
    m_indexes=gsp.m_indexes;
    m_motionModel=gsp.m_motionModel;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_matcher=gsp.m_matcher;
   // m_map = gsp.m_map;
    m_count=gsp.m_count;
    m_readingCount=gsp.m_readingCount;
    m_lastPartPose=gsp.m_lastPartPose;
    m_pose=gsp.m_pose;
    m_odoPose=gsp.m_odoPose;
    m_linearDistance=gsp.m_linearDistance;
    m_angularDistance=gsp.m_angularDistance;
    m_neff=gsp.m_neff;	
    m_xmin=gsp.m_xmin;
    m_ymin=gsp.m_ymin;
    m_xmax=gsp.m_xmax;
    m_ymax=gsp.m_ymax;
    m_delta=gsp.m_delta;
    
    m_regScore=gsp.m_regScore;
    m_critScore=gsp.m_critScore;
    m_maxMove=gsp.m_maxMove;
    
    m_linearThresholdDistance=gsp.m_linearThresholdDistance;
    m_angularThresholdDistance=gsp.m_angularThresholdDistance;
    m_obsSigmaGain=gsp.m_obsSigmaGain;
    
    TNodeVector v=gsp.getTrajectories();
    for (unsigned int i=0; i<v.size(); i++)
    {
		  m_particles[i].node=v[i];
    }
    updateTreeWeights(false);
    
  };
  
  GridSlamProcessor* GridSlamProcessor::clone() const 
  {
    /*
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __func__ << ": performing preclone_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      if (a1.m_reference){
		PointerMap::iterator f=pmap.find(a1.m_reference);
		if (f==pmap.end())
		  pmap.insert(make_pair(a1.m_reference, 1));
		else
		  f->second++;
	      }
	    }
	  }
	}
	cerr << __func__ <<  ": Number of allocated chunks" << pmap.size() << endl;
	for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
	  assert(it->first->shares==(unsigned int)it->second);

	cerr << __func__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
*/
	GridSlamProcessor* cloned=new GridSlamProcessor(*this);
	/*
# ifdef MAP_CONSISTENCY_CHECK
	cerr << __func__ <<  ": trajectories end" << endl;
	cerr << __func__ << ": performing afterclone_fit_test" << endl;
	ParticleVector::const_iterator jt=cloned->m_particles.begin();
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const ScanMatcherMap& m2(jt->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
	  const HierarchicalArray2D<PointAccumulator>& h2(m2.storage());
	  jt++;
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      const autoptr< Array2D<PointAccumulator> >& a2(h2.m_cells[x][y]);
	      assert(a1.m_reference==a2.m_reference);
	      assert((!a1.m_reference) || !(a1.m_reference->shares%2));
	    }
	  }
	}
	cerr << __func__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
*/
	return cloned;
}
  
  GridSlamProcessor::~GridSlamProcessor()
  {
    for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
    {
      /*
#ifdef TREE_CONSISTENCY_CHECK		
      TNode* node=it->node;
      while(node)
	node=node->parent;
      cerr << "@" << endl;
#endif
*/
      if (it->node)delete it->node;
    }
    /*
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __func__ << ": performing predestruction_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
    for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      const ScanMatcherMap& m1(it->map);
      const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
      for (int x=0; x<h1.getXSize(); x++){
	for (int y=0; y<h1.getYSize(); y++){
	  const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	  if (a1.m_reference){
	    PointerMap::iterator f=pmap.find(a1.m_reference);
	    if (f==pmap.end())
	      pmap.insert(make_pair(a1.m_reference, 1));
	    else
	      f->second++;
	  }
	}
      }
    }
    cerr << __func__ << ": Number of allocated chunks" << pmap.size() << endl;
    for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
      assert(it->first->shares>=(unsigned int)it->second);
    cerr << __func__ << ": SUCCESS, the error is somewhere else" << endl;
# endif
*/
  }


		
  void GridSlamProcessor::setMatchingParameters (double urange, double range, double sigma, int kernsize, double lopt, double aopt, 
						 int iterations, double likelihoodSigma, double likelihoodGain, unsigned int likelihoodSkip)
  {
    m_obsSigmaGain=likelihoodGain;
    m_matcher.setMatchingParameters(urange, range, sigma, kernsize, lopt, aopt, iterations, likelihoodSigma, likelihoodSkip);
       
  }
  
void GridSlamProcessor::setMotionModelParameters(double srr, double srt, double str, double stt)
{
  m_motionModel.srr=srr;
  m_motionModel.srt=srt;
  m_motionModel.str=str;
  m_motionModel.stt=stt;	  
  
}

  void GridSlamProcessor::setUpdateDistances(double linear, double angular, double resampleThreshold)
  {
    m_linearThresholdDistance=linear; 
    m_angularThresholdDistance=angular;
    m_resampleThreshold=resampleThreshold;	
    
  }  

  GridSlamProcessor::Particle::Particle(const ScanMatcherMap& map): pose(0,0,0), weight(0), weightSum(0), gweight(0), previousIndex(0)
  {
    node=0;
  }
  
  
  void GridSlamProcessor::setSensorMap(const SensorMap& smap)
  {    
    
    SensorMap::const_iterator laser_it=smap.find(std::string("FLASER"));
    
    if (laser_it==smap.end())
    {
      laser_it=smap.find(std::string("ROBOTLASER1"));
      assert(laser_it!=smap.end());
    }
    const RangeSensor* rangeSensor=dynamic_cast<const RangeSensor*>((laser_it->second));
    assert(rangeSensor && rangeSensor->beams().size());
    
    m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
    double* angles=new double[rangeSensor->beams().size()];
    for (unsigned int i=0; i<m_beams; i++)
    {
      angles[i]=rangeSensor->beams()[i].pose.theta;
    }
    m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
    delete [] angles;
    
  }
  
  void GridSlamProcessor::init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, OrientedPoint initialPose)
  {
    
    m_xmin=xmin;
    m_ymin=ymin;
    m_xmax=xmax;
    m_ymax=ymax;
    m_delta=delta;    

    m_particles.clear();
    TNode* node=new TNode(initialPose, 0, 0, 0);
   
    m_map =new ScanMatcherMap(Point(xmin+xmax, ymin+ymax)*.5, xmax-xmin, ymax-ymin, delta);  

    for (unsigned int i=0; i<size; i++)
    {
      m_particles.push_back(Particle(*m_map));
      m_particles.back().pose=initialPose;
      m_particles.back().previousPose=initialPose;
      m_particles.back().setWeight(0);
      m_particles.back().previousIndex=0;
      
		  m_particles.back().node= node;
    }
    m_neff=(double)size;
    m_count=0;
    m_readingCount=0;
    m_linearDistance=m_angularDistance=0;
    
  }
  bool GridSlamProcessor::processScan(const RangeReading & reading, int adaptParticles)
  {
    
    OrientedPoint relPose=reading.getPose();
    if (!m_count)
    {
      m_lastPartPose=m_odoPose=  relPose;
    }    
    
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
    {
      OrientedPoint& pose(it->pose);
      pose=m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose);
    }    

    OrientedPoint move=relPose - m_lastPartPose;
    move.theta=atan2(sin(move.theta), cos(move.theta));
    m_linearDistance+=sqrt(move*move);
    m_angularDistance+=fabs(move.theta);

    m_odoPose=relPose;

    bool processed=false;
/*
    std::cout << " ........................................................................................................................." << std::endl;
    std::cout << " start                                         process scan                      "    << std::endl;

    std::cout << "linear_distance: "<<(float)m_linearDistance << std::endl;
    std::cout << "angular_distance: "<<(float)m_angularDistance << std::endl;
    std::cout << "process_scan_time: "<<reading.getTime() - last_update_time_ << std::endl;
    std::cout << "m_count :" <<m_count << std::endl;
    std::cout << "period_ :" <<period_ << std::endl;
*/
    if (! m_count 	|| m_linearDistance>=m_linearThresholdDistance 	|| m_angularDistance>=m_angularThresholdDistance    || (period_ >= 0.0 && (reading.getTime() - last_update_time_) > period_))
    {

      double * plainReading = new double[m_beams];
      for(unsigned int i=0; i<m_beams; i++)
      {
	      plainReading[i]=reading[i];
      }

      RangeReading* reading_copy =  new RangeReading(reading.size(),   &(reading[0]), static_cast<const RangeSensor*>(reading.getSensor()),  reading.getTime());

      if (m_count > 0)
      {
	      scanMatch(plainReading);	
	      updateTreeWeights(false);				
 	      resample(plainReading, adaptParticles, reading_copy);
      } 
      else 
      {
	      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
        {	

	        m_matcher.invalidateActiveArea();
	        m_matcher.computeActiveArea(*m_map, it->pose, plainReading);
	        m_matcher.registerScan(*m_map, it->pose, plainReading);
	  
	        TNode* node=new	TNode(it->pose, 0., it->node,  0);
          node->reading = reading_copy;
	        it->node=node;	  
	      }
        
        Particle best =    getParticles()[getBestParticleIndex()];  
        for(TNode* n = best.node;  n;  n = n->parent)
        {
          if(!n->reading)
          {
            continue;
          }
          m_matcher.invalidateActiveArea();
          m_matcher.computeActiveArea(*m_map, n->pose, &((*n->reading)[0]));
          m_matcher.registerScan(*m_map, n->pose, &((*n->reading)[0]));
        } 
        
      }
     updateTreeWeights(false);
      
      delete [] plainReading;

      m_lastPartPose=m_odoPose; 
      m_linearDistance=0;
      m_angularDistance=0;
      m_count = 1;
      processed=true;
      
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++)
      {
	      it->previousPose=it->pose;
      }

      last_update_time_ = reading.getTime();    

    }
    m_readingCount = 1;
    return processed;
    
  }
  
  int GridSlamProcessor::getBestParticleIndex() const
  {
    unsigned int bi=0;
    double bw=-std::numeric_limits<double>::max();
    for (unsigned int i=0; i<m_particles.size(); i++)
    {
      if (bw<m_particles[i].weightSum)
      {
	      bw=m_particles[i].weightSum;
	      bi=i;
      }
    }
    return (int) bi;
  }
  
};



