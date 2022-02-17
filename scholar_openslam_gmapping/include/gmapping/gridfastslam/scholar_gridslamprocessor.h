#ifndef GRIDSLAMPROCESSOR_H
#define GRIDSLAMPROCESSOR_H

#include <climits>
#include <limits>
#include <fstream>
#include <vector>
#include <deque>
#include <gmapping/particlefilter/scholar_particlefilter.h>
#include <gmapping/utils/point.h>
#include <gmapping/utils/macro_params.h>

#include <gmapping/log/scholar_sensorlog.h>

#include <gmapping/sensor/sensor_range/rangesensor.h>
#include <gmapping/sensor/sensor_range/rangereading.h>
#include <gmapping/scanmatcher/scholar_scanmatcher.h>
#include "gmapping/gridfastslam/scholar_motionmodel.h"
#include <gmapping/gridfastslam/gridfastslam_export.h>

#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
namespace Scholar_GMapping {


  class GRIDFASTSLAM_EXPORT GridSlamProcessor
  {
    public:

      GridSlamProcessor();

      ScanMatcherMap* m_map;

     // ScanMatcherMap llmap(Point(0,0)*.5, 20, 20, 0.05);
      //Point(const double x,const double y);
      struct TNode
      {
     
        TNode(const OrientedPoint& pose, double weight, TNode* parent=0, unsigned int childs=0);
        ~TNode();
        OrientedPoint pose; 
      
        double weight;

        double accWeight;

        double gweight;

        TNode* parent;

        const RangeReading* reading;

        unsigned int childs;

        mutable unsigned int visitCounter;

        mutable bool flag;
      };
    
      typedef std::vector<GridSlamProcessor::TNode*> TNodeVector;
      typedef std::deque<GridSlamProcessor::TNode*> TNodeDeque;

      struct Particle
      {
    
        Particle(const ScanMatcherMap& map);

        inline operator double() const {return weight;}
        inline operator OrientedPoint() const {return pose;}
     
        inline void setWeight(double w) {weight=w;}
        OrientedPoint pose;

        OrientedPoint previousPose;

        double weight;

        double weightSum;

        double gweight;

        int previousIndex;

        TNode* node; 
      };
	
    
      typedef std::vector<Particle> ParticleVector;
    
      GridSlamProcessor* clone() const;    
      virtual ~GridSlamProcessor();
    
      void setSensorMap(const SensorMap& smap);
      void init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, OrientedPoint initialPose=OrientedPoint(0,0,0));
      void setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt, double aopt, 	int iterations, double likelihoodSigma=1, double likelihoodGain=1, unsigned int likelihoodSkip=0);
      void setMotionModelParameters(double srr, double srt, double str, double stt);

      void setUpdateDistances(double linear, double angular, double resampleThreshold);
      void setUpdatePeriod(double p) {period_=p;}
    
      bool processScan(const RangeReading & reading, int adaptParticles=0);
    
      TNodeVector getTrajectories() const;
      void integrateScanSequence(TNode* node);
    
      ScanMatcher m_matcher;
    
      inline const ParticleVector& getParticles() const {return m_particles; }
    
      inline const std::vector<unsigned int>& getIndexes() const{return m_indexes; }
      int getBestParticleIndex() const;
    
      MEMBER_PARAM_SET_GET(m_matcher, double, laserMaxRange, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, usableRange, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher,double, gaussianSigma, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher,double, likelihoodSigma, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, int,    kernelSize, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, optAngularDelta, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, optLinearDelta, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, unsigned int, optRecursiveIterations, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, unsigned int, likelihoodSkip, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, llsamplerange, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, lasamplerange, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, llsamplestep, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, double, lasamplestep, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, bool, generateMap, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, bool, enlargeStep, protected, public, public);

      MEMBER_PARAM_SET_GET(m_matcher, OrientedPoint, laserPose, protected, public, public);


      STRUCT_PARAM_SET_GET(m_motionModel, double, srr, protected, public, public);

      STRUCT_PARAM_SET_GET(m_motionModel, double, srt, protected, public, public);

      STRUCT_PARAM_SET_GET(m_motionModel, double, str, protected, public, public);

      STRUCT_PARAM_SET_GET(m_motionModel, double, stt, protected, public, public);
		
      PARAM_SET_GET(double, minimumScore, protected, public, public);

  protected:

      GridSlamProcessor(const GridSlamProcessor& gsp);
 
      unsigned int m_beams;
      double last_update_time_;
      double period_;
    
      ParticleVector m_particles;

      std::vector<unsigned int> m_indexes;

      std::vector<double> m_weights;
    
      MotionModel m_motionModel;

      PARAM_SET_GET(double, resampleThreshold, protected, public, public);
      
      int  m_count, m_readingCount;
      OrientedPoint m_lastPartPose;
      OrientedPoint m_odoPose;
      OrientedPoint m_pose;

      double m_linearDistance, m_angularDistance;
      PARAM_GET(double, neff, protected, public);
      
      PARAM_GET(double, xmin, protected, public);
      PARAM_GET(double, ymin, protected, public);
      PARAM_GET(double, xmax, protected, public);
      PARAM_GET(double, ymax, protected, public);
      PARAM_GET(double, delta, protected, public);
	
      PARAM_SET_GET(double, regScore, protected, public, public);
      PARAM_SET_GET(double, critScore, protected, public, public);
      PARAM_SET_GET(double, maxMove, protected, public, public);
	
      PARAM_SET_GET(double, linearThresholdDistance, protected, public, public);

      PARAM_SET_GET(double, angularThresholdDistance, protected, public, public);
    
      PARAM_SET_GET(double, obsSigmaGain, protected, public, public);
	
    
  private:
    
      inline void scanMatch(const double *plainReading);
      inline void normalize();
    
      inline bool resample(const double* plainReading, int adaptParticles, 
	  	const RangeReading* rr=0);
    
    
      void updateTreeWeights(bool weightsAlreadyNormalized = false);
      void resetTree();
      double propagateWeights();
    
  };

typedef std::multimap<const GridSlamProcessor::TNode*, GridSlamProcessor::TNode*> TNodeMultimap;


#include "gmapping/gridfastslam/scholar_gridslamprocessor.hxx"

};

#endif
