#include <cstring>
#include <limits>
#include <list>
#include <iostream>

#include "gmapping/scanmatcher/scholar_scanmatcher.h"
#include "gmapping/scanmatcher/scholar_gridlinetraversal.h"

namespace Scholar_GMapping {

using namespace std;

const double ScanMatcher::nullLikelihood=-.5;

ScanMatcher::ScanMatcher(): m_laserPose(0,0,0)
{
	m_laserBeams=0;
	m_optRecursiveIterations=3;
	m_activeAreaComputed=false;

	m_llsamplerange=0.01;
	m_llsamplestep=0.01;
	m_lasamplerange=0.005;
	m_lasamplestep=0.005;
	m_enlargeStep=10.;
	m_fullnessThreshold=0.1;
	m_angularOdometryReliability=0.0;
	m_linearOdometryReliability=0.0;
	m_freeCellRatio=sqrt(2.);
	m_initialBeamsSkip=0;
	
   m_linePoints = new IntPoint[20000];
}

ScanMatcher::~ScanMatcher()
{
	delete [] m_linePoints;
}

void ScanMatcher::invalidateActiveArea()
{
	m_activeAreaComputed=false;
}

void ScanMatcher::computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings)
{
	if (m_activeAreaComputed)return;
	OrientedPoint lp=p;
	lp.x+=(cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y);
	lp.y+=(sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y);
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);

	Point min(map.map2world(0,0));
	Point max(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
	     
	if (lp.x<min.x) min.x=lp.x;
	if (lp.y<min.y) min.y=lp.y;
	if (lp.x>max.x) max.x=lp.x;
	if (lp.y>max.y) max.y=lp.y;
	

	const double * angle=m_laserAngles+m_initialBeamsSkip;
	
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{

		if (*r > m_laserMaxRange || *r==0.0 || isnan(*r)) continue;

		double d=*r>m_usableRange?m_usableRange:*r;
		Point phit=lp;
		phit.x+=d*cos(lp.theta+*angle);
		phit.y+=d*sin(lp.theta+*angle);
		
		if (phit.x<min.x) min.x=phit.x;
		if (phit.y<min.y) min.y=phit.y;
		if (phit.x>max.x) max.x=phit.x;
		if (phit.y>max.y) max.y=phit.y;
		
	}
	
	if ( !map.isInside(min)	|| !map.isInside(max))
	{
		Point lmin(map.map2world(0,0));
		Point lmax(map.map2world(map.getMapSizeX()-1,map.getMapSizeY()-1));
		
		min.x=( min.x >= lmin.x )? lmin.x: min.x-m_enlargeStep;
		max.x=( max.x <= lmax.x )? lmax.x: max.x+m_enlargeStep;
		min.y=( min.y >= lmin.y )? lmin.y: min.y-m_enlargeStep;
		max.y=( max.y <= lmax.y )? lmax.y: max.y+m_enlargeStep;
		map.resize(min.x, min.y, max.x, max.y);
	}
	
	HierarchicalArray2D<PointAccumulator>::PointSet activeArea;
	angle=m_laserAngles+m_initialBeamsSkip;

	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		if (m_generateMap)
		{
			double d=*r;
			if (d>m_laserMaxRange || d==0.0 || isnan(d))continue;
			if (d>m_usableRange)d=m_usableRange;
			
			  Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			  IntPoint p0=map.world2map(lp);
			  IntPoint p1=map.world2map(phit);
			
			  GridLineTraversalLine line;
			  line.points=m_linePoints;
			  GridLineTraversal::gridLine(p0, p1, &line);
			  for (int i=0; i<line.num_points-1; i++)
			  {
				  assert(map.isInside(m_linePoints[i]));
				  activeArea.insert(map.storage().patchIndexes(m_linePoints[i]));
				  assert(m_linePoints[i].x>=0 && m_linePoints[i].y>=0);
			  }
			  if (d<m_usableRange)
			  {
				  IntPoint cp=map.storage().patchIndexes(p1);
				  assert(cp.x>=0 && cp.y>=0);
				  activeArea.insert(cp);
			  }
			
		} 
		else 
		{
			   if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;

			   Point phit=lp;
			   phit.x+=*r*cos(lp.theta+*angle);
			   phit.y+=*r*sin(lp.theta+*angle);
			   IntPoint p1=map.world2map(phit);
			   assert(p1.x>=0 && p1.y>=0);
			   IntPoint cp=map.storage().patchIndexes(p1);
			   assert(cp.x>=0 && cp.y>=0);
			   activeArea.insert(cp);
			
		}
	}
	
	map.storage().setActiveArea(activeArea, true);
	m_activeAreaComputed=true;

}

double ScanMatcher::registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings)
{
	if (!m_activeAreaComputed)computeActiveArea(map, p, readings);

	map.storage().allocActiveArea();

	OrientedPoint lp=p;
	lp.x+=(cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y);
	lp.y+=(sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y);
	lp.theta+=m_laserPose.theta;
	IntPoint p0=map.world2map(lp);
		
	const double * angle=m_laserAngles+m_initialBeamsSkip;
	double esum=0;

	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++)
	{
		if (m_generateMap)
		{

			double d=*r;
			if (d>m_laserMaxRange||d==0.0||isnan(d))continue;
			if (d>m_usableRange)d=m_usableRange;			

			   Point phit=lp+Point(d*cos(lp.theta+*angle),d*sin(lp.theta+*angle));
			   IntPoint p1=map.world2map(phit);
			   GridLineTraversalLine line;
			   line.points=m_linePoints;
			   GridLineTraversal::gridLine(p0, p1, &line);
			   for (int i=0; i<line.num_points-1; i++)
			  {
				  PointAccumulator& cell=map.cell(line.points[i]);
				  double e=-cell.entropy();
				  cell.update(false, Point(0,0));
				  e+=cell.entropy();
				  esum+=e;
			  }
			  if (d<m_usableRange)
			  {
				double e=-map.cell(p1).entropy();
				map.cell(p1).update(true, phit);
				e+=map.cell(p1).entropy();
				esum+=e;
			  }
			
		} 
		else 
		{
			if (*r>m_laserMaxRange||*r>m_usableRange||*r==0.0||isnan(*r)) continue;

			Point phit=lp;
			phit.x+=*r*cos(lp.theta+*angle);
			phit.y+=*r*sin(lp.theta+*angle);
			IntPoint p1=map.world2map(phit);
			assert(p1.x>=0 && p1.y>=0);
			map.cell(p1).update(true,phit);
			
		}
	}

	return esum;
}
double ScanMatcher::optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& init, const double* readings) const
{
	double bestScore=-1;
	OrientedPoint currentPose=init;
	double currentScore=score(map, currentPose, readings);//对当前位置进行评价
	double adelta=m_optAngularDelta, ldelta=m_optLinearDelta;
	unsigned int refinement=0;
	enum Move{Front, Back, Left, Right, TurnLeft, TurnRight, Done};

	do{
		if (bestScore>=currentScore)//当前评分小于最优解
		{
			refinement++;//循环次数
			adelta*=.5;
			ldelta*=.5;
		}
		bestScore=currentScore;

		OrientedPoint bestLocalPose=currentPose;
		OrientedPoint localPose=currentPose;

		Move move=Front;
		do 
		{
			localPose=currentPose;
			switch(move)
			{
				case Front:
					localPose.x+=ldelta;
					move=Back;
					break;
				case Back:
					localPose.x-=ldelta;
					move=Left;
					break;
				case Left:
					localPose.y-=ldelta;
					move=Right;
					break;
				case Right:
					localPose.y+=ldelta;
					move=TurnLeft;
					break;
				case TurnLeft:
					localPose.theta+=adelta;
					move=TurnRight;
					break;
				case TurnRight:
					localPose.theta-=adelta;
					move=Done;
					break;
				default:;
			}
			
			double odo_gain=1;
			if (m_angularOdometryReliability>0.)
			{
				double dth=init.theta-localPose.theta; 	
				dth=atan2(sin(dth), cos(dth)); 	
				dth*=dth;
				odo_gain*=exp(-m_angularOdometryReliability*dth);
			}
			if (m_linearOdometryReliability>0.)
			{
				double dx=init.x-localPose.x;
				double dy=init.y-localPose.y;
				double drho=dx*dx+dy*dy;
				odo_gain*=exp(-m_linearOdometryReliability*drho);
			}
			double localScore=odo_gain*score(map, localPose, readings);
			
			if (localScore>currentScore)
			{
				currentScore=localScore;
				bestLocalPose=localPose;
			}
		} while(move!=Done);
		currentPose=bestLocalPose;
	}while (currentScore > bestScore || refinement < m_optRecursiveIterations);

	pnew=currentPose;
	return bestScore;
}

void ScanMatcher::setLaserParameters(unsigned int beams, double* angles, const OrientedPoint& lpose)
{
	assert(beams<LASER_MAXBEAMS);
	m_laserPose=lpose;
	m_laserBeams=beams;
	memcpy(m_laserAngles, angles, sizeof(double)*m_laserBeams);	
}
	
void ScanMatcher::setMatchingParameters	(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations,  double likelihoodSigma, unsigned int likelihoodSkip)
{	
	m_usableRange=urange;
	m_laserMaxRange=range;
	m_kernelSize=kernsize;
	m_optLinearDelta=lopt;
	m_optAngularDelta=aopt;
	m_optRecursiveIterations=iterations;
	m_gaussianSigma=sigma;
	m_likelihoodSigma=likelihoodSigma;
	m_likelihoodSkip=likelihoodSkip;
}

};

