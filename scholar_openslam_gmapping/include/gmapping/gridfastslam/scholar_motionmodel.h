#ifndef MOTIONMODEL_H
#define MOTIONMODEL_H

#include <gmapping/utils/point.h>
#include <gmapping/utils/stat.h>
#include <gmapping/utils/macro_params.h>

namespace  Scholar_GMapping { 

struct MotionModel
{
	OrientedPoint drawFromMotion(const OrientedPoint& p, const OrientedPoint& pnew, const OrientedPoint& pold) const;

	double srr, str, srt, stt;
	double xx,xy,xth,yy,yth,thth;
};

};

#endif
