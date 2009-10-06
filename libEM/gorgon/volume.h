// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju (taoju@cse.wustl.edu)
// Description:   Volumetric data definition

#ifndef VOLUME_H
#define VOLUME_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "priority_queue.h"
#include <vector>
#include "emdata.h"

namespace EMAN {

class Volume
{
private:
	EMData* emdata;
	PriorityQueue<int, int>* que;
public:
	//I commented out the overloaded versions that aren't used by PerformPureJuSkeletonization().	

	Volume(int que_max);

	Volume()
	{	emdata = new EMData;
	}
	Volume(EMData* wrap_this)
	{	emdata = wrap_this;
	}	
	//~Volume();
	EMData* get_emdata()
	{	return emdata;
	}
	void set_emdata(EMData* wrap_this)
	{	emdata = wrap_this;
	}
	void delete_emdata()
	{	delete emdata;
	}
	int getSizeX()
	{	return emdata->get_xsize();
	}
	int getSizeY()
	{	return emdata->get_ysize();
	}
	int getSizeZ()
	{	return emdata->get_zsize();
	}
	//double getDataAt(int index);
	double getDataAt(int x, int y, int z)
	{	return emdata->get_value_at(x,y,z);
	}	
	//void setDataAt(int index, double d);
	void setDataAt(int x, int y, int z, double d )
	{	emdata->set_value_at(x,y,z, static_cast<float>(d));
	}
		
	
	void curveSkeleton(float thr, Volume* svol);
	//void curveSkeleton(Volume* grayvol, float lowthr, float highthr, Volume* svol);
	void curveSkeleton2D(float thr, Volume * svol);
	//void erodeHelix();
	void erodeHelix(int disthr);
	//int erodeSheet();
	int erodeSheet(int disthr);
	void pad(int padBy, double padValue);
	//void skeleton(float thr, int off);
	void skeleton(float thr, Volume* svol, Volume* hvol);
	//void surfaceSkeletonPres(float thr);
	void surfaceSkeletonPres(float thr, Volume* svol);
};

}
#endif
