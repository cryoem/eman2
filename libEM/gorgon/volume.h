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

#define MAX_SHEETS 100000
#define MAX_QUEUELEN 5000000

#define MAX_ERODE 1000

using namespace std;

namespace EMAN {

const int neighbor6[6][3]={{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}} ;
const int neighbor4[4][2]={{0,1},{0,-1},{1,0},{-1,0}} ;
const int neighbor64[6][4][3] = { 
	{{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}},
	{{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}},
	{{0,0,1},{0,0,-1},{1,0,0},{-1,0,0}},
	{{0,0,1},{0,0,-1},{1,0,0},{-1,0,0}},
	{{0,0,1},{0,0,-1},{0,1,0},{0,-1,0}},
	{{0,0,1},{0,0,-1},{0,1,0},{0,-1,0}}} ;
	
const int sheetNeighbor[12][4][3] = {
	{{0,-1,-1},{0,-1,0},{0,0,-1},{0,0,0}},
	{{0,-1,0},{0,-1,1},{0,0,0},{0,0,1}},
	{{0,0,-1},{0,0,0},{0,1,-1},{0,1,0}},
	{{0,0,0},{0,0,1},{0,1,0},{0,1,1}},

	{{-1,0,-1},{-1,0,0},{0,0,-1},{0,0,0}},
	{{-1,0,0},{-1,0,1},{0,0,0},{0,0,1}},
	{{0,0,-1},{0,0,0},{1,0,-1},{1,0,0}},
	{{0,0,0},{0,0,1},{1,0,0},{1,0,1}},

	{{-1,-1,0},{-1,0,0},{0,-1,0},{0,0,0}},
	{{-1,0,0},{-1,1,0},{0,0,0},{0,1,0}},
	{{0,-1,0},{0,0,0},{1,-1,0},{1,0,0}},
	{{0,0,0},{0,1,0},{1,0,0},{1,1,0}}
	};

const int faceCells[12][2]={{0,4},{1,5},{2,6},{3,7},{0,2},{1,3},{4,6},{5,7},{0,1},{2,3},{4,5},{6,7}};

const int cubeFaces[6][4] = 
{ {1,5,7,3},{0,2,6,4},{2,3,7,6},{0,4,5,1},{5,4,6,7},{0,1,3,2}};


//const int faceEdges[12][2] = {{2,4},{2,5},{3,4},{3,5},
//{0,4},{0,5},{1,4},{1,5},
//{0,2},{0,3},{1,2},{1,3}};

const int faceEdges[12][2] = {{3,1},{3,0},{2,1},{2,0},
{5,1},{5,0},{4,1},{4,0},
{5,3},{5,2},{4,3},{4,2}};

const int edgeFaces[6][4] = {{1,3,5,7},{0,2,4,6},{2,3,9,11},{0,1,8,10},{6,7,10,11},{4,5,8,9}} ;

struct gridPoint
{
	int x, y, z;
};

struct gridQueueEle
{
	int x, y, z;
	int score ;
	gridQueueEle* next ;
};

class gridQueue2;



//**********************************************************************
// Volume
//**********************************************************************
class Volume
{
private:
	/*EMData has a float array named data, and so does this class.
	 *However, the indexing works differently.
	 * Volume: index = x*sizey*sizez + y*sizez + z;
	 * EMData: index = x + y*nx + z*nx*ny; //nx means the same thing as sizex
	 * Because of this,  we'll define new coordinates (x', y', z') for
	 * this Volume class such that
	 * (x', y', z') = (z, y, x)
	 * where (x,y,z) is the EMData coordinates. We can also think of
	 * this as moving the density value every point from (x, y, z) to (z, y, x).
	 * Thus, we are reflecing our density map over the plane x = z.
	 */
	EMData* emdata;

	/* Sizes */
	int sizex, sizey, sizez ;
	//float spacingX, spacingY, spacingZ;
	//float originX, originY, originZ;

	/* Data array */
	float * data ;
public:
	//If the Volume object created an EMData object, it owns that object
	//so it will be deleted when the Volume object is deleted
	bool owns_emdata; 
	
	//I commented out the overloaded versions that aren't used by PerformPureJuSkeletonization().	

	Volume(EMData* wrap_this)
	{
		set_emdata(wrap_this);
	}

	Volume(int nx, int ny, int nz)
	{
		emdata = new EMData(nx,ny,nz);
		owns_emdata = true;
		sizex = nz;
		sizey = ny;
		sizez = nx;
		data = emdata->get_data();
	}

	~Volume()
	{
		if (owns_emdata)
			delete_emdata();
	}

	void update_emdata_dimensions()
	{
		emdata->set_size(sizez, sizey, sizex);
	}

	EMData* get_emdata()
	{	return emdata;
	}
	void set_emdata(EMData* wrap_this)
	{	emdata = wrap_this;
		owns_emdata = false;
		sizex = emdata->get_zsize();
		sizey = emdata->get_ysize();
		sizez = emdata->get_xsize();
		data = emdata->get_data();//Returns a pointer to a 1D float array
		
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
	void setDataAt(int index, double d)
	{	emdata->set_value_at(index, float(d));
	}
	void setDataAt(int x, int y, int z, double d )
	{	emdata->set_value_at(x,y,z, static_cast<float>(d));
	}

	void pad(int padBy, double padValue);



	void curveSkeleton(float thr, Volume* svol);
	//void curveSkeleton(Volume* grayvol, float lowthr, float highthr, Volume* svol);
	void curveSkeleton2D(float thr, Volume * svol);
	//void erodeHelix();
	void erodeHelix(int disthr);
	//int erodeSheet();
	int erodeSheet(int disthr);
	void skeleton(float thr, int off);
	void skeleton(float thr, Volume* svol, Volume* hvol);
	//void surfaceSkeletonPres(float thr);
	void surfaceSkeletonPres(float thr, Volume* svol);



	int countIntEuler( int ox, int oy, int oz );
	int components6( int vox[3][3][3] );
	int components26( int vox[3][3][3] );
	int countInt( double vox[3][3][3] );
	int countExt( double vox[3][3][3] );
	int hasCell( int ox, int oy, int oz );
	int hasCompleteSheet( int ox, int oy, int oz );
	int getNumNeighbor6( int ox, int oy, int oz );
	int getNumPotComplex( int ox, int oy, int oz );
	int getNumPotComplex2( int ox, int oy, int oz );
	int hasCompleteHelix( int ox, int oy, int oz );
	//int hasCompleteHelix( int ox, int oy, int oz, Volume* fvol );
	int isFeatureFace( int ox, int oy, int oz );
	int isHelixEnd( int ox, int oy, int oz, Volume* nvol );
	int isHelixEnd( int ox, int oy, int oz );
	int isSheetEnd( int ox, int oy, int oz );
	//int isSheetEnd( int ox, int oy, int oz, Volume* nvol );
	int isSimple( int ox, int oy, int oz );
	Volume* markCellFace( );
	/**
	 * Normalize to a given range 
	 */
	void threshold( double thr, int out, int in, int boundary, bool markBoundary);
	void threshold( double thr )
	{	threshold( thr, 0, 1, 0, true) ;
	}
	void threshold( double thr, int out, int in )
	{	threshold( thr, out, in, out, true) ;
	}
	void threshold( double thr, int out, int in, int boundary)
	{	threshold(thr, out, in, boundary, true);
	}
};

//**********************************************************************
// GridQueue2
//**********************************************************************
class GridQueue2
{
	gridQueueEle* head ;
	gridQueueEle* pre ;
	gridQueueEle* prepre ;
	gridQueueEle* cur ;
	int numEles ;

public:

	GridQueue2( )
	{
		head = NULL ;
		cur = NULL ;
		pre = NULL ;
		prepre = NULL ;
		numEles = 0 ;
	}

	~GridQueue2()
	{
		gridQueueEle* ele;
		reset();
		ele=getNext();
		while ( (ele=remove()) != NULL ){};
	}
	gridQueueEle* getNext( )
	{
		if ( cur == NULL )
		{
			prepre = NULL ;
			pre = NULL ;
			cur = head ;
		}
		else
		{
			prepre = pre ;
			pre = cur ;
			cur = cur->next ;
		}

		return cur ;
	}

	void reset( )
	{
		prepre = NULL ;
		pre = NULL ;
		cur = NULL ;
	}

	int getNumElements( )
	{
		return numEles ;
	}

	void prepend( int xx, int yy, int zz )
	{
		gridQueueEle* ele = new gridQueueEle ;
		ele->x = xx ;
		ele->y = yy ;
		ele->z = zz ;
		ele->score = 0 ;
		ele->next = head;
		head = ele ;
		numEles ++ ;

		reset() ;
	}

	/* Remove current element pointed by cur */
	gridQueueEle* remove( )
	{
		gridQueueEle* temp = cur ;
		if ( cur != NULL )
		{
			cur = cur->next ;
			delete temp ;

			if ( pre != NULL )
			{
				pre->next = cur ;
			}
			else
			{
				head = cur ;
			}
			numEles -- ;
		}


		return cur ;
	}

	/* Switching pre and cur */
	gridQueueEle* swap( )
	{
		if ( prepre != NULL )
		{
			pre->next = cur->next ;
			cur->next = prepre->next ;
			prepre->next = cur ;

		}
		else
		{
			pre->next = cur->next ;
			cur->next = pre ;
			head = cur ;
		}
		
		gridQueueEle* temp = pre ;
		pre = cur ;
		cur = temp ;

		return cur ;
	}
};

}
#endif
