// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju (taoju@cse.wustl.edu), Refactored by Sasakthi Abeysinghe (sasakthi.abeysinghe@wustl.edu)
// Description:   Volumetric data definition

#include "volume_data.h"
//#include "ThinningTemplate.h"
#include "grid_queue.h"
#include "grid_queue2.h"
//#include <cstdio>
//#include <cstdlib>
//#include <cmath>
#include "priority_queue.h"
//#include <vector>

#ifndef SKELETON_MAKER_VOLUME_H
#define SKELETON_MAKER_VOLUME_H

#define MAX_SHEETS 100000
#define MAX_QUEUELEN 5000000
#define MAX_ERODE 1000

using namespace std;

namespace wustl_mm {
	namespace SkeletonMaker {

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

		const int faceEdges[12][2] = {{3,1},{3,0},{2,1},{2,0},
									  {5,1},{5,0},{4,1},{4,0},
									  {5,3},{5,2},{4,3},{4,2}};

		const int edgeFaces[6][4] = {{1,3,5,7},{0,2,4,6},{2,3,9,11},{0,1,8,10},{6,7,10,11},{4,5,8,9}} ;

		struct gridPoint
		{
			int x, y, z;
		};

		class Volume {
		public:
			Volume(EMData* em);//eman2
			Volume(int x, int y, int z);
			Volume(int x, int y, int z, float val);
			Volume(int x, int y, int z, int offx, int offy, int offz, Volume * vol);
			~Volume( );

			EMData* get_emdata(); //eman2
			float getSpacingX();
			float getSpacingY();
			float getSpacingZ();
			float getOriginX();
			float getOriginY();
			float getOriginZ();
			int getSizeX();
			int getSizeY();
			int getSizeZ();
			int getIndex(int x, int y, int z);
			double getDataAt( int x, int y, int z );
			double getDataAt( int index );
			void setSpacing(float spx, float spy, float spz );
			void setOrigin(float orgX, float orgY, float orgZ);
			void setDataAt( int x, int y, int z, double d );
			void setDataAt( int index, double d );


			//Volume * getPseudoDensity();
			//Volume * getDistanceField(int rad, float randf);
			//int getNonZeroVoxelCount();
			//void print();
			//void subtract(Volume * vol);
			void pad (int padBy, double padValue);
			//void applyMask(Volume * maskVol, double maskValue, bool keepMaskValue);
			//double getMin();
			//double getMax();
			//double getMaxValuePosition(int& maxX, int& maxY, int& maxZ);
			//double getLocalMax(int x, int y, int z, int radius);
			//double getLocalMin(int x, int y, int z, int radius);
			//void fill(double val);
			//int isBertrandBorder(int ox, int oy, int oz, int dir);
			//int isBertrandEndPoint(int ox, int oy, int oz);
			//int isHelix(int ox, int oy, int oz);
			//int isSheet(int ox, int oy, int oz);
			//Volume * getSheets(int minSize);
			//Volume * getHelices(int minSize);
			//int isEndPoint(int ox, int oy, int oz);
			int getNumNeighbor6(int ox, int oy, int oz);
			//int testIsSheetEnd(int ox, int oy, int oz);
			//int isNoiseSheetEnd(int ox, int oy, int oz);
			//int isInternal(int ox, int oy, int oz);
			//int isInternal2(int ox, int oy, int oz);
			//int hasIsolatedFace(int ox, int oy, int oz);
			//int hasIsolatedEdge(int ox, int oy, int oz);
			//int countFace(int ox, int oy, int oz, int m);
			int hasCell(int ox, int oy, int oz);
			Volume * markCellFace();
			//Volume * markFaceEdge();
			int hasCompleteSheet(int ox, int oy, int oz, Volume * fvol);
			int hasCompleteSheet(int ox, int oy, int oz);
			//int hasCompleteSheetSlow(int ox, int oy, int oz);
			int hasCompleteHelix(int ox, int oy, int oz);
			int hasCompleteHelix(int ox, int oy, int oz, Volume * fvol);
			int isHelixEnd(int ox, int oy, int oz, Volume * nvol);
			//int isFeature18(int ox, int oy, int oz);
			//int isEdgeEnd(int ox, int oy, int oz);
			//int isFaceEnd(int ox, int oy, int oz);
			//int isNoise(int ox, int oy, int oz, int noise);
			//int isNoiseHelixEnd(int ox, int oy, int oz);
			int isHelixEnd(int ox, int oy, int oz);
			int isSheetEnd(int ox, int oy, int oz, Volume * nvol);
			//int getNumFaces( int ox, int oy, int oz );
			//int getNumCells( int ox, int oy, int oz );
			//int getNumIsolatedEdges( int ox, int oy, int oz );
			//int getNumIsolatedFaces( int ox, int oy, int oz );
			//int isFeatureFace2( int ox, int oy, int oz );
			int isFeatureFace( int ox, int oy, int oz );
			//int hasFeatureFace( int ox, int oy, int oz );
			int isSheetEnd( int ox, int oy, int oz );
			int isSimple( int ox, int oy, int oz );
			int isPiercable( int ox, int oy, int oz );
			//int isSimple2( int v[3][3][3] );
			//int getNumPotComplex3( int ox, int oy, int oz );
			//int getNumPotComplex4( int ox, int oy, int oz );
			int getNumPotComplex( int ox, int oy, int oz );
			int getNumPotComplex2( int ox, int oy, int oz );
			//int getNumNeighbor( int ox, int oy, int oz );
			//void setScoreNeighbor( GridQueue* queue );
			int components6( int vox[3][3][3] );
			int components26( int vox[3][3][3] );
			int countExt( double vox[3][3][3] );
			int countInt( double vox[3][3][3] );
			int countIntEuler( int ox, int oy, int oz );
			//void erodeNoTopo( float thr, int wid );
			//void erodeTopo( float thr, int wid );
			//void erode2( float thr, int wid );
			//void erodeShapeTopo( float thr, int wid );
			//void erodeAtom( float thr, int wid, Volume* avol );
			void curveSkeleton( Volume* grayvol, float lowthr, float highthr, Volume* svol );
			void curveSkeleton( float thr, Volume* svol );
			void curveSkeleton2D( float thr, Volume* svol );
			void skeleton( float thr, int off );
			//void skeleton2( float thr, int off );
			//void pointSkeleton( Volume* grayvol, float lowthr, float highthr, Volume* svol, Volume* hvol );
			void skeleton( float thr, Volume* svol, Volume* hvol );
			void erodeHelix( );
			void erodeHelix( int disthr );
			int erodeSheet( );
			int erodeSheet( int disthr );
			//void erodeSheetOld( int disthr );
			//void addNoise( float thr, float pos );
			//void sequentialSkeleton( float thr, int type, int noise );
			//void dumbsurfaceSkeleton( float thr );
			//void surfaceSkeleton( Volume* grayvol, float lowthr, float highthr );
			//void surfaceSkeleton( float thr );
			//void surfaceSkeleton( float thr, Volume* svol );
			//void surfaceSkeletonOld( float thr );
			void surfaceSkeletonPres( float thr, Volume * preserve );
			//void bertrandSurfaceSkeleton2( float thr );
			//void bertrandSurfaceSkeleton( float thr );
			//void palagyiSurfaceSkeleton( float thr );
			void threshold( double thr );
			void threshold( double thr, int out, int in );
			void threshold( double thr, int out, int in, int boundary);
			void threshold( double thr, int out, int in, int boundary, bool markBoundary);
			//void threshold2( double thr, int out, int in );
			//void smooth( float alpha );
			//void normalize( double min, double max );
			//void normalize( double min, double max, double thresh, double ithresh );
			//Volume * getDataRange(int x, int y, int z, int radius);
			//double getInterpDataAt( double x, double y, double z );
			//void rotateX ( double a );
			//void toMathematicaFile( char* fname );
			//void toMathematicaFile( char* fname, int lx, int hx, int ly, int hy, int lz, int hz );
			//void toOFFCells( char* fname );
			//void toOFFCells2( char* fname );
			//void toOFFCells2( char* fname, float thr );
			//void toOFFCells( char* fname, float thr );
			//void segment( float threshold, Volume* lowvol, Volume* highvol, char* mrcfile );
			//void segment( float threshold, Volume* vol, int maxDis, char* mrcfile );
			//void writeSegmentation( float threshold, Volume* segvol, char* txtfile, char* mrcfile );
			//void floodFill( float thr );
			//void reduceComponent( int size );
			//void reduceComponent2( int num );
			//void floodFillPQR( int offset );
			//void writeDistances( char* fname, int maxDis );
			//void toPQRFile( char* fname, float spc, float minx, float miny, float minz, int padding );
			//void toMRCFile( char* fname );
		//private:
			VolumeData * getVolumeData();

		private:
			VolumeData * volData;
		};



	}
}

#endif
