// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju (taoju@cse.wustl.edu), Refactored by Sasakthi Abeysinghe (sasakthi.abeysinghe@wustl.edu)
// Description:   Volumetric data definition

#include "volume.h"

using namespace wustl_mm::SkeletonMaker;

		Volume::Volume(EMData* em) //eman2
		{
			this->volData = new VolumeData(em);
		}
		Volume::Volume(int x, int y, int z) {
			volData = new VolumeData(x, y, z);
		}

		Volume::Volume(int x, int y, int z, float val) {
			volData = new VolumeData(x, y, z, val);
		}

		Volume::Volume(int x, int y, int z, int offx, int offy, int offz, Volume * vol) {
			volData = new VolumeData(x, y, z, offx, offy, offz, vol->getVolumeData());
		}


		Volume::~Volume( )
		{
			delete volData;
		}

		EMData* Volume::get_emdata() //eman2
		{
			return this->getVolumeData()->get_emdata();
		}

		int Volume::getSizeX() {
			return volData->GetSizeX();
		}

		int Volume::getSizeY() {
			return volData->GetSizeY();
		}

		int Volume::getSizeZ() {
			return volData->GetSizeZ();
		}

		int Volume::getIndex(int x, int y, int z) {
			return volData->GetIndex(x, y, z);
		}

		void Volume::setDataAt( int x, int y, int z, double d ) {
			volData->SetDataAt(x, y, z, (float)d);
		}

		void Volume::setDataAt( int index, double d ) {
			volData->SetDataAt(index, (float)d);
		}

		double Volume::getDataAt( int x, int y, int z )
		{
			return volData->GetDataAt(x, y, z);
		}

		double Volume::getDataAt( int index ) {
			return volData->GetDataAt(index);
		}

		VolumeData * Volume::getVolumeData() {
			return volData;
		}

		void Volume::setSpacing(float spx, float spy, float spz ) {
			volData->SetSpacing(spx, spy, spz);
		}

		void Volume::setOrigin(float orgX, float orgY, float orgZ) {
			volData->SetOrigin(orgX, orgY, orgZ);
		}

		float Volume::getSpacingX() {
			return volData->GetSpacingX();
		}

		float Volume::getSpacingY() {
			return volData->GetSpacingY();
		}

		float Volume::getSpacingZ() {
			return volData->GetSpacingZ();
		}

		float Volume::getOriginX() {
			return volData->GetOriginX();
		}

		float Volume::getOriginY() {
			return volData->GetOriginY();
		}

		float Volume::getOriginZ() {
			return volData->GetOriginZ();
		}

		//Volume * Volume::getPseudoDensity( ) {
			//// This function assumes the volume is binary (1/0), and builds a pseudo density volume from the 1 voxels
			//// First: assign a density value at each point
			//int i, j, k ;
			//Volume * res = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this);
			//int size = getSizeX() * getSizeY() * getSizeZ() ;
			//srand(123) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//if ( res->getDataAt( i, j, k ) > 0 ) {
							//int ct = 0 ;
							//for ( int m = 0 ; m < 6 ; m ++ ) {
								//if ( res->getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) > 0 ) {
									//ct ++ ;
								//}
							//}
							//res->setDataAt( i,j,k, (k/(float)getSizeZ())*(k/(float)getSizeZ()) ) ;
							//if ( ct > 2 ) {
								////res->setDataAt( i,j,k, rand() / (double) RAND_MAX / 2.0f ) ;
							//} else {
								////res->setDataAt( i,j,k, rand() / (double) RAND_MAX ) ;
							//}
						//}
					//}

			///* Next, smooth
			//for ( i = 0 ; i < 20 ; i ++ )
			//{
				//printf("Smoothing round %d\n", i) ;
				//res->smooth( 0.5f ) ;
			//}
			//*/

			//Volume * tvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, res);
			//float d, ad, ct, temp;
			//for ( int it = 0 ; it < 3 ; it ++ )
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//if ( (d = (float)tvol->getDataAt( i, j, k )) > 0 ) {
							//ad = 0 ; ct = 0 ;
							//for ( int m = 0 ; m < 6 ; m ++ ) {
								//if ( (temp = (float)tvol->getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] )) > 0 ) {
									//ad += temp;
									//ct ++ ;
								//}
							//}
							//if ( ct > 0 ) {
								//res->setDataAt( i, j, k, ( d + ad/ct ) / 2 ) ;
							//}
						//}
					//}

			//delete tvol;
			//tvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, res ) ;
			//for ( i = 0 ; i < 40 ; i ++ )
			//{
				//printf("Smoothing round %d\n", i) ;
				//res->smooth( 0.5f ) ;
				//continue ;

				////for ( j = 0 ; j < size ; j ++ )
				////{
				////	if ( tvol->getDataAt( j ) > 0 )
				////	{
				////		res->setDataAt( j, tvol->getDataAt( j ) ) ;
				////	}
				////
				////}

			//}


			//return res ;
		//}


		//Volume * Volume::getDistanceField(int rad, float randf) {
			//// This function assumes the volume is binary (1/0), and builds a pseudo density volume from the 1 voxels
			//// rad is the radius of each distance function (e.g., center pixel gets 1, pixels at rad from the center gets 0)
			//// randf is how much noise you want to add. this means the center pixel will maximally have value 1+randf.

			//// First: assign a density value at each point
			//int i, j, k ;
			//Volume * res = new Volume(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this);
			//srand( 123 ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//if ( getDataAt(i, j, k) > 0 ) {
							//float mag = 1 + randf * (float) rand() / (float) RAND_MAX ;
							//int lx = max(0,i-rad) ;
							//int ly = max(0,j-rad) ;
							//int lz = max(0,k-rad) ;
							//int hx = min(getSizeX()-1,i+rad) ;
							//int hy = min(getSizeY()-1,j+rad) ;
							//int hz = min(getSizeZ()-1,k+rad) ;
							//int x,y,z;
							//for ( x = lx ; x <= hx ; x ++ )
								//for ( y = ly ; y <= hy ; y ++ )
									//for ( z = lz ; z <= hz ; z ++ ) {
										//float val = 1 - (float) sqrt((double)((x-i)*(x-i) + (y-j)*(y-j) + (z-k)*(z-k))) / (float) rad ;
										//val *= mag ;
										//if ( res->getDataAt( x, y, z ) < val ) {
											//res->setDataAt( x, y, z, val ) ;
										//}
									//}
						//}
					//}

			///* Next, smooth */
			//for ( i = 0 ; i < 2 ; i ++ )
			//{
				//printf("Smoothing round %d\n", i) ;
				//res->smooth( 0.5f ) ;
			//}


			//return res ;
		//}


		//int Volume::getNonZeroVoxelCount() {
			//int count = 0;
			//for(int x = 0; x < getSizeX(); x++){
				//for(int y = 0; y < getSizeY(); y++){
					//for(int z = 0; z < getSizeZ(); z++){
						//if(this->getDataAt(x, y, z) > 0.0) {
							//count++;
						//}
					//}
				//}
			//}
			//return count;
		//}
		//void Volume::print() {
			//for(int x = 0; x < getSizeX(); x++) {
				//printf("{ ");
				//for(int y = 0; y < getSizeY(); y++) {
					//printf("{ ");
					//for(int z = 0; z < getSizeZ(); z++) {
						//printf("%f, ", getDataAt(x, y, z));
					//}
					//printf("} ");
				//}
				//printf("} ");
			//}
			//printf("\n");
		//}


		//void Volume::subtract( Volume* vol ) {
			//int i, j, k ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//if ( getDataAt( i, j, k ) > 0 ) {
							//if ( vol->getDataAt(i,j,k) > 0 ) {
								//setDataAt( i, j, k, 0 ) ;
							//}
						//}
					//}

		//}

		void Volume::pad (int padBy, double padValue) {
			volData->Pad(padBy, padValue);
		}

		//void Volume::applyMask(Volume * maskVol, double maskValue, bool keepMaskValue) {
			//for(int x = 0; x < maskVol->getSizeX(); x++) {
				//for(int y = 0; y < maskVol->getSizeY(); y++) {
					//for(int z = 0; z < maskVol->getSizeZ(); z++) {
						//if(((maskVol->getDataAt(x, y, z) == maskValue) && !keepMaskValue) ||
							//((maskVol->getDataAt(x, y, z) != maskValue) && keepMaskValue)) {
							//setDataAt(x, y, z, 0);
						//}
					//}
				//}
			//}
		//}

		//double Volume::getMin() {
			//int size = volData->GetMaxIndex();
			//double rvalue = volData->GetDataAt(0);
			//for (int i=1; i < size; i++) {
				//float val = volData->GetDataAt(i);
				//if ( rvalue > val) {
					//rvalue = val;
				//}
			//}
			//return rvalue;
		//}

		//double Volume::getMax() {
			//int size = volData->GetMaxIndex();
			//double rvalue = volData->GetDataAt(0);
			//for (int i=1; i < size; i++) {
				//float val = volData->GetDataAt(i);
				//if ( rvalue < val) {
					//rvalue = val;
				//}
			//}
			//return rvalue ;
		//}

		//double Volume::getMaxValuePosition(int& maxX, int& maxY, int& maxZ) {
			//double maxVal = getDataAt(0,0,0);
			//maxX = 0; maxY = 0; maxZ = 0;
			//double data;

			//for(int x = 0; x < getSizeX(); x++) {
				//for(int y = 0; y < getSizeY(); y++) {
					//for(int z = 0; z < getSizeZ(); z++) {
						//data = getDataAt(x, y, z);
						//if(data > maxVal) {
							//maxVal = data;
							//maxX = x; maxY = y; maxZ = z;
						//}
					//}
				//}
			//}
			//return maxVal;
		//}

		//double Volume::getLocalMax(int x, int y, int z, int radius) {
			//double mx = getDataAt(x, y, z);
			//for(int xx = x - radius; xx <= x + radius; xx++) {
				//for(int yy = y - radius; yy <= y + radius; yy++) {
					//for(int zz = z - radius; zz <= z + radius; zz++) {
						//mx = max(mx, getDataAt(xx, yy, zz));
					//}
				//}
			//}
			//return mx;
		//}

		//double Volume::getLocalMin(int x, int y, int z, int radius) {
			//double mn = getDataAt(x, y, z);
			//for(int xx = x - radius; xx <= x + radius; xx++) {
				//for(int yy = y - radius; yy <= y + radius; yy++) {
					//for(int zz = z - radius; zz <= z + radius; zz++) {
						//mn = min(mn, getDataAt(xx, yy, zz));
					//}
				//}
			//}
			//return mn;
		//}

		//void Volume::fill( double val )
		//{
			//for(int x = 0; x < getSizeX(); x++) {
				//for(int y = 0; y < getSizeY(); y++) {
					//for(int z = 0; z < getSizeZ(); z++) {
						//setDataAt(x, y, z, val);
					//}
				//}
			//}
		//}

		//int Volume::isBertrandBorder( int ox, int oy, int oz, int dir ) {
			//int nx = ox + neighbor6[dir][0] ;
			//int ny = oy + neighbor6[dir][1] ;
			//int nz = oz + neighbor6[dir][2] ;
			//if ( getDataAt( nx, ny, nz) < 0 ) {
				//return 1 ;
			//}
			//return 0 ;
		//}

		//int Volume::isBertrandEndPoint( int ox, int oy, int oz ) {
			//double vox[3][3][3] ;

			//int i, j, k ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ ) {
						//double tval = getDataAt( ox + i, oy + j, oz + k ) ;
							//vox[ i + 1 ][ j + 1 ][ k + 1 ] = tval ;
					//}

			//// Count #X^x_6
			//int xx6 = 0 ;
			//for ( i = 0 ; i < 6 ; i ++ ) {
				//if ( vox[ neighbor6[i][0] + 1 ][ neighbor6[i][1] + 1 ][ neighbor6[i][2] + 1 ] >= 0 ) {
					//xx6 ++ ;
				//}
			//}

			//// Now, count A26 and /A26
			//int a26 = 0, na26 = 0 ;
			//for ( i = 0 ; i < 3 ; i += 2 )
				//for ( j = 0 ; j < 3 ; j += 2 )
					//for ( k = 0 ; k < 3 ; k += 2 ) {
						//if ( vox[1][j][k] < 0 ) {
							//continue ;
						//}
						//if ( vox[i][1][k] < 0 ) {
							//continue ;
						//}
						//if ( vox[i][j][1] < 0 ) {
							//continue ;
						//}
						//if ( vox[i][1][1] < 0 ) {
							//continue ;
						//}
						//if ( vox[1][j][1] < 0 ) {
							//continue ;
						//}
						//if ( vox[1][1][k] < 0 ) {
							//continue ;
						//}
						//if ( vox[i][j][k] >= 0 ) {
							//a26 ++ ;
						//} else {
							//na26 ++ ;
						//}
					//}

			//if (( na26 == 0 ) && ( a26 != 0 ) && ( xx6 <= a26 + 2 ))  {
				//return 0 ;
			//}
			//if ( isFeatureFace(ox,oy,oz) ) {
				////return 0 ;
			//}
			//return 1 ;
		//}

		//int Volume::isHelix( int ox, int oy, int oz ) {
			//int cn = 12 ;
			//int nx, ny, nz ;
			//int i, j, k ;

			//double vox[3][3][3] ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ ) {
						//vox[i+1][j+1][k+1] = getDataAt( ox + i, oy + j, oz + k ) ;
					//}

			//for ( i = 0 ; i < 12 ; i ++ ) {
				//for ( j = 0 ; j < 4 ; j ++ ) {
					//nx = sheetNeighbor[i][j][0] + 1;
					//ny = sheetNeighbor[i][j][1] + 1;
					//nz = sheetNeighbor[i][j][2] + 1;

					//if ( vox[nx][ny][nz] <= 0 ) {
						//cn -- ;
						//break ;
					//}
				//}
			//}

			//if ( cn >= 1 ) {
				//return 0 ;
			//} else {
				//return 1 ;
			//}
		//}

		//int Volume::isSheet( int ox, int oy, int oz ) {
			//int cn = 12 ;
			//int nx, ny, nz ;

			//for ( int i = 0 ; i < 12 ; i ++ ) {
				//for ( int j = 0 ; j < 4 ; j ++ ) {
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) <= 0 ) {
						//cn -- ;
						//break ;
					//}
				//}
			//}
			//return ( cn >= 3 ) ;
		//}

		//Volume * Volume::getSheets( int minSize ) {
			//int i, j, k ;

			////Initialize volume
			//printf("Initialize volume at %d %d %d\n", getSizeX(), getSizeY(), getSizeZ() ) ;
			//Volume* svol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;

			////Initialize cluster counters
			//int sheets[MAX_SHEETS] ;
			//for ( i = 0 ; i < MAX_SHEETS ; i ++ ) {
				//sheets[ i ] = 0 ;
			//}
			//int totSheets = 1 ;

			////Start clustering
			//printf("Start clustering...\n" ) ;
			//int ox, oy, oz ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//if ( getDataAt(i,j,k) <= 0 || svol->getDataAt(i,j,k) != 0 ) {
							//// Not a data point or has been visited
							//continue ;
						//}
						//if ( ! isSheet( i, j, k ) ) {
							//// Not a sheet point
							//continue ;
						//}

						////Initialize queue
						//int numNodes = 1 ;
						//svol->setDataAt( i, j, k, totSheets ) ;
						//GridQueue* queue = new GridQueue() ;
						//queue->pushQueue( i, j, k ) ;
						//while ( queue->popQueue(ox, oy, oz) ) {
							//// Test if neighbors satisfy sheet condition
							//if ( isSheet( ox, oy, oz ) ) {
								//for ( int m = 0 ; m < 6 ; m ++ ) {
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;

									//if ( getDataAt(nx,ny,nz) > 0 && svol->getDataAt(nx,ny,nz) == 0 ) {
										//svol->setDataAt(nx,ny,nz,totSheets);
										//queue->pushQueue(nx,ny,nz) ;
										//numNodes ++ ;
									//}
								//}
							//}
						//}

						//delete queue ;
						//if ( numNodes > 0 ) {
						////	printf("Sheet %d contain %d nodes.\n", totSheets, numNodes) ;
							//sheets[ totSheets ] = numNodes ;
							//totSheets ++ ;
						//}
					//}

			//// Removing clusters less than minSize
			//printf("Removing small clusters.\n") ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//int cnt = (int) svol->getDataAt(i,j,k) ;
						//if ( cnt > 0 && sheets[ cnt ] < minSize ) {
							//svol->setDataAt(i,j,k,-1) ;
						//}
					//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//svol->threshold( 0.1, 0, 1 ) ;

			//return svol ;
		//}

		//Volume * Volume::getHelices( int minSize ) {
			//printf("Segmenting helices from eroded volume.\n") ;
			//int i, j, k ;

			////Initialize volume
			//printf("Initialize volume at %d %d %d\n", getSizeX(), getSizeY(), getSizeZ() ) ;
			//Volume* svol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;

			////Initialize cluster counters
			//int helices[MAX_SHEETS] ;
			//for ( i = 0 ; i < MAX_SHEETS ; i ++ ) {
				//helices[ i ] = 0 ;
			//}
			//int totHelices = 1 ;

			////Start clustering
			//printf("Start clustering...\n" ) ;
			//int ox, oy, oz ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//if ( getDataAt(i,j,k) <= 0 || svol->getDataAt(i,j,k) != 0 ) {
							//// Not a data point or has been visited
							//continue ;
						//}
						//if ( ! isHelix( i, j, k ) ) {
							//// Not a helix point
							//continue ;
						//}

						////Initialize queue
						//int numNodes = 1 ;
						//svol->setDataAt( i, j, k, totHelices ) ;
						//GridQueue* queue = new GridQueue() ;
						//queue->pushQueue( i, j, k ) ;
						//while ( queue->popQueue(ox, oy, oz) )
						//{
							//// Test if neighbors satisfy helix condition
							//if ( isHelix( ox, oy, oz ) ) {
								//for ( int m = 0 ; m < 6 ; m ++ ) {
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;

									//if ( getDataAt(nx,ny,nz) > 0 && svol->getDataAt(nx,ny,nz) == 0 ) {
										//svol->setDataAt(nx,ny,nz,totHelices);
										//queue->pushQueue(nx,ny,nz) ;
										//numNodes ++ ;
									//}
								//}
							//}
						//}

						//delete queue ;
						//if ( numNodes > 0 ) {
						////	printf("Helix %d contain %d nodes.\n", totHelices, numNodes) ;
							//helices[ totHelices ] = numNodes ;
							//totHelices ++ ;
						//}
					//}

			//// Removing clusters less than minSize
			//printf("Removing small clusters.\n") ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ ) {
						//int cnt = (int) svol->getDataAt(i,j,k) ;
						//if ( cnt > 0 && helices[ cnt ] < minSize ) {
							//svol->setDataAt(i,j,k,-1) ;
						//}
					//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//svol->threshold( 0.1, 0, 1 ) ;

			//return svol ;

		//}

		//int Volume::isEndPoint( int ox, int oy, int oz ) {
			//if ( getDataAt( ox - 1, oy, oz ) < 0 && getDataAt( ox + 1, oy, oz ) < 0 ) {
				//return 1 ;
			//}
			//if ( getDataAt( ox, oy - 1, oz ) < 0 && getDataAt( ox, oy + 1, oz ) < 0 ) {
				//return 1 ;
			//}
			//if ( getDataAt( ox, oy, oz - 1 ) < 0 && getDataAt( ox, oy, oz + 1 ) < 0 ) {
				//return 1 ;
			//}
			//return 0 ;
		//}

		int Volume::getNumNeighbor6( int ox, int oy, int oz ) {
			int rvalue = 0 ;
			for ( int i = 0 ; i < 6 ; i ++ ) {
				int nx = ox + neighbor6[i][0] ;
				int ny = oy + neighbor6[i][1] ;
				int nz = oz + neighbor6[i][2] ;
				if ( getDataAt( nx, ny, nz ) >= 0 ) {
					rvalue ++ ;
				}
			}

			return rvalue ;
		}
		//int Volume::testIsSheetEnd( int ox, int oy, int oz ) {
			//// Returns 1 if it lies on the edge of a non-noise sheet
			//int i, j ;
			//int nx, ny, nz ;

			//int edge[6] = { 0,0,0,0,0,0 } ;
			//int faceflag[ 12 ] ;
			//int hasNoiseFace = 0 ;
			//int tot = 0 ;

			//for ( i = 0 ; i < 12 ; i ++ ) {
				//faceflag[ i ] = 1 ;
				//int hasNoise = 0 ;

				//for ( j = 0 ; j < 4 ; j ++ ) {
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) == 0 ) {
						//hasNoise = 1 ;
					//}
					//if ( getDataAt( nx, ny, nz ) < 0 ) {
						//faceflag[ i ] = 0 ;
						//break ;
					//}
				//}
				//if ( faceflag[ i ] == 1 && hasNoise ) {
					//hasNoiseFace ++ ;
					//// return 0 ;
				//}

				//if ( faceflag[ i ] ) {
					//int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					//edge[ e0 ] ++ ;
					//edge[ e1 ] ++ ;
					//tot ++ ;
				//}
			//}

			//if ( hasNoiseFace == tot ) {
				//return 0 ;
			//}

			//if ( tot == 0 ) {
				//return 0 ;
			//}

			//// Removing 1s
			//int numones = 0 ;
			//for ( i = 0 ; i < 6 ; i ++ ) {
				//if ( edge[ i ] == 1 ) {
					//numones ++ ;
				//}
			//}
			//while( numones > 0 ) {
				//int e ;
				//for ( i = 0 ; i < 6 ; i ++ ) {
					//if ( edge[ i ] == 1 ) {
						//e = i ;
						//break ;
					//}
				//}

				//int f, e2 ;
				//for ( j = 0 ; j < 4 ; j ++ ) {
					//f = edgeFaces[ e ][ j ] ;
					//if ( faceflag[ f ] ) {
						//break ;
					//}
				//}

				//if ( faceEdges[ f ][ 0 ] == e ) {
					//e2 = faceEdges[ f ][ 1 ] ;
				//}  else {
					//e2 = faceEdges[ f ][ 0 ] ;
				//}

				//edge[ e ] -- ;
				//numones -- ;
				//edge[ e2 ] -- ;
				//faceflag[ f ] = 0 ;
				//tot -- ;

				//if ( edge[ e2 ] == 1 ) {
					//numones ++ ;
				//} else if ( edge[ e2 ] == 0 ) {
					//numones -- ;
				//}
			//}

			//if ( tot > 0 ) {
				//return 0 ;
			//}
			//return 1 ;
		//}

		//int Volume::isNoiseSheetEnd( int ox, int oy, int oz ) {
			//int i, j ;
			//int nx, ny, nz ;

			//int edge[6] = { 0,0,0,0,0,0 } ;
			//int faceflag[ 12 ] ;
			//int hasNoiseFace = 0 ;
			//int tot = 0 ;

			//for ( i = 0 ; i < 12 ; i ++ ) {
				//faceflag[ i ] = 1 ;
				//int hasNoise = 0 ;

				//for ( j = 0 ; j < 4 ; j ++ ) {
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) == 0 ) {
						//hasNoise = 1 ;
					//}
					//if ( getDataAt( nx, ny, nz ) < 0 ) {
						//faceflag[ i ] = 0 ;
						//break ;
					//}
				//}
				//if ( faceflag[ i ] == 1 && hasNoise ) {
					//hasNoiseFace ++ ;
					//// return 0 ;
				//}

				//if ( faceflag[ i ] ) {
					//int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					//edge[ e0 ] ++ ;
					//edge[ e1 ] ++ ;
					//tot ++ ;
				//}
			//}

			//if ( hasNoiseFace < tot ) {
				//return 0 ;
			//}

			//if ( tot == 0 ) {
				//return 0 ;
			//}

			//// Removing 1s
			//int numones = 0 ;
			//for ( i = 0 ; i < 6 ; i ++ ) {
				//if ( edge[ i ] == 1 )
				//{
					//numones ++ ;
				//}
			//}
			//while( numones > 0 ) {
				//int e ;
				//for ( i = 0 ; i < 6 ; i ++ ) {
					//if ( edge[ i ] == 1 ) {
						//e = i ;
						//break ;
					//}
				//}

				//int f, e2 ;
				//for ( j = 0 ; j < 4 ; j ++ ) {
					//f = edgeFaces[ e ][ j ] ;
					//if ( faceflag[ f ] ) {
						//break ;
					//}
				//}

				//if ( faceEdges[ f ][ 0 ] == e ) {
					//e2 = faceEdges[ f ][ 1 ] ;
				//} else {
					//e2 = faceEdges[ f ][ 0 ] ;
				//}

				//edge[ e ] -- ;
				//numones -- ;
				//edge[ e2 ] -- ;
				//faceflag[ f ] = 0 ;
				//tot -- ;

				//if ( edge[ e2 ] == 1 ) {
					//numones ++ ;
				//} else if ( edge[ e2 ] == 0 ) {
					//numones -- ;
				//}
			//}

			//if ( tot > 0 ) {
				//return 0 ;
			//}
			//return 1 ;
		//}

		//int Volume::isInternal( int ox, int oy, int oz ) {
			//// assuming it's 0/1 volume
			//int i, j, k ;

			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ ) {
						//if ( getDataAt( ox + i, oy + j, oz + k ) <= 0 ) {
							//return 0 ;
						//}
					//}
			//return 1 ;
		//}

		//int Volume::isInternal2( int ox, int oy, int oz ) {
			//// assuming it's -1/0 volume
			//int i, j, k ;

			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ ) {
						//if ( getDataAt( ox + i, oy + j, oz + k ) < 0 ) {
							//return 0 ;
						//}
					//}

			//return 1 ;
		//}

		//int Volume::hasIsolatedFace( int ox, int oy, int oz ) {
			//int i, j, k ;
			//int nx, ny, nz ;

			//double vox[3][3][3] ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ ) {
						//vox[ i + 1 ][ j + 1 ][ k + 1 ] = getDataAt( ox + i, oy + j, oz + k ) ;
					//}

			//int cells[8] = { 1, 1, 1, 1, 1, 1, 1, 1 } ;
			//for ( i = 0 ; i < 8 ; i ++ ) {
				//int x = ( ( i >> 2 ) & 1 ) ;
				//int y = ( ( i >> 1 ) & 1 ) ;
				//int z = ( i & 1 ) ;
				//for ( j = 0 ; j < 8 ; j ++ ) {
					//nx = x + ( ( j >> 2 ) & 1 ) ;
					//ny = y + ( ( j >> 1 ) & 1 ) ;
					//nz = z + ( j & 1 ) ;

					//if ( vox[nx][ny][nz] < 0 ) {
						//cells[i] = 0 ;
						//break;
					//}
				//}
			//}

			//for ( i = 0 ; i < 12 ; i ++ ) {
				//if ( cells[ faceCells[i][0] ] == 1 || cells[ faceCells[i][1] ] == 1 ) {
					//continue ;
				//}
				//int flag = 1 ;
				//for ( j = 0 ; j < 4 ; j ++ ) {
					//nx = 1 + sheetNeighbor[i][j][0] ;
					//ny = 1 + sheetNeighbor[i][j][1] ;
					//nz = 1 + sheetNeighbor[i][j][2] ;

					//if ( vox[nx][ny][nz] < 0 ) {
						//flag = 0 ;
						//break;
					//}
				//}
				//if ( flag ) {
					//return 1 ;
				//}
			//}

			//return 0 ;
		//}

		//int Volume::hasIsolatedEdge( int ox, int oy, int oz ) {
			//int i, j, k ;
			//int nx, ny, nz ;

			//double vox[3][3][3] ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ ) {
						//vox[ i + 1 ][ j + 1 ][ k + 1 ] = getDataAt( ox + i, oy + j, oz + k ) ;
					//}

			//int edge[6] = { 0,0,0,0,0,0 } ;

			//for ( i = 0 ; i < 12 ; i ++ ) {
				//int flag = 1 ;
				//for ( j = 0 ; j < 4 ; j ++ ) {
					//nx = 1 + sheetNeighbor[i][j][0] ;
					//ny = 1 + sheetNeighbor[i][j][1] ;
					//nz = 1 + sheetNeighbor[i][j][2] ;

					//if ( vox[nx][ny][nz] < 0 ) {
						//flag = 0 ;
						//break ;
					//}
				//}

				//if ( flag ) {
					//int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					//edge[ e0 ] ++ ;
					//edge[ e1 ] ++ ;
				//}
			//}

			//for ( i = 0 ; i < 6 ; i ++ ) {
				//if ( edge[i] ) {
					//continue ;
				//}

				//nx = 1 + neighbor6[i][0] ;
				//ny = 1 + neighbor6[i][1] ;
				//nz = 1 + neighbor6[i][2] ;

				//if ( vox[nx][ny][nz] >= 0 ) {
					//return 1 ;
				//}

			//}

			//return 0 ;
		//}

		//int Volume::countFace( int ox, int oy, int oz, int m ) {
			//int facenum = 4 ;
			//for ( int i = 0 ; i < 4 ; i ++ ) {
				//for ( int j = 0 ; j < 4 ; j ++ ) {
					//int nx = ox + sheetNeighbor[edgeFaces[m][i]][j][0] ;
					//int ny = oy + sheetNeighbor[edgeFaces[m][i]][j][1] ;
					//int nz = oz + sheetNeighbor[edgeFaces[m][i]][j][2] ;

					//if ( getDataAt( nx, ny, nz ) < 0 ) {
						//facenum -- ;
						//break ;
					//}
				//}
			//}

			//return facenum;
		//}
		int Volume::hasCell( int ox, int oy, int oz ) {
			for ( int i = 0 ; i < 2 ; i ++ )
				for ( int j = 0 ; j < 2 ; j ++ )
					for ( int k = 0 ; k < 2 ; k ++ ) {
						if ( getDataAt( ox + i, oy + j, oz + k ) < 0 ) {
							return 0 ;
						}
					}
			return 1 ;
		}

		Volume * Volume::markCellFace( ) {
			int i, j, k ;
			Volume* fvol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;

			//return fvol ;

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if( getDataAt(i,j,k) >= 0 )
						{
							if ( hasCell(i,j,k) )
							{
								for ( int m = 0 ; m < 6 ; m ++ )
								{
									int nx = i + neighbor6[m][0] ;
									int ny = j + neighbor6[m][1] ;
									int nz = k + neighbor6[m][2] ;
									if ( ! hasCell(nx,ny,nz) )
									{
										fvol->setDataAt(i,j,k,(double)(1<<m)) ;
										break ;
									}
								}
							}
						}
					}


			return fvol ;
		}

		//Volume * Volume::markFaceEdge( ) {
			//int x,y,z, i,j ;
			//Volume* fvol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;

			////return fvol ;

			//for ( x = 0 ; x < getSizeX() ; x ++ )
				//for ( y = 0 ; y < getSizeY() ; y ++ )
					//for ( z = 0 ; z < getSizeZ() ; z ++ )
					//{
						//if( getDataAt(x,y,z) >= 0 )
						//{

							//for ( i = 0 ; i < 3 ; i ++ )
							//{
								//int hasFace = 1 ;
								//for ( j = 0 ; j < 4 ; j ++ )
								//{
									//int nx = x + sheetNeighbor[4 * i + 3][j][0] ;
									//int ny = y + sheetNeighbor[4 * i + 3][j][1] ;
									//int nz = z + sheetNeighbor[4 * i + 3][j][2] ;

									//if ( getDataAt( nx, ny, nz ) < 0 )
									//{
										//hasFace = 0 ;
										//break ;
									//}
								//}

								//if ( hasFace )
								//{
									//// Look for open edges
									//switch( i )
									//{
									//case 0:
										//if ( countFace( x, y, z, 0 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z, (double)(1<<0)) ;
											//break ;
										//}
										//if ( countFace( x, y, z, 2 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z, (double)(1<<2)) ;
											//break ;
										//}
										//if ( countFace( x, y + 1, z, 0 ) == 1 )
										//{
											//fvol->setDataAt(x, y + 1, z, (double)(1<<0)) ;
											//break ;
										//}
										//if ( countFace( x, y, z + 1, 2 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z + 1, (double)(1<<2)) ;
											//break ;
										//}
										//printf("Hmmm... a face with no open edges.\n");
										//break ;
									//case 1:
										//if ( countFace( x, y, z, 0 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z, (double)(1<<0)) ;
											//break ;
										//}
										//if ( countFace( x, y, z, 4 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z, (double)(1<<4)) ;
											//break ;
										//}
										//if ( countFace( x + 1, y, z, 0 ) == 1 )
										//{
											//fvol->setDataAt(x + 1, y, z, (double)(1<<0)) ;
											//break ;
										//}
										//if ( countFace( x, y, z + 1, 4 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z + 1, (double)(1<<4)) ;
											//break ;
										//}
										//printf("Hmmm... a face with no open edges.\n");
										//break ;
									//case 2:
										//if ( countFace( x, y, z, 2 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z, (double)(1<<2)) ;
											//break ;
										//}
										//if ( countFace( x, y, z, 4 ) == 1 )
										//{
											//fvol->setDataAt(x, y, z, (double)(1<<4)) ;
											//break ;
										//}
										//if ( countFace( x + 1, y, z, 2 ) == 1 )
										//{
											//fvol->setDataAt(x + 1, y, z, (double)(1<<2)) ;
											//break ;
										//}
										//if ( countFace( x, y + 1, z, 4 ) == 1 )
										//{
											//fvol->setDataAt(x, y + 1, z, (double)(1<<4)) ;
											//break ;
										//}
										//printf("Hmmm... a face with no open edges.\n");
										//break ;
									//}
								//}
							//}

						//}
					//}


			//return fvol ;
		//}

		int Volume::hasCompleteSheet( int ox, int oy, int oz, Volume* fvol ) {
			int i, j, k ;
			int nx, ny, nz ;

			int edge[6] = { 0,0,0,0,0,0 } ;
			int faceflag[ 12 ] ;
			int tot = 0 ;
			int cellflag[ 8 ] ;

			int ct = 0 ;
			for (  i = -1 ; i < 1 ; i ++ )
				for (  j = -1 ; j < 1 ; j ++ )
					for (  k = -1 ; k < 1 ; k ++ )
					{
						if ( hasCell( ox + i, oy + j, oz + k ) )
						{
							cellflag[ ct ] = 1 ;
						}
						else
						{
							cellflag[ ct ] = 0 ;
						}
						ct ++ ;
					}

			for ( i = 0 ; i < 12 ; i ++ )
			{
				faceflag[ i ] = 1 ;
				for ( j = 0 ; j < 4 ; j ++ )
				{
					nx = ox + sheetNeighbor[i][j][0] ;
					ny = oy + sheetNeighbor[i][j][1] ;
					nz = oz + sheetNeighbor[i][j][2] ;

					if ( getDataAt( nx, ny, nz ) < 0 )
					{
						faceflag[ i ] = 0 ;
						break ;
					}
				}

				if ( faceflag[ i ] )
				{
					if ( cellflag[ faceCells[i][0] ] ^ cellflag[ faceCells[i][1] ] )
					{
						int v1 = (int)( fvol->getDataAt(
							ox - 1 + (( faceCells[i][0] >> 2 ) & 1 ),
							oy - 1 + (( faceCells[i][0] >> 1 ) & 1 ),
							oz - 1 + (( faceCells[i][0] ) & 1)) ) ;
						int v2 = (int)( fvol->getDataAt(
							ox - 1 + (( faceCells[i][1] >> 2 ) & 1 ),
							oy - 1 + (( faceCells[i][1] >> 1 ) & 1 ),
							oz - 1 + (( faceCells[i][1] ) & 1)) ) ;
						if ( ((v1 >> (2 * (2 - i/4))) & 1) ||
							 ((v2 >> (2 * (2 - i/4) + 1 )) & 1) )
						{
							faceflag[ i ] = 0 ;
						}
					}
				}

				if ( faceflag[ i ] )
				{
					int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					edge[ e0 ] ++ ;
					edge[ e1 ] ++ ;
					tot ++ ;
				}
			}

			// Removing 1s
			int numones = 0 ;
			for ( i = 0 ; i < 6 ; i ++ )
			{
				if ( edge[ i ] == 1 )
				{
					numones ++ ;
				}
			}
			while( numones > 0 )
			{
				int e = 0;
				for ( i = 0 ; i < 6 ; i ++ )
				{
					if ( edge[ i ] == 1 )
					{
						e = i ;
						break ;
					}
				}
				/*
				if ( edge[ e ] != 1 )
				{
					printf("Wrong Again!********\n") ;
				}
				*/

				int f, e2 ;
				for ( j = 0 ; j < 4 ; j ++ )
				{
					f = edgeFaces[ e ][ j ] ;
					if ( faceflag[ f ] )
					{
						break ;
					}
				}

				/*
				if ( faceflag[ f ] == 0 )
				{
					printf("Wrong!********\n") ;
				}
				*/

				if ( faceEdges[ f ][ 0 ] == e )
				{
					e2 = faceEdges[ f ][ 1 ] ;
				}
				else
				{
					e2 = faceEdges[ f ][ 0 ] ;
				}


				edge[ e ] -- ;
				numones -- ;
				edge[ e2 ] -- ;
				faceflag[ f ] = 0 ;
				tot -- ;

				if ( edge[ e2 ] == 1 )
				{
					numones ++ ;
				}
				else if ( edge[ e2 ] == 0 )
				{
					numones -- ;
				}
			}

			if ( tot > 0 )
			{
				return 1 ;
			}

			return 0 ;
		}

		int Volume::hasCompleteSheet( int ox, int oy, int oz ) {
			// Returns 1 if it lies in the middle of a sheet
			int temp = countIntEuler( ox, oy, oz ) ;
			if ( temp > 0 )
			{
				return 1 ;
			}
			else
			{
				return 0 ;
			}
		}

		//int Volume::hasCompleteSheetSlow( int ox, int oy, int oz ) {
			//int i, j ;
			//int nx, ny, nz ;

			//int edge[6] = { 0,0,0,0,0,0 } ;
			//int faceflag[ 12 ] ;
			//int tot = 0 ;

			//for ( i = 0 ; i < 12 ; i ++ )
			//{
				//faceflag[ i ] = 1 ;
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) < 0 )
					//{
						//faceflag[ i ] = 0 ;
						//break ;
					//}
				//}

				//if ( faceflag[ i ] )
				//{
					//int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					//edge[ e0 ] ++ ;
					//edge[ e1 ] ++ ;
					//tot ++ ;
				//}
			//}

			//// Removing 1s
			//int numones = 0 ;
			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//if ( edge[ i ] == 1 )
				//{
					//numones ++ ;
				//}
			//}
			//while( numones > 0 )
			//{
				//int e ;
				//for ( i = 0 ; i < 6 ; i ++ )
				//{
					//if ( edge[ i ] == 1 )
					//{
						//e = i ;
						//break ;
					//}
				//}
				///*
				//if ( edge[ e ] != 1 )
				//{
					//printf("Wrong Again!********\n") ;
				//}
				//*/

				//int f, e2 ;
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//f = edgeFaces[ e ][ j ] ;
					//if ( faceflag[ f ] )
					//{
						//break ;
					//}
				//}

				///*
				//if ( faceflag[ f ] == 0 )
				//{
					//printf("Wrong!********\n") ;
				//}
				//*/

				//if ( faceEdges[ f ][ 0 ] == e )
				//{
					//e2 = faceEdges[ f ][ 1 ] ;
				//}
				//else
				//{
					//e2 = faceEdges[ f ][ 0 ] ;
				//}


				//edge[ e ] -- ;
				//numones -- ;
				//edge[ e2 ] -- ;
				//faceflag[ f ] = 0 ;
				//tot -- ;

				//if ( edge[ e2 ] == 1 )
				//{
					//numones ++ ;
				//}
				//else if ( edge[ e2 ] == 0 )
				//{
					//numones -- ;
				//}
			//}

	////		int temp = countIntEuler( ox, oy, oz ) ;
			//if ( tot > 0 )
			//{
	////			if ( temp <= 0 )
				//{
	////				printf("Counting wrong: %d\n", temp ) ;
				//}
				//return 1 ;
			//}
			//else
			//{
	////			if ( temp > 0 )
				//{
	////				printf("Counting wrong: %d\n", temp ) ;
				//}
				//return 0 ;
			//}
		//}
		int Volume::hasCompleteHelix( int ox, int oy, int oz )
		{
			// Returns 1 if it has a complete helix
			int i ;
			int c1 = 0;
			int nx, ny, nz ;
			int j ;

			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;
				if ( getDataAt( nx, ny, nz ) >= 0 )
				{
					c1 ++ ;
					j = i ;
				}

			}

			if ( c1 > 1 ) // || c1 == 0 )
			{
				return 1 ;
			}

			return 0 ;

			/*
			ox = ox + neighbor6[j][0] ;
			oy = oy + neighbor6[j][1] ;
			oz = oz + neighbor6[j][2] ;
			c1 = 0 ;
			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;
				if ( getDataAt( nx, ny, nz ) >= 0 )
				{
					c1 ++ ;
				}

			}

			if ( c1 > 1 )
			{
				return 0 ;
			}
			else
			{
				return 1 ;
			}
			*/
		}

		int Volume::hasCompleteHelix( int ox, int oy, int oz, Volume* fvol )
		{

			int i ;
			int c1 = 0;
			int nx, ny, nz ;
			int j ;

			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;
				if ( getDataAt( nx, ny, nz ) >= 0 )
				{
					if ( i % 2 == 0 )
					{
						nx = ox ;
						ny = oy ;
						nz = oz ;
					}

					int val = (int)fvol->getDataAt( nx, ny, nz) ;
					if ( (( val >> ( 2 * ( i / 2 ) ) ) & 1) == 0 )
					{
						c1 ++ ;
						j = i ;
					}
				}

			}

			if ( c1 > 1 )
			{
				return 1 ;
			}

			return 0 ;

			/*
			ox = ox + neighbor6[j][0] ;
			oy = oy + neighbor6[j][1] ;
			oz = oz + neighbor6[j][2] ;
			c1 = 0 ;
			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;
				if ( getDataAt( nx, ny, nz ) >= 0 )
				{
					c1 ++ ;
				}

			}

			if ( c1 > 1 )
			{
				return 0 ;
			}
			else
			{
				return 1 ;
			}
			*/
		}

		int Volume::isHelixEnd( int ox, int oy, int oz, Volume* nvol ) {
			// Returns 1 if it is a curve endpoint
			int i ;
			int c1 = 0 , c2 = 0 ;
			int nx, ny, nz ;

			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;

				double val = getDataAt( nx, ny, nz ) ;

				if ( val >= 0 )
				{
					c1 ++ ;
					if ( val > 0 && val < MAX_ERODE && nvol->getDataAt( nx, ny, nz ) == 0 )
					{
						c2 ++ ;
					}
				}

			}

			if ( c1 == 1 && c2 == 1 )
			{
				return 1 ;
			}

			return 0 ;
		}

		//int Volume::isFeature18( int ox, int oy, int oz )
		//{
			//// Border: > 0
			//// Interior: == 0
			//// Outside: < 0
			//if ( getDataAt( ox, oy, oz ) <= 0 )
			//{
				//return 0 ;
			//}

			//for ( int i = 0 ; i < 3 ; i ++ )
				//for ( int j = 0 ; j < 3 ; j ++ )
					//for ( int k = 0 ; k < 3 ; k ++ )
					//{
						//if ( i == 1 || j == 1 || k == 1 )
						//{
							//if( getDataAt( ox - 1 + i, oy - 1 + j, oz - 1 + k ) == 0 )
							//{
								//return 0 ;
							//}
						//}
					//}
			//return 1 ;
		//}

		//int Volume::isEdgeEnd( int ox, int oy, int oz )
		//{
			//int i ;
			//int c1 = 0 ;
			//int nx, ny, nz ;

			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//nx = ox + neighbor6[i][0] ;
				//ny = oy + neighbor6[i][1] ;
				//nz = oz + neighbor6[i][2] ;

				//double val = getDataAt( nx, ny, nz ) ;

				//if ( val >= 0 )
				//{
					//c1 ++ ;
				//}

			//}

			//if ( c1 == 1 )
			//{
				//return 1 ;
			//}

			//return 0 ;
		//}

		//int Volume::isFaceEnd( int ox, int oy, int oz )
		//{
			//// return isSheetEnd(ox,oy,oz) ;

			//int i, j, k ;
			//int nx, ny, nz ;

			//double vox[3][3][3] ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ )
					//{
						//vox[ i + 1 ][ j + 1 ][ k + 1 ] = getDataAt( ox + i, oy + j, oz + k ) ;
					//}

			//int edge[6] = { 4,4,4,4,4,4 } ;
			//int edge2[6] = { 4,4,4,4,4,4 } ;

			//for ( i = 0 ; i < 12 ; i ++ )
			//{
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = 1 + sheetNeighbor[i][j][0] ;
					//ny = 1 + sheetNeighbor[i][j][1] ;
					//nz = 1 + sheetNeighbor[i][j][2] ;

					//if ( vox[nx][ny][nz] < 0 )
					//{
						//edge[ faceEdges[ i ][ 0 ] ] -- ;
						//edge[ faceEdges[ i ][ 1 ] ] -- ;
						//break ;
					//}
				//}
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = 1 + sheetNeighbor[i][j][0] ;
					//ny = 1 + sheetNeighbor[i][j][1] ;
					//nz = 1 + sheetNeighbor[i][j][2] ;

					//if ( vox[nx][ny][nz] < 2 )
					//{
						//edge2[ faceEdges[ i ][ 0 ] ] -- ;
						//edge2[ faceEdges[ i ][ 1 ] ] -- ;
						//break ;
					//}
				//}
			//}

			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//if ( edge[ i ] == 1 && edge2[ i ] == 1 )
				//{
					//return 1 ;
				//}
			//}

			//return 0 ;
		//}

		//int Volume::isNoise( int ox, int oy, int oz, int noise )
		//{
			//if ( getDataAt(ox,oy,oz) == 1 )
			//{
				//return 1 ;
			//}

			//if ( noise > 0 )
			//{
				//for ( int i = 0 ; i < 6 ; i ++ )
				//{
					//int nx = ox + neighbor6[i][0] ;
					//int ny = oy + neighbor6[i][1] ;
					//int nz = oz + neighbor6[i][2] ;

					//if ( getDataAt( nx, ny, nz ) > 0 )
					//{
						//if ( isNoise( nx, ny, nz, noise - 1 ) )
						//{
							//return 1 ;
						//}
					//}
				//}
			//}

			//return 0 ;

		//}

		//int Volume::isNoiseHelixEnd( int ox, int oy, int oz )
		//{

			//int i ;
			//int c1 = 0 , c2 = 0 ;
			//int nx, ny, nz ;

			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//nx = ox + neighbor6[i][0] ;
				//ny = oy + neighbor6[i][1] ;
				//nz = oz + neighbor6[i][2] ;

				//double val = getDataAt( nx, ny, nz ) ;

				//if ( val >= 0 )
				//{
					//c1 ++ ;
					//if ( val > 0 && val < MAX_ERODE )
					//{
						//c2 ++ ;
					//}
				//}

			//}

			//if ( c1 == 1 && c2 == 0 )
			//{
				//return 1 ;
			//}

			//return 0 ;
		//}


		int Volume::isHelixEnd( int ox, int oy, int oz ) {

			int i ;
			int c1 = 0 , c2 = 0 ;
			int nx, ny, nz ;

			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;

				double val = getDataAt( nx, ny, nz ) ;

				if ( val >= 0 )
				{
					c1 ++ ;
					if ( getNumNeighbor6(nx,ny,nz) < 6 ) // if ( val > 0 && val < MAX_ERODE )
					{
						c2 ++ ;
					}
				}

			}

			if ( c1 == 1 && c2 == 1 )
			{
				return 1 ;
			}

			return 0 ;
		}
		int Volume::isSheetEnd( int ox, int oy, int oz, Volume* nvol )
		{
			// Returns 1 if it contains a sheet boundary. Noise-resistant
			int i, j ;
			int nx, ny, nz ;

			int edge[6] = { 0,0,0,0,0,0 } ;
			int faceflag[ 12 ] ;
			int hasFeatureFace = 0 ;
			int tot = 0 ;

			for ( i = 0 ; i < 12 ; i ++ )
			{
				faceflag[ i ] = 1 ;
				int hasFeature = 1 ;

				for ( j = 0 ; j < 4 ; j ++ )
				{
					nx = ox + sheetNeighbor[i][j][0] ;
					ny = oy + sheetNeighbor[i][j][1] ;
					nz = oz + sheetNeighbor[i][j][2] ;

					if ( getDataAt( nx, ny, nz ) == 0 || nvol->getDataAt( nx, ny, nz ) == 1 )
					{
						hasFeature = 0 ;
					}
					if ( getDataAt( nx, ny, nz ) < 0 )
					{
						faceflag[ i ] = 0 ;
						break ;
					}
				}
				if ( faceflag[ i ] == 1 && hasFeature )
				{
					hasFeatureFace ++ ;
					// return 0 ;
				}

				if ( faceflag[ i ] )
				{
					int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					edge[ e0 ] ++ ;
					edge[ e1 ] ++ ;
					tot ++ ;
				}
			}

			if ( tot == 0 || hasFeatureFace == 0 )
			{
				return 0 ;
			}

			// Removing 1s
			int numones = 0 ;
			for ( i = 0 ; i < 6 ; i ++ )
			{
				if ( edge[ i ] == 1 )
				{
					numones ++ ;
				}
			}
			while( numones > 0 )
			{
				int e = 0;
				for ( i = 0 ; i < 6 ; i ++ )
				{
					if ( edge[ i ] == 1 )
					{
						e = i ;
						break ;
					}
				}
				/*
				if ( edge[ e ] != 1 )
				{
					printf("Wrong Again!********\n") ;
				}
				*/

				int f, e2 ;
				for ( j = 0 ; j < 4 ; j ++ )
				{
					f = edgeFaces[ e ][ j ] ;
					if ( faceflag[ f ] )
					{
						break ;
					}
				}

				/*
				if ( faceflag[ f ] == 0 )
				{
					printf("Wrong!********\n") ;
				}
				*/

				if ( faceEdges[ f ][ 0 ] == e )
				{
					e2 = faceEdges[ f ][ 1 ] ;
				}
				else
				{
					e2 = faceEdges[ f ][ 0 ] ;
				}


				edge[ e ] -- ;
				numones -- ;
				edge[ e2 ] -- ;
				faceflag[ f ] = 0 ;
				tot -- ;

				if ( edge[ e2 ] == 1 )
				{
					numones ++ ;
				}
				else if ( edge[ e2 ] == 0 )
				{
					numones -- ;
				}
			}

			if ( tot > 0 )
			{
				return 0 ;
			}

			return 1 ;
		}

		//int Volume::getNumFaces( int ox, int oy, int oz )
		//{
			//int i, j ;
			//int nx, ny, nz ;

			//int faces = 12 ;
			//for ( i = 0 ; i < 12 ; i ++ )
			//{
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) < 0 )
					//{
						//faces -- ;
						//break ;
					//}
				//}
			//}

			//return faces ;
		//}

		//int Volume::getNumCells( int ox, int oy, int oz )
		//{
			//int i, j, k ;

			//int cells = 0 ;

			//for (  i = -1 ; i < 1 ; i ++ )
				//for (  j = -1 ; j < 1 ; j ++ )
					//for (  k = -1 ; k < 1 ; k ++ )
					//{
						//if ( hasCell( ox + i, oy + j, oz + k ) )
						//{
							//cells ++ ;
						//}
					//}

			//// printf("Faces: %d\n", faces);
			//return cells ;
		//}


		//int Volume::getNumIsolatedEdges( int ox, int oy, int oz )
		//{
			//int i, j, k ;
			//int nx, ny, nz ;

			//double vox[3][3][3] ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ )
					//{
						//vox[ i + 1 ][ j + 1 ][ k + 1 ] = getDataAt( ox + i, oy + j, oz + k ) ;
					//}


			//int edge[6] = { 0,0,0,0,0,0 } ;

			//for ( i = 0 ; i < 12 ; i ++ )
			//{
				//int flag = 1 ;
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = 1 + sheetNeighbor[i][j][0] ;
					//ny = 1 + sheetNeighbor[i][j][1] ;
					//nz = 1 + sheetNeighbor[i][j][2] ;

					//if ( vox[nx][ny][nz] < 0 )
					//{
						//flag = 0 ;
						//break ;
					//}
				//}

				//if ( flag )
				//{
					//int e0 = faceEdges[ i ][ 0 ], e1 = faceEdges[ i ][ 1 ] ;
					//edge[ e0 ] ++ ;
					//edge[ e1 ] ++ ;
				//}
			//}

			//int edges = 0 ;
			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//if ( edge[i] )
				//{
					//continue ;
				//}

				//nx = 1 + neighbor6[i][0] ;
				//ny = 1 + neighbor6[i][1] ;
				//nz = 1 + neighbor6[i][2] ;

				//if ( vox[nx][ny][nz] >= 0 )
				//{
					//edges ++ ;
				//}

			//}

			//return edges ;
		//}

		//int Volume::getNumIsolatedFaces( int ox, int oy, int oz )
		//{
			//int i, j, k ;
			//int nx, ny, nz ;

			//int faces = 0 ;
			//int cellflag[ 8 ] ;

			//int ct = 0 ;
			//for (  i = -1 ; i < 1 ; i ++ )
				//for (  j = -1 ; j < 1 ; j ++ )
					//for (  k = -1 ; k < 1 ; k ++ )
					//{
						//if ( hasCell( ox + i, oy + j, oz + k ) )
						//{
							//cellflag[ ct ] = 1 ;
						//}
						//else
						//{
							//cellflag[ ct ] = 0 ;
						//}
						//ct ++ ;
					//}

			//for ( i = 0 ; i < 12 ; i ++ )
			//{
				//int hasFace = 1 ;
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) < 0 )
					//{
						//hasFace = 0 ;
						//break ;
					//}
				//}

				//if ( hasFace )
				//{
					//if (cellflag[ faceCells[i][0] ] == 0 && cellflag[ faceCells[i][1] ] == 0 )
					//{
						//faces ++ ;
					//}
				//}
			//}

			//// printf("Faces: %d\n", faces);
			//return faces ;
		//}

		//int Volume::isFeatureFace2( int ox, int oy, int oz )
		//{
			//int i ;
			//int nx, ny, nz ;

			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//nx = ox + neighbor6[i][0] ;
				//ny = oy + neighbor6[i][1] ;
				//nz = oz + neighbor6[i][2] ;

				//double val = getDataAt( nx, ny, nz ) ;

				//if ( val == 0 )
				//{
					//return 0 ;
				//}

			//}

			//return 1 ;
		//}

		int Volume::isFeatureFace( int ox, int oy, int oz )
		{
			// return 1 ;

			int i, j ;
			int nx, ny, nz ;

			int faces = 12 ;
			for ( i = 0 ; i < 12 ; i ++ )
			{
				int ct = 0 ;
				for ( j = 0 ; j < 4 ; j ++ )
				{
					nx = ox + sheetNeighbor[i][j][0] ;
					ny = oy + sheetNeighbor[i][j][1] ;
					nz = oz + sheetNeighbor[i][j][2] ;

					if ( getDataAt( nx, ny, nz ) < 0 )
					{
						ct = -1 ;
						break ;
					}
					else if ( getNumNeighbor6( nx, ny, nz ) == 6 )
					{
						ct = -1 ;
						break ;
					}
	//				else if ( getDataAt( nx, ny, nz ) == 0 )
	//				{
	//					ct ++ ;
	//				}


				}
				if ( ct == -1 || ct >= 1 )
				{
					faces -- ;
				}
			}

			if ( faces > 0 )
			{
				return 1 ;
			}
			return 0 ;

		}

		//int Volume::hasFeatureFace( int ox, int oy, int oz )
		//{
			//int i, j, k ;
			//int nx, ny, nz ;

			//int faceflag ;
			//int cellflag[ 8 ] ;

			//// Get cells
			//int ct = 0 ;
			//for (  i = -1 ; i < 1 ; i ++ )
				//for (  j = -1 ; j < 1 ; j ++ )
					//for (  k = -1 ; k < 1 ; k ++ )
					//{
						//if ( hasCell( ox + i, oy + j, oz + k ) )
						//{
							//cellflag[ ct ] = 1 ;
						//}
						//else
						//{
							//cellflag[ ct ] = 0 ;
						//}
						//ct ++ ;
					//}

			//// Find isolated and boundary faces
			//for ( i = 0 ; i < 12 ; i ++ )
			//{
				//faceflag = 1 ;
				//for ( j = 0 ; j < 4 ; j ++ )
				//{
					//nx = ox + sheetNeighbor[i][j][0] ;
					//ny = oy + sheetNeighbor[i][j][1] ;
					//nz = oz + sheetNeighbor[i][j][2] ;

					//if ( getDataAt( nx, ny, nz ) < 0 )
					//{
						//faceflag = 0 ;
						//break ;
					//}
				//}

				//if ( faceflag )
				//{
					//if ( cellflag[ faceCells[i][0] ] == 0 && cellflag[ faceCells[i][1] ] == 0 )
					//{
						//return 1 ;
					//}
				//}
			//}

			//return 0 ;

		//}

		int Volume::isSheetEnd( int ox, int oy, int oz )
		{
			return ( ( hasCompleteSheet(ox,oy,oz) == 0 ) && isFeatureFace(ox,oy,oz) ) ;

			//// return testIsSheetEnd( ox, oy, oz ) ;
			//
			//int i, j, k ;
			//int nx, ny, nz ;

			//double vox[3][3][3] ;
			//for ( i = -1 ; i < 2 ; i ++ )
			//	for ( j = -1 ; j < 2 ; j ++ )
			//		for ( k = -1 ; k < 2 ; k ++ )
			//		{
			//			vox[ i + 1 ][ j + 1 ][ k + 1 ] = getDataAt( ox + i, oy + j, oz + k ) ;
			//		}

			//int edge[6] = { 4,4,4,4,4,4 } ;
			//int edge2[6] = { 4,4,4,4,4,4 } ;

			//for ( i = 0 ; i < 12 ; i ++ )
			//{
			//	for ( j = 0 ; j < 4 ; j ++ )
			//	{
			//		nx = 1 + sheetNeighbor[i][j][0] ;
			//		ny = 1 + sheetNeighbor[i][j][1] ;
			//		nz = 1 + sheetNeighbor[i][j][2] ;

			//		if ( vox[nx][ny][nz] < 0 )
			//		{
			//			edge[ faceEdges[ i ][ 0 ] ] -- ;
			//			edge[ faceEdges[ i ][ 1 ] ] -- ;
			//			break ;
			//		}
			//	}

			//	for ( j = 0 ; j < 4 ; j ++ )
			//	{
			//		nx = 1 + sheetNeighbor[i][j][0] ;
			//		ny = 1 + sheetNeighbor[i][j][1] ;
			//		nz = 1 + sheetNeighbor[i][j][2] ;

			//		if ( vox[nx][ny][nz] <= 0 )
			//		{
			//			edge2[ faceEdges[ i ][ 0 ] ] -- ;
			//			edge2[ faceEdges[ i ][ 1 ] ] -- ;
			//			break ;
			//		}
			//	}
			//}

			//
			///*
			//for ( i = 0 ; i < 6 ; i ++ )
			//{
			//	nx = 1 + neighbor6[i][0] ;
			//	ny = 1 + neighbor6[i][1] ;
			//	nz = 1 + neighbor6[i][2] ;
			//	if ( edge[i] == 0 && vox[nx][ny][nz] >= 0 )
			//	{
			//		return 0 ;
			//	}
			//}
			//*/
			//


			//for ( i = 0 ; i < 6 ; i ++ )
			//{
			//	if ( edge[ i ] == 1 && edge2[ i ] == 1 )
			//	{
			//		return 1 ;
			//	}
			//}

			//return 0 ;
		}

		int Volume::isSimple( int ox, int oy, int oz )
		{
			/* Test if this is a simple voxel */
			// int flag = 0 ;
			double vox[3][3][3] ;

			int i, j, k ;
			for ( i = -1 ; i < 2 ; i ++ )
				for ( j = -1 ; j < 2 ; j ++ )
					for ( k = -1 ; k < 2 ; k ++ )
					{
						double tval = getDataAt( ox + i, oy + j, oz + k ) ;

						/*
						if ( tval == 0 || tval > (va )
						{
							flag = 1 ;
						}
						*/
						/*
						if ( tval < 0 && tval == - curwid )
						{
						printf("Here %d", (int)-tval) ;
						vox[ i + 1 ][ j + 1 ][ k + 1 ] = - tval ;
						}
						else
						*/
						{
							vox[ i + 1 ][ j + 1 ][ k + 1 ] = tval ;
						}
					}

				/* Debugging
				printf("{") ;
				for ( i = 0 ; i < 3 ; i ++ )
				{
					if ( i ) printf(",") ;
					printf("{") ;
					for ( j = 0 ; j < 3 ; j ++ )
					{
						if ( j ) printf(",") ;
						printf("{") ;
						for ( k = 0 ; k < 3 ; k ++ )
						{
							if ( k ) printf(",") ;
							printf("%d", (vox[i][j][k] >=0 ? 1: 0));
						}
						printf("}") ;
					}
					printf("}") ;
				}
				printf("} Int: %d, Ext: %d\n", countInt( vox ), countExt( vox )) ;
				*/

			if ( countInt( vox ) == 1 && countExt( vox ) == 1 )
			{
				return 1 ;
			}
			else
			{
				return 0 ;
			}
		}


		int Volume::isPiercable( int ox, int oy, int oz )
		{
			/* Test if this is a simple voxel */
			// int flag = 0 ;
			double vox[3][3][3] ;

			int i, j, k ;
			for ( i = -1 ; i < 2 ; i ++ )
				for ( j = -1 ; j < 2 ; j ++ )
					for ( k = -1 ; k < 2 ; k ++ )
					{
						double tval = getDataAt( ox + i, oy + j, oz + k ) ;

						/*
						if ( tval == 0 || tval > (va )
						{
							flag = 1 ;
						}
						*/
						/*
						if ( tval < 0 && tval == - curwid )
						{
						printf("Here %d", (int)-tval) ;
						vox[ i + 1 ][ j + 1 ][ k + 1 ] = - tval ;
						}
						else
						*/
						{
							vox[ i + 1 ][ j + 1 ][ k + 1 ] = tval ;
						}
					}

				/* Debugging
				printf("{") ;
				for ( i = 0 ; i < 3 ; i ++ )
				{
					if ( i ) printf(",") ;
					printf("{") ;
					for ( j = 0 ; j < 3 ; j ++ )
					{
						if ( j ) printf(",") ;
						printf("{") ;
						for ( k = 0 ; k < 3 ; k ++ )
						{
							if ( k ) printf(",") ;
							printf("%d", (vox[i][j][k] >=0 ? 1: 0));
						}
						printf("}") ;
					}
					printf("}") ;
				}
				printf("} Int: %d, Ext: %d\n", countInt( vox ), countExt( vox )) ;
				*/

			if ( countInt( vox ) == 1 && countExt( vox ) != 1 )
			{
				return 1 ;
			}
			else
			{
				return 0 ;
			}
		}


		//int Volume::isSimple2( int v[3][3][3] )
		//{
			//// int flag = 0 ;
			//double vox[3][3][3] ;

			//int i, j, k ;
			//for ( i = 0 ; i < 3 ; i ++ )
				//for ( j = 0 ; j < 3 ; j ++ )
					//for ( k = 0 ; k < 3 ; k ++ )
					//{
						//if ( v[i][j][k] == 0 )
						//{
							//vox[ i ][ j ][ k ] = 1 ;
						//}
						//else
						//{
							//vox[i][j][k] = -1 ;
						//}
					//}
			//if ( countInt( vox ) == 1 && countExt( vox ) == 1 )
			//{
				//return 1 ;
			//}
			//else
			//{
				//printf("Int: %d Ext: %d\n",countInt( vox ),countExt( vox ) );
				//return 0 ;
			//}
		//}

		//int Volume::getNumPotComplex3( int ox, int oy, int oz )
		//{
			//// return 0 ;


			//int i, j, k ;
			//if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz ) )
			//{
				//return 60 ;
			//}
			//double val = getDataAt( ox, oy, oz ) ;

			//int nx, ny, nz ;


			//int numSimple = 0 ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ )
					//{
						//nx = ox + i ;
						//ny = oy + j ;
						//nz = oz + k ;
						//if ( getDataAt( nx, ny, nz ) > 0 )
						//{
							//if ( isSimple( nx, ny, nz ) && ! isSheetEnd( nx, ny, nz ) )
							//{
								//numSimple ++ ;
							//}
						//}
					//}
			//numSimple -- ;

			//setDataAt( ox, oy, oz, -val ) ;

			//int numPotSimple = 0 ;
			//for ( i = -1 ; i < 2 ; i ++ )
				//for ( j = -1 ; j < 2 ; j ++ )
					//for ( k = -1 ; k < 2 ; k ++ )
					//{
						//nx = ox + i ;
						//ny = oy + j ;
						//nz = oz + k ;
						//if ( getDataAt( nx, ny, nz ) > 0 )
						//{
							//if ( isSimple( nx, ny, nz ) && ! isSheetEnd( nx, ny, nz ) )
							//{
								//numPotSimple ++ ;
							//}
						//}
					//}


			//setDataAt( ox, oy, oz, val ) ;




			//return 30 - ( numPotSimple  - numSimple ) ;
		//}

		//int Volume::getNumPotComplex4( int ox, int oy, int oz )
		//{

			//int cells = getNumCells(ox,oy,oz) ;
			////int ifaces = getNumIsolatedFaces(ox,oy,oz) ;
			//int faces = getNumFaces(ox,oy,oz) ;
			////int iedges = getNumIsolatedEdges(ox,oy,oz) ;

			//return ( cells * 6 - 2 * faces ) ; // + 2 * ( faces * 4 - 4 * iedges ) ;
		//}

		int Volume::getNumPotComplex( int ox, int oy, int oz )
		{
			//return 0 ;

			int i;
			double val = getDataAt( ox, oy, oz ) ;
			if ( val <= 0 )
			{
		//		return 70 ;
			}

			// return getNumNeighbor6( ox, oy, oz ) ;

			// int v = ((getNumNeighbor6( ox, oy, oz ) & 255) << 24) ;
			//int v = 0  ;

			int rvalue = 0, nx, ny, nz ;
			setDataAt( ox, oy, oz, -val ) ;

			/*
			for ( i = -1 ; i < 2 ; i ++ )
				for ( j = -1 ; j < 2 ; j ++ )
					for ( k = -1 ; k < 2 ; k ++ )
					{
						nx = ox + i ;
						ny = oy + j ;
						nz = oz + k ;
						if ( getDataAt( nx, ny, nz ) == val )
						{
							if ( isSheetEnd( nx, ny, nz) || ! isSimple ( nx, ny, nz ) )
							{
								rvalue ++ ;
							}
						}
					}
			*/

			for ( i = 0 ; i < 6 ; i ++ )
			{
				nx = ox + neighbor6[i][0] ;
				ny = oy + neighbor6[i][1] ;
				nz = oz + neighbor6[i][2] ;

				if ( getDataAt(nx,ny,nz) >= 0 )
				{
					int num = getNumNeighbor6( nx, ny, nz ) ;
					if ( num > rvalue )
					{
						rvalue = num ;
					}
				}
			}


			setDataAt( ox, oy, oz, val ) ;

			return rvalue + getNumNeighbor6( ox, oy, oz ) * 10 ;
			/*
			int v = (((rvalue + getNumNeighbor6( ox, oy, oz ) * 10) & 255) << 24) ;
			v |= ( ( ox & 255 ) << 16 )  ;
			v |= ( ( oy & 255 ) << 8 ) ;
			v |= ( ( oz & 255 ) ) ;
			return v ;
			*/

		}

		int Volume::getNumPotComplex2( int ox, int oy, int oz )
		{
			return getNumPotComplex( ox, oy, oz ) ;

			//int i, j, k ;
			//double val = getDataAt( ox, oy, oz ) ;
			//if ( val <= 0 )
			//{
			//	return 0 ;
			//}

			//int rvalue = 0, nx, ny, nz ;
			//setDataAt( ox, oy, oz, -val ) ;

			//for ( i = -1 ; i < 2 ; i ++ )
			//	for ( j = -1 ; j < 2 ; j ++ )
			//		for ( k = -1 ; k < 2 ; k ++ )
			//		{
			//			nx = ox + i ;
			//			ny = oy + j ;
			//			nz = oz + k ;
			//			if ( getDataAt( nx, ny, nz ) == val )
			//			{
			//				if ( isHelixEnd( nx, ny, nz) || ! isSimple ( nx, ny, nz ) )
			//				{
			//					rvalue ++ ;
			//				}
			//			}
			//		}

			//setDataAt( ox, oy, oz, val ) ;

			//return rvalue + getNumNeighbor6( ox, oy, oz ) * 30 ;
		}

		//int Volume::getNumNeighbor( int ox, int oy, int oz )
		//{
			//int rvalue = 0 ;
			//double val = getDataAt( ox, oy, oz ) ;
			//for ( int i = 0 ; i < 6 ; i ++ )
			//{
				//int nx = ox + neighbor6[i][0] ;
				//int ny = oy + neighbor6[i][1] ;
				//int nz = oz + neighbor6[i][2] ;

				//if ( getDataAt( nx, ny, nz ) == val )
				//{
					//rvalue ++ ;
				//}
			//}
			///*
			//for ( int i = -1 ; i < 2 ; i ++ )
				//for ( int j = -1 ; j < 2 ; j ++ )
					//for ( int k = -1 ; k < 2 ; k ++ )
					//{
						//int nx = ox + i ;
						//int ny = oy + j ;
						//int nz = oz + k ;

						//if ( getDataAt( nx, ny, nz ) == val )
						//{
							//rvalue ++ ;
						//}
					//}
			//*/
			//return rvalue ;
		//}


		//void Volume::setScoreNeighbor( GridQueue* queue )
		//{
			////printf("Scoring each node with number of neighbors...\n") ;
			//gridQueueEle* ele = queue->getHead() ;
			//while ( ele != NULL )
			//{
				//ele->score = getNumNeighbor( ele->x, ele->y, ele->z ) ;
				//ele = ele->next ;
			//}

			//queue->sort( queue->getNumElements() ) ;
		//}


		int Volume::components6( int vox[3][3][3] )
		{
			// Stupid flood fill
			int tot = 0 ;
			int queue[27][3] ;
			int vis[3][3][3] ;
			int head = 0, tail = 1 ;
			int i, j, k ;
			for ( i = 0 ; i < 3 ; i ++ )
				for ( j = 0 ; j < 3 ; j ++ )
					for ( k = 0 ; k < 3 ; k ++ )
					{
						vis[i][j][k] = 0 ;
						if ( vox[i][j][k] )
						{
							if ( tot == 0 )
							{
								queue[0][0] = i ;
								queue[0][1] = j ;
								queue[0][2] = k ;
								vis[i][j][k] = 1 ;
							}
							tot ++ ;
						}
					}
			if ( tot == 0 )
			{
				return 0 ;
			}
			// printf("total: %d\n", tot) ;

			int ct = 1 ;
			while ( head != tail )
			{
				int x = queue[head][0] ;
				int y = queue[head][1] ;
				int z = queue[head][2] ;
				head ++ ;

				for ( i = 0 ; i < 6 ; i ++ )
				{
					int nx = x + neighbor6[i][0] ;
					int ny = y + neighbor6[i][1] ;
					int nz = z + neighbor6[i][2] ;
					if ( nx >=0 && nx < 3 && ny >=0 && ny < 3 && nz >=0 && nz < 3 )
					{
						if ( vox[nx][ny][nz] && ! vis[nx][ny][nz] )
						{
							queue[tail][0] = nx ;
							queue[tail][1] = ny ;
							queue[tail][2] = nz ;
							tail ++ ;
							vis[nx][ny][nz] = 1 ;
							ct ++ ;
						}
					}
				}
			}

			if ( ct == tot )
			{
				return 1 ;
			}
			else
			{
				return 2 ;
			}

		}
		int Volume::components26( int vox[3][3][3] )
		{
			// Stupid flood fill
			int tot = 0 ;
			int queue[27][3] ;
			int vis[3][3][3] ;
			int head = 0, tail = 1 ;
			int i, j, k ;
			for ( i = 0 ; i < 3 ; i ++ )
				for ( j = 0 ; j < 3 ; j ++ )
					for ( k = 0 ; k < 3 ; k ++ )
					{
						vis[i][j][k] = 0 ;
						if ( vox[i][j][k] )
						{
							if ( tot == 0 )
							{
								queue[0][0] = i ;
								queue[0][1] = j ;
								queue[0][2] = k ;
								vis[i][j][k] = 1 ;
							}
							tot ++ ;
						}
					}
			if ( tot == 0 )
			{
				return 0 ;
			}

			int ct = 1 ;
			while ( head != tail )
			{
				int x = queue[head][0] ;
				int y = queue[head][1] ;
				int z = queue[head][2] ;
				head ++ ;

				for ( i = -1 ; i < 2 ; i ++ )
					for ( j = -1 ; j < 2 ; j ++ )
						for ( k = -1 ; k < 2 ; k ++ )
						{
							int nx = x + i ;
							int ny = y + j ;
							int nz = z + k ;
							if ( nx >=0 && nx < 3 && ny >=0 && ny < 3 && nz >=0 && nz < 3 )
							{
								if ( vox[nx][ny][nz] && ! vis[nx][ny][nz] )
								{
									queue[tail][0] = nx ;
									queue[tail][1] = ny ;
									queue[tail][2] = nz ;
									tail ++ ;
									vis[nx][ny][nz] = 1 ;
									ct ++ ;
								}
							}
						}
			}

			if ( ct == tot )
			{
				return 1 ;
			}
			else
			{
				return 2 ;
			}

		}

		int Volume::countExt( double vox[3][3][3] )
		{
			int tvox[3][3][3] ;

			for ( int i = 0 ; i < 3 ; i ++ )
				for ( int j = 0 ; j < 3 ; j ++ )
					for ( int k = 0 ; k < 3 ; k ++ )
					{
						if ( vox[i][j][k] < 0 )
						{
							tvox[i][j][k] = 1 ;
						}
						else
						{
							tvox[i][j][k] = 0 ;
						}
					}

			return components26( tvox ) ;
		}

		int Volume::countInt( double vox[3][3][3] )
		{
			int i, j, k ;
			int tvox[3][3][3] ;

			for ( i = 0 ; i < 3 ; i ++ )
				for ( j = 0 ; j < 3 ; j ++ )
					for ( k = 0 ; k < 3 ; k ++ )
					{
						tvox[i][j][k] = 0 ;
					}

			for ( i = 0 ; i < 6 ; i ++ )
			{
				int nx = 1 + neighbor6[i][0] ;
				int ny = 1 + neighbor6[i][1] ;
				int nz = 1 + neighbor6[i][2] ;
				if ( vox[ nx ][ ny ][ nz ] >= 0 )
				{
					tvox[ nx ][ ny ][ nz ] = 1 ;
					for ( j = 0 ; j < 4 ; j ++ )
					{
						int nnx = neighbor64[i][j][0] + nx ;
						int nny = neighbor64[i][j][1] + ny ;
						int nnz = neighbor64[i][j][2] + nz ;
						if ( vox[ nnx ][ nny ][ nnz ] >= 0 )
						{
							tvox[ nnx ][ nny ][ nnz ] = 1 ;
						}
					}
				}
			}

			return components6( tvox ) ;
		}

		int Volume::countIntEuler( int ox, int oy, int oz )
		{
			int nv = 0 , ne = 0, nc = 0 ;

			int i, j, k ;
			int tvox[3][3][3] ;
			double vox[3][3][3] ;

			for ( i = 0 ; i < 3 ; i ++ )
				for ( j = 0 ; j < 3 ; j ++ )
					for ( k = 0 ; k < 3 ; k ++ )
					{
						vox[i][j][k] = getDataAt( ox - 1 + i, oy - 1 + j, oz - 1 + k ) ;
						tvox[i][j][k] = 0 ;
					}

			for ( i = 0 ; i < 6 ; i ++ )
			{
				int nx = 1 + neighbor6[i][0] ;
				int ny = 1 + neighbor6[i][1] ;
				int nz = 1 + neighbor6[i][2] ;
				if ( vox[nx][ny][nz] >= 0 )
				{
					tvox[ nx ][ ny ][ nz ] = 1 ;

					nv ++ ;

					for ( j = 0 ; j < 4 ; j ++ )
					{
						int nnx = neighbor64[i][j][0] + nx ;
						int nny = neighbor64[i][j][1] + ny ;
						int nnz = neighbor64[i][j][2] + nz ;
						if ( vox[nnx][nny][nnz] >= 0 )
						{
							if ( tvox[nnx][nny][nnz] == 0 )
							{
								tvox[nnx][nny][nnz] = 1 ;
								nv ++ ;
							}

							ne ++ ;
						}
					}
				}
			}

			nc = components6( tvox ) ;

			return ( nc - ( nv - ne ) ) ;
		}

		//void Volume::erodeNoTopo( float thr, int wid )
		//{
		///* No topology check */
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue* queue = new GridQueue( ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//setDataAt( i, j, k, 1 ) ;
									//queue->pushQueue( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue->getNumElements() ) ;

			//// Perform erosion
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//double val = 0;
			//int ox, oy, oz ;
			//int curwid = 0 ;
			//int total = 0, ct = 0 ;
			//while ( 1 )
			//{
				//if ( ct == total )
				//{
					//#ifdef VERBOSE
					//printf("Layer %d has %d nodes.\n", (int) curwid, total ) ;
					//#endif
					//curwid ++ ;
					//ct = 0 ;
					//total = queue->getNumElements() ;
					//if ( total == 0 )
					//{
						//break ;
					//}
				//}

				//queue->popQueue(ox, oy, oz) ;
				//val = getDataAt( ox, oy, oz ) ;
				//if ( val > wid )
				//{
					//break ;
				//}
				//ct ++ ;

				//setDataAt( ox, oy, oz, -val ) ;


				//for ( int m = 0 ; m < 6 ; m ++ )
				//{
					//int nx = ox + neighbor6[m][0] ;
					//int ny = oy + neighbor6[m][1] ;
					//int nz = oz + neighbor6[m][2] ;
					//if ( getDataAt( nx, ny, nz ) == 0 )
					//{
						//setDataAt( nx, ny, nz, val + 1 ) ;
						//queue->pushQueue( nx, ny, nz ) ;
					//}
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::erodeTopo( float thr, int wid )
		//{
			///* Minimal topology check */
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue* queue = new GridQueue( ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//setDataAt( i, j, k, 1 ) ;
									//queue->pushQueue( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue->getNumElements() ) ;

			//// Perform erosion
			//printf("Start erosion to %d...\n", wid) ;
			//#endif

			//double val = 0;
			//int ox, oy, oz ;
			//int curwid = 0 ;
			//int total = 0, ct = 0 ;
			//while ( 1 )
			//{
				//if ( ct == total )
				//{
					//#ifdef VERBOSE
					//printf("Layer %d has %d nodes.\n", (int) curwid, total ) ;
					//#endif
					//curwid ++ ;
					//ct = 0 ;
					//total = queue->getNumElements() ;
					//if ( total == 0 )
					//{
						//break ;
					//}
				//}

				//queue->popQueue(ox, oy, oz) ;
				//val = getDataAt( ox, oy, oz ) ;
				//if ( val > wid )
				//{
					//break ;
				//}
				//ct ++ ;

				//if (  isSimple( ox, oy, oz ) )
				//{
					//// Simple node, remove
					//setDataAt( ox, oy, oz, -val ) ;
				//}
				//else
				//{
					//// Preserve for topology
					//setDataAt( ox, oy, oz, val + 1 ) ;
					//queue->pushQueue( ox, oy, oz ) ;
				//}


				//for ( int m = 0 ; m < 6 ; m ++ )
				//{
					//int nx = ox + neighbor6[m][0] ;
					//int ny = oy + neighbor6[m][1] ;
					//int nz = oz + neighbor6[m][2] ;
					//if ( getDataAt( nx, ny, nz ) == 0 )
					//{
						//setDataAt( nx, ny, nz, val + 1 ) ;
						//queue->pushQueue( nx, ny, nz ) ;
					//}
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::erode2( float thr, int wid )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue = new GridQueue2( ) ;
			//GridQueue2* queue2 = new GridQueue2( ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//setDataAt( i, j, k, 1 ) ;
									//queue->prepend( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue->getNumElements() ) ;

			//// Perform erosion
			//// wid = MAX_ERODE ;
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//gridQueueEle* ele ;
			//int ox, oy, oz ;

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue->getNumElements(), curwid) ;
				//#endif

				///* set nodes for next layer
				//while ( ( ele = queue->getNext() ) != NULL )
				//{
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ele->x + neighbor6[m][0] ;
						//int ny = ele->y + neighbor6[m][1] ;
						//int nz = ele->z + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

				//}
				//*/

				//// erosion
				//int seed[3] = {-1,-1,-1};
				//queue->reset() ;
				//while ( queue->getNumElements() > 0 )
				//{
					//if ( seed[0] < 0 ) printf("After initial scoring...\n");
					//queue->reset() ;
					//ele = queue->getNext() ;

					//// Secure complex nodes
					//while ( ele != NULL )
					//{
						//ox = ele->x ;
						//oy = ele->y ;
						//oz = ele->z ;

						//// Check simple only if within the last modified range
						//if ( seed[0] < 0 ||
							//( ox < seed[0] + 2 && ox > seed[0] - 2 &&
							//oy < seed[1] + 2 && oy > seed[1] - 2 &&
							//oz < seed[2] + 2 && oz > seed[2] - 2 ) )
						//{
							//if ( ! isSimple( ox, oy, oz ) )
							//{
								//// Complex, set to next layer
								//setDataAt( ox, oy, oz, curwid + 1 ) ;
								//queue2->prepend( ox, oy, oz ) ;
								//ele = queue->remove() ;

								//numComplex ++ ;
							//}
							//else
							//{
								//ele = queue->getNext() ;
							//}
						//}
						//else
						//{
							//ele = queue->getNext() ;
						//}
					//}

					//// Remove the simple node with the most potential neighbors
					//queue->reset() ;
					//ele = queue->getNext() ;
					//int preScore = -1;
					//while ( ele != NULL )
					//{
						//ox = ele->x ;
						//oy = ele->y ;
						//oz = ele->z ;

						//// Update score only if within the last modified range
						//if ( seed[0] < 0 ||
							//( ox < seed[0] + 3 && ox > seed[0] - 3 &&
							  //oy < seed[1] + 3 && oy > seed[1] - 3 &&
							  //oz < seed[2] + 3 && oz > seed[2] - 3 ) )
						//{
							//ele->score = getNumPotComplex( ox, oy, oz ) ;
						//}


						//if ( ele->score < preScore )
						//{
							//// Swap
							//ele = queue->swap() ;
						//}
						//else
						//{
							//preScore = ele->score ;
						//}

						//// At the end of the queue, remove this simple node
						//if ( ele->next == NULL )
						//{
							//ox = ele->x ;
							//oy = ele->y ;
							//oz = ele->z ;
							//setDataAt( ox, oy, oz, -1 ) ;
	////						printf("%d %d %d\n", ox, oy, oz) ;
							//seed[0] = ox ;
							//seed[1] = oy ;
							//seed[2] = oz ;
							//queue->remove() ;
							//// printf("Highest score: %d\n", preScore) ;


							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//int nx = ox + neighbor6[m][0] ;
								//int ny = oy + neighbor6[m][1] ;
								//int nz = oz + neighbor6[m][2] ;
								//if ( getDataAt( nx, ny, nz ) == 0 )
								//{
									//setDataAt( nx, ny, nz, curwid + 1 ) ;
									//queue2->prepend( nx, ny, nz ) ;
								//}
							//}


							//numSimple ++ ;
							//ele = NULL ;
						//}
						//else
						//{
							//ele = queue->getNext() ;
						//}
					//}

				//}

				//delete queue ;
				//queue = queue2 ;
				//queue2 = new GridQueue2() ;
				//#ifdef VERBOSE
				//printf("%d complex, %d simple\n", numComplex, numSimple) ;
				//#endif

				//if ( numSimple == 0 )
				//{
					//break ;
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::erodeShapeTopo( float thr, int wid )
		//{
			///* Faster version of erode2 using priority queue */
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//setDataAt( i, j, k, 1 ) ;
									//queue2->prepend( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;


			//// Perform erosion
			//// wid = MAX_ERODE ;
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				//#endif

				//// First,
				//// check for complex nodes in queue2
				//// move them from queue2 to queue3
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Check simple
					//if ( ! isSimple( ox, oy, oz ) )
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue3->prepend( ox, oy, oz ) ;
						//ele = queue2->remove() ;

						//numComplex ++ ;
					//}
					//else
					//{
						//ele = queue2->getNext() ;
					//}
				//}

				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//queue->add( gp, -score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
					//score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					//// Remove this simple node
					//setDataAt( ox, oy, oz, -1 ) ;
					//numSimple ++ ;
					//// printf("Highest score: %d\n", score) ;

					//// Move its neighboring unvisited node to queue2
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

					//// Find complex nodes in its 3x3 neighborhood
					//// move them to queue2
					//for ( i = -1 ; i < 2 ; i ++ )
						//for ( j = -1 ; j < 2 ; j ++ )
							//for ( k = -1 ; k < 2 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//// Check simple
								//if ( getDataAt( nx, ny, nz ) == curwid && ! isSimple( nx, ny, nz ) )
								//{
									//// Complex, set to next layer
									//setDataAt( nx, ny, nz, curwid + 1 ) ;
									//queue2->prepend( nx, ny, nz ) ;
									//numComplex ++ ;
								//}
							//}

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue
							///*
					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//queue->add( gp, -score ) ;
									//}
								//}
							//}
							//*/

				//}

				//#ifdef VERBOSE
				//printf("%d complex, %d simple\n", numComplex, numSimple) ;
				//#endif
				//if ( numSimple == 0 )
				//{
					//break ;
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
			//delete queue;

		//}


		//void Volume::erodeAtom( float thr, int wid, Volume* avol )
		//{
			///* Erode to atoms */
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//if ( avol->getDataAt(i,j,k) > 0 )
							//{
								//setDataAt( i, j, k, MAX_ERODE ) ;
							//}
							//else
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//setDataAt( i, j, k, 1 ) ;
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;


			//// Perform erosion
			//// wid = MAX_ERODE ;
			//printf("Start erosion to %d...\n", wid) ;
			//#endif

			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				//#endif

				//// First,
				//// check for complex nodes in queue2
				//// move them from queue2 to queue3
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Check simple
					//if ( ! isSimple( ox, oy, oz ) )
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue3->prepend( ox, oy, oz ) ;
						//ele = queue2->remove() ;

						//numComplex ++ ;
					//}
					//else
					//{
						//ele = queue2->getNext() ;
					//}
				//}

				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//queue->add( gp, -score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
					//score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					//// Remove this simple node
					//setDataAt( ox, oy, oz, -1 ) ;
					//numSimple ++ ;
					//// printf("Highest score: %d\n", score) ;

					//// Move its neighboring unvisited node to queue2
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

					//// Find complex nodes in its 3x3 neighborhood
					//// move them to queue2
					//for ( i = -1 ; i < 2 ; i ++ )
						//for ( j = -1 ; j < 2 ; j ++ )
							//for ( k = -1 ; k < 2 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//// Check simple
								//if ( getDataAt( nx, ny, nz ) == curwid && ! isSimple( nx, ny, nz ) )

								//{
									//// Complex, set to next layer
									//setDataAt( nx, ny, nz, curwid + 1 ) ;
									//queue2->prepend( nx, ny, nz ) ;
									//numComplex ++ ;
								//}
							//}

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue
							///*
					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//queue->add( gp, -score ) ;
									//}
								//}
							//}
							//*/

				//}
				//#ifdef VERBOSE
				//printf("%d complex, %d simple\n", numComplex, numSimple) ;
				//#endif

				//if ( numSimple == 0 )
				//{
						//break ;
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
			//delete queue;

		//}

		///* Thin the current volume while preserving voxels with values > highthr or <= lowthr in grayvol
		//*  Assuming the current volume has already been thresholded to 0/1
		//*/
		void Volume::curveSkeleton( Volume* grayvol, float lowthr, float highthr, Volume* svol )
		{
			int i, j, k ;
			// First, threshold the volume
			#ifdef VERBOSE
			printf("Thresholding the volume to -1/0...\n") ;
			#endif
			threshold( 0.5f, -1, 0 ) ;

			// Next, apply convergent erosion
			// by preserving: complex nodes, curve end-points, and sheet points

			// Next, initialize the linked queue
			#ifdef VERBOSE
			printf("Initializing queue...\n") ;
			#endif
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			GridQueue2* queue4 = new GridQueue2( ) ;
			PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							float v = (float)grayvol->getDataAt(i,j,k) ;
							if ( v <= lowthr || v > highthr || svol->getDataAt(i,j,k) > 0 )
							{
								setDataAt( i, j, k, MAX_ERODE ) ;
							}
							else
							{
								for ( int m = 0 ; m < 6 ; m ++ )
								{
									if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									{
										// setDataAt( i, j, k, 1 ) ;
										queue2->prepend( i, j, k ) ;
										break ;
									}
								}
							}
						}
					}
			int wid = MAX_ERODE ;
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;
			printf("Start erosion to %d...\n", wid) ;
			#endif


			// Perform erosion
			gridQueueEle* ele ;
			gridPoint* gp ;
			int ox, oy, oz ;
			int score ;
			Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			{
				scrvol->setDataAt( i, -1 ) ;
			}

	#ifdef  NOISE_DIS_HELIX
			Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	#endif

			for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			{
				// At the start of each iteration,
				// queue2 holds all the nodes for this layer
				// queue3 and queue are empty

				int numComplex = 0, numSimple = 0 ;
				#ifdef VERBOSE
				printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				#endif

				/*
				We first need to assign curwid + 1 to every node in this layer
				*/
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					if ( getDataAt(ox,oy,oz) == curwid )
					{
						ele = queue2->remove() ;
					}
					else
					{
						setDataAt(ox,oy,oz, curwid) ;
						ele = queue2->getNext() ;
					}
				}
				queue4->reset() ;
				ele = queue4->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					queue2->prepend(ox,oy,oz) ;
					ele = queue4->remove() ;
				}

				// Now queue2 holds all the nodes for this layer

	#ifdef NOISE_DIS_HELIX
				/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				queue2->reset() ;

				// First run
				int flag = 0 ;
				while ( ( ele = queue2->getNext() ) != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;
					if ( NOISE_DIS_HELIX <= 1 )
					{
						noisevol->setDataAt( ox, oy, oz, 0 ) ;
					}
					else
					{
						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) == 0 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								flag = 1 ;
								break ;
							}
						}
						if ( ! flag )
						{
							noisevol->setDataAt( ox, oy, oz, 0 ) ;
						}
					}
				}

				int cur, visited ;
				for ( cur = 1 ; cur < NOISE_DIS_HELIX ; cur ++ )
				{
					queue2->reset() ;
					int count = 0 ;
					visited = 0 ;

					while ( ( ele = queue2->getNext() ) != NULL )
					{
						ox = ele->x ;
						oy = ele->y ;
						oz = ele->z ;

						if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						{
							visited ++ ;
							continue ;
						}

						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								visited ++ ;
								count ++ ;
								break ;
							}
						}
					}

					if ( count == 0 )
					{
						break ;
					}
				}
				printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited ) ;


	#endif
				/* Commented out for debugging

				// First,
				// check for complex nodes in queue2
				// move them from queue2 to queue3
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Check simple
	#ifndef NOISE_DIS_HELIX
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#else
					if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue3->prepend( ox, oy, oz ) ;
						ele = queue2->remove() ;

						numComplex ++ ;
					}
					else
					{
						ele = queue2->getNext() ;
					}
				}
				*/

				// Next,
				// Compute score for each node left in queue2
				// move them into priority queue
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Compute score
					score = getNumPotComplex2( ox, oy, oz ) ;
					scrvol->setDataAt( ox, oy, oz, score ) ;

					// Push to queue
					gp = new gridPoint ;
					gp->x = ox ;
					gp->y = oy ;
					gp->z = oz ;
					// queue->add( gp, -score ) ;
					queue->add( gp, score ) ;

					ele = queue2->remove() ;
				}

				// Rename queue3 to be queue2,
				// Clear queue3
				// From now on, queue2 holds nodes of next level
				delete queue2 ;
				queue2 = queue3 ;
				queue3 = new GridQueue2( ) ;

				// Next, start priority queue iteration
				while ( ! queue->isEmpty() )
				{
					// Retrieve the node with the highest score
					queue->remove( gp, score ) ;
					ox = gp->x ;
					oy = gp->y ;
					oz = gp->z ;
					delete gp ;
	//				score = -score ;

					// Ignore the node
					// if it has been processed before
					// or it has an updated score
					if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					{
						continue ;
					}

					/* Commented out for debugging

					// Remove this simple node
					setDataAt( ox, oy, oz, -1 ) ;
					numSimple ++ ;
					// printf("Highest score: %d\n", score) ;
					*/

					/* Added for debugging */
					// Check simple
	#ifndef NOISE_DIS_HELIX
					// if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) )
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#else
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue4->prepend( ox, oy, oz ) ;
						numComplex ++ ;
					}
					else
					{
						setDataAt( ox, oy, oz, -1 ) ;
						numSimple ++ ;


					/* Adding ends */
						// Move its neighboring unvisited node to queue2
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) == 0 )
							{
								// setDataAt( nx, ny, nz, curwid + 1 ) ;
								queue2->prepend( nx, ny, nz ) ;
							}
						}

					}


					/* Commented out for debugging

					// Find complex nodes in its 3x3 neighborhood
					// move them to queue2
					for ( i = -1 ; i < 2 ; i ++ )
						for ( j = -1 ; j < 2 ; j ++ )
							for ( k = -1 ; k < 2 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								// Check simple
								if ( getDataAt( nx, ny, nz ) == curwid &&
									// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	#ifndef NOISE_DIS_HELIX
									( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	#else
									( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	#endif

								{
									// Complex, set to next layer
									setDataAt( nx, ny, nz, curwid + 1 ) ;
									queue2->prepend( nx, ny, nz ) ;
									numComplex ++ ;
								}
							}
					*/

					// Update scores for nodes in its 5x5 neighborhood
					// insert them back into priority queue

					for ( i = -2 ; i < 3 ;i ++ )
						for ( j = -2 ; j < 3 ; j ++ )
							for ( k = -2 ; k < 3 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								if ( getDataAt( nx, ny, nz ) == curwid )
								{
									// Compute score
									score = getNumPotComplex2( nx, ny, nz ) ;

									if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									{
										// printf("Update\n") ;
										scrvol->setDataAt( nx, ny, nz, score ) ;
										// Push to queue
										gp = new gridPoint ;
										gp->x = nx ;
										gp->y = ny ;
										gp->z = nz ;
										// queue->add( gp, -score ) ;
										queue->add( gp, score ) ;
									}
								}
							}


				}

				#ifdef VERBOSE
				printf("%d complex, %d simple\n", numComplex, numSimple) ;
				#endif

				if ( numSimple == 0 )
				{
					if ( queue2->getNumElements() > 0 )
					{
						printf("*************************wierd here*************************\n");
					}
						break ;
				}
			}

			// Remove all internal voxels (contained in manifold surfaces)
			queue2->reset() ;
			queue4->reset() ;
			ele = queue4->getNext() ;
			while ( ele != NULL )
			{
				ox = ele->x ;
				oy = ele->y ;
				oz = ele->z ;

				if ( isPiercable(ox,oy,oz) == 1 )  // hasCompleteSheet( ox, oy, oz ) == 1 ) //
				{
					queue2->prepend(ox,oy,oz) ;
				//	setDataAt( ox, oy, oz, -1 ) ;
				}
				ele = queue4->remove() ;
			}

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) == 0 && isPiercable(i,j,k) ) //hasCompleteSheet(i,j,k) == 1) //
						{
							queue2->prepend( i, j, k ) ;
						}
					}
			queue2->reset() ;
			ele = queue2->getNext() ;
			while ( ele != NULL )
			{
				ox = ele->x ;
				oy = ele->y ;
				oz = ele->z ;
				setDataAt( ox, oy, oz, -1 ) ;
				ele = queue2->remove() ;
			}


			// Finally, clean up
			delete scrvol;
			delete queue;
			delete queue2;
			delete queue3;
			delete queue4;
			#ifdef VERBOSE
			printf("Thresholding the volume to 0/1...\n") ;
			#endif
			threshold( 0, 0, 1 ) ;
		}

		// Compute curve skeleton
		void Volume::curveSkeleton( float thr, Volume* svol )
		{
			int i, j, k ;
			// First, threshold the volume
			#ifdef VERBOSE
			printf("Thresholding the volume to -1/0...\n") ;
			#endif
			threshold( thr, -1, 0 ) ;

			// Next, apply convergent erosion
			// by preserving: complex nodes, curve end-points, and sheet points

			// Next, initialize the linked queue
			#ifdef VERBOSE
			printf("Initializing queue...\n") ;
			#endif
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			GridQueue2* queue4 = new GridQueue2( ) ;
			PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							if ( svol->getDataAt(i,j,k) > 0 )
							{
								setDataAt( i, j, k, MAX_ERODE ) ;
							}
							else
							{
								for ( int m = 0 ; m < 6 ; m ++ )
								{
									if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									{
										// setDataAt( i, j, k, 1 ) ;
										queue2->prepend( i, j, k ) ;
										break ;
									}
								}
							}
						}
					}

			int wid = MAX_ERODE ;
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;
			printf("Start erosion to %d...\n", wid) ;
			#endif


			// Perform erosion
			gridQueueEle* ele ;
			gridPoint* gp ;
			int ox, oy, oz ;
			int score ;
			Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			{
				scrvol->setDataAt( i, -1 ) ;
			}

	#ifdef  NOISE_DIS_HELIX
			Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	#endif

			for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			{
				// At the start of each iteration,
				// queue2 holds all the nodes for this layer
				// queue3 and queue are empty

				int numComplex = 0, numSimple = 0 ;
				#ifdef VERBOSE
				printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				#endif

				/*
				We first need to assign curwid + 1 to every node in this layer
				*/
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					if ( getDataAt(ox,oy,oz) == curwid )
					{
						ele = queue2->remove() ;
					}
					else
					{
						setDataAt(ox,oy,oz, curwid) ;
						ele = queue2->getNext() ;
					}
				}
				queue4->reset() ;
				ele = queue4->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					queue2->prepend(ox,oy,oz) ;
					ele = queue4->remove() ;
				}

				// Now queue2 holds all the nodes for this layer

	#ifdef NOISE_DIS_HELIX
				/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				queue2->reset() ;

				// First run
				int flag = 0 ;
				while ( ( ele = queue2->getNext() ) != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;
					if ( NOISE_DIS_HELIX <= 1 )
					{
						noisevol->setDataAt( ox, oy, oz, 0 ) ;
					}
					else
					{
						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) == 0 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								flag = 1 ;
								break ;
							}
						}
						if ( ! flag )
						{
							noisevol->setDataAt( ox, oy, oz, 0 ) ;
						}
					}
				}

				int cur, visited ;
				for ( cur = 1 ; cur < NOISE_DIS_HELIX ; cur ++ )
				{
					queue2->reset() ;
					int count = 0 ;
					visited = 0 ;

					while ( ( ele = queue2->getNext() ) != NULL )
					{
						ox = ele->x ;
						oy = ele->y ;
						oz = ele->z ;

						if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						{
							visited ++ ;
							continue ;
						}

						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								visited ++ ;
								count ++ ;
								break ;
							}
						}
					}

					if ( count == 0 )
					{
						break ;
					}
				}
				printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited ) ;


	#endif
				/* Commented out for debugging

				// First,
				// check for complex nodes in queue2
				// move them from queue2 to queue3
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Check simple
	#ifndef NOISE_DIS_HELIX
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#else
					if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue3->prepend( ox, oy, oz ) ;
						ele = queue2->remove() ;

						numComplex ++ ;
					}
					else
					{
						ele = queue2->getNext() ;
					}
				}
				*/

				// Next,
				// Compute score for each node left in queue2
				// move them into priority queue
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Compute score
					score = getNumPotComplex2( ox, oy, oz ) ;
					scrvol->setDataAt( ox, oy, oz, score ) ;

					// Push to queue
					gp = new gridPoint ;
					gp->x = ox ;
					gp->y = oy ;
					gp->z = oz ;
					// queue->add( gp, -score ) ;
					queue->add( gp, score ) ;

					ele = queue2->remove() ;
				}

				// Rename queue3 to be queue2,
				// Clear queue3
				// From now on, queue2 holds nodes of next level
				delete queue2 ;
				queue2 = queue3 ;
				queue3 = new GridQueue2( ) ;

				// Next, start priority queue iteration
				while ( ! queue->isEmpty() )
				{
					// Retrieve the node with the highest score
					queue->remove( gp, score ) ;
					ox = gp->x ;
					oy = gp->y ;
					oz = gp->z ;
					delete gp ;
	//				score = -score ;

					// Ignore the node
					// if it has been processed before
					// or it has an updated score
					if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					{
						continue ;
					}

					/* Commented out for debugging

					// Remove this simple node
					setDataAt( ox, oy, oz, -1 ) ;
					numSimple ++ ;
					// printf("Highest score: %d\n", score) ;
					*/

					/* Added for debugging */
					// Check simple
	#ifndef NOISE_DIS_HELIX
					// if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) )
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#else
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue4->prepend( ox, oy, oz ) ;
						numComplex ++ ;
					}
					else
					{
						setDataAt( ox, oy, oz, -1 ) ;
						numSimple ++ ;
					}
					/* Adding ends */

					// Move its neighboring unvisited node to queue2
					for ( int m = 0 ; m < 6 ; m ++ )
					{
						int nx = ox + neighbor6[m][0] ;
						int ny = oy + neighbor6[m][1] ;
						int nz = oz + neighbor6[m][2] ;
						if ( getDataAt( nx, ny, nz ) == 0 )
						{
							// setDataAt( nx, ny, nz, curwid + 1 ) ;
							queue2->prepend( nx, ny, nz ) ;
						}
					}

					/* Commented out for debugging

					// Find complex nodes in its 3x3 neighborhood
					// move them to queue2
					for ( i = -1 ; i < 2 ; i ++ )
						for ( j = -1 ; j < 2 ; j ++ )
							for ( k = -1 ; k < 2 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								// Check simple
								if ( getDataAt( nx, ny, nz ) == curwid &&
									// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	#ifndef NOISE_DIS_HELIX
									( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	#else
									( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	#endif

								{
									// Complex, set to next layer
									setDataAt( nx, ny, nz, curwid + 1 ) ;
									queue2->prepend( nx, ny, nz ) ;
									numComplex ++ ;
								}
							}
					*/

					// Update scores for nodes in its 5x5 neighborhood
					// insert them back into priority queue

					for ( i = -2 ; i < 3 ;i ++ )
						for ( j = -2 ; j < 3 ; j ++ )
							for ( k = -2 ; k < 3 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								if ( getDataAt( nx, ny, nz ) == curwid )
								{
									// Compute score
									score = getNumPotComplex2( nx, ny, nz ) ;

									if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									{
										// printf("Update\n") ;
										scrvol->setDataAt( nx, ny, nz, score ) ;
										// Push to queue
										gp = new gridPoint ;
										gp->x = nx ;
										gp->y = ny ;
										gp->z = nz ;
										// queue->add( gp, -score ) ;
										queue->add( gp, score ) ;
									}
								}
							}


				}

				#ifdef VERBOSE
				printf("%d complex, %d simple\n", numComplex, numSimple) ;
				#endif

				if ( numSimple == 0 )
				{
						break ;
				}
			}

			// Finally, clean up
			delete scrvol;
			delete queue;
			delete queue2;
			delete queue3;
			delete queue4;
			#ifdef VERBOSE
			printf("Thresholding the volume to 0/1...\n") ;
			#endif
			threshold( 0, 0, 1 ) ;
		}

		// Compute curve skeleton in 2D
		void Volume::curveSkeleton2D( float thr, Volume* svol )
		{
			int i, j, k ;
			// First, threshold the volume
			#ifdef VERBOSE
			printf("Thresholding the volume to -1/0...\n") ;
			#endif
			threshold( thr, -1, 0 ) ;

			// Next, apply convergent erosion
			// by preserving: complex nodes, curve end-points, and sheet points

			// Next, initialize the linked queue
			#ifdef VERBOSE
			printf("Initializing queue...\n") ;
			#endif
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			GridQueue2* queue4 = new GridQueue2( ) ;
			PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							if ( svol->getDataAt(i,j,k) > 0 )
							{
								setDataAt( i, j, k, MAX_ERODE ) ;
							}
							else
							{
								for ( int m = 0 ; m < 4 ; m ++ )
								{
									if ( getDataAt( i + neighbor4[m][0], j + neighbor4[m][1], k ) < 0 )
									{
										// setDataAt( i, j, k, 1 ) ;
										queue2->prepend( i, j, k ) ;
										break ;
									}
								}
							}
						}
					}
			int wid = MAX_ERODE ;
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;
			printf("Start erosion to %d...\n", wid) ;
			#endif


			// Perform erosion
			gridQueueEle* ele ;
			gridPoint* gp ;
			int ox, oy, oz ;
			int score ;
			Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			{
				scrvol->setDataAt( i, -1 ) ;
			}

	#ifdef  NOISE_DIS_HELIX
			Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	#endif

			for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			{
				// At the start of each iteration,
				// queue2 holds all the nodes for this layer
				// queue3 and queue are empty

				int numComplex = 0, numSimple = 0 ;
				#ifdef VERBOSE
				printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				#endif

				/*
				We first need to assign curwid + 1 to every node in this layer
				*/
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					if ( getDataAt(ox,oy,oz) == curwid )
					{
						ele = queue2->remove() ;
					}
					else
					{
						setDataAt(ox,oy,oz, curwid) ;
						ele = queue2->getNext() ;
					}
				}
				queue4->reset() ;
				ele = queue4->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					queue2->prepend(ox,oy,oz) ;
					ele = queue4->remove() ;
				}

				// Now queue2 holds all the nodes for this layer

				#ifdef NOISE_DIS_HELIX
				/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				queue2->reset() ;

				// First run
				int flag = 0 ;
				while ( ( ele = queue2->getNext() ) != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;
					if ( NOISE_DIS_HELIX <= 1 )
					{
						noisevol->setDataAt( ox, oy, oz, 0 ) ;
					}
					else
					{
						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) == 0 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								flag = 1 ;
								break ;
							}
						}
						if ( ! flag )
						{
							noisevol->setDataAt( ox, oy, oz, 0 ) ;
						}
					}
				}

				int cur, visited ;
				for ( cur = 1 ; cur < NOISE_DIS_HELIX ; cur ++ )
				{
					queue2->reset() ;
					int count = 0 ;
					visited = 0 ;

					while ( ( ele = queue2->getNext() ) != NULL )
					{
						ox = ele->x ;
						oy = ele->y ;
						oz = ele->z ;

						if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						{
							visited ++ ;
							continue ;
						}

						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								visited ++ ;
								count ++ ;
								break ;
							}
						}
					}

					if ( count == 0 )
					{
						break ;
					}
				}
				printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited ) ;


				#endif
				/* Commented out for debugging

				// First,
				// check for complex nodes in queue2
				// move them from queue2 to queue3
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Check simple
	#ifndef NOISE_DIS_HELIX
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#else
					if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue3->prepend( ox, oy, oz ) ;
						ele = queue2->remove() ;

						numComplex ++ ;
					}
					else
					{
						ele = queue2->getNext() ;
					}
				}
				*/

				// Next,
				// Compute score for each node left in queue2
				// move them into priority queue
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Compute score
					score = getNumPotComplex2( ox, oy, oz ) ;
					//score = getNumNeighbor6( ox, oy, oz ) ;
					scrvol->setDataAt( ox, oy, oz, score ) ;

					// Push to queue
					gp = new gridPoint ;
					gp->x = ox ;
					gp->y = oy ;
					gp->z = oz ;
					// queue->add( gp, -score ) ;
					queue->add( gp, score ) ;

					ele = queue2->remove() ;
				}

				// Rename queue3 to be queue2,
				// Clear queue3
				// From now on, queue2 holds nodes of next level
				delete queue2 ;
				queue2 = queue3 ;
				queue3 = new GridQueue2( ) ;

				// Next, start priority queue iteration
				while ( ! queue->isEmpty() )
				{
					// Retrieve the node with the highest score
					queue->remove( gp, score ) ;
					ox = gp->x ;
					oy = gp->y ;
					oz = gp->z ;
					delete gp ;
	//				score = -score ;

					// Ignore the node
					// if it has been processed before
					// or it has an updated score
					if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					{
						continue ;
					}

					/* Commented out for debugging

					// Remove this simple node
					setDataAt( ox, oy, oz, -1 ) ;
					numSimple ++ ;
					// printf("Highest score: %d\n", score) ;
					*/

					/* Added for debugging */
					// Check simple
					#ifndef NOISE_DIS_HELIX
					// if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) )
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
					#else
					if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
					#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue4->prepend( ox, oy, oz ) ;
						numComplex ++ ;
					}
					else
					{
						setDataAt( ox, oy, oz, -1 ) ;
						numSimple ++ ;
					}
					/* Adding ends */

					// Move its neighboring unvisited node to queue2
					for ( int m = 0 ; m < 4 ; m ++ )
					{
						int nx = ox + neighbor4[m][0] ;
						int ny = oy + neighbor4[m][1] ;
						int nz = oz ;
						if ( getDataAt( nx, ny, nz ) == 0 )
						{
							// setDataAt( nx, ny, nz, curwid + 1 ) ;
							queue2->prepend( nx, ny, nz ) ;
						}
					}

					/* Commented out for debugging

					// Find complex nodes in its 3x3 neighborhood
					// move them to queue2
					for ( i = -1 ; i < 2 ; i ++ )
						for ( j = -1 ; j < 2 ; j ++ )
							for ( k = -1 ; k < 2 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								// Check simple
								if ( getDataAt( nx, ny, nz ) == curwid &&
									// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	#ifndef NOISE_DIS_HELIX
									( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	#else
									( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	#endif

								{
									// Complex, set to next layer
									setDataAt( nx, ny, nz, curwid + 1 ) ;
									queue2->prepend( nx, ny, nz ) ;
									numComplex ++ ;
								}
							}
					*/

					// Update scores for nodes in its 5x5 neighborhood
					// insert them back into priority queue

					for ( i = -2 ; i < 3 ;i ++ )
						for ( j = -2 ; j < 3 ; j ++ )
							for ( k = -2 ; k < 3 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								if ( getDataAt( nx, ny, nz ) == curwid )
								{
									// Compute score
									score = getNumPotComplex2( nx, ny, nz ) ;
									//score = getNumNeighbor6( nx, ny, nz ) ;

									if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									{
										// printf("Update\n") ;
										scrvol->setDataAt( nx, ny, nz, score ) ;
										// Push to queue
										gp = new gridPoint ;
										gp->x = nx ;
										gp->y = ny ;
										gp->z = nz ;
										// queue->add( gp, -score ) ;
										queue->add( gp, score ) ;
									}
								}
							}


				}

				#ifdef VERBOSE
				printf("%d complex, %d simple\n", numComplex, numSimple) ;
				#endif

				if ( numSimple == 0 )
				{
						break ;
				}
			}

			// Finally, clean up
			#ifdef VERBOSE
			printf("Thresholding the volume to 0/1...\n") ;
			#endif
			threshold( 0, 0, 1 ) ;
			delete scrvol;
			delete queue;
			delete queue2;
			delete queue3;
			delete queue4;
		}

		// Compute minimal skeleton
		void Volume::skeleton( float thr, int )
		{
			int i, j, k ;
			// First, threshold the volume
			#ifdef VERBOSE
			printf("Thresholding the volume to -1/0...\n") ;
			#endif
			threshold( thr, -1, 0 ) ;

			// Next, apply convergent erosion
			// by preserving: complex nodes, curve end-points, and sheet points

			// Next, initialize the linked queue
			#ifdef VERBOSE
			printf("Initializing queue...\n") ;
			#endif
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue = new GridQueue2( ) ;

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							{
								for ( int m = 0 ; m < 6 ; m ++ )
								{
									if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									{
										setDataAt( i, j, k, 1 ) ;
										queue2->prepend( i, j, k ) ;
										break ;
									}
								}
							}
						}
					}
			int wid = 0 ;
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;
			printf("Start erosion to %d...\n", wid) ;
			#endif

			// Perform erosion
			gridQueueEle* ele ;
			int ox, oy, oz ;

			while( 1 )
			{
				wid ++ ;

				// At the start of each iteration,
				// queue2 holds all the nodes for this layer
				// queue is empty

				int numComplex = 0, numSimple = 0 ;
				#ifdef VERBOSE
				printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), wid) ;
				#endif

				// Rename queue2 to be queue,
				// Clear queue2
				// From now on, queue2 holds nodes of next level
				delete queue ;
				queue = queue2 ;
				queue2 = new GridQueue2( ) ;

				// Next, start queue iteration
				queue->reset() ;
				ele = queue->getNext();
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;
	//				delete ele ;

					// Check simple
					if ( ! isSimple( ox, oy, oz ) )
					{
						// Complex, set to next layer
						queue2->prepend( ox, oy, oz ) ;
						numComplex ++ ;
					}
					/*
					else if ( ox == off || oy == off || oz == off ||
						ox == getSizeX() - off - 1 || oy == getSizeY() - off - 1 || oz == getSizeZ() - off - 1 )
					{
						// Wall, don't erode, set to next layer
						queue2->prepend( ox, oy, oz ) ;
						numComplex ++ ;
					}
					*/
					else
					{
						setDataAt( ox, oy, oz, -1 ) ;
						numSimple ++ ;

						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) == 0 )
							{
								setDataAt( nx, ny, nz, 1 ) ;
								queue2->prepend( nx, ny, nz ) ;
							}
						}

					}

					ele = queue->remove() ;
				}
				#ifdef VERBOSE
				printf("Level %d: %d complex, %d simple\n", wid, numComplex, numSimple) ;
				#endif

				if ( numSimple == 0 )
				{
					break ;
				}
			}

			// Finally, clean up
			delete queue;
			delete queue2;
			#ifdef VERBOSE
			printf("Thresholding the volume to 0/1...\n") ;
			#endif
			threshold( 0, 0, 1 ) ;
		}

		//// Compute minimal skeleton
		//void Volume::skeleton2( float thr, int off )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -1/0...\n") ;
			//#endif
			//threshold( thr, -1, 0 ) ;

			//// Next, apply convergent erosion
			//// by preserving: complex nodes, curve end-points, and sheet points

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//if ( i == off || j == off || k == off ||
								 //i == getSizeX() - off - 1 || j == getSizeY() - off - 1 || k == getSizeZ() - off - 1 )
							//{
								//setDataAt( i, j, k, MAX_ERODE ) ;
							//}
							//else
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//setDataAt( i, j, k, 1 ) ;
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}
			//int wid = MAX_ERODE ;
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;


			//// Perform erosion
			//printf("Start erosion to %d...\n", wid) ;
			//#endif

			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;

			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}


			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				//#endif

				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex2( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//// queue->add( gp, -score ) ;
					//queue->add( gp, score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//// From now on, queue2 holds nodes of next level
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
	////				score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					///* Added for debugging */
					//// Check simple
					//if ( ! isSimple( ox, oy, oz ) )
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue2->prepend( ox, oy, oz ) ;
						//numComplex ++ ;
					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//numSimple ++ ;
					//}
					///* Adding ends */

					//// Move its neighboring unvisited node to queue2
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue
						///*
					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex2( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//// queue->add( gp, -score ) ;
										//queue->add( gp, score ) ;
									//}
								//}
							//}

						//*/
				//}
				//#ifdef VERBOSE
				//printf("%d complex, %d simple\n", numComplex, numSimple) ;
				//#endif

				//if ( numSimple == 0 )
				//{
						//break ;
				//}
			//}

			//// Finally, clean up
			//delete scrvol ;
			//delete queue;
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
		//}

		///* Thin the current volume while preserving voxels with values > highthr or <= lowthr in grayvol
		//*  Assuming the current volume has already been thresholded to 0/1
		//*/
		//void Volume::pointSkeleton( Volume* grayvol, float lowthr, float highthr, Volume* svol, Volume* hvol )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -1/0...\n") ;
			//#endif
			//threshold( 0.5f, -1, 0 ) ;

			//// Next, apply convergent erosion
			//// by preserving: complex nodes, curve end-points, and sheet points

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//float v = (float)grayvol->getDataAt( i, j, k ) ;
							//if ( v <= lowthr || v > highthr || svol->getDataAt(i,j,k) > 0 || hvol->getDataAt(i,j,k) > 0 )
							//{
								//setDataAt( i, j, k, MAX_ERODE ) ;
							//}
							//else
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//setDataAt( i, j, k, 1 ) ;
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;
			//#endif


			//// Perform erosion
			//int wid = MAX_ERODE ;
			//#ifdef VERBOSE
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}


			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				//#endif


				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex2( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//// queue->add( gp, -score ) ;
					//queue->add( gp, score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//// From now on, queue2 holds nodes of next level
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
	////				score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					///* Added for debugging */
					//// Check simple
					//if ( ! isSimple( ox, oy, oz ) )
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue2->prepend( ox, oy, oz ) ;
						//numComplex ++ ;
					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//numSimple ++ ;
					///* Adding ends */

						//// Move its neighboring unvisited node to queue2
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//setDataAt( nx, ny, nz, curwid + 1 ) ;
								//queue2->prepend( nx, ny, nz ) ;
							//}
						//}

					//}

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue

					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex2( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//// queue->add( gp, -score ) ;
										//queue->add( gp, score ) ;
									//}
								//}
							//}


				//}
				//#ifdef VERBOSE
				//printf("%d complex, %d simple\n", numComplex, numSimple) ;
				//#endif

				//if ( numSimple == 0 )
				//{
						//break ;
				//}
			//}

			//// Remove all internal voxels (contained in manifold curves)
			//queue2->reset() ;
			//ele = queue2->getNext() ;
			//while ( ele != NULL )
			//{
				//ox = ele->x ;
				//oy = ele->y ;
				//oz = ele->z ;

				//if ( hasCompleteHelix( ox, oy, oz ) == 1 )
				//{
					//ele = queue2->getNext() ;
				//}
				//else
				//{
					//ele = queue2->remove() ;
				//}
			//}

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) == 0 && hasCompleteHelix(i,j,k) == 1)
						//{
							//queue2->prepend( i, j, k ) ;
						//}
					//}
			//queue2->reset() ;
			//ele = queue2->getNext() ;
			//while ( ele != NULL )
			//{
				//ox = ele->x ;
				//oy = ele->y ;
				//oz = ele->z ;
				//setDataAt( ox, oy, oz, -1 ) ;
				//ele = queue2->remove() ;
			//}

			//// Finally, clean up
			//delete scrvol;
			//delete queue;
			//delete queue2;
			//delete queue3;
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
		//}


		// Compute minimal skeleton
		void Volume::skeleton( float thr, Volume* svol, Volume* hvol )
		{
			int i, j, k ;
			// First, threshold the volume
			#ifdef VERBOSE
			printf("Thresholding the volume to -1/0...\n") ;
			#endif
			threshold( thr, -1, 0 ) ;

			// Next, apply convergent erosion
			// by preserving: complex nodes, curve end-points, and sheet points

			// Next, initialize the linked queue
			#ifdef VERBOSE
			printf("Initializing queue...\n") ;
			#endif
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							if ( svol->getDataAt(i,j,k) > 0 || hvol->getDataAt(i,j,k) > 0 )
							{
								setDataAt( i, j, k, MAX_ERODE ) ;
							}
							else
							{
								for ( int m = 0 ; m < 6 ; m ++ )
								{
									if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									{
										setDataAt( i, j, k, 1 ) ;
										queue2->prepend( i, j, k ) ;
										break ;
									}
								}
							}
						}
					}
			int wid = MAX_ERODE ;
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;


			// Perform erosion
			printf("Start erosion to %d...\n", wid) ;
			#endif
			gridQueueEle* ele ;
			gridPoint* gp ;
			int ox, oy, oz ;
			int score ;
			Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			{
				scrvol->setDataAt( i, -1 ) ;
			}


			for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			{
				// At the start of each iteration,
				// queue2 holds all the nodes for this layer
				// queue3 and queue are empty

				int numComplex = 0, numSimple = 0 ;
				#ifdef VERBOSE
				printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				#endif


				// Next,
				// Compute score for each node left in queue2
				// move them into priority queue
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Compute score
					score = getNumPotComplex2( ox, oy, oz ) ;
					scrvol->setDataAt( ox, oy, oz, score ) ;

					// Push to queue
					gp = new gridPoint ;
					gp->x = ox ;
					gp->y = oy ;
					gp->z = oz ;
					// queue->add( gp, -score ) ;
					queue->add( gp, score ) ;

					ele = queue2->remove() ;
				}

				// Rename queue3 to be queue2,
				// Clear queue3
				// From now on, queue2 holds nodes of next level
				delete queue2 ;
				queue2 = queue3 ;
				queue3 = new GridQueue2( ) ;

				// Next, start priority queue iteration
				while ( ! queue->isEmpty() )
				{
					// Retrieve the node with the highest score
					queue->remove( gp, score ) ;
					ox = gp->x ;
					oy = gp->y ;
					oz = gp->z ;
					delete gp ;
	//				score = -score ;

					// Ignore the node
					// if it has been processed before
					// or it has an updated score
					if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					{
						continue ;
					}

					/* Added for debugging */
					// Check simple
					if ( ! isSimple( ox, oy, oz ) )
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue2->prepend( ox, oy, oz ) ;
						numComplex ++ ;
					}
					else
					{
						setDataAt( ox, oy, oz, -1 ) ;
						numSimple ++ ;
					}
					/* Adding ends */

					// Move its neighboring unvisited node to queue2
					for ( int m = 0 ; m < 6 ; m ++ )
					{
						int nx = ox + neighbor6[m][0] ;
						int ny = oy + neighbor6[m][1] ;
						int nz = oz + neighbor6[m][2] ;
						if ( getDataAt( nx, ny, nz ) == 0 )
						{
							setDataAt( nx, ny, nz, curwid + 1 ) ;
							queue2->prepend( nx, ny, nz ) ;
						}
					}

					// Update scores for nodes in its 5x5 neighborhood
					// insert them back into priority queue

					for ( i = -2 ; i < 3 ;i ++ )
						for ( j = -2 ; j < 3 ; j ++ )
							for ( k = -2 ; k < 3 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								if ( getDataAt( nx, ny, nz ) == curwid )
								{
									// Compute score
									score = getNumPotComplex2( nx, ny, nz ) ;

									if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									{
										// printf("Update\n") ;
										scrvol->setDataAt( nx, ny, nz, score ) ;
										// Push to queue
										gp = new gridPoint ;
										gp->x = nx ;
										gp->y = ny ;
										gp->z = nz ;
										// queue->add( gp, -score ) ;
										queue->add( gp, score ) ;
									}
								}
							}


				}

				#ifdef VERBOSE
				printf("%d complex, %d simple\n", numComplex, numSimple) ;
				#endif

				if ( numSimple == 0 )
				{
					delete queue2 ;
					break ;
				}
			}

			// Finally, clean up
			#ifdef VERBOSE
			printf("Thresholding the volume to 0/1...\n") ;
			#endif
			threshold( 0, 0, 1 ) ;
			delete scrvol;
			delete queue;
			delete queue3;
		}


		// Apply helix erosion
		void Volume::erodeHelix( )
		{
			erodeHelix( 3 ) ;
		}


		void Volume::erodeHelix( int disthr )
		{
			int i, j, k ;
			// First, threshold the volume
			//printf("Thresholding the volume to -1/0...\n") ;
			threshold( 0.1f, -1, 0 ) ;

			/* Debug: remove faces */
			//Volume* facevol = markFaceEdge() ;
			/* End debugging */

			// Next, initialize the linked queue
			//printf("Initializing queue...\n") ;
			Volume * fvol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			GridQueue2** queues = new GridQueue2* [ 10000 ] ;

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							if ( ! hasCompleteHelix( i, j, k ) )
							// if ( ! hasCompleteHelix( i, j, k, facevol ) )
							{
								queue2->prepend( i, j, k ) ;
								fvol->setDataAt( i, j, k, -1 ) ;
							}
						}
					}
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;

			// Start erosion
			gridQueueEle* ele ;
			int dis = -1 ;
			while ( queue2->getNumElements() > 0 )
			{
				// First, set distance
				dis -- ;
				queues[ - dis ] = new GridQueue2( ) ;
				//printf("Distance transform to %d...", dis) ;
				queue2->reset() ;
				while( ( ele = queue2->getNext() ) != NULL )
				{
					setDataAt( ele->x, ele->y, ele->z, dis ) ;
					queues[ -dis ]->prepend( ele->x, ele->y, ele->z ) ;
				}
				//printf("%d nodes\n", queues[-dis]->getNumElements()) ;

				// Next, find next layer
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					for ( int m = 0 ; m < 6 ; m ++ )
					{
						int nx = ele->x + neighbor6[m][0] ;
						int ny = ele->y + neighbor6[m][1] ;
						int nz = ele->z + neighbor6[m][2] ;
						if ( getDataAt( nx, ny, nz ) == 0 )
						{
							fvol->setDataAt( nx, ny, nz, dis ) ;

							if ( ! hasCompleteHelix( nx, ny, nz ) )
							// if ( ! hasCompleteHelix( nx, ny, nz, facevol ) )
							{
								setDataAt( nx, ny, nz, 1 ) ;
								queue3->prepend( nx, ny, nz ) ;
							}
						}
					}

					ele = queue2->remove() ;
				}

				// Next, swap queues
				GridQueue2* temp = queue2 ;
				queue2 = queue3 ;
				queue3 = temp ;
			}

			// Deal with closed rings
			dis -- ;
			queues[ - dis ] = new GridQueue2( ) ;
			#ifdef VERBOSE
			printf("Detecting closed rings %d...", dis) ;
			#endif
			int ftot = 0 ;
			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) == 0 )
						{
							setDataAt( i, j, k, dis ) ;
							queues[ -dis ]->prepend( i, j, k ) ;
	/*
							int fval = (int) fvol->getDataAt( i, j, k ) ;
							if ( fval == 0)
							{
								// queues[ -dis ]->prepend( i, j, k ) ;
							}
							else
							{
								setDataAt( i, j, k, fval - 1 ) ;
								// queues[ -fval + 1 ]->prepend( i, j, k ) ;
							}
	*/
							ftot ++ ;
						}
					}
			#ifdef VERBOSE
			printf("%d nodes\n", ftot) ;
			#endif


			// return ;

			/* Find local minimum: to help determine erosion level
			int cts[ 64 ] ;
			for ( i = 0 ; i <= - dis ; i ++ )
			{
				cts[ i ] = 0 ;
			}
			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						double val = getDataAt( i, j, k ) ;
						if ( val < -1 )
						{
							int min = 1 ;
							for ( int m = 0 ; m < 6 ; m ++ )
							{
								int nx = i + neighbor6[m][0] ;
								int ny = j + neighbor6[m][1] ;
								int nz = k + neighbor6[m][2] ;
								if ( getDataAt( nx, ny, nz ) < val )
								{
									min = 0 ;
									break ;
								}
							}

							if ( min )
							{
								cts[ (int) - val ] ++ ;
							}
						}
					}

			for ( i = 2 ; i <= - dis ; i ++ )
			{
				printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
			}
			*/

			// Dilate back
			// Starting from nodes with distance - 2 - disthr

			if ( disthr + 2 > - dis )
			{
				disthr = - dis - 2 ;
			}
			int d;

			for ( d = - dis ; d > disthr + 1 ; d -- )
			{
				queues[ d ]->reset() ;
				while ( (ele = queues[ d ]->getNext() ) != NULL )
				{
					setDataAt( ele->x, ele->y, ele->z, d ) ;
				}
			}


			for ( d = disthr + 1 ; d >= 2 ; d -- )
			{
				//delete queue3 ;
				//queue3 = new GridQueue2( ) ;
				queues[ d ]->reset() ;
				ele = queues[ d ]->getNext() ;
				while ( ele != NULL )
				{
					int dilatable = 0 ;
					for ( int m = 0 ; m < 6 ; m ++ )
							{
								int nx = ele->x + neighbor6[m][0] ;
								int ny = ele->y + neighbor6[m][1] ;
								int nz = ele->z + neighbor6[m][2] ;
								if ( getDataAt( nx, ny, nz ) == d + 1 )
								{
									dilatable = 1 ;
									break ;
								}
							}


					if ( ! dilatable )
					{
						/*
						setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
						*/

						setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
						if ( d > 2 )
						{
							queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
						}
						ele = queues[ d ]->remove() ;
					}
					else
					{
						setDataAt( ele->x, ele->y, ele->z, d ) ;
						ele = queues[ d ]->getNext() ;
					}

				}

				/* Debug: extract minimal set */
				while ( 1 )
				{
					int num = 0 ;
					queues[ d ]->reset() ;
					ele = queues[ d ]->getNext() ;
					while ( ele != NULL )
					{
						int critical = 0 ;
						setDataAt( ele->x, ele->y, ele->z, -1 ) ;

						for ( i = -1 ; i < 2 ; i ++ )
						{
							for ( j = -1 ; j < 2 ; j ++ )
							{
								for ( k = -1 ; k < 2 ; k ++ )
								{
									if ( i != 0 && j != 0 && k != 0 )
									{
										continue ;
									}
									int nx = ele->x + i ;
									int ny = ele->y + j ;
									int nz = ele->z + k ;
									if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteHelix( nx,ny,nz ) ) //, facevol ) )
									{
										critical = 1 ;
										break ;
									}
								}
								if ( critical )
								{
									break ;
								}
							}
							if ( critical )
							{
								break ;
							}
						}

						if ( critical )
						{
							setDataAt( ele->x, ele->y, ele->z, d ) ;
							ele = queues[ d ]->getNext() ;
						}
						else
						{
							setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
							if ( d > 2 )
							{
								queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
							}
							ele = queues[ d ]->remove() ;
							num ++ ;
						}

					}

					#ifdef VERBOSE
					printf("Non-minimal: %d\n", num) ;
					#endif

					if ( num == 0 )
					{
						break ;
					}
				}


				/* End of debugging */

				/*
				queue3->reset() ;
				ele = queue3->getNext() ;
				while ( ele != NULL )
				{
					setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
					ele = queue3->remove() ;
				}
				*/
			}

			// Finally, threshold the volume
			#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			#endif
			//threshold( -1, 1, 0, 0 ) ;
			threshold( 0, 0, 1 ) ;
			delete fvol ;
			delete queue2;
			delete queue3;
			for ( d = -dis ; d >= 2 ; d -- ) {
				delete queues[d];
			}
			delete [] queues;

		}



		// Apply sheet erosion
		int Volume::erodeSheet( )
		{
			return erodeSheet( 3 ) ;
		}


		int Volume::erodeSheet( int disthr )
		{
			int i, j, k ;
			// First, threshold the volume
			//printf("Thresholding the volume to -1/0...\n") ;
			threshold( 0.1f, -1, 0 ) ;

			/* Debug: remove cells */
			Volume* facevol = markCellFace() ;
			/* End debugging */

			// Next, initialize the linked queue
			//printf("Initializing queue...\n") ;
			Volume * fvol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			GridQueue2** queues = new GridQueue2* [ 10000 ] ;

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 )
						{
							if ( ! hasCompleteSheet( i, j, k ) )
							//if ( ! hasCompleteSheet( i, j, k, facevol ) )
							{
								queue2->prepend( i, j, k ) ;
								fvol->setDataAt( i, j, k, -1 ) ;
							}
						}
					}
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;
			#endif

			// Start erosion
			gridQueueEle* ele ;
			int dis = -1 ;
			while ( queue2->getNumElements() > 0 )
			{
				// First, set distance
				dis -- ;
				queues[ - dis ] = new GridQueue2( ) ;
				//printf("Distance transform to %d...", dis) ;
				queue2->reset() ;
				while( ( ele = queue2->getNext() ) != NULL )
				{
					setDataAt( ele->x, ele->y, ele->z, dis ) ;
					queues[ -dis ]->prepend( ele->x, ele->y, ele->z ) ;
				}
				//printf("%d nodes\n", queues[-dis]->getNumElements()) ;

				// Next, find next layer
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					// for ( int m = 0 ; m < 6 ; m ++ )
					for ( int mx = -1 ; mx < 2 ; mx ++ )
						for ( int my = -1 ; my < 2 ; my ++ )
							for ( int mz = -1 ; mz < 2 ; mz ++ )
							{
								if ( mx != 0 && my != 0 && mz != 0 )
								{
									continue ;
								}
								//int nx = ele->x + neighbor6[m][0] ;
								//int ny = ele->y + neighbor6[m][1] ;
								//int nz = ele->z + neighbor6[m][2] ;
								int nx = ele->x + mx ;
								int ny = ele->y + my ;
								int nz = ele->z + mz ;

								if ( getDataAt( nx, ny, nz ) == 0 )
								{
									fvol->setDataAt( nx, ny, nz, dis ) ;

									if  ( ! hasCompleteSheet( nx, ny, nz ) )
									// if  ( ! hasCompleteSheet( nx, ny, nz, facevol ) )
									{
										setDataAt( nx, ny, nz, 1 ) ;
										queue3->prepend( nx, ny, nz ) ;
									}
								}
							}

					ele = queue2->remove() ;
				}

				// Next, swap queues
				GridQueue2* temp = queue2 ;
				queue2 = queue3 ;
				queue3 = temp ;
			}

			/* Deal with closed rings */

			dis -- ;
			queues[ - dis ] = new GridQueue2( ) ;
			#ifdef VERBOSE
			printf("Detecting closed rings %d...", dis) ;
			#endif
			int ftot = 0 ;
			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) == 0 )
						{
							/*
							int fval = (int) fvol->getDataAt( i, j, k ) ;
							if ( fval == 0)
							{
								setDataAt( i, j, k, dis - 2 ) ;
								// queues[ -dis ]->prepend( i, j, k ) ;
							}
							else
							{
								setDataAt( i, j, k, fval - 1 ) ;
								queues[ -fval + 1 ]->prepend( i, j, k ) ;
							}
							*/
							setDataAt( i, j, k, dis ) ;
							queues[ -dis ]->prepend( i, j, k ) ;

							ftot ++ ;
						}
					}
			#ifdef VERBOSE
			printf("%d nodes\n", ftot) ;
			#endif


			/* Find local minimum: to help determine erosion level
			int cts[ 64 ] ;
			for ( i = 0 ; i <= - dis ; i ++ )
			{
				cts[ i ] = 0 ;
			}
			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						double val = getDataAt( i, j, k ) ;
						if ( val < -1 )
						{
							int min = 1 ;
							for ( int m = 0 ; m < 6 ; m ++ )
							{
								int nx = i + neighbor6[m][0] ;
								int ny = j + neighbor6[m][1] ;
								int nz = k + neighbor6[m][2] ;
								if ( getDataAt( nx, ny, nz ) < val )
								{
									min = 0 ;
									break ;
								}
							}

							if ( min )
							{
								cts[ (int) - val ] ++ ;
							}
						}
					}

			for ( i = 2 ; i <= - dis ; i ++ )
			{
				printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
			}
			*/

			// return ;

			// Dilate back
			// Starting from nodes with distance - 2 - disthr
			int d ;
			if ( disthr + 2 > - dis )
			{
				disthr = - dis - 2 ;

			}
			for ( d = - dis ; d > disthr + 1 ; d -- )
			{
				queues[ d ]->reset() ;
				while ( (ele = queues[ d ]->getNext() ) != NULL )
				{
					setDataAt( ele->x, ele->y, ele->z, d ) ;
				}
			}

			for (d = disthr + 1 ; d >= 2 ; d -- )
			{

				//delete queue3 ;
				//queue3 = new GridQueue2( ) ;
				queues[ d ]->reset() ;
				ele = queues[ d ]->getNext() ;
				while ( ele != NULL )
				{
					int dilatable = 0 ;
					// for ( int m = 0 ; m < 6 ; m ++ )
					/*
					for ( int mx = -1 ; mx < 2 ; mx ++ )
					for ( int my = -1 ; my < 2 ; my ++ )
					for ( int mz = -1 ; mz < 2 ; mz ++ )
					{
					if ( mx == 0 || my == 0 || mz == 0 )
					{
					int nx = ele->x + mx ; // neighbor6[m][0] ;
					int ny = ele->y + my ; // neighbor6[m][1] ;
					int nz = ele->z + mz ; // neighbor6[m][2] ;
					if ( getDataAt( nx, ny, nz ) == - d - 1 )
					{
					dilatable = 1 ;
					break ;
					}
					}
					}
					*/
					for ( i = 0 ; i < 12 ; i ++ )
					{
						int flag = 1, flag2 = 0 ;
						for ( j = 0 ; j < 4 ; j ++ )
						{
							int nx = ele->x + sheetNeighbor[i][j][0] ;
							int ny = ele->y + sheetNeighbor[i][j][1] ;
							int nz = ele->z + sheetNeighbor[i][j][2] ;

							double val = getDataAt( nx, ny, nz ) ;

							if ( val > - d && val < 0)
							{
								flag = 0 ;
								break ;
							}
							if ( val == d + 1 )
							{
								flag2 ++ ;
							}
						}

						if ( flag && flag2 )
						{
							dilatable = 1 ;
							break ;
						}
					}

					if ( ! dilatable )
					{
						// setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
						// queue3->prepend( ele->x, ele->y, ele->z ) ;

						setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
						if ( d > 2 )
						{
							queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
						}
						ele = queues[ d ]->remove() ;
					}
					else
					{
						setDataAt( ele->x, ele->y, ele->z, d ) ;
						ele = queues[ d ]->getNext() ;
					}
				}

				/* Debug: extract minimal set */
				while ( 1 )
				{
					int num = 0 ;
					queues[ d ]->reset() ;
					ele = queues[ d ]->getNext() ;
					while ( ele != NULL )
					{
						int critical = 0 ;
						setDataAt( ele->x, ele->y, ele->z, -1 ) ;

						for ( i = -1 ; i < 2 ; i ++ )
						{
							for ( j = -1 ; j < 2 ; j ++ )
							{
								for ( k = -1 ; k < 2 ; k ++ )
								{
									if ( i != 0 && j != 0 && k != 0 )
									{
										continue ;
									}
									int nx = ele->x + i ;
									int ny = ele->y + j ;
									int nz = ele->z + k ;
									// if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteSheet( nx,ny,nz, facevol ) )
									if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteSheet( nx,ny,nz ) )
									{
										critical = 1 ;
										break ;
									}
								}
								if ( critical )
								{
									break ;
								}
							}
							if ( critical )
							{
								break ;
							}
						}

						if ( critical )
						{
							setDataAt( ele->x, ele->y, ele->z, d ) ;
							ele = queues[ d ]->getNext() ;
						}
						else
						{
							setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
							if ( d > 2 )
							{
								queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
							}
							ele = queues[ d ]->remove() ;
							num ++ ;
						}

					}
					#ifdef VERBOSE
					printf("Non-minimal: %d\n", num) ;
					#endif

					if ( num == 0 )
					{
						break ;
					}
				}


				/* End of debugging */

				/*
				queue3->reset() ;
				ele = queue3->getNext() ;
				while ( ele != NULL )
				{
					setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
					ele = queue3->remove() ;
				}
				*/
			}


			// Finally, threshold the volume
			#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			#endif
			//threshold( -1, 1, 0, 0 ) ;
			threshold( 0, 0, 1 ) ;

			delete facevol ;
			delete fvol ;
			delete queue2;
			delete queue3;
			for (d = -dis ; d >= 2 ; d -- ) {
				delete queues[d];
			}
			delete [] queues;

			return - dis - 1 ;
		}

		//void Volume::erodeSheetOld( int disthr )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -1/0...\n") ;
			//#endif
			//threshold( 0.1f, -1, 0 ) ;

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2** queues = new GridQueue2* [ 64 ] ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//if ( ! hasCompleteSheet( i, j, k ) )
							//{
								//queue2->prepend( i, j, k ) ;
							//}
						//}
					//}
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;

			//// Start erosion
			//gridQueueEle* ele ;
			//int dis = -1 ;
			//while ( queue2->getNumElements() > 0 )
			//{
				//// First, set distance
				//dis -- ;
				//queues[ - dis ] = new GridQueue2( ) ;
				//printf("Distance transform to %d...", dis) ;
				//queue2->reset() ;
				//while( ( ele = queue2->getNext() ) != NULL )
				//{
					//setDataAt( ele->x, ele->y, ele->z, dis ) ;
					//queues[ -dis ]->prepend( ele->x, ele->y, ele->z ) ;
				//}
				//printf("%d nodes\n", queues[-dis]->getNumElements()) ;

				//// Next, find next layer
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ele->x + neighbor6[m][0] ;
						//int ny = ele->y + neighbor6[m][1] ;
						//int nz = ele->z + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 && ! hasCompleteSheet( nx, ny, nz ) )
						//{
							//setDataAt( nx, ny, nz, 1 ) ;
							//queue3->prepend( nx, ny, nz ) ;
						//}
					//}

					//ele = queue2->remove() ;
				//}

				//// Next, swap queues
				//GridQueue2* temp = queue2 ;
				//queue2 = queue3 ;
				//queue3 = temp ;
			//}

			///* Deal with closed rings
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) == 0 )
						//{
							//setDataAt( i, j, k, dis - 1 ) ;
						//}
					//}
			//*/

			///* Find local minimum: to help determine erosion level */
			//int cts[ 64 ] ;
			//for ( i = 0 ; i <= - dis ; i ++ )
			//{
				//cts[ i ] = 0 ;
			//}
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//double val = getDataAt( i, j, k ) ;
						//if ( val < -1 )
						//{
							//int min = 1 ;
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//int nx = i + neighbor6[m][0] ;
								//int ny = j + neighbor6[m][1] ;
								//int nz = k + neighbor6[m][2] ;
								//if ( getDataAt( nx, ny, nz ) < val )
								//{
									//min = 0 ;
									//break ;
								//}
							//}

							//if ( min )
							//{
								//cts[ (int) - val ] ++ ;
							//}
						//}
					//}

			//for ( i = 2 ; i <= - dis ; i ++ )
			//{
				//printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
			//}


			//// return ;

			//// Dilate back
			//// Starting from nodes with distance - 2 - disthr

			//for ( int d = disthr + 1 ; d >= 2 ; d -- )
			//{
				//delete queue3 ;
				//queue3 = new GridQueue2( ) ;
				//queues[ d ]->reset() ;
				//while ( (ele = queues[ d ]->getNext() ) != NULL )
				//{
					//int dilatable = 0 ;
					//// for ( int m = 0 ; m < 6 ; m ++ )
					///*
					//for ( int mx = -1 ; mx < 2 ; mx ++ )
						//for ( int my = -1 ; my < 2 ; my ++ )
							//for ( int mz = -1 ; mz < 2 ; mz ++ )
							//{
								//if ( mx == 0 || my == 0 || mz == 0 )
								//{
									//int nx = ele->x + mx ; // neighbor6[m][0] ;
									//int ny = ele->y + my ; // neighbor6[m][1] ;
									//int nz = ele->z + mz ; // neighbor6[m][2] ;
									//if ( getDataAt( nx, ny, nz ) == - d - 1 )
									//{
										//dilatable = 1 ;
										//break ;
									//}
								//}
							//}
					//*/
					//for ( i = 0 ; i < 12 ; i ++ )
					//{
						//int flag = 1, flag2 = 0 ;
						//for ( j = 0 ; j < 4 ; j ++ )
						//{
							//int nx = ele->x + sheetNeighbor[i][j][0] ;
							//int ny = ele->y + sheetNeighbor[i][j][1] ;
							//int nz = ele->z + sheetNeighbor[i][j][2] ;

							//double val = getDataAt( nx, ny, nz ) ;

							//if ( val > - d )
							//{
								//flag = 0 ;
								//break ;
							//}
							//if ( val == - d - 1 )
							//{
								//flag2 ++ ;
							//}
						//}

						//if ( flag && flag2 )
						//{
							//dilatable = 1 ;
							//break ;
						//}
					//}

					//if ( ! dilatable )
					//{
						//// setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
						//queue3->prepend( ele->x, ele->y, ele->z ) ;
						///*
						//setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
						//if ( d > 2 )
						//{
							//queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
						//}
						//*/
					//}
				//}

				//queue3->reset() ;
				//ele = queue3->getNext() ;
				//while ( ele != NULL )
				//{
					//setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
					//ele = queue3->remove() ;
				//}
			//}

			//// Finally, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( -1, 1, 0, 0 ) ;
		//}



		//void Volume::addNoise( float thr, float pos )
		//{
			//int i, j, k ;
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;
			//Volume* tvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( tvol->getDataAt( i, j, k ) >= 0  && isSimple( i, j, k ) )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( tvol->getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//if ( rand() < RAND_MAX * pos )
									//{
										//setDataAt( i, j, k, thr - 1 ) ;
									//}
									//break ;
								//}
							//}
						//}
					//}

		//}

		///************************************************************************/
		//// A sequential thinning method for extracting 6-connected skeletons
		//// Properties: medial, topology preserving, shape preserving, noise resistance!
		//// @param thr: threshold
		//// @param type: 0 for curve preserving, 1 for surface preserving
		//// @param noise: 0 for no-noise-reduction, n for level-n noise reduction
		///************************************************************************/
		//void Volume::sequentialSkeleton( float thr, int type, int noise )
		//{
			//int i, j, k, m ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -1/0...\n") ;
			//#endif
			//threshold( thr, -1, 0 ) ;

			//// Values in the volume:
			//// -1:	outside
			//// 0:	unvisited
			//// 1:	next layer
			//// 2:	this layer - non complex nodes
			//// 3:	complex nodes

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//gridQueueEle* ele ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//for ( m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//setDataAt( i, j, k, 1 ) ;
									//queue2->prepend( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//// printf("Total %d nodes\n", queue2->getNumElements() ) ;


			//// Perform erosion
			//int dis = 0 ;
			//int ox, oy, oz ;
			//int nx, ny, nz ;
			//while ( queue2->getNumElements() > 0 )
			//{
				//dis ++ ;
				//printf("At distance %d, there are %d nodes.\n", dis, queue2->getNumElements()) ;

				//// At the beginning of each iteration:
				//// queue2 holds all nodes dis away from the background in city block metric
				//// queue3 is empty

				//// First, find out next layer and put into queue3
				//queue2->reset() ;
				//while ( (ele = queue2->getNext()) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//setDataAt( ox, oy, oz, 2 ) ;
					//for ( m = 0 ; m < 6 ; m ++ )
					//{
							//nx = ox + neighbor6[m][0] ;
							//ny = oy + neighbor6[m][1] ;
							//nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//setDataAt( nx, ny, nz, 1 ) ;
								//queue3->prepend( nx, ny, nz ) ;
							//}
					//}
				//}

				//// Next, find out complex nodes in this layer and remove from queue2
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//if ( (! isSimple( ox, oy, oz ))  ||
						 //( ( ( type == 0 && isEdgeEnd( ox, oy, oz ) ) ||
							 //( type == 1 && isFaceEnd( ox, oy, oz ) ) ) && ! isNoise( ox, oy, oz, noise ) ) )
					//{
						//setDataAt( ox, oy, oz, 3 ) ;
						//ele = queue2->remove() ;
					//}
					//else
					//{
						//ele = queue2->getNext() ;
					//}
				//}

				//// Now, the main loop: delete non-complex nodes until none
				//queue2->reset() ;
				//while ( ( ele = queue2->getNext() ) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//queue2->remove() ;
					//queue2->reset() ;

					//if ( getDataAt(ox,oy,oz) != 2 )
					//{
						//continue ;
					//}

					//if ( ( ( ( type == 0 && isEdgeEnd( ox, oy, oz ) ) ||
									 //( type == 1 && isFaceEnd( ox, oy, oz ) ) ) && ! isNoise( ox, oy, oz, noise ) ) )
					//{
						//setDataAt( ox, oy, oz, 3 ) ;
					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
					//}

					//for ( i = -1 ; i < 2 ; i ++ )
						//for ( j = -1 ; j < 2 ; j ++ )
							//for ( k = -1 ; k < 2 ; k ++ )
							//{
								//nx = ox + i ;
								//ny = oy + j ;
								//nz = oz + k ;
								//int val = (int)( getDataAt( nx, ny, nz ) ) ;
								//if ( val > 1 )
								//{
									//int complex = 0 ;
									//if ( (! isSimple( nx, ny, nz ))  ||
										//( ( ( type == 0 && isEdgeEnd( nx, ny, nz ) ) ||
										//( type == 1 && isFaceEnd( nx, ny, nz ) ) ) && ! isNoise( nx, ny, nz, noise ) ) )
									//{
										//complex = 1 ;
									//}

									//if ( val == 2 && complex )
									//{
										//// A non-complex node becomes complex
										//setDataAt( nx, ny, nz, 3 ) ;
									//}
									//else if ( val == 3 && ! complex )
									//{
										//// A complex node turns non-complex
										//setDataAt( nx, ny, nz, 2 ) ;
										//queue2->prepend( nx, ny, nz ) ;
									//}

								//}
							//}

					//queue2->reset() ;
				//}

				//// Finally, swap queue3 and queue2
				//GridQueue2* temp = queue3 ;
				//queue3 = queue2 ;
				//queue2 = temp ;
			//}


			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::dumbsurfaceSkeleton( float thr )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//gridQueueEle* ele ;


			//while ( 1 )
			//{
				//int n = 0 ;

				//queue2->reset() ;
				//for ( i = 0 ; i < getSizeX() ; i ++ )
					//for ( j = 0 ; j < getSizeY() ; j ++ )
						//for ( k = 0 ; k < getSizeZ() ; k ++ )
						//{
							//if ( getDataAt( i, j, k ) == 0 )
							//{
								//{
									//for ( int m = 0 ; m < 6 ; m ++ )
									//{
										//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
										//{
											//queue2->prepend( i, j, k ) ;
											//break ;
										//}
									//}
								//}
							//}
						//}

				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//int ox = ele->x ;
					//int oy = ele->y ;
					//int oz = ele->z ;
					//if ( isSimple( ox, oy, oz ) && hasCompleteSheet( ox, oy, oz ) == 1 )
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//n ++ ;
					//}

					//ele = queue2->remove() ;
				//}

				//if ( n == 0 )
				//{
					//break ;
				//}

				//printf("%d simple nodes found.\n", n);
			//}


			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
		//}

		///* Thin the current volume while preserving voxels with values > highthr or <= lowthr in grayvol
		//*  Assuming the current volume has already been thresholded to 0/1
		//*/
		//void Volume::surfaceSkeleton( Volume* grayvol, float lowthr, float highthr ) {
			//int i, j, k ;
			//threshold( 0.5f, -MAX_ERODE, 0 ) ;

			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2* queue4 = new GridQueue2( ) ;

			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );
			//int ct = 0 ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//float v = (float)grayvol->getDataAt(i,j,k) ;
							//if ( v > highthr || v <= lowthr )
							//{
								//setDataAt( i, j, k, MAX_ERODE ) ;
							//}
							//else
							//{
								//ct ++ ;

								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}


			//// Perform erosion
			//int wid = MAX_ERODE ;
			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}

			//#ifdef  NOISE_DIS_SHEET
			//Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
			//#endif

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{


				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				//#endif

				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//if ( getDataAt(ox,oy,oz) == curwid )
					//{
						//ele = queue2->remove() ;
					//}
					//else
					//{
						//setDataAt(ox,oy,oz, curwid) ;
						//ele = queue2->getNext() ;
					//}
				//}
				//queue4->reset() ;
				//ele = queue4->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//queue2->prepend(ox,oy,oz) ;
					//ele = queue4->remove() ;
				//}

				//// Now queue2 holds all the nodes for this layer

				//#ifdef NOISE_DIS_SHEET
				///* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				//queue2->reset() ;

				//// First run
				//int flag = 0 ;
				//while ( ( ele = queue2->getNext() ) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//if ( NOISE_DIS_SHEET <= 1 )
					//{
						//noisevol->setDataAt( ox, oy, oz, 0 ) ;
					//}
					//else
					//{
						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//flag = 1 ;
								//break ;
							//}
						//}
						//if ( ! flag )
						//{
							//noisevol->setDataAt( ox, oy, oz, 0 ) ;
						//}
					//}
				//}

				//for ( int cur = 1 ; cur < NOISE_DIS_SHEET ; cur ++ )
				//{
					//queue2->reset() ;
					//int count = 0 ;

					//while ( ( ele = queue2->getNext() ) != NULL )
					//{
						//ox = ele->x ;
						//oy = ele->y ;
						//oz = ele->z ;

						//if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						//{
							//continue ;
						//}

						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//count ++ ;
								//break ;
							//}
						//}
					//}

					//if ( count == 0 )
					//{
						//break ;
					//}
				//}


				//#endif

				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//// queue->add( gp, -score ) ;
					//queue->add( gp, score ) ;

					//ele = queue2->remove() ;
				//}


				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//int nowComplex = 0 ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;

					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}


					//#ifndef NOISE_DIS_SHEET
					//// if ( hasFeatureFace( ox, oy, oz ) )
					//if ( (! isSimple( ox, oy, oz )) || isSheetEnd( ox, oy, oz ) )
					//// if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz)))
					//#else
					//// if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) )
					//if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) || isHelixEnd( ox, oy, oz, noisevol ))
					//// if ( isBertrandEndPoint( ox, oy, oz ) )
					//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue4->prepend( ox, oy, oz ) ;
						//numComplex ++ ;

						//nowComplex = 1 ;
					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//numSimple ++ ;

						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//// setDataAt( nx, ny, nz, curwid + 1 ) ;
								//queue2->prepend( nx, ny, nz ) ;
							//}
						//}

						//if ( nowComplex )
						//{

							//// printf("Error: %d\n", score);
						//}
					//}

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue

					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//// queue->add( gp, -score ) ;
										//queue->add( gp, score ) ;
									//}
								//}
							//}


				//}

				//if ( numSimple == 0 )
				//{
					//if ( queue2->getNumElements() > 0 )
					//{
						//printf("*************************wierd**********************\n");
					//}
						//break ;
				//}
			//}

			//// Remove all internal voxels (contained in cells)

			//queue4->reset() ;
			//ele = queue4->getNext() ;
			//while ( ele != NULL )
			//{
				//ele = queue4->remove() ;
			//}

			//queue2->reset() ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) == 0 && isInternal2( i,j,k ) == 1 )
						//{
							//queue2->prepend( i, j, k ) ;
						//}
					//}
			//queue2->reset() ;
			//ele = queue2->getNext() ;
			//while ( ele != NULL )
			//{
				//ox = ele->x ;
				//oy = ele->y ;
				//oz = ele->z ;
				//setDataAt( ox, oy, oz, -1 ) ;
				//ele = queue2->remove() ;
			//}



			//// Finally, clean up
			//delete scrvol;
			//delete queue;
			//delete queue2;
			//delete queue3;
			//delete queue4;
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::surfaceSkeleton( float thr )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the linked queue
			//printf("Initializing queue...\n") ;
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2* queue4 = new GridQueue2( ) ;

			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//// setDataAt( i, j, k, 1 ) ;
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}
			//int wid = MAX_ERODE ;
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;


			//// Perform erosion
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}

	//#ifdef  NOISE_DIS_SHEET
			//Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	//#endif

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 and queue4 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;

				///*
				//We first need to assign curwid + 1 to every node in this layer
				//*/
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//if ( getDataAt(ox,oy,oz) == curwid )
					//{
						//ele = queue2->remove() ;
					//}
					//else
					//{
						//setDataAt(ox,oy,oz, curwid) ;
						//ele = queue2->getNext() ;
					//}
				//}
				//queue4->reset() ;
				//ele = queue4->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//queue2->prepend(ox,oy,oz) ;
					//ele = queue4->remove() ;
				//}

				//// Now queue2 holds all the nodes for this layer

	//#ifdef NOISE_DIS_SHEET
				///* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				//queue2->reset() ;

				//// First run
				//int flag = 0 ;
				//while ( ( ele = queue2->getNext() ) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//if ( NOISE_DIS_SHEET <= 1 )
					//{
						//noisevol->setDataAt( ox, oy, oz, 0 ) ;
					//}
					//else
					//{
						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//flag = 1 ;
								//break ;
							//}
						//}
						//if ( ! flag )
						//{
							//noisevol->setDataAt( ox, oy, oz, 0 ) ;
						//}
					//}
				//}

				//for ( int cur = 1 ; cur < NOISE_DIS_SHEET ; cur ++ )
				//{
					//queue2->reset() ;
					//int count = 0 ;

					//while ( ( ele = queue2->getNext() ) != NULL )
					//{
						//ox = ele->x ;
						//oy = ele->y ;
						//oz = ele->z ;

						//if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						//{
							//continue ;
						//}

						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//count ++ ;
								//break ;
							//}
						//}
					//}

					//if ( count == 0 )
					//{
						//break ;
					//}
				//}


	//#endif

				///* Commented for debugging

				//// First,
				//// check for complex nodes in queue2
				//// move them from queue2 to queue3
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Check simple
	//#ifndef NOISE_DIS_SHEET
					//if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	//#else
					//if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue3->prepend( ox, oy, oz ) ;
						//ele = queue2->remove() ;

						//numComplex ++ ;
					//}
					//else
					//{
						//ele = queue2->getNext() ;
					//}
				//}
				//*/


				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//// queue->add( gp, -score ) ;
					//queue->add( gp, score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//// From now on, queue2 holds nodes of next level
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//int nowComplex = 0 ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
					//// printf("%d\n", score);
	////				score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					///* Commented for debugging

					//// Remove this simple node
					//setDataAt( ox, oy, oz, -1 ) ;
					//numSimple ++ ;
					//// printf("Highest score: %d\n", score) ;
					//*/

					///* Added for debugging */
					//// Check simple
	//#ifndef NOISE_DIS_SHEET
					//// if ( hasFeatureFace( ox, oy, oz ) )
					//if ( (! isSimple( ox, oy, oz )) || isSheetEnd( ox, oy, oz ) )
					//// if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz)))
	//#else
					//// if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) )
					//if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) || isHelixEnd( ox, oy, oz, noisevol ))
					//// if ( isBertrandEndPoint( ox, oy, oz ) )
	//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue4->prepend( ox, oy, oz ) ;
						//numComplex ++ ;

						//nowComplex = 1 ;
					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//numSimple ++ ;

						//if ( nowComplex )
						//{

							//// printf("Error: %d\n", score);
						//}
					//}
					///* Adding ends */

					//// Move its neighboring unvisited node to queue2
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//// setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

					///* Commented for debugging

					//// Find complex nodes in its 3x3 neighborhood
					//// move them to queue2
					//for ( i = -1 ; i < 2 ; i ++ )
						//for ( j = -1 ; j < 2 ; j ++ )
							//for ( k = -1 ; k < 2 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//// Check simple
								//if ( getDataAt( nx, ny, nz ) == curwid &&
									//// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	//#ifndef NOISE_DIS_SHEET
									//( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	//#else
									//( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	//#endif

								//{
									//// Complex, set to next layer
									//setDataAt( nx, ny, nz, curwid + 1 ) ;
									//queue2->prepend( nx, ny, nz ) ;
									//numComplex ++ ;
								//}
							//}
					//*/

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue

					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//// queue->add( gp, -score ) ;
										//queue->add( gp, score ) ;
									//}
								//}
							//}


				//}

				//printf("%d complex, %d simple\n", numComplex, numSimple) ;

				//if ( numSimple == 0 )
				//{
						//break ;
				//}
			//}

			//// Finally, clean up
			//printf("Thresholding the volume to 0/1...\n") ;
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::surfaceSkeleton( float thr, Volume* svol )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the linked queue
			//printf("Initializing queue...\n") ;
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2* queue4 = new GridQueue2( ) ;

			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//if ( svol->getDataAt(i,j,k) > 0 )
							//{
								//setDataAt( i, j, k, MAX_ERODE ) ;
							//}
							//else
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//// setDataAt( i, j, k, 1 ) ;
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}
			//int wid = MAX_ERODE ;
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;


			//// Perform erosion
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}

	//#ifdef  NOISE_DIS_SHEET
			//Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	//#endif

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 and queue4 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;

				///*
				//We first need to assign curwid + 1 to every node in this layer
				//*/
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//if ( getDataAt(ox,oy,oz) == curwid )
					//{
						//ele = queue2->remove() ;
					//}
					//else
					//{
						//setDataAt(ox,oy,oz, curwid) ;
						//ele = queue2->getNext() ;
					//}
				//}
				//queue4->reset() ;
				//ele = queue4->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//queue2->prepend(ox,oy,oz) ;
					//ele = queue4->remove() ;
				//}

				//// Now queue2 holds all the nodes for this layer

	//#ifdef NOISE_DIS_SHEET
				///* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				//queue2->reset() ;

				//// First run
				//int flag = 0 ;
				//while ( ( ele = queue2->getNext() ) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//if ( NOISE_DIS_SHEET <= 1 )
					//{
						//noisevol->setDataAt( ox, oy, oz, 0 ) ;
					//}
					//else
					//{
						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//flag = 1 ;
								//break ;
							//}
						//}
						//if ( ! flag )
						//{
							//noisevol->setDataAt( ox, oy, oz, 0 ) ;
						//}
					//}
				//}

				//for ( int cur = 1 ; cur < NOISE_DIS_SHEET ; cur ++ )
				//{
					//queue2->reset() ;
					//int count = 0 ;

					//while ( ( ele = queue2->getNext() ) != NULL )
					//{
						//ox = ele->x ;
						//oy = ele->y ;
						//oz = ele->z ;

						//if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						//{
							//continue ;
						//}

						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//count ++ ;
								//break ;
							//}
						//}
					//}

					//if ( count == 0 )
					//{
						//break ;
					//}
				//}


	//#endif

				///* Commented for debugging

				//// First,
				//// check for complex nodes in queue2
				//// move them from queue2 to queue3
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Check simple
	//#ifndef NOISE_DIS_SHEET
					//if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	//#else
					//if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue3->prepend( ox, oy, oz ) ;
						//ele = queue2->remove() ;

						//numComplex ++ ;
					//}
					//else
					//{
						//ele = queue2->getNext() ;
					//}
				//}
				//*/


				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//// queue->add( gp, -score ) ;
					//queue->add( gp, score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//// From now on, queue2 holds nodes of next level
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//int nowComplex = 0 ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
					//// printf("%d\n", score);
	////				score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					///* Commented for debugging

					//// Remove this simple node
					//setDataAt( ox, oy, oz, -1 ) ;
					//numSimple ++ ;
					//// printf("Highest score: %d\n", score) ;
					//*/

					///* Added for debugging */
					//// Check simple
	//#ifndef NOISE_DIS_SHEET
					//// if ( hasFeatureFace( ox, oy, oz ) )
					//if ( (! isSimple( ox, oy, oz )) || isSheetEnd( ox, oy, oz ) )
					//// if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz)))
	//#else
					//// if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) )
					//if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) || isHelixEnd( ox, oy, oz, noisevol ))
					//// if ( isBertrandEndPoint( ox, oy, oz ) )
	//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue4->prepend( ox, oy, oz ) ;
						//numComplex ++ ;

						//nowComplex = 1 ;
					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//numSimple ++ ;

						//if ( nowComplex )
						//{

							//// printf("Error: %d\n", score);
						//}
					//}
					///* Adding ends */

					//// Move its neighboring unvisited node to queue2
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//// setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

					///* Commented for debugging

					//// Find complex nodes in its 3x3 neighborhood
					//// move them to queue2
					//for ( i = -1 ; i < 2 ; i ++ )
						//for ( j = -1 ; j < 2 ; j ++ )
							//for ( k = -1 ; k < 2 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//// Check simple
								//if ( getDataAt( nx, ny, nz ) == curwid &&
									//// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	//#ifndef NOISE_DIS_SHEET
									//( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	//#else
									//( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	//#endif

								//{
									//// Complex, set to next layer
									//setDataAt( nx, ny, nz, curwid + 1 ) ;
									//queue2->prepend( nx, ny, nz ) ;
									//numComplex ++ ;
								//}
							//}
					//*/

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue

					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//// queue->add( gp, -score ) ;
										//queue->add( gp, score ) ;
									//}
								//}
							//}


				//}

				//printf("%d complex, %d simple\n", numComplex, numSimple) ;

				//if ( numSimple == 0 )
				//{
						//break ;
				//}
			//}

			//// Finally, clean up
			//printf("Thresholding the volume to 0/1...\n") ;
			//threshold( 0, 0, 1 ) ;

		//}

		//void Volume::surfaceSkeletonOld( float thr )
		//{
			//int i, j, k ;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Next, initialize the linked queue
			//#ifdef VERBOSE
			//printf("Initializing queue...\n") ;
			//#endif
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2* queue4 = new GridQueue2( ) ;

			//PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= 0 )
						//{
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									//{
										//// setDataAt( i, j, k, 1 ) ;
										//queue2->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}
					//}
			//#ifdef VERBOSE
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;
			//#endif

			//// Perform erosion
			//int wid = MAX_ERODE ;
			//#ifdef VERBOSE
			//printf("Start erosion to %d...\n", wid) ;
			//#endif
			//gridQueueEle* ele ;
			//gridPoint* gp ;
			//int ox, oy, oz ;
			//int score ;
			//Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//scrvol->setDataAt( i, -1 ) ;
			//}

	//#ifdef  NOISE_DIS_SHEET
			//Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	//#endif

			//for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			//{
				//// At the start of each iteration,
				//// queue2 and queue4 holds all the nodes for this layer
				//// queue3 and queue are empty

				//int numComplex = 0, numSimple = 0 ;
				//#ifdef VERBOSE
				//printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				//#endif

				///*
				//We first need to assign curwid + 1 to every node in this layer
				//*/
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//if ( getDataAt(ox,oy,oz) == curwid )
					//{
						//ele = queue2->remove() ;
					//}
					//else
					//{
						//setDataAt(ox,oy,oz, curwid) ;
						//ele = queue2->getNext() ;
					//}
				//}
				//queue4->reset() ;
				//ele = queue4->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//queue2->prepend(ox,oy,oz) ;
					//ele = queue4->remove() ;
				//}

				//// Now queue2 holds all the nodes for this layer

	//#ifdef NOISE_DIS_SHEET
				///* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				//queue2->reset() ;

				//// First run
				//int flag = 0 ;
				//while ( ( ele = queue2->getNext() ) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//if ( NOISE_DIS_SHEET <= 1 )
					//{
						//noisevol->setDataAt( ox, oy, oz, 0 ) ;
					//}
					//else
					//{
						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) == 0 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//flag = 1 ;
								//break ;
							//}
						//}
						//if ( ! flag )
						//{
							//noisevol->setDataAt( ox, oy, oz, 0 ) ;
						//}
					//}
				//}

				//for ( int cur = 1 ; cur < NOISE_DIS_SHEET ; cur ++ )
				//{
					//queue2->reset() ;
					//int count = 0 ;

					//while ( ( ele = queue2->getNext() ) != NULL )
					//{
						//ox = ele->x ;
						//oy = ele->y ;
						//oz = ele->z ;

						//if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						//{
							//continue ;
						//}

						//flag = 0 ;
						//for ( int m = 0 ; m < 6 ; m ++ )
						//{
							//int nx = ox + neighbor6[m][0] ;
							//int ny = oy + neighbor6[m][1] ;
							//int nz = oz + neighbor6[m][2] ;
							//if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							//{
								//noisevol->setDataAt( ox, oy, oz, 1 ) ;
								//count ++ ;
								//break ;
							//}
						//}
					//}

					//if ( count == 0 )
					//{
						//break ;
					//}
				//}


	//#endif

				///* Commented for debugging

				//// First,
				//// check for complex nodes in queue2
				//// move them from queue2 to queue3
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Check simple
	//#ifndef NOISE_DIS_SHEET
					//if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	//#else
					//if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue3->prepend( ox, oy, oz ) ;
						//ele = queue2->remove() ;

						//numComplex ++ ;
					//}
					//else
					//{
						//ele = queue2->getNext() ;
					//}
				//}
				//*/


				//// Next,
				//// Compute score for each node left in queue2
				//// move them into priority queue
				//queue2->reset() ;
				//ele = queue2->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//// Compute score
					//score = getNumPotComplex( ox, oy, oz ) ;
					//scrvol->setDataAt( ox, oy, oz, score ) ;

					//// Push to queue
					//gp = new gridPoint ;
					//gp->x = ox ;
					//gp->y = oy ;
					//gp->z = oz ;
					//// queue->add( gp, -score ) ;
					//queue->add( gp, score ) ;

					//ele = queue2->remove() ;
				//}

				//// Rename queue3 to be queue2,
				//// Clear queue3
				//// From now on, queue2 holds nodes of next level
				//delete queue2 ;
				//queue2 = queue3 ;
				//queue3 = new GridQueue2( ) ;

				//// Next, start priority queue iteration
				//while ( ! queue->isEmpty() )
				//{
					//// Retrieve the node with the highest score
					//queue->remove( gp, score ) ;
					//ox = gp->x ;
					//oy = gp->y ;
					//oz = gp->z ;
					//delete gp ;
					//// printf("%d\n", score);
	////				score = -score ;

					//// Ignore the node
					//// if it has been processed before
					//// or it has an updated score
					//if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					//{
						//continue ;
					//}

					///* Commented for debugging

					//// Remove this simple node
					//setDataAt( ox, oy, oz, -1 ) ;
					//numSimple ++ ;
					//// printf("Highest score: %d\n", score) ;
					//*/

					///* Added for debugging */
					//// Check simple
	//#ifndef NOISE_DIS_SHEET
					//// if ( hasFeatureFace( ox, oy, oz ) )
					//if ( (! isSimple( ox, oy, oz )) || isSheetEnd( ox, oy, oz ) )
					//// if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz)))
	//#else
					//// if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) )
					//if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) || isHelixEnd( ox, oy, oz, noisevol ))
					//// if ( isBertrandEndPoint( ox, oy, oz ) )
	//#endif
					//{
						//// Complex, set to next layer
						//setDataAt( ox, oy, oz, curwid + 1 ) ;
						//queue4->prepend( ox, oy, oz ) ;
						//numComplex ++ ;

					//}
					//else
					//{
						//setDataAt( ox, oy, oz, -1 ) ;
						//numSimple ++ ;

					//}
					///* Adding ends */

					//// Move its neighboring unvisited node to queue2
					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) == 0 )
						//{
							//// setDataAt( nx, ny, nz, curwid + 1 ) ;
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}

					///* Commented for debugging

					//// Find complex nodes in its 3x3 neighborhood
					//// move them to queue2
					//for ( i = -1 ; i < 2 ; i ++ )
						//for ( j = -1 ; j < 2 ; j ++ )
							//for ( k = -1 ; k < 2 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//// Check simple
								//if ( getDataAt( nx, ny, nz ) == curwid &&
									//// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	//#ifndef NOISE_DIS_SHEET
									//( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	//#else
									//( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	//#endif

								//{
									//// Complex, set to next layer
									//setDataAt( nx, ny, nz, curwid + 1 ) ;
									//queue2->prepend( nx, ny, nz ) ;
									//numComplex ++ ;
								//}
							//}
					//*/

					//// Update scores for nodes in its 5x5 neighborhood
					//// insert them back into priority queue

					//for ( i = -2 ; i < 3 ;i ++ )
						//for ( j = -2 ; j < 3 ; j ++ )
							//for ( k = -2 ; k < 3 ; k ++ )
							//{
								//int nx = ox + i ;
								//int ny = oy + j ;
								//int nz = oz + k ;

								//if ( getDataAt( nx, ny, nz ) == curwid )
								//{
									//// Compute score
									//score = getNumPotComplex( nx, ny, nz ) ;

									//if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									//{
										//// printf("Update\n") ;
										//scrvol->setDataAt( nx, ny, nz, score ) ;
										//// Push to queue
										//gp = new gridPoint ;
										//gp->x = nx ;
										//gp->y = ny ;
										//gp->z = nz ;
										//// queue->add( gp, -score ) ;
										//queue->add( gp, score ) ;
									//}
								//}
							//}


				//}
				//#ifdef VERBOSE
				//printf("%d complex, %d simple\n", numComplex, numSimple) ;
				//#endif

				//if ( numSimple == 0 )
				//{
						//break ;
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
			//delete queue;

		//}

		void Volume::surfaceSkeletonPres( float thr, Volume * preserve )
		{
			int i, j, k ;
			// First, threshold the volume
			#ifdef VERBOSE
			printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			#endif
			threshold( thr, -MAX_ERODE, 0 ) ;

			// Next, initialize the linked queue
			#ifdef VERBOSE
			printf("Initializing queue...\n") ;
			#endif
			GridQueue2* queue2 = new GridQueue2( ) ;
			GridQueue2* queue3 = new GridQueue2( ) ;
			GridQueue2* queue4 = new GridQueue2( ) ;

			PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );

			for ( i = 0 ; i < getSizeX() ; i ++ )
				for ( j = 0 ; j < getSizeY() ; j ++ )
					for ( k = 0 ; k < getSizeZ() ; k ++ )
					{
						if ( getDataAt( i, j, k ) >= 0 ) {
							if(preserve->getDataAt(i, j, k) > 0) {
								setDataAt(i, j, k, MAX_ERODE);
							} else {
								for ( int m = 0 ; m < 6 ; m ++ )
								{
									if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
									{
										// setDataAt( i, j, k, 1 ) ;
										queue2->prepend( i, j, k ) ;
										break ;
									}
								}
							}
						}
					}
			int wid = MAX_ERODE ;
			#ifdef VERBOSE
			printf("Total %d nodes\n", queue2->getNumElements() ) ;
			printf("Start erosion to %d...\n", wid) ;
			#endif


			// Perform erosion
			gridQueueEle* ele ;
			gridPoint* gp ;
			int ox, oy, oz ;
			int score ;
			Volume* scrvol = new Volume( this->getSizeX() , this->getSizeY(), this->getSizeZ() ) ;
			for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			{
				scrvol->setDataAt( i, -1 ) ;
			}

	#ifdef  NOISE_DIS_SHEET
			Volume* noisevol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
	#endif

			for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
			{
				// At the start of each iteration,
				// queue2 and queue4 holds all the nodes for this layer
				// queue3 and queue are empty

				int numComplex = 0, numSimple = 0 ;
				#ifdef VERBOSE
				printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
				#endif

				/*
				We first need to assign curwid + 1 to every node in this layer
				*/
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					if ( getDataAt(ox,oy,oz) == curwid )
					{
						ele = queue2->remove() ;
					}
					else
					{
						setDataAt(ox,oy,oz, curwid) ;
						ele = queue2->getNext() ;
					}
				}
				queue4->reset() ;
				ele = queue4->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					queue2->prepend(ox,oy,oz) ;
					ele = queue4->remove() ;
				}

				// Now queue2 holds all the nodes for this layer

	#ifdef NOISE_DIS_SHEET
				/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
				queue2->reset() ;

				// First run
				int flag = 0 ;
				while ( ( ele = queue2->getNext() ) != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;
					if ( NOISE_DIS_SHEET <= 1 )
					{
						noisevol->setDataAt( ox, oy, oz, 0 ) ;
					}
					else
					{
						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) == 0 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								flag = 1 ;
								break ;
							}
						}
						if ( ! flag )
						{
							noisevol->setDataAt( ox, oy, oz, 0 ) ;
						}
					}
				}

				for ( int cur = 1 ; cur < NOISE_DIS_SHEET ; cur ++ )
				{
					queue2->reset() ;
					int count = 0 ;

					while ( ( ele = queue2->getNext() ) != NULL )
					{
						ox = ele->x ;
						oy = ele->y ;
						oz = ele->z ;

						if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
						{
							continue ;
						}

						flag = 0 ;
						for ( int m = 0 ; m < 6 ; m ++ )
						{
							int nx = ox + neighbor6[m][0] ;
							int ny = oy + neighbor6[m][1] ;
							int nz = oz + neighbor6[m][2] ;
							if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
							{
								noisevol->setDataAt( ox, oy, oz, 1 ) ;
								count ++ ;
								break ;
							}
						}
					}

					if ( count == 0 )
					{
						break ;
					}
				}


	#endif

				/* Commented for debugging

				// First,
				// check for complex nodes in queue2
				// move them from queue2 to queue3
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Check simple
	#ifndef NOISE_DIS_SHEET
					if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) )
	#else
					if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue3->prepend( ox, oy, oz ) ;
						ele = queue2->remove() ;

						numComplex ++ ;
					}
					else
					{
						ele = queue2->getNext() ;
					}
				}
				*/


				// Next,
				// Compute score for each node left in queue2
				// move them into priority queue
				queue2->reset() ;
				ele = queue2->getNext() ;
				while ( ele != NULL )
				{
					ox = ele->x ;
					oy = ele->y ;
					oz = ele->z ;

					// Compute score
					score = getNumPotComplex( ox, oy, oz ) ;
					scrvol->setDataAt( ox, oy, oz, score ) ;

					// Push to queue
					gp = new gridPoint ;
					gp->x = ox ;
					gp->y = oy ;
					gp->z = oz ;
					// queue->add( gp, -score ) ;
					queue->add( gp, score ) ;

					ele = queue2->remove() ;
				}

				// Rename queue3 to be queue2,
				// Clear queue3
				// From now on, queue2 holds nodes of next level
				delete queue2 ;
				queue2 = queue3 ;
				queue3 = new GridQueue2( ) ;


				// Next, start priority queue iteration
				while ( ! queue->isEmpty() )
				{
					// Retrieve the node with the highest score
					queue->remove( gp, score ) ;
					ox = gp->x ;
					oy = gp->y ;
					oz = gp->z ;
					delete gp ;
					// printf("%d\n", score);
	//				score = -score ;

					// Ignore the node
					// if it has been processed before
					// or it has an updated score
					if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
					{
						continue ;
					}

					/* Commented for debugging

					// Remove this simple node
					setDataAt( ox, oy, oz, -1 ) ;
					numSimple ++ ;
					// printf("Highest score: %d\n", score) ;
					*/

					/* Added for debugging */
					// Check simple
	#ifndef NOISE_DIS_SHEET
					// if ( hasFeatureFace( ox, oy, oz ) )
					if ( (! isSimple( ox, oy, oz )) || isSheetEnd( ox, oy, oz ) )
					// if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz)))
	#else
					// if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) )
					if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) || isHelixEnd( ox, oy, oz, noisevol ))
					// if ( isBertrandEndPoint( ox, oy, oz ) )
	#endif
					{
						// Complex, set to next layer
						setDataAt( ox, oy, oz, curwid + 1 ) ;
						queue4->prepend( ox, oy, oz ) ;
						numComplex ++ ;

					}
					else
					{
						setDataAt( ox, oy, oz, -1 ) ;
						numSimple ++ ;

					}
					/* Adding ends */

					// Move its neighboring unvisited node to queue2
					for ( int m = 0 ; m < 6 ; m ++ )
					{
						int nx = ox + neighbor6[m][0] ;
						int ny = oy + neighbor6[m][1] ;
						int nz = oz + neighbor6[m][2] ;
						if ( getDataAt( nx, ny, nz ) == 0 )
						{
							// setDataAt( nx, ny, nz, curwid + 1 ) ;
							queue2->prepend( nx, ny, nz ) ;
						}
					}

					/* Commented for debugging

					// Find complex nodes in its 3x3 neighborhood
					// move them to queue2
					for ( i = -1 ; i < 2 ; i ++ )
						for ( j = -1 ; j < 2 ; j ++ )
							for ( k = -1 ; k < 2 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								// Check simple
								if ( getDataAt( nx, ny, nz ) == curwid &&
									// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
	#ifndef NOISE_DIS_SHEET
									( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
	#else
									( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
	#endif

								{
									// Complex, set to next layer
									setDataAt( nx, ny, nz, curwid + 1 ) ;
									queue2->prepend( nx, ny, nz ) ;
									numComplex ++ ;
								}
							}
					*/

					// Update scores for nodes in its 5x5 neighborhood
					// insert them back into priority queue

					for ( i = -2 ; i < 3 ;i ++ )
						for ( j = -2 ; j < 3 ; j ++ )
							for ( k = -2 ; k < 3 ; k ++ )
							{
								int nx = ox + i ;
								int ny = oy + j ;
								int nz = oz + k ;

								if ( getDataAt( nx, ny, nz ) == curwid )
								{
									// Compute score
									score = getNumPotComplex( nx, ny, nz ) ;

									if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
									{
										// printf("Update\n") ;
										scrvol->setDataAt( nx, ny, nz, score ) ;
										// Push to queue
										gp = new gridPoint ;
										gp->x = nx ;
										gp->y = ny ;
										gp->z = nz ;
										// queue->add( gp, -score ) ;
										queue->add( gp, score ) ;
									}
								}
							}


				}

				#ifdef VERBOSE
				printf("%d complex, %d simple\n", numComplex, numSimple) ;
				#endif

				if ( numSimple == 0 )
				{
						break ;
				}
			}

			// Finally, clean up
			#ifdef VERBOSE
			printf("Thresholding the volume to 0/1...\n") ;
			#endif
			threshold( 0, 0, 1 ) ;

			delete scrvol;
			delete queue;
			delete queue2;
			delete queue3;
			delete queue4;

		}

		///* Bertrand's parallel 6-connected surface thinning */
		//void Volume::bertrandSurfaceSkeleton2( float thr )
		//{
			//int i, j, k ;

			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;


			//// Initialize queues
			//printf("Initializing queues...\n") ;
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2* queue4 = new GridQueue2( ) ;
			//Volume* fvol = new Volume( this->getSizeX(), this->getSizeY(), this->getSizeZ(), 0 ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) == 0 )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
								//{
									//fvol->setDataAt( i, j, k, 1 ) ;
									//queue2->prepend( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;

			//// Start iteration
			//int it = 0 ;
			//gridQueueEle* ele ;
			//while ( ++ it )
			//{
				//printf("Iteration %d... nodes in queue: %d...", it, queue2->getNumElements()) ;

				//// queue2 holds all boundary nodes

				//int deleted = 0 ;
				//for ( i = 0 ; i < 6 ; i ++ )
				//{
					//// At the beginning of each iteration,
					//// queue2 holds all remaining boundary nodes
					//// queue3 is a deletable array, starting empty
					//// queue4 holds the candidates for next layer

					//queue2->reset() ;
					//ele = queue2->getNext() ;

					//// For each sub-iteration, go through queue2 first and find points to delete
					//while ( ele != NULL )
					//{
						//int ox = ele->x ;
						//int oy = ele->y ;
						//int oz = ele->z ;

						//if ( isBertrandBorder( ox,oy,oz, i ) )
						//{
							//if ( ! isBertrandEndPoint( ox,oy,oz ) )
							//{
								//// Move this node to queue3
								//ele = queue2->remove() ;
								//queue3->prepend( ox, oy, oz ) ;

								//// Move its neighboring unvisited node to queue4
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;
									//if ( fvol->getDataAt( nx, ny, nz ) == 0 )
									//{
										//fvol->setDataAt( nx, ny, nz, 1 ) ;
										//queue4->prepend( nx, ny, nz ) ;
									//}
								//}
							//}
							//else
							//{
								//ele = queue2->getNext() ;
							//}
						//}
						//else
						//{
							//ele = queue2->getNext() ;
						//}
					//}

					//// Now, queue2 holds all remaining nodes,
					//// queue3 holds all deletable nodes,
					//// and queue4 holds nodes to be added to the next layer

					//// Simultaneous deletion
					//if ( queue3->getNumElements() == 0 )
					//{
						//// break ;
					//}
					//else
					//{
						//queue3->reset() ;
						//ele = queue3->getNext() ;
						//while ( ele != NULL )
						//{
							//setDataAt( ele->x, ele->y, ele->z, -1 ) ;
							//ele = queue3->remove() ;
							//deleted ++ ;
						//}
					//}

					//// return ;
				//}

				//// After all sub-iterations
				//// Move all queue4 nodes to queue2
				//queue4->reset() ;
				//ele = queue4->getNext() ;
				//while ( ele != NULL )
				//{
					//queue2->prepend( ele->x, ele->y, ele->z) ;
					//ele = queue4->remove() ;
				//}

				//if ( deleted == 0 )
				//{
					//printf("No more deletable nodes.\n");
					//break ;
				//}
				//else
				//{
					//printf("Deleted: %d\n", deleted) ;
				//}
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;

		//}



		//void Volume::bertrandSurfaceSkeleton( float thr )
		//{
			//int i, j, k ;

			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
			//#endif
			//threshold( thr, -MAX_ERODE, 0 ) ;

			//// Start iteration
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//int dir = 0, ct = 0 ;
			//while( ++ct )
			//{
				//printf("Round %d, direction %d...", ct, dir) ;
				//queue2->reset() ;

				//for ( i = 0 ; i < getSizeX() ; i ++ )
					//for ( j = 0 ; j < getSizeY() ; j ++ )
						//for ( k = 0 ; k < getSizeZ() ; k ++ )
						//{
							//if ( getDataAt(i,j,k) >= 0 )
							//{
								//if ( isBertrandBorder( i,j,k, dir ) )
								//{
									//if ( ! isBertrandEndPoint( i,j,k ) )
									//{
										//queue2->prepend( i,j,k ) ;
									//}
								//}
							//}

						//}


				//if ( queue2->getNumElements() == 0 )
				//{
					//printf("Done.\n");
					//break ;
				//}
				//else
				//{
					//queue2->reset() ;
					//printf("%d nodes deleted.\n", queue2->getNumElements()) ;
					//gridQueueEle* ele = queue2->getNext() ;
					//while ( ele != NULL )
					//{
						//setDataAt( ele->x, ele->y, ele->z, -1 ) ;
						//ele = queue2->remove() ;
					//}
				//}
				//dir = ( dir + 1 ) % 6 ;
			//}

			//// Finally, clean up
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( 0, 0, 1 ) ;
		//}


		///* Palagyi's parallel surface thinning */
		//void Volume::palagyiSurfaceSkeleton( float thr )
		//{
			//int i, j, k ;

			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( thr, 0, 1 ) ;

			//// Next, initialize the templates
			//printf("Initializing surface endpoints templates...\n") ;
			//ThinningTemplate* US[6] ;

			//int b0[] = {12,13} ;
			//int w0[] = {2, 5, 8, 11, 14, 17, 20, 23, 26} ;
			//int ob0[] = {10, 16} ;
			//int ob20[] = {4, 22} ;
			//int nu[] = {0} ;
			//US[0] = new ThinningTemplate( b0, 2, w0, 9, ob0, 2, ob20, 2, nu, 0, nu, 0 ) ;
			//US[1] = new ThinningTemplate( US[0], 0, 1 ) ;
			///*
			//int b01[] = {13,16} ;
			//int w01[] = {0,1,2,9,10,11,18,19,20} ;
			//int ob01[] = {12, 14} ;
			//int ob201[] = {4, 22} ;
			//US[1] = new ThinningTemplate( b01, 2, w01, 9, ob01, 2, ob201, 2, nu, 0, nu, 0 ) ;
			//*/

			//int b1[] = {12,13,16,22} ;
			//int w1[] = {2,10,11,14} ;
			//int ow[] = {1,5} ;
			//US[2] = new ThinningTemplate( b1, 4, w1, 4, nu, 0, nu, 0, ow, 2, nu, 0 ) ;
			//US[3] = new ThinningTemplate( US[2], 0 ) ;

			//int b2[] = {2,12,13,16,22} ;
			//int w2[] = {10,11,14} ;
			//int op[] = {1,5} ;
			//US[4] = new ThinningTemplate( b2, 5, w2, 3, nu, 0, nu, 0, nu, 0, op, 2 ) ;
			//US[5] = new ThinningTemplate( US[4], 0 ) ;

			//ThinningTemplate * NE[6], * WD[6], * ES[6], * UW[6], * ND[6], * SW[6], * UN[6], * ED[6], * NW[6], * UE[6], * SD[6] ;

			//for ( i = 0 ; i < 6 ; i ++ )
			//{
				//SD[i] = new ThinningTemplate( US[i], 0, 1 ) ;
				//ND[i] = new ThinningTemplate( SD[i], 0, 1 ) ;
				//UN[i] = new ThinningTemplate( ND[i], 0, 1 ) ;

				//ES[i] = new ThinningTemplate( US[i], 1, 1 ) ;
				//NE[i] = new ThinningTemplate( ES[i], 2, 1 ) ;
				//NW[i] = new ThinningTemplate( NE[i], 2, 1 ) ;
				//SW[i] = new ThinningTemplate( NW[i], 2, 1 ) ;

				//UE[i] = new ThinningTemplate( US[i], 2, 1 ) ;
				//ED[i] = new ThinningTemplate( UE[i], 1, 1 ) ;
				//WD[i] = new ThinningTemplate( ED[i], 1, 1 ) ;
				//UW[i] = new ThinningTemplate( WD[i], 1, 1 ) ;
			//}

			//ThinningTemplate** alltemps[12] = { US, NE, WD, ES, UW, ND, SW, UN, ED, NW, UE, SD } ;

			//// Initialize queues
			//printf("Initializing queues...\n") ;
			//GridQueue2* queue2 = new GridQueue2( ) ;
			//GridQueue2* queue3 = new GridQueue2( ) ;
			//GridQueue2* queue4 = new GridQueue2( ) ;
			//Volume* fvol = new Volume( this->getSizeX(), this->getSizeY(), this->getSizeZ(), 0 ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) == 1 )
						//{
							//for ( int m = 0 ; m < 6 ; m ++ )
							//{
								//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) == 0 )
								//{
									//fvol->setDataAt( i, j, k, 1 ) ;
									//queue2->prepend( i, j, k ) ;
									//break ;
								//}
							//}
						//}
					//}
			//printf("Total %d nodes\n", queue2->getNumElements() ) ;

			//// Start iteration
			//int it = 0 ;
			//int vox[3][3][3] ;
			//gridQueueEle* ele ;
			//while ( queue2->getNumElements() > 0 )
			//{
				//printf("Iteration %d... nodes in queue: %d...", it, queue2->getNumElements()) ;

				//// queue2 holds all boundary nodes

				//int deleted = 0 ;
				//for ( i = 0 ; i < 12 ; i ++ )
				//{
					//// At the beginning of each iteration,
					//// queue2 holds all remaining boundary nodes
					//// queue3 is a deletable array, starting empty
					//// queue4 holds the candidates for next layer

					//queue2->reset() ;
					//ele = queue2->getNext()		;

					//// For each sub-iteration, go through queue2 first and find points to delete
					//while ( ele != NULL )
					//{
						//int ox = ele->x ;
						//int oy = ele->y ;
						//int oz = ele->z ;

						//// Check with templates
						//int match = 0 ;
						//for ( int ci = -1 ; ci < 2 ; ci ++ )
							//for ( int cj = -1 ; cj < 2 ; cj ++ )
								//for ( int ck = -1 ; ck < 2 ; ck ++ )
								//{
									//vox[ ci + 1 ][cj + 1][ck + 1] = (int)getDataAt( ox + ci, oy + cj, oz + ck ) ;
								//}

						//for ( j = 0 ; j < 6 ; j ++ )
						//// j = 1 ;
						//{
							//if ( alltemps[i][j]->isMatch( vox ) )
							//{
								///* Debug */
								//if ( ! isSimple2( vox ) )
								//{
									//printf("Wrong! %d %d\n", i, j) ;
									//for ( int cci = 0 ; cci < 3 ; cci ++ )
									//{
										//for ( int ccj = 0 ; ccj < 3 ; ccj ++ )
										//{
											//for ( int cck = 0 ; cck < 3 ; cck ++ )
											//{
												//printf("%d ",vox[ cci ][ccj][cck]);
											//}
											//printf(" , ") ;
										//}
										//printf("\n") ;
									//}
									//exit(0) ;
								//}

								//// Move this node to queue3
								//ele = queue2->remove() ;
								//queue3->prepend( ox, oy, oz ) ;

								//// Move its neighboring unvisited node to queue4
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;
									//if ( fvol->getDataAt( nx, ny, nz ) == 0 )
									//{
										//fvol->setDataAt( nx, ny, nz, 1 ) ;
										//queue4->prepend( nx, ny, nz ) ;
									//}
								//}

								//match = 1 ;
								//break ;
							//}
						//}

						//if ( match == 0 )
						//{
							//ele = queue2->getNext() ;
						//}
					//}

					//// Now, queue2 holds all remaining nodes,
					//// queue3 holds all deletable nodes,
					//// and queue4 holds nodes to be added to the next layer

					//// Simultaneous deletion
					//queue3->reset() ;
					//ele = queue3->getNext() ;
					//while ( ele != NULL )
					//{
						//setDataAt( ele->x, ele->y, ele->z, 0 ) ;
						//ele = queue3->remove() ;
						//deleted ++ ;
					//}

					//// return ;
				//}

				//// After all sub-iterations
				//// Move all queue4 nodes to queue2
				//queue4->reset() ;
				//ele = queue4->getNext() ;
				//while ( ele != NULL )
				//{
					//queue2->prepend( ele->x, ele->y, ele->z) ;
					//ele = queue4->remove() ;
				//}

				//if ( deleted == 0 )
				//{
					//printf("No more deletable nodes.\n");
					//break ;
				//}
				//else
				//{
					//printf("Deleted: %d\n", deleted) ;
				//}
			//}
		//}

		/**
		 * Normalize to a given range
		 */
		void Volume::threshold( double thr )
		{
			threshold( thr, 0, 1, 0, true) ;
		}

		void Volume::threshold( double thr, int out, int in )
		{
			threshold( thr, out, in, out, true) ;
		}

		void Volume::threshold( double thr, int out, int in, int boundary)
		{
			threshold(thr, out, in, boundary, true);
		}

		void Volume::threshold( double thr, int out, int in, int boundary, bool markBoundary)
		{
			float val;
			for ( int i = 0 ; i < getSizeX() ; i ++ )
				for ( int j = 0 ; j < getSizeY() ; j ++ )
					for ( int k = 0 ; k < getSizeZ() ; k ++ )
					{
						val = (float)getDataAt(i, j, k);
						if(markBoundary) {
							if ( i > 1 && i < getSizeX() - 2 && j > 1 && j < getSizeY() - 2 && k > 1 && k < getSizeZ() - 2 ) {
								if(val < thr) {
									setDataAt(i, j, k, out);
								} else {
									setDataAt(i, j, k, in);
								}
							}
							else
							{
								setDataAt(i, j, k, boundary);
							}
						} else {
							if(val < thr) {
								setDataAt(i, j, k, out);
							} else {
								setDataAt(i, j, k, in);
							}
						}
					}
		}

		//void Volume::threshold2( double thr, int out, int in )
		//{
			//for ( int i = 0 ; i < getSizeX() ; i ++ )
				//for ( int j = 0 ; j < getSizeY() ; j ++ )
					//for ( int k = 0 ; k < getSizeZ() ; k ++ ) {
						//double val = getDataAt(i, j, k);
						//if(val <= thr) {
							//setDataAt(i, j, k, out);
						//} else {
							//setDataAt(i, j, k, in);
						//}
					//}
		//}

		//void Volume::smooth( float alpha )
		//{
			//VolumeData * smoothedData = new VolumeData(getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, volData);

			//for (int i = 1; i < getSizeX() - 1; i++)
				//for (int j = 1; j < getSizeY() - 1; j++)
					//for (int k = 1; k < getSizeZ() - 1; k++) {
						//float v = (float)getDataAt( i - 1, j, k ) +
							//(float)getDataAt( i + 1, j, k ) +
							//(float)getDataAt( i, j - 1, k ) +
							//(float)getDataAt( i, j + 1, k ) +
							//(float)getDataAt( i, j, k - 1 ) +
							//(float)getDataAt( i, j, k + 1 ) ;
						//smoothedData->SetDataAt(i, j, k, smoothedData->GetDataAt(i, j, k) * alpha + ( 1 - alpha ) * v / 6);
					//}
			//delete volData;
			//volData = smoothedData;
		//}

		//void Volume::normalize( double min, double max )
		//{
			//double imin = getMin() ;
			//double imax = getMax() ;
			//double irange = imax - imin ;
			//double range = max - min ;

			//int size = volData->GetMaxIndex();
			//for(int i = 0 ; i < size ; i ++) {
				//setDataAt(i, ((getDataAt(i) - (float)imin ) / (float)irange) * (float)range + (float)min);
			//}
		//}

		//void Volume::normalize( double min, double max, double thresh, double ithresh )
		//{
			//double imin = getMin() ;
			//double imax = getMax() ;
			//double irange1 = ithresh - imin ;
			//double irange2 = imax - ithresh ;
			//double range1 = thresh - min;
			//double range2 = max - thresh ;

			//int size = volData->GetMaxIndex();
			//for (int i = 0; i < size; i++) {
				//if (getDataAt(i) < ithresh) {
					//setDataAt(i, ((getDataAt(i) - (float)imin ) / (float)irange1) * (float)range1 + (float)min);
				//}
				//else
				//{
					//setDataAt(i, (float)max - (( (float)imax - getDataAt(i)) / (float)irange2) * (float)range2);
				//}
			//}
		//}

		///* Set data at a pixel */

		//Volume * Volume::getDataRange(int x, int y, int z, int radius) {
			//Volume * range = new Volume(radius*2+1, radius*2+1, radius*2+1);
			//for(int xx = x-radius; xx <= x+radius; xx++) {
				//for(int yy = y-radius; yy <= y+radius; yy++) {
					//for(int zz = z-radius; zz <= z+radius; zz++) {
						//range->setDataAt(xx-x+radius, yy-y+radius, zz-z+radius, getDataAt(xx, yy, zz));
					//}
				//}
			//}
			//return range;
		//}

		///* Get data at an interpolated voxel */
		//double Volume::getInterpDataAt( double x, double y, double z )
		//{
			///*
			//double rad = getSizeX() / 4.0 ;
			//double cent = ( getSizeX() - 1 ) / 2.0 ;

			//double ox = x - cent ;
			//double oy = y - cent ;
			//double oz = z - cent ;

			//double a = -0.3 ;
			//double nx = ox ;
			//double ny = cos( a ) * oy + sin( a ) * oz ;
			//double nz = - sin( a ) * oy + cos( a ) * oz ;

			//double b = 1.4 ;
			//double nnx = cos( b ) * nx + sin( b ) * ny - 2;
			//double nny = -sin( b ) * nx + cos ( b ) * ny - 1;
			//double nnz = nz + 1;

			//double dis = nnx * nnx + nny * nny ;
			//return 10 - 10 * dis / ( rad * rad ) ;
			//*/

			//double rvalue ;
			//int hx = (int) ceil( x ) ;
			//int lx = (int) floor( x ) ;
			//int hy = (int) ceil( y ) ;
			//int ly = (int) floor( y ) ;
			//int hz = (int) ceil( z ) ;
			//int lz = (int) floor( z ) ;

			//double x1 = x - lx, x2 = 1 - x1 ;
			//double r1 = x2 * getDataAt( lx, ly, lz) + x1 * getDataAt( hx, ly, lz ) ;
			//double r2 = x2 * getDataAt( lx, ly, hz) + x1 * getDataAt( hx, ly, hz ) ;
			//double r3 = x2 * getDataAt( lx, hy, lz) + x1 * getDataAt( hx, hy, lz ) ;
			//double r4 = x2 * getDataAt( lx, hy, hz) + x1 * getDataAt( hx, hy, hz ) ;

			//double y1 = y - ly, y2 = 1 - y1 ;
			//double r5 = y2 * r1 + y1 * r3 ;
			//double r6 = y2 * r2 + y1 * r4 ;

			//double z1 = z - lz, z2 = 1 - z1 ;
			//rvalue = z2 * r5 + z1 * r6 ;

			//return rvalue ;
		//}

		///* Rotation routine */
		//void Volume::rotateX ( double a )
		//{
			//int i ;
			//int sizeX = getSizeX(), sizeY = getSizeY(), sizeZ = getSizeZ();

			//if ( sizeX != sizeY || sizeX != sizeZ) {
				//return ;
			//}

			//VolumeData * newData = new VolumeData(sizeX, sizeY, sizeZ, 0, 0, 0, volData);

			//double cent = ( sizeX - 1 ) / 2.0 ;
			//for ( i = 0 ; i < sizeX ; i ++ )
				//for ( int j = 0 ; j < sizeY ; j ++ )
					//for ( int k = 0 ; k < sizeZ ; k ++ )
					//{
						//double x = i - cent ;
						//double y = j - cent ;
						//double z = k - cent ;

						//double nx = x + cent ;
						//double ny = cos( a ) * y + sin( a ) * z + cent ;
						//double nz = - sin( a ) * y + cos( a ) * z + cent ;

						//if ( nx < 0 ) {
							//nx = 0 ;
						//} else if ( nx > sizeX - 1 ) {
							//nx = sizeX - 1 ;
						//}

						//if ( ny < 0 ) {
							//ny = 0 ;
						//} else if ( ny > sizeY - 1 ) {
							//ny = sizeY - 1 ;
						//}

						//if ( nz < 0 ) {
							//nz = 0 ;
						//} else if ( nz > sizeZ - 1 ) {
							//nz = sizeZ - 1 ;
						//}

						//newData->SetDataAt(i, j, k, (float)getInterpDataAt( nx, ny, nz ));
					//}

				//delete volData;
				//volData = newData;
		//}


		///* Write to file */
		//void Volume::toMathematicaFile( char* fname )
		//{
			//FILE* fout = fopen( fname, "w" ) ;

			//fprintf( fout, "{" ) ;
			//for ( int i = 0 ; i < getSizeX() ; i ++ )
			//{
				//fprintf( fout, "{" ) ;
				//for ( int j = 0 ; j < getSizeY() ; j ++ )
				//{
					//fprintf( fout, "{" ) ;
					//for ( int k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//fprintf( fout, "%.15f", getDataAt( i, j, k ) ) ;
						//if ( k < getSizeZ() - 1 )
						//{
							//fprintf( fout, "," ) ;
						//}
					//}
					//fprintf( fout, "}" ) ;
					//if ( j < getSizeY() - 1 )
					//{
						//fprintf( fout, ",\n" ) ;
					//} else {
						//fprintf( fout, "\n" ) ;
					//}
				//}
				//fprintf( fout, "}" ) ;
				//if ( i < getSizeX() - 1 )
				//{
					//fprintf( fout, ",\n\n\n" ) ;
				//} else {
					//fprintf( fout, "\n\n\n" ) ;
				//}
			//}
			//fprintf(fout,"}") ;

			//fclose( fout ) ;

		//}

		///* Write to file */
		//void Volume::toMathematicaFile( char* fname, int lx, int hx, int ly, int hy, int lz, int hz )
		//{
			//FILE* fout = fopen( fname, "w" ) ;

			//fprintf( fout, "{" ) ;
			//for ( int i = lx ; i < hx ; i ++ )
			//{
				//fprintf( fout, "{" ) ;
				//for ( int j = ly ; j < hy ; j ++ )
				//{
					//fprintf( fout, "{" ) ;
					//for ( int k = lz ; k < hz ; k ++ )
					//{
						//fprintf( fout, "%.15f", getDataAt( i, j, k ) ) ;
						//if ( k < hz - 1 )
						//{
							//fprintf( fout, "," ) ;
						//}
					//}
					//fprintf( fout, "}" ) ;
					//if ( j < hy - 1 )
					//{
						//fprintf( fout, "," ) ;
					//}
				//}
				//fprintf( fout, "}" ) ;
				//if ( i < hx - 1 )
				//{
					//fprintf( fout, "," ) ;
				//}
			//}
			//fprintf(fout,"}") ;

			//fclose( fout ) ;

		//}

		//void Volume::toOFFCells( char* fname )
		//{
			//toOFFCells( fname, 0.0001f ) ;
		//}
		//void Volume::toOFFCells2( char* fname )
		//{
			//toOFFCells2( fname, 0.0001f ) ;
		//}

		//void Volume::toOFFCells2( char* fname, float thr )
		//{
			//int i, j, k ;
			//Volume* indvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), -1 ) ;

			//// Get number of cells to write
			//int numverts = 0, numfaces = 0 ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= thr )
						//{
							//indvol->setDataAt( i,j,k, numverts ) ;
							//numverts ++ ;

							//for ( int mi = 0 ; mi < 3 ; mi ++ )
							//{
								//int find = mi * 4 + 3 ;
								//int isFace = 1 ;
								//for ( int mj = 0 ; mj < 4 ; mj ++ )
								//{
									//int nx = i + sheetNeighbor[find][mj][0] ;
									//int ny = j + sheetNeighbor[find][mj][1] ;
									//int nz = k + sheetNeighbor[find][mj][2] ;

									//if ( getDataAt( nx, ny, nz ) < thr )
									//{
										//isFace = 0 ;
										//break ;
									//}
								//}
								//if ( isFace )
								//{
									//numfaces ++ ;
								//}

								//int eind = mi * 2 + 1 ;
								//if ( getDataAt( i + neighbor6[eind][0], j + neighbor6[eind][1], k + neighbor6[eind][2]) >= thr )
								//{
									//numfaces ++ ;
								//}
							//}
						//}
					//}

			//FILE* fin = fopen( fname, "w" ) ;
			//fprintf( fin, "OFF\n" ) ;
			//fprintf( fin, "%d %d 0\n", numverts, numfaces ) ;

			//// Write vertices
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i,j,k ) >= thr )
						//{
							//fprintf( fin, "%d %d %d\n", i, j, k ) ;
						//}
					//}

			//// Write faces
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i,j,k ) >= thr )
						//{
							//int thisvt = (int)( indvol->getDataAt( i,j,k ) );
							//for ( int mi = 0 ; mi < 3 ; mi ++ )
							//{
								//int find = mi * 4 + 3 ;
								//int isFace = 1 ;
								//int vts[4] ;
								//for ( int mj = 0 ; mj < 4 ; mj ++ )
								//{
									//int nx = i + sheetNeighbor[find][mj][0] ;
									//int ny = j + sheetNeighbor[find][mj][1] ;
									//int nz = k + sheetNeighbor[find][mj][2] ;

									//vts[ mj ] = (int)( indvol->getDataAt( nx, ny, nz ) );

									//if ( getDataAt( nx, ny, nz ) < thr )
									//{
										//isFace = 0 ;
										//break ;
									//}
								//}
								//if ( isFace )
								//{
									//fprintf( fin, "4 %d %d %d %d\n", vts[0], vts[1], vts[3], vts[2] ) ;
								//}

								//int eind = mi * 2 + 1 ;
								//int mx = i + neighbor6[eind][0] ;
								//int my = j + neighbor6[eind][1] ;
								//int mz = k + neighbor6[eind][2] ;
								//int vt = (int)( indvol->getDataAt( mx, my, mz ) );
								//if ( getDataAt( mx, my, mz ) >= thr )
								//{
									//fprintf( fin, "4 %d %d %d %d\n", thisvt, thisvt, vt, vt ) ;
								//}
							//}
						//}
					//}

			//fclose( fin ) ;
			//delete indvol ;
		//}

		//void Volume::toOFFCells( char* fname, float thr )
		//{
			//int i, j, k ;
			//// Volume* indvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), -1 ) ;

			//// Get number of cells to write
			//int numverts = 0, numfaces = 0 ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i, j, k ) >= thr )
						//{
							//numverts += 8 ;
							//for ( int mi = 0 ; mi < 6 ; mi ++ )
							//{
								//if ( getDataAt( i + neighbor6[mi][0], j + neighbor6[mi][1], k + neighbor6[mi][2] ) < thr )
								//{
									//numfaces ++ ;
								//}
							//}
						//}
					//}

			//FILE* fin = fopen( fname, "w" ) ;
			//fprintf( fin, "OFF\n" ) ;
			//fprintf( fin, "%d %d 0\n", numverts, numfaces ) ;

			//// Write vertices
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i,j,k ) >= thr )
						//{
							//float ox = i - 0.5f ;
							//float oy = j - 0.5f ;
							//float oz = k - 0.5f ;

							//for ( int mi = 0 ; mi < 2 ; mi ++ )
								//for ( int mj = 0 ; mj < 2 ; mj ++ )
									//for ( int mk = 0 ; mk < 2 ; mk ++ )
									//{
										//fprintf( fin, "%f %f %f\n", ox + mi, oy + mj, oz + mk ) ;
									//}
						//}
					//}

			//// Write faces
			//int ct = 0 ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i,j,k ) >= thr )
						//{
							//for ( int mi = 0 ; mi < 6 ; mi ++ )
							//{
								//if ( getDataAt( i + neighbor6[mi][0], j + neighbor6[mi][1], k + neighbor6[mi][2] ) < thr )
								//{
									//fprintf( fin, "4 %d %d %d %d\n", cubeFaces[mi][0] + ct, cubeFaces[mi][1] + ct, cubeFaces[mi][2] + ct, cubeFaces[mi][3] + ct ) ;
								//}
							//}

							//ct += 8 ;
						//}
					//}

			//fclose( fin ) ;
			////delete indvol ;
		//}

		//void Volume::segment( float threshold, Volume* lowvol, Volume* highvol, char* mrcfile )
		//{
			//int i,j,k ;
			//Volume* segvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( lowvol->getDataAt(i,j,k) > 0 )
						//{
							//segvol->setDataAt( i,j,k, 1 ) ;
						//}
						//else if ( highvol->getDataAt(i,j,k) > 0 )
						//{
							//segvol->setDataAt( i,j,k, 2 ) ;
						//}
						//else
						//{
							//segvol->setDataAt( i,j,k, 0 ) ;
						//}
					//}

			//writeSegmentation( threshold, segvol, NULL, mrcfile ) ;
		//}

		//void Volume::segment( float threshold, Volume* vol, int maxDis, char* mrcfile )
		//{
			//int i,j;
			//Volume* testvol = NULL ;
			//Volume* disvol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
			//printf("Writing distance transform to %d levels.\n", maxDis);

			//int totNodes = 0 ;
			//int size = getSizeX() * getSizeY() * getSizeZ() ;
			//for ( i = maxDis ; i >= 0 ; i -- )
			//{
				//if ( i == 1 ) continue ;

				//int nodes = 0 ;
				//testvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, vol ) ;
				//testvol->erodeSheet( i ) ;

				//for ( j = 0 ; j < size ; j ++ )
				//{
					//if ( disvol->getDataAt(j) == 0 && testvol->getDataAt(j) > 0 )
					//{
						//disvol->setDataAt( j, i + 1 ) ;
						//nodes ++ ;
					//}
				//}
				//printf("Level %d has %d nodes.\n", i, nodes );
				//totNodes += nodes ;
				//delete testvol ;
			//}
			//printf("Totally %d nodes.\n", totNodes );

			//writeSegmentation( threshold, disvol, NULL, mrcfile ) ;
		//}

		//// Segment the volume using segvol
		//// background voxels have values 0, others have values > 0
		//void Volume::writeSegmentation( float threshold, Volume* segvol, char* txtfile, char* mrcfile )
		//{
			//printf("Start segmentation.\n") ;
			//int i,j,k ;
			//Volume* vvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;
			//Volume* tvol1 = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;
			//Volume* tvol2 = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;

			//// Initialize
			//GridQueue2* queue = new GridQueue2() ;
			//GridQueue2* queue2 = new GridQueue2() ;
			//GridQueue2* queue3 = new GridQueue2() ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt( i,j,k ) < threshold || segvol->getDataAt(i,j,k) <= 0 )
						//{
							//continue ;
						//}

						//vvol->setDataAt(i,j,k,1) ;
						//queue->prepend( i,j,k ) ;
					//}

			//// Dilation
			//printf("Dilation...") ;
			//int ox, oy, oz ;
			//gridQueueEle* ele ;
			//while ( queue->getNumElements() > 0 )
			//{
				//// At the beginning
				//// queue holds all nodes of previous layer
				//// queue2 is empty, ready to fill with the current layer

				//// First, fill current layer and assign values
				//queue->reset() ;
				//ele = queue->getNext() ;
				//while ( ele != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//double seg = segvol->getDataAt(ox,oy,oz) ;
					//int isBorder = 0 ;

					//for ( int m = 0 ; m < 6 ; m ++ )
					//{
						//int nx = ox + neighbor6[m][0] ;
						//int ny = oy + neighbor6[m][1] ;
						//int nz = oz + neighbor6[m][2] ;
						//if ( getDataAt( nx, ny, nz ) >= threshold && vvol->getDataAt( nx, ny, nz ) == 0 )
						//{
							//double ct = (int) tvol1->getDataAt(nx,ny,nz) ;
							//double val = tvol2->getDataAt(nx,ny,nz) ;

							//tvol1->setDataAt(nx,ny,nz, ct + 1 ) ;
							//tvol2->setDataAt(nx,ny,nz, val + seg ) ;

							//if (ct == 0)
							//{
								//queue2->prepend( nx, ny, nz ) ;
							//}
						//}
						//else if ( getDataAt( nx, ny, nz ) < threshold )
						//{
							//isBorder = 1 ;
						//}
					//}

					//if ( isBorder )
					//{
						//queue3->prepend( ox, oy, oz ) ;
					//}

					//ele = queue->remove() ;
				//}

				//// Normalize values of nodes in queue2
				//queue2->reset() ;
				//while ( ( ele = queue2->getNext() ) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;

					//double ct = (int) tvol1->getDataAt(ox,oy,oz) ;
					//double val = tvol2->getDataAt(ox,oy,oz) ;

					//if ( ct == 0 )
					//{
						//printf("Wrong! %f\n", ct) ;
					//}
					//segvol->setDataAt(ox,oy,oz, val / ct ) ;
					//vvol->setDataAt(ox,oy,oz, 1 ) ;
				//}

				//// Finally, swap queues
				//GridQueue2* temp = queue2 ;
				//queue2 = queue ;
				//queue = temp ;
			//}
			//printf("Done.\n") ;

			///* Debug
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//if ( getDataAt(i) >= threshold && segvol->getDataAt(i) == 0 )
				//{
					//segvol->setDataAt(i, 1) ;
				//}
				//else
				//{
					//segvol->setDataAt(i, 0) ;
				//}
			//}
			//*/

			//// Write to MRC
			////printf("Writing to %s...", mrcfile) ;
			////segvol->toMRCFile( mrcfile ) ;


			//// Growing out one layer for display
			//printf("Growing out one layer...\n") ;
			//// First, fill current layer and assign values
			//queue2->reset() ;
			//queue3->reset() ;
			//ele = queue3->getNext() ;
			//while ( ele != NULL )
			//{
				//ox = ele->x ;
				//oy = ele->y ;
				//oz = ele->z ;
				//double seg = segvol->getDataAt(ox,oy,oz) ;

				//for ( int mx = -1 ; mx < 2 ; mx ++ )
					//for ( int my = -1 ; my < 2 ; my ++ )
						//for ( int mz = -1 ; mz < 2 ; mz ++ )
				//{
					//int nx = ox + mx ; // neighbor6[m][0] ;
					//int ny = oy + my ; // neighbor6[m][1] ;
					//int nz = oz + mz ; // neighbor6[m][2] ;
					//if ( vvol->getDataAt( nx, ny, nz ) == 0 )
					//{
						//double ct = (int) tvol1->getDataAt(nx,ny,nz) ;
						//double val = tvol2->getDataAt(nx,ny,nz) ;

						//tvol1->setDataAt(nx,ny,nz, ct + 1 ) ;
						//tvol2->setDataAt(nx,ny,nz, val + seg ) ;

						//if (ct == 0)
						//{
							//queue2->prepend( nx, ny, nz ) ;
						//}
					//}
				//}
				//ele = queue3->remove() ;
			//}

			//// Normalize values of nodes in queue2
			//queue2->reset() ;
			//while ( ( ele = queue2->getNext() ) != NULL )
			//{
				//ox = ele->x ;
				//oy = ele->y ;
				//oz = ele->z ;

				//double ct = tvol1->getDataAt(ox,oy,oz) ;
				//double val = tvol2->getDataAt(ox,oy,oz) ;

				//if ( ct == 0 )
				//{
					//printf("Wrong! %f\n", ct) ;
				//}
				//segvol->setDataAt(ox,oy,oz, val / ct ) ;
			//}

			//printf("Writing...") ;
			//segvol->toMRCFile( mrcfile ) ;
			//segvol->toMRCFile( "../colors.mrc" ) ;
			//printf("Done.\n") ;

			//printf("Segmentation...") ;
			//for ( i = 0 ; i < getSizeX() * getSizeY() * getSizeZ() ; i ++ )
			//{
				//float segval = (float)segvol->getDataAt(i) ;
				//// High values
				//if ( segval > 1.5f )
				//{
					//tvol1->setDataAt(i, getDataAt(i) ) ;
				//}
				//else
				//{
					//tvol1->setDataAt(i, -1) ;
				//}

				//// Low values
				//if ( segval < 1.5f && segval >= 1 )
				//{
					//tvol2->setDataAt(i, getDataAt(i) ) ;
				//}
				//else
				//{
					//tvol2->setDataAt(i, -1) ;
				//}
			//}
			//char nname[1024] ;
			//sprintf( nname, "%s_sheet.mrc", mrcfile ) ;
			//tvol1->toMRCFile( nname ) ;
			//sprintf( nname, "%s_helix.mrc", mrcfile ) ;
			//tvol2->toMRCFile( nname ) ;
			//printf("Done.\n") ;
			//return ;

			///* Write to text
			//if ( txtfile != NULL )
			//{
				//printf("Writing to %s...", txtfile) ;
				//// Count border points
				//queue->reset() ;
				//for ( i = 0 ; i < getSizeX() ; i ++ )
					//for ( j = 0 ; j < getSizeY() ; j ++ )
						//for ( k = 0 ; k < getSizeZ() ; k ++ )
						//{
							//if ( getDataAt( i, j, k ) >= threshold )
							//{
								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < threshold )
									//{
										//queue->prepend( i, j, k ) ;
										//break ;
									//}
								//}
							//}
						//}

				//FILE* fout = fopen( txtfile, "w" ) ;
				//fprintf( fout, "%d\n", queue->getNumElements() ) ;
				//queue->reset() ;
				//while( (ele=queue->getNext()) != NULL )
				//{
					//ox = ele->x ;
					//oy = ele->y ;
					//oz = ele->z ;
					//fprintf( fout, "%d %d %d %f\n", ox, oy, oz, segvol->getDataAt(ox,oy,oz) ) ;
				//}
				//fclose( fout ) ;
				//printf("Done.\n") ;
			//}
			//*/
		//}

		//void Volume::floodFill( float thr )
		//{
			//int i;
			//// First, threshold the volume
			//#ifdef VERBOSE
			//printf("Thresholding the volume to 0/1...\n") ;
			//#endif
			//threshold( thr, 0, 1 ) ;

			//// Next, initialize queue
			//Volume* tvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;
			//GridQueue2* queue = new GridQueue2() ;
			//gridQueueEle* ele ;
			//queue->prepend( 0, 0, 0 ) ;
			//tvol->setDataAt( 0,0,0, -1 ) ;

			//// Iteration
			//printf("Flood filling...\n") ;
			//queue->reset() ;
			//ele = queue->getNext() ;
			//int ct = 1 ;
			//while( ele != NULL )
			//{
				//int ox = ele->x ;
				//int oy = ele->y ;
				//int oz = ele->z ;
				//queue->remove() ;

				//for ( int m = 0 ; m < 6 ; m ++ )
				//{
					//int nx = ox + neighbor6[m][0] ;
					//int ny = oy + neighbor6[m][1] ;
					//int nz = oz + neighbor6[m][2] ;
					//if ( nx < 0 || nx >= getSizeX() || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
					//{
						//continue ;
					//}
					//if ( tvol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 0 )
					//{
						//queue->prepend( nx, ny, nz ) ;
						//tvol->setDataAt(nx,ny,nz, -1) ;
						//ct ++ ;

						//if ( ct % 100000 == 0 )
						//{
							//printf("%d nodes processed.\n", ct);
						//}
					//}
				//}

				//queue->reset() ;
				//ele = queue->getNext() ;
			//}
			//printf("Done.\n") ;

			//// Done
			//for ( i = 0 ; i < getSizeX()*getSizeY()*getSizeZ() ; i ++ )
			//{
				//this->setDataAt(i, tvol->getDataAt(i) ) ;
			//}
		//}

		//void Volume::reduceComponent( int size )
		//{
			//int i, j, k ;
			//Volume* tvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;
			//GridQueue2* queue = new GridQueue2() ;
			//GridQueue2* queue2 = new GridQueue2() ;
			//gridQueueEle* ele ;
			//int numC = 0, numBC = 0 ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt(i,j,k) == 1 && tvol->getDataAt(i,j,k) == 0 )
						//{
							//numC ++ ;
							//// Start searching for this component
							//queue->prepend( i,j,k ) ;
							//queue2->prepend( i,j,k ) ;
							//tvol->setDataAt( i,j,k, 1 ) ;
							//int ct = 1 ;

							//queue->reset() ;
							//ele = queue->getNext() ;
							//while( ele != NULL )
							//{
								//int ox = ele->x ;
								//int oy = ele->y ;
								//int oz = ele->z ;
								//queue->remove() ;

								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;
									//if ( nx < 0 || nx >= getSizeX()  || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
									//{
										//continue ;
									//}
									//if ( tvol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 1 )
									//{
										//queue->prepend( nx, ny, nz ) ;
										//queue2->prepend( nx, ny, nz ) ;
										//tvol->setDataAt( nx, ny, nz, 1 ) ;
										//ct ++ ;
									//}
								//}

								//queue->reset() ;
								//ele = queue->getNext() ;
							//}

							//if ( ct < size )
							//{
								//// remove component

								//queue2->reset() ;
								//ele = queue2->getNext() ;
								//while ( ele != NULL )
								//{
									//setDataAt( ele->x, ele->y, ele->z, 0 ) ;
									//ele = queue2->remove() ;
								//}
								//queue2->reset() ;

							//}
							//else
							//{
								//queue2 = new GridQueue2() ;
								//numBC ++ ;
							//}
						//}
					//}

			//printf("Total number of components: %d Remained: %d\n", numC, numBC ) ;

		//}

		//void Volume::reduceComponent2( int num )
		//{
			//int i, j, k ;
			//Volume* tvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;
			//GridQueue2* queue = new GridQueue2() ;
			//GridQueue2* queue2 = new GridQueue2() ;
			//gridQueueEle* ele ;
			//int numC = 0, numBC = 0 ;
			//int* tops = new int[ num ] ;
			//int* topinds = new int[ num ] ;
			//for ( i = 0 ; i < num ; i ++ )
			//{
				//tops[i] = 0 ;
				//topinds[i] = -1 ;
			//}

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt(i,j,k) == 1 && tvol->getDataAt(i,j,k) == 0 )
						//{
							//numC ++ ;
							//// Start searching for this component
							//queue->prepend( i,j,k ) ;
							//queue2->prepend( i,j,k ) ;
							//tvol->setDataAt( i,j,k, 1 ) ;
							//int ct = 1 ;

							//queue->reset() ;
							//ele = queue->getNext() ;
							//while( ele != NULL )
							//{
								//int ox = ele->x ;
								//int oy = ele->y ;
								//int oz = ele->z ;
								//queue->remove() ;

								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;
									//if ( nx < 0 || nx >= getSizeX()  || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
									//{
										//continue ;
									//}
									//if ( tvol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 1 )
									//{
										//queue->prepend( nx, ny, nz ) ;
										//queue2->prepend( nx, ny, nz ) ;
										//tvol->setDataAt( nx, ny, nz, 1 ) ;
										//ct ++ ;
									//}
								//}

								//queue->reset() ;
								//ele = queue->getNext() ;
							//}

							//queue2->reset() ;
							//ele = queue2->getNext() ;
							//while ( ele != NULL )
							//{
								//ele = queue2->remove() ;
							//}
							//queue2->reset() ;

							//for ( int ind = 0 ; ind < num ; ind ++ )
							//{
								//if ( ct > tops[ ind ] )
								//{
									//for ( int nind = num - 1 ; nind > ind ; nind -- )
									//{
										//tops[nind] = tops[ nind - 1 ] ;
										//topinds[nind] = topinds[ nind - 1 ] ;
									//}
									//tops[ ind ] = ct ;
									//topinds[ ind ] = numC ;
									//break ;
								//}
							//}
						//}
					//}

			//printf("%d components total, %d selected components:\n", numC, num) ;
			//for ( i = 0 ; i < num ; i ++ )
			//{
				//printf("Id %d: size %d\n", topinds[i], tops[i] ) ;
			//}

			//tvol->fill( 0 ) ;

			//numC = 0 ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( getDataAt(i,j,k) == 1 && tvol->getDataAt(i,j,k) == 0 )
						//{
							//numC ++ ;
							//// Start searching for this component
							//queue->prepend( i,j,k ) ;
							//queue2->prepend( i,j,k ) ;
							//tvol->setDataAt( i,j,k, 1 ) ;
							//int ct = 1 ;

							//queue->reset() ;
							//ele = queue->getNext() ;
							//while( ele != NULL )
							//{
								//int ox = ele->x ;
								//int oy = ele->y ;
								//int oz = ele->z ;
								//queue->remove() ;

								//for ( int m = 0 ; m < 6 ; m ++ )
								//{
									//int nx = ox + neighbor6[m][0] ;
									//int ny = oy + neighbor6[m][1] ;
									//int nz = oz + neighbor6[m][2] ;
									//if ( nx < 0 || nx >= getSizeX()  || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
									//{
										//continue ;
									//}
									//if ( tvol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 1 )
									//{
										//queue->prepend( nx, ny, nz ) ;
										//queue2->prepend( nx, ny, nz ) ;
										//tvol->setDataAt( nx, ny, nz, 1 ) ;
										//ct ++ ;
									//}
								//}

								//queue->reset() ;
								//ele = queue->getNext() ;
							//}

							//int removing = 1 ;
							//for ( int ind = 0 ; ind < num ; ind ++ )
							//{
								//if ( topinds[ ind ] == numC )
								//{
									//removing = 0 ;
									//break ;
								//}
							//}

							//if ( removing )
							//{
								//// remove component

								//queue2->reset() ;
								//ele = queue2->getNext() ;
								//while ( ele != NULL )
								//{
									//setDataAt( ele->x, ele->y, ele->z, 0 ) ;
									//ele = queue2->remove() ;
								//}
								//queue2->reset() ;

							//}
							//else
							//{
								//queue2 = new GridQueue2() ;
								//numBC ++ ;
							//}
						//}
					//}

			//printf("Total number of components: %d Remained: %d\n", numC, numBC ) ;
			//delete tvol ;
		//}

		//void Volume::floodFillPQR( int offset )
		//{
			//int i;

			//// Next, initialize queue
			//GridQueue2* queue = new GridQueue2() ;
			//gridQueueEle* ele ;
			//queue->prepend( 0, 0, 0 ) ;

			//// Iteration
			//printf("Flood filling outside from  (0,0,0)...\n") ;
			//int ct = 1 ;
			///*
			//queue->reset() ;
			//ele = queue->getNext() ;
			//while( ele != NULL )
			//{
				//int ox = ele->x ;
				//int oy = ele->y ;
				//int oz = ele->z ;
				//queue->remove() ;

				//int boundary = 0 ;
				//for ( int m = 0 ; m < 6 ; m ++ )
				//{
					//int nx = ox + neighbor6[m][0] ;
					//int ny = oy + neighbor6[m][1] ;
					//int nz = oz + neighbor6[m][2] ;
					//if ( nx < 0 || nx >= getSizeX()  || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
					//{
						//continue ;
					//}
					//if ( tvol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 0 )
					//{
						//queue->prepend( nx, ny, nz ) ;
						//tvol->setDataAt( nx, ny, nz, 1 ) ;
						//ct ++ ;
						//if ( ct % 100000 == 0 )
						//{
							//printf("%d nodes processed.\n", ct);
						//}
					//}
					//else if ( getDataAt( nx, ny, nz ) == 1 )
					//{
						//boundary = 1 ;
					//}
				//}

				//if ( boundary )
				//{
					//tvol->setDataAt( ox, oy, oz, 2 ) ;
				//}

				//queue->reset() ;
				//ele = queue->getNext() ;
			//}
			//printf("Done.\n") ;
			//for ( i = 0 ; i < getSizeX()*getSizeY()*getSizeZ() ; i ++ )
			//{

				//if ( tvol->getDataAt(i) == 2 )
				//{
					//setDataAt( i, 1 ) ;
				//}
			//}
			//*/


			//// Find inside seed point
			//int maxRounds = 1 ;
			//for ( int rounds = 0 ; rounds < maxRounds ; rounds ++ )
			//{

			//int isSolid = 0 ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
			//{
				//if ( getDataAt( i, getSizeY() / 2, getSizeZ() / 2 ) == 1 )
				//{
					//isSolid = 1 ;
				//}
				//else if ( isSolid == 1 )
				//{
					//break ;
				//}
			//}

			//Volume* invol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0 ) ;
			//queue->prepend( i, getSizeY()/2, getSizeZ()/2 ) ;
			//invol->setDataAt( i, getSizeY()/2, getSizeZ()/2, 1 ) ;
			//printf("Flood filling inside from  (%d,%d,%d)...\n",i, getSizeY()/2, getSizeZ()/2) ;

			//// Iteration
			//queue->reset() ;
			//ele = queue->getNext() ;
			//ct = 1 ;
			//while( ele != NULL )
			//{
				//int ox = ele->x ;
				//int oy = ele->y ;
				//int oz = ele->z ;
				//queue->remove() ;

				//int boundary = 0 ;
				//for ( int m = 0 ; m < 6 ; m ++ )
				//{
					//int nx = ox + neighbor6[m][0] ;
					//int ny = oy + neighbor6[m][1] ;
					//int nz = oz + neighbor6[m][2] ;
					//if ( nx < 0 || nx >= getSizeX()  || ny < 0 || ny >= getSizeY() || nz < 0 || nz >= getSizeZ() )
					//{
						//continue ;
					//}
					//if ( invol->getDataAt( nx, ny, nz ) == 0 && getDataAt( nx, ny, nz ) == 0 )
					//{
						//queue->prepend( nx, ny, nz ) ;
						//invol->setDataAt( nx, ny, nz, 1 ) ;
						//ct ++ ;
						//if ( ct % 100000 == 0 )
						//{
							//printf("%d nodes processed.\n", ct);
						//}
					//}
					//else if ( getDataAt( nx, ny, nz ) == 1 )
					//{
						//boundary = 1 ;
					//}
				//}

				//if ( boundary )
				//{
					//invol->setDataAt( ox, oy, oz, 2 ) ;
				//}

				//queue->reset() ;
				//ele = queue->getNext() ;
			//}
			//printf("Done.\n") ;

			//// Done
			//for ( i = 0 ; i < getSizeX()*getSizeY()*getSizeZ() ; i ++ )
			//{
	///*
				//if ( tvol->getDataAt(i) == 2 )
				//{
					//setDataAt( i, 1 ) ;
				//}
				//*/
				//if ( invol->getDataAt(i) == 2 )
				//{
					//setDataAt( i, 1 ) ;
				//}

	///*
				//else if ( tvol->getDataAt(i) == 0 && invol->getDataAt(i) == 0 )
				//{
					//setDataAt( i, 1 ) ;
				//}
	//*/
			//}

					//delete invol ;

			//}

	////		delete tvol ;

	//}


		//void Volume::writeDistances( char* fname, int maxDis )
		//{
			//int i, j, k ;
			//Volume* testvol = NULL ;
			//Volume* disvol = new Volume( getSizeX(), getSizeY(), getSizeZ() ) ;
			//printf("Writing distance transform to %d levels.\n", maxDis);

			//int totNodes = 0 ;
			//int size = getSizeX() * getSizeY() * getSizeZ() ;
			//float score = 10 ;
			//for ( i = maxDis ; i >= 0 ; i -- )
			//{
				//int nodes = 0 ;
				//testvol = new Volume( getSizeX(), getSizeY(), getSizeZ(), 0, 0, 0, this ) ;
				//testvol->erodeSheet( i ) ;

				//for ( j = 0 ; j < size ; j ++ )
				//{
					//if ( disvol->getDataAt(j) == 0 && testvol->getDataAt(j) > 0 )
					//{
						//disvol->setDataAt( j, i + 1 ) ;
						//nodes ++ ;
					//}
				//}
				//printf("Level %d has %d nodes.\n", i, nodes );
				//totNodes += nodes ;
				//score -= 0.5f ;
				//delete testvol ;
			//}
			//printf("Totally %d nodes.\n", totNodes );

			////disvol->toMRCFile( "..\distance.mrc" ) ;
			//FILE* fout = fopen( fname, "w" ) ;
			//fprintf( fout, "%d\n", totNodes ) ;
			//int ct = 0 ;
			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//float val = (float)disvol->getDataAt(i,j,k) ;
						//if ( val > 0 )
						//{
							//fprintf( fout, "%d %d %d %f\n", i, j, k, val ) ;
							//ct ++ ;
						//}
					//}

			//if ( ct != totNodes )
			//{
				//printf("Counting wrong! %d %d\n", totNodes, ct) ;
			//}
			//fclose( fout ) ;

		//}

		//void Volume::toPQRFile( char* fname, float spc, float minx, float miny, float minz, int padding )
		//{
			//FILE* fout = fopen( fname, "w" ) ;
			//int i, j, k ;

			//for ( i = 0 ; i < getSizeX() ; i ++ )
				//for ( j = 0 ; j < getSizeY() ; j ++ )
					//for ( k = 0 ; k < getSizeZ() ; k ++ )
					//{
						//if ( i < padding || i >= getSizeX() - padding || j < padding || j >= getSizeY() - padding || k < padding || k >= getSizeZ() - padding )
						//{
							//continue ;
						//}
						//float val = (float)this->getDataAt(i,j,k) ;
						//if ( val > 0 )
						//{
							//float x = (i - padding) * spc + minx ;
							//float y = (j - padding) * spc + miny ;
							//float z = (k - padding) * spc + minz ;
							//fprintf( fout, "ATOM      1  X   DUM     1     %4.3f %4.3f %4.3f 0.000 0.000\n", x, y, z ) ;
						//}
					//}

			//fclose( fout ) ;
		//}

		//void Volume::toMRCFile( char* fname )
		//{
			//FILE* fout = fopen( fname, "wb" ) ;

			//// Write header
			//int sizeX = getSizeX();
			//int sizeY = getSizeY();
			//int sizeZ = getSizeZ();
			//fwrite( &sizeX, sizeof( int ), 1, fout ) ;
			//fwrite( &sizeY, sizeof( int ), 1, fout ) ;
			//fwrite( &sizeZ, sizeof( int ), 1, fout ) ;

			//int mode = 2 ;
			//fwrite( &mode, sizeof ( int ), 1, fout ) ;

			//int off[3] = {0,0,0} ;
			//int intv[3] = { getSizeX() - 1, getSizeY() - 1, getSizeZ() - 1 } ;
			//fwrite( off, sizeof( int ), 3, fout ) ;
			//fwrite( intv, sizeof( int ), 3, fout ) ;

			//float cella[3] = {getSpacingX() * (float)(getSizeX()), getSpacingY() * (float)(getSizeY()), getSpacingZ() * (float)(getSizeZ())};
			//float cellb[3] = {90,90,90} ;
			//fwrite( cella, sizeof( float ), 3, fout ) ;
			//fwrite( cellb, sizeof( float ), 3, fout ) ;

			//int cols[3] = {1,2,3} ;
			//fwrite( cols, sizeof( int ), 3, fout ) ;

			//double dmin = 100000, dmax = -100000 ;
			//int i ;
			//int size = volData->GetMaxIndex();
			//for (i = 0 ; i < size; i++) {
				//float val = (float)getDataAt(i);
				//if (val < dmin) {
					//dmin = val;
				//}
				//if (val > dmax) {
					//dmax = val;
				//}
			//}
			//float ds[3] = {(float)dmin, (float)dmax, (float)0} ;
			//fwrite( ds, sizeof( float ), 3, fout ) ;

			//int zero = 0 ;
			//for (i = 23 ; i <= 49 ; i ++ )
			//{
				//fwrite( &zero, sizeof( int ), 1, fout );
			//}

			//float origins[3];
			//origins[0] = getOriginX() / getSpacingX() + 0.5f * (float)getSizeX();
			//origins[1] = getOriginY() / getSpacingY() + 0.5f * (float)getSizeY();
			//origins[2] = getOriginZ() / getSpacingZ() + 0.5f * (float)getSizeZ();

			//fwrite( origins, sizeof( float ), 3, fout) ;

			//for (i = 53 ; i <= 256 ; i ++ )
			//{
				//fwrite( &zero, sizeof( int ), 1, fout ) ;
			//}

			//// Write contents
			//for ( int z = 0 ; z < getSizeZ() ; z ++ )
				//for ( int y = 0 ; y < getSizeY() ; y ++ )
					//for ( int x = 0 ; x < getSizeX() ; x ++ )
					//{
						//float d = (float)getDataAt(x,y,z) ;
						//fwrite( &d, sizeof( float ), 1, fout ) ;
					//}

			//fclose( fout ) ;
		//}
