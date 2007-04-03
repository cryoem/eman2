/**
 * Volumetric data definition
 *
 * Author: Tao Ju
 * Date: 02/16/2005
 */


#ifndef VOLUME_H
#define VOLUME_H

#include <cstdio>
#include <cstdlib>

namespace EMAN
{

	class Volume
	{
	public:
		/* Initialization */
		Volume( int x, int y, int z ) 
		{
			sizex = x ;
			sizey = y ; 
			sizez = z ;
	
			spcx = 1 ;
			spcy = 1 ;
			spcz = 1 ;
	
			data = new double [ x * y * z ] ;
			for ( int i = 0 ; i < x * y * z ; i ++ )
			{
				data[ i ] = 0 ;
			}
		};
	
		Volume( int x, int y, int z, float val ) 
		{
			sizex = x ;
			sizey = y ; 
			sizez = z ;
	
			spcx = 1 ;
			spcy = 1 ;
			spcz = 1 ;
	
			data = new double [ x * y * z ] ;
			for ( int i = 0 ; i < x * y * z ; i ++ )
			{
				data[ i ] = val ;
			}
		};
	
		Volume( int x, int y, int z, int offx, int offy, int offz, Volume* vol )
		{
			sizex = x ;
			sizey = y ; 
			sizez = z ;
	
			spcx = vol->getSpacingX() ;
			spcy = vol->getSpacingY() ;
			spcz = vol->getSpacingZ() ;
	
			data = new double [ x * y * z ] ;
	
			int ct = 0 ;
			for ( int i = offx ; i < x + offx; i ++ )
				for ( int j = offy ; j < y + offy; j ++ )
					for ( int k = offz ; k < z + offz; k ++ )
					{
						data[ ct ] = vol->getDataAt( i, j, k ) ;
						ct ++ ;
					}
		}
	
		/* Statistics function */
		int getSizeX( )
		{
			return sizex ;
		}
		int getSizeY( )
		{
			return sizey ;
		}
		int getSizeZ( )
		{
			return sizez ;
		}
	
		void setSpacing( float sx, float sy, float sz )
		{
			spcx = sx ;
			spcy = sy ;
			spcz = sz ;
		}
	
		float getSpacingX( )
		{
			return spcx ;
		}
		float getSpacingY( )
		{
			return spcy ;
		}
		float getSpacingZ( )
		{
			return spcz ;
		}
	
		double getMin()
		{
			int size = sizex * sizey * sizez ;
			double rvalue = data[0] ;
			for ( int i = 1 ; i < size ; i ++ )
			{
				if ( rvalue > data[ i ] )
				{
					rvalue = data[ i ] ;
				}
			}
			return rvalue ;
		}
	
		double getMax()
		{
			int size = sizex * sizey * sizez ;
			double rvalue = data[0] ;
			for ( int i = 1 ; i < size ; i ++ )
			{
				if ( rvalue < data[ i ] )
				{
					rvalue = data[ i ] ;
				}
			}
			return rvalue ;
		}
	
		/**
		 * Normalize to a given range 
		 */
		void threshold( double threshold )
		{
			int size = sizex * sizey * sizez ;
			for ( int i = 0 ; i < size ; i ++ )
			{
				data[ i ] = data[ i ] > threshold ? 1 : 0 ;
			}
		}
	
		void smooth( float alpha )
		{
			int size = sizex * sizey * sizez ;
			float* ndata = new float[ size ] ;
			for ( int i = 0 ; i < size ; i ++ )
			{
				ndata[ i ] = data[ i ] ;
			}
	
			for ( int i = 1 ; i < sizex - 1 ; i ++ )
				for ( int j = 1 ; j < sizey - 1 ; j ++ )
					for ( int k = 1 ; k < sizez - 1 ; k ++ )
					{
						int ct =  i * sizey * sizez + j * sizez + k ;
	
						float v = getDataAt( i - 1, j, k ) + 
							getDataAt( i + 1, j, k ) +
							getDataAt( i, j - 1, k ) +
							getDataAt( i, j + 1, k ) +
							getDataAt( i, j, k - 1 ) +
							getDataAt( i, j, k + 1 ) ;
						ndata[ ct ] = ndata[ ct ] * alpha + ( 1 - alpha ) * v / 6 ;
					}
	
			for ( int i = 0 ; i < size ; i ++ )
			{
				data[ i ] = ndata[ i ] ;
			}
	
			delete ndata ;
		}
	
		void expand( float thresh )
		{
			printf("Expanding... %d %d %d\n", sizex, sizey, sizez) ;
			int size = sizex * sizey * sizez ;
			char* flag = new char[ size ] ;
			int num = 0 ;
			for ( int i = 0 ; i < size ; i ++ )
			{
				flag[ i ] = ( data[ i ] > thresh? 1 : 0 ) ;
				if ( flag[i] )
				{
					num ++ ;
				}
			}
			printf("Number: %d\n", num) ;
	
			num = 0 ;
			int num2 = 0 ;
			int ct = 0 ;
			for ( int i = 0 ; i < sizex ; i ++ )
				for ( int j = 0 ; j < sizey ; j ++ )
					for ( int k = 0 ; k < sizez ; k ++ )
					{
	
						if ( flag[ ct ] )
						{
							num2 ++ ;
							float v = getDataAt( i, j, k ) ; // data[ ct ] ;
							if ( v <= thresh )
							{
								printf("%f\n", v) ;
							}
							for ( int x = -1 ; x < 2 ; x ++ )
								for ( int y = -1 ; y < 2 ; y ++ )
									for ( int z = -1 ; z < 2 ; z ++ )
									{
										if ( getDataAt( i + x, j + y, k + z ) <= thresh )
										{
											setDataAt( i + x, j + y, k + z, v ) ;
											num ++ ;
										}
									}
						}
	
						ct ++ ;
					}
			printf("Number: %d Number2: %d\n", num, num2) ;
	
			delete flag ;
		}
	
		void normalize( double min, double max )
		{
			double imin = getMin() ;
			double imax = getMax() ;
			double irange = imax - imin ;
			double range = max - min ;
	
			int size = sizex * sizey * sizez ;
			for ( int i = 0 ; i < size ; i ++ )
			{
				data[ i ] = (( data[ i ] - imin ) / irange) * range + min ;
			}
		}
	
		void normalize( double min, double max, double thresh, double ithresh )
		{
			double imin = getMin() ;
			double imax = getMax() ;
			double irange1 = ithresh - imin ;
			double irange2 = imax - ithresh ;
			double range1 = thresh - min;
			double range2 = max - thresh ;
	
			int size = sizex * sizey * sizez ;
			for ( int i = 0 ; i < size ; i ++ )
			{
				if ( data[ i ] < ithresh )
				{
					data[ i ] = (( data[ i ] - imin ) / irange1) * range1 + min ;
				}
				else
				{
					data[ i ] = max - (( imax - data[ i ] ) / irange2) * range2 ;
				}
			}
		}
	
		/* Set data at a pixel */
		void setDataAt( int x, int y, int z, double d )
		{
			data[ x * sizey * sizez + y * sizez + z ] = d ;
		}
	
		/* Get data at a single voxel */
		double getDataAt( int x, int y, int z ) 
		{
			return data[ x * sizey * sizez + y * sizez + z ] ;
		}
		double getDataAt( int index ) 
		{
			return data[ index ] ;
		}
		
		/* Get data at an interpolated voxel */
		double getInterpDataAt( double x, double y, double z ) 
		{
			/*
			double rad = sizex / 4.0 ;
			double cent = ( sizex - 1 ) / 2.0 ;
			
			double ox = x - cent ;
			double oy = y - cent ;
			double oz = z - cent ;
			
			double a = -0.3 ;
			double nx = ox ;
			double ny = cos( a ) * oy + sin( a ) * oz ;
			double nz = - sin( a ) * oy + cos( a ) * oz ;
	
			double b = 1.4 ;
			double nnx = cos( b ) * nx + sin( b ) * ny - 2;
			double nny = -sin( b ) * nx + cos ( b ) * ny - 1;
			double nnz = nz + 1;
			
			double dis = nnx * nnx + nny * nny ;
			return 10 - 10 * dis / ( rad * rad ) ;
			*/
	
			double rvalue ;
			int hx = (int) ceil( x ) ;
			int lx = (int) floor( x ) ;
			int hy = (int) ceil( y ) ;
			int ly = (int) floor( y ) ;
			int hz = (int) ceil( z ) ;
			int lz = (int) floor( z ) ;
			
			double x1 = x - lx, x2 = 1 - x1 ;
			double r1 = x2 * getDataAt( lx, ly, lz) + x1 * getDataAt( hx, ly, lz ) ; 
			double r2 = x2 * getDataAt( lx, ly, hz) + x1 * getDataAt( hx, ly, hz ) ; 
			double r3 = x2 * getDataAt( lx, hy, lz) + x1 * getDataAt( hx, hy, lz ) ; 
			double r4 = x2 * getDataAt( lx, hy, hz) + x1 * getDataAt( hx, hy, hz ) ; 
			
			double y1 = y - ly, y2 = 1 - y1 ;
			double r5 = y2 * r1 + y1 * r3 ;
			double r6 = y2 * r2 + y1 * r4 ;
	
			double z1 = z - lz, z2 = 1 - z1 ;
			rvalue = z2 * r5 + z1 * r6 ;
	
			return rvalue ;
		}
	
		/* Rotation routine */
		void rotateX ( double a )
		{
			double * ndata = new double[ sizex * sizey * sizez ] ;
			if ( sizex != sizey || sizex != sizez )
			{
				return ;
			}
	
			int ct = 0 ;
			double cent = ( sizex - 1 ) / 2.0 ;
			for ( int i = 0 ; i < sizex ; i ++ )
				for ( int j = 0 ; j < sizey ; j ++ )
					for ( int k = 0 ; k < sizez ; k ++ )
					{
						double x = i - cent ;
						double y = j - cent ;
						double z = k - cent ;
	
						double nx = x + cent ;
						double ny = cos( a ) * y + sin( a ) * z + cent ;
						double nz = - sin( a ) * y + cos( a ) * z + cent ;
	
						if ( nx < 0 )
						{
							nx = 0 ;
						}
						else if ( nx > sizex - 1 )
						{
							nx = sizex - 1 ;
						}
	
						if ( ny < 0 )
						{
							ny = 0 ;
						}
						else if ( ny > sizey - 1 )
						{
							ny = sizey - 1 ;
						}
	
						if ( nz < 0 )
						{
							nz = 0 ;
						}
						else if ( nz > sizez - 1 )
						{
							nz = sizez - 1 ;
						}
						
						ndata[ ct ] = getInterpDataAt( nx, ny, nz );
	
						ct ++ ;
					}
	
				for ( int i = 0 ; i < sizex * sizey * sizez ; i ++ )
				{
					data[ct] = ndata[ct] ;
				}
	
				delete ndata ;
		}
	
	
		/* Destruction */
		~Volume( )
		{
			delete data ;
		}
	
		/* Write to file */
		void toMathematicaFile( char* fname )
		{
			FILE* fout = fopen( fname, "w" ) ;
	
			fprintf( fout, "{" ) ;
			for ( int i = 0 ; i < sizex ; i ++ )
			{
				fprintf( fout, "{" ) ;
				for ( int j = 0 ; j < sizey ; j ++ )
				{
					fprintf( fout, "{" ) ;
					for ( int k = 0 ; k < sizez ; k ++ )
					{
						fprintf( fout, "%.15f", getDataAt( i, j, k ) ) ;
						if ( k < sizez - 1 )
						{
							fprintf( fout, "," ) ;
						}
					}
					fprintf( fout, "}" ) ;
					if ( j < sizey - 1 )
					{
						fprintf( fout, "," ) ;
					}
				}
				fprintf( fout, "}" ) ;
				if ( i < sizex - 1 )
				{
					fprintf( fout, "," ) ;
				}
			}
			fprintf(fout,"}") ;
	
			fclose( fout ) ;
	
		}
	
		/* Write to file */
		void toMathematicaFile( char* fname, int lx, int hx, int ly, int hy, int lz, int hz )
		{
			FILE* fout = fopen( fname, "w" ) ;
	
			fprintf( fout, "{" ) ;
			for ( int i = lx ; i < hx ; i ++ )
			{
				fprintf( fout, "{" ) ;
				for ( int j = ly ; j < hy ; j ++ )
				{
					fprintf( fout, "{" ) ;
					for ( int k = lz ; k < hz ; k ++ )
					{
						fprintf( fout, "%.15f", getDataAt( i, j, k ) ) ;
						if ( k < hz - 1 )
						{
							fprintf( fout, "," ) ;
						}
					}
					fprintf( fout, "}" ) ;
					if ( j < hy - 1 )
					{
						fprintf( fout, "," ) ;
					}
				}
				fprintf( fout, "}" ) ;
				if ( i < hx - 1 )
				{
					fprintf( fout, "," ) ;
				}
			}
			fprintf(fout,"}") ;
	
			fclose( fout ) ;
	
		}
	
		void toMRCFile( char* fname )
		{
			FILE* fout = fopen( fname, "wb" ) ;
	
			// Write header
			fwrite( &sizex, sizeof( int ), 1, fout ) ;
			fwrite( &sizey, sizeof( int ), 1, fout ) ;
			fwrite( &sizez, sizeof( int ), 1, fout ) ;
	
			int mode = 2 ;
			fwrite( &mode, sizeof ( int ), 1, fout ) ;
			
			int off[3] = {0,0,0} ;
			int intv[3] = { sizex - 1, sizey - 1, sizez - 1 } ;
			fwrite( off, sizeof( int ), 3, fout ) ;
			fwrite( intv, sizeof( int ), 3, fout ) ;
	
			float cella[3] = {2,2,2} ;
			float cellb[3] = {90,90,90} ;
			fwrite( cella, sizeof( float ), 3, fout ) ;
			fwrite( cellb, sizeof( float ), 3, fout ) ;
	
			int cols[3] = {1,2,3} ;
			fwrite( cols, sizeof( int ), 3, fout ) ;
	
			double dmin = 100000, dmax = -100000 ;
			for ( int i = 0 ; i < sizex * sizey * sizez ; i ++ )
			{
				if ( data[ i ] < dmin )
				{
					dmin = data[ i ] ;
				}
				if ( data[i] > dmax )
				{
					dmax = data[ i ] ;
				}
			}
			float ds[3] = {dmin, dmax, 0} ;
			fwrite( ds, sizeof( float ), 3, fout ) ;
	
			int zero = 0 ;
			for ( int i = 22 ; i < 256 ; i ++ )
			{
				fwrite( &zero, sizeof( int ), 1, fout ) ;
			}
	
			// Write contents
			for ( int z = 0 ; z < sizez ; z ++ )
				for ( int y = 0 ; y < sizey ; y ++ )
					for ( int x = 0 ; x < sizex ; x ++ )
					{
						float d = (float)getDataAt(x,y,z) ;
						fwrite( &d, sizeof( float ), 1, fout ) ;
					}
	
			fclose( fout ) ;
		}
	
	private:
	
		/* Sizes */
		int sizex, sizey, sizez ;
	
		/* Data array */
		double * data ;
	
		/* Spacing */
		float spcx, spcy, spcz ;
	};

}

#endif
