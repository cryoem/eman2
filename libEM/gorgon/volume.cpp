// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju (taoju@cse.wustl.edu)
// Description:   Volumetric data definition

#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include "priority_queue.h"
#include <vector>
#include "volume.h"

using namespace std;
using namespace EMAN;

// void Volume::curveSkeleton(float thr, Volume* svol)
// 	{
// 		int i, j, k ;
// 		// First, threshold the volume
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to -1/0...\n") ;
// 		#endif
// 		threshold( thr, -1, 0 ) ;
// 
// 		// Next, apply convergent erosion 
// 		// by preserving: complex nodes, curve end-points, and sheet points
// 
// 		// Next, initialize the linked queue
// 		#ifdef VERBOSE
// 		printf("Initializing queue...\n") ;
// 		#endif
// 		GridQueue2* queue2 = new GridQueue2( ) ;
// 		GridQueue2* queue3 = new GridQueue2( ) ;
// 		GridQueue2* queue4 = new GridQueue2( ) ;
// 		PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );
// 
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) >= 0 )
// 					{
// 						if ( svol->getDataAt(i,j,k) > 0 )
// 						{
// 							setDataAt( i, j, k, MAX_ERODE ) ;
// 						}
// 						else
// 						{
// 							for ( int m = 0 ; m < 6 ; m ++ )
// 							{
// 								if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
// 								{
// 									// setDataAt( i, j, k, 1 ) ;
// 									queue2->prepend( i, j, k ) ;
// 									break ;
// 								}
// 							}
// 						}
// 					}
// 				}
// 		
// 		int wid = MAX_ERODE ;
// 		#ifdef VERBOSE
// 		printf("Total %d nodes\n", queue2->getNumElements() ) ;
// 		printf("Start erosion to %d...\n", wid) ;
// 		#endif
// 
// 
// 		// Perform erosion 
// 		gridQueueEle* ele ;
// 		gridPoint* gp ;
// 		int ox, oy, oz ;
// 		int score ;
// 		Volume* scrvol = new Volume( this->sizex , this->sizey, this->sizez ) ;
// 		for ( i = 0 ; i < sizex * sizey * sizez ; i ++ )
// 		{
// 			scrvol->setDataAt( i, -1 ) ;
// 		}
// 
// #ifdef  NOISE_DIS_HELIX
// 		Volume* noisevol = new Volume( sizex, sizey, sizez ) ;
// #endif
// 
// 		for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
// 		{
// 			// At the start of each iteration, 
// 			// queue2 holds all the nodes for this layer
// 			// queue3 and queue are empty
// 
// 			int numComplex = 0, numSimple = 0 ;
// 			#ifdef VERBOSE
// 			printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
// 			#endif
// 			
// 			/*
// 			We first need to assign curwid + 1 to every node in this layer
// 			*/
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				if ( getDataAt(ox,oy,oz) == curwid )
// 				{
// 					ele = queue2->remove() ;
// 				}
// 				else
// 				{
// 					setDataAt(ox,oy,oz, curwid) ;
// 					ele = queue2->getNext() ;
// 				}
// 			}
// 			queue4->reset() ;
// 			ele = queue4->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 
// 				queue2->prepend(ox,oy,oz) ;
// 				ele = queue4->remove() ;
// 			}
// 
// 			// Now queue2 holds all the nodes for this layer
// 
// #ifdef NOISE_DIS_HELIX
// 			/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
// 			queue2->reset() ;
// 			
// 			// First run
// 			int flag = 0 ;
// 			while ( ( ele = queue2->getNext() ) != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				if ( NOISE_DIS_HELIX <= 1 )
// 				{
// 					noisevol->setDataAt( ox, oy, oz, 0 ) ;
// 				}
// 				else
// 				{
// 					flag = 0 ;
// 					for ( int m = 0 ; m < 6 ; m ++ )
// 					{
// 						int nx = ox + neighbor6[m][0] ;
// 						int ny = oy + neighbor6[m][1] ;
// 						int nz = oz + neighbor6[m][2] ;
// 						if ( getDataAt( nx, ny, nz ) == 0 )
// 						{
// 							noisevol->setDataAt( ox, oy, oz, 1 ) ;
// 							flag = 1 ;
// 							break ;
// 						}
// 					}
// 					if ( ! flag )
// 					{
// 						noisevol->setDataAt( ox, oy, oz, 0 ) ;
// 					}
// 				}
// 			}
// 
// 			int cur, visited ;
// 			for ( cur = 1 ; cur < NOISE_DIS_HELIX ; cur ++ )
// 			{
// 				queue2->reset() ;
// 				int count = 0 ;
// 				visited = 0 ;
// 
// 				while ( ( ele = queue2->getNext() ) != NULL )
// 				{
// 					ox = ele->x ;
// 					oy = ele->y ;
// 					oz = ele->z ;
// 					
// 					if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
// 					{
// 						visited ++ ;
// 						continue ;
// 					}
// 
// 					flag = 0 ;
// 					for ( int m = 0 ; m < 6 ; m ++ )
// 					{
// 						int nx = ox + neighbor6[m][0] ;
// 						int ny = oy + neighbor6[m][1] ;
// 						int nz = oz + neighbor6[m][2] ;
// 						if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
// 						{
// 							noisevol->setDataAt( ox, oy, oz, 1 ) ;
// 							visited ++ ;
// 							count ++ ;
// 							break ;
// 						}
// 					}
// 				}
// 
// 				if ( count == 0 )
// 				{
// 					break ;
// 				}
// 			}
// 			printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited ) ;
// 
// 
// #endif
// 			/* Commented out for debugging
// 			
// 			// First, 
// 			// check for complex nodes in queue2 
// 			// move them from queue2 to queue3
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Check simple 
// #ifndef NOISE_DIS_HELIX
// 				if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #else
// 				if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) ) 
// #endif
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue3->prepend( ox, oy, oz ) ;
// 					ele = queue2->remove() ;
// 					
// 					numComplex ++ ;
// 				}
// 				else
// 				{
// 					ele = queue2->getNext() ;
// 				}
// 			}
// 			*/
// 
// 			// Next,
// 			// Compute score for each node left in queue2
// 			// move them into priority queue
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Compute score
// 				score = getNumPotComplex2( ox, oy, oz ) ;
// 				scrvol->setDataAt( ox, oy, oz, score ) ;
// 
// 				// Push to queue
// 				gp = new gridPoint ;
// 				gp->x = ox ;
// 				gp->y = oy ;
// 				gp->z = oz ;
// 				// queue->add( gp, -score ) ;
// 				queue->add( gp, score ) ;
// 				
// 				ele = queue2->remove() ;
// 			}
// 
// 			// Rename queue3 to be queue2, 
// 			// Clear queue3
// 			// From now on, queue2 holds nodes of next level
// 			delete queue2 ;
// 			queue2 = queue3 ;
// 			queue3 = new GridQueue2( ) ;
// 
// 			// Next, start priority queue iteration
// 			while ( ! queue->isEmpty() )
// 			{
// 				// Retrieve the node with the highest score
// 				queue->remove( gp, score ) ;
// 				ox = gp->x ;
// 				oy = gp->y ;
// 				oz = gp->z ;
// 				delete gp ;
// //				score = -score ;
// 
// 				// Ignore the node 
// 				// if it has been processed before
// 				// or it has an updated score
// 				if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
// 				{
// 					continue ;
// 				}
// 
// 				/* Commented out for debugging
// 
// 				// Remove this simple node
// 				setDataAt( ox, oy, oz, -1 ) ;
// 				numSimple ++ ;
// 				// printf("Highest score: %d\n", score) ;
// 				*/
// 
// 				/* Added for debugging */
// 				// Check simple 
// #ifndef NOISE_DIS_HELIX
// 				// if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) ) 
// 				if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #else
// 				if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #endif
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue4->prepend( ox, oy, oz ) ;
// 					numComplex ++ ;
// 				}
// 				else
// 				{
// 					setDataAt( ox, oy, oz, -1 ) ;
// 					numSimple ++ ;
// 				}
// 				/* Adding ends */
// 
// 				// Move its neighboring unvisited node to queue2
// 				for ( int m = 0 ; m < 6 ; m ++ )
// 				{
// 					int nx = ox + neighbor6[m][0] ;
// 					int ny = oy + neighbor6[m][1] ;
// 					int nz = oz + neighbor6[m][2] ;
// 					if ( getDataAt( nx, ny, nz ) == 0 )
// 					{
// 						// setDataAt( nx, ny, nz, curwid + 1 ) ;
// 						queue2->prepend( nx, ny, nz ) ;
// 					}
// 				}
// 				
// 				/* Commented out for debugging
// 			
// 				// Find complex nodes in its 3x3 neighborhood
// 				// move them to queue2
// 				for ( i = -1 ; i < 2 ; i ++ )
// 					for ( j = -1 ; j < 2 ; j ++ )
// 						for ( k = -1 ; k < 2 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							// Check simple 
// 							if ( getDataAt( nx, ny, nz ) == curwid && 
// 								// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
// #ifndef NOISE_DIS_HELIX
// 								( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
// #else
// 								( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
// #endif
// 
// 							{
// 								// Complex, set to next layer
// 								setDataAt( nx, ny, nz, curwid + 1 ) ;
// 								queue2->prepend( nx, ny, nz ) ;
// 								numComplex ++ ;
// 							}
// 						}
// 				*/
// 
// 				// Update scores for nodes in its 5x5 neighborhood
// 				// insert them back into priority queue
// 						
// 				for ( i = -2 ; i < 3 ;i ++ )
// 					for ( j = -2 ; j < 3 ; j ++ )
// 						for ( k = -2 ; k < 3 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							if ( getDataAt( nx, ny, nz ) == curwid )
// 							{
// 								// Compute score
// 								score = getNumPotComplex2( nx, ny, nz ) ;
// 								
// 								if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
// 								{
// 									// printf("Update\n") ;
// 									scrvol->setDataAt( nx, ny, nz, score ) ;
// 									// Push to queue
// 									gp = new gridPoint ;
// 									gp->x = nx ;
// 									gp->y = ny ;
// 									gp->z = nz ;
// 									// queue->add( gp, -score ) ;
// 									queue->add( gp, score ) ;
// 								}
// 							}
// 						}
// 						
// 
// 			}
// 
// 			#ifdef VERBOSE
// 			printf("%d complex, %d simple\n", numComplex, numSimple) ;
// 			#endif
// 			
// 			if ( numSimple == 0 )
// 			{
// 					break ;
// 			}
// 		}
// 
// 		// Finally, clean up
// 		delete scrvol;
// 		delete queue;
// 		delete queue2;
// 		delete queue3;
// 		delete queue4;
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to 0/1...\n") ;
// 		#endif
// 		threshold( 0, 0, 1 ) ;		
// 	}
// 
// void Volume::curveSkeleton2D(float thr, Volume * svol)
// 	{
// 		int i, j, k ;
// 		// First, threshold the volume
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to -1/0...\n") ;
// 		#endif
// 		threshold( thr, -1, 0 ) ;
// 
// 		// Next, apply convergent erosion 
// 		// by preserving: complex nodes, curve end-points, and sheet points
// 
// 		// Next, initialize the linked queue
// 		#ifdef VERBOSE
// 		printf("Initializing queue...\n") ;
// 		#endif
// 		GridQueue2* queue2 = new GridQueue2( ) ;
// 		GridQueue2* queue3 = new GridQueue2( ) ;
// 		GridQueue2* queue4 = new GridQueue2( ) ;
// 		PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );
// 
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) >= 0 )
// 					{
// 						if ( svol->getDataAt(i,j,k) > 0 )
// 						{
// 							setDataAt( i, j, k, MAX_ERODE ) ;
// 						}
// 						else
// 						{
// 							for ( int m = 0 ; m < 4 ; m ++ )
// 							{
// 								if ( getDataAt( i + neighbor4[m][0], j + neighbor4[m][1], k ) < 0 )
// 								{
// 									// setDataAt( i, j, k, 1 ) ;
// 									queue2->prepend( i, j, k ) ;
// 									break ;
// 								}
// 							}
// 						}
// 					}
// 				}
// 		int wid = MAX_ERODE ;
// 		#ifdef VERBOSE
// 		printf("Total %d nodes\n", queue2->getNumElements() ) ;
// 		printf("Start erosion to %d...\n", wid) ;
// 		#endif
// 
// 
// 		// Perform erosion 
// 		gridQueueEle* ele ;
// 		gridPoint* gp ;
// 		int ox, oy, oz ;
// 		int score ;
// 		Volume* scrvol = new Volume( this->sizex , this->sizey, this->sizez ) ;
// 		for ( i = 0 ; i < sizex * sizey * sizez ; i ++ )
// 		{
// 			scrvol->setDataAt( i, -1 ) ;
// 		}
// 
// #ifdef  NOISE_DIS_HELIX
// 		Volume* noisevol = new Volume( sizex, sizey, sizez ) ;
// #endif
// 
// 		for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
// 		{
// 			// At the start of each iteration, 
// 			// queue2 holds all the nodes for this layer
// 			// queue3 and queue are empty
// 
// 			int numComplex = 0, numSimple = 0 ;
// 			#ifdef VERBOSE
// 			printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
// 			#endif
// 			
// 			/*
// 			We first need to assign curwid + 1 to every node in this layer
// 			*/
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				if ( getDataAt(ox,oy,oz) == curwid )
// 				{
// 					ele = queue2->remove() ;
// 				}
// 				else
// 				{
// 					setDataAt(ox,oy,oz, curwid) ;
// 					ele = queue2->getNext() ;
// 				}
// 			}
// 			queue4->reset() ;
// 			ele = queue4->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 
// 				queue2->prepend(ox,oy,oz) ;
// 				ele = queue4->remove() ;
// 			}
// 
// 			// Now queue2 holds all the nodes for this layer
// 
// #ifdef NOISE_DIS_HELIX
// 			/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
// 			queue2->reset() ;
// 			
// 			// First run
// 			int flag = 0 ;
// 			while ( ( ele = queue2->getNext() ) != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				if ( NOISE_DIS_HELIX <= 1 )
// 				{
// 					noisevol->setDataAt( ox, oy, oz, 0 ) ;
// 				}
// 				else
// 				{
// 					flag = 0 ;
// 					for ( int m = 0 ; m < 6 ; m ++ )
// 					{
// 						int nx = ox + neighbor6[m][0] ;
// 						int ny = oy + neighbor6[m][1] ;
// 						int nz = oz + neighbor6[m][2] ;
// 						if ( getDataAt( nx, ny, nz ) == 0 )
// 						{
// 							noisevol->setDataAt( ox, oy, oz, 1 ) ;
// 							flag = 1 ;
// 							break ;
// 						}
// 					}
// 					if ( ! flag )
// 					{
// 						noisevol->setDataAt( ox, oy, oz, 0 ) ;
// 					}
// 				}
// 			}
// 
// 			int cur, visited ;
// 			for ( cur = 1 ; cur < NOISE_DIS_HELIX ; cur ++ )
// 			{
// 				queue2->reset() ;
// 				int count = 0 ;
// 				visited = 0 ;
// 
// 				while ( ( ele = queue2->getNext() ) != NULL )
// 				{
// 					ox = ele->x ;
// 					oy = ele->y ;
// 					oz = ele->z ;
// 					
// 					if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
// 					{
// 						visited ++ ;
// 						continue ;
// 					}
// 
// 					flag = 0 ;
// 					for ( int m = 0 ; m < 6 ; m ++ )
// 					{
// 						int nx = ox + neighbor6[m][0] ;
// 						int ny = oy + neighbor6[m][1] ;
// 						int nz = oz + neighbor6[m][2] ;
// 						if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
// 						{
// 							noisevol->setDataAt( ox, oy, oz, 1 ) ;
// 							visited ++ ;
// 							count ++ ;
// 							break ;
// 						}
// 					}
// 				}
// 
// 				if ( count == 0 )
// 				{
// 					break ;
// 				}
// 			}
// 			printf("Maximum feature distance: %d Un-touched: %d\n", cur, queue2->getNumElements() - visited ) ;
// 
// 
// #endif
// 			/* Commented out for debugging
// 			
// 			// First, 
// 			// check for complex nodes in queue2 
// 			// move them from queue2 to queue3
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Check simple 
// #ifndef NOISE_DIS_HELIX
// 				if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #else
// 				if ( isHelixEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) ) 
// #endif
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue3->prepend( ox, oy, oz ) ;
// 					ele = queue2->remove() ;
// 					
// 					numComplex ++ ;
// 				}
// 				else
// 				{
// 					ele = queue2->getNext() ;
// 				}
// 			}
// 			*/
// 
// 			// Next,
// 			// Compute score for each node left in queue2
// 			// move them into priority queue
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Compute score
// 				score = getNumPotComplex2( ox, oy, oz ) ;
// 				//score = getNumNeighbor6( ox, oy, oz ) ;
// 				scrvol->setDataAt( ox, oy, oz, score ) ;
// 
// 				// Push to queue
// 				gp = new gridPoint ;
// 				gp->x = ox ;
// 				gp->y = oy ;
// 				gp->z = oz ;
// 				// queue->add( gp, -score ) ;
// 				queue->add( gp, score ) ;
// 				
// 				ele = queue2->remove() ;
// 			}
// 
// 			// Rename queue3 to be queue2, 
// 			// Clear queue3
// 			// From now on, queue2 holds nodes of next level
// 			delete queue2 ;
// 			queue2 = queue3 ;
// 			queue3 = new GridQueue2( ) ;
// 
// 			// Next, start priority queue iteration
// 			while ( ! queue->isEmpty() )
// 			{
// 				// Retrieve the node with the highest score
// 				queue->remove( gp, score ) ;
// 				ox = gp->x ;
// 				oy = gp->y ;
// 				oz = gp->z ;
// 				delete gp ;
// //				score = -score ;
// 
// 				// Ignore the node 
// 				// if it has been processed before
// 				// or it has an updated score
// 				if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
// 				{
// 					continue ;
// 				}
// 
// 				/* Commented out for debugging
// 
// 				// Remove this simple node
// 				setDataAt( ox, oy, oz, -1 ) ;
// 				numSimple ++ ;
// 				// printf("Highest score: %d\n", score) ;
// 				*/
// 
// 				/* Added for debugging */
// 				// Check simple 
// #ifndef NOISE_DIS_HELIX
// 				// if ( hasIsolatedEdge( ox, oy, oz ) && ! isNoiseHelixEnd( ox, oy, oz ) ) 
// 				if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #else
// 				if ( isHelixEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #endif
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue4->prepend( ox, oy, oz ) ;
// 					numComplex ++ ;
// 				}
// 				else
// 				{
// 					setDataAt( ox, oy, oz, -1 ) ;
// 					numSimple ++ ;
// 				}
// 				/* Adding ends */
// 
// 				// Move its neighboring unvisited node to queue2
// 				for ( int m = 0 ; m < 4 ; m ++ )
// 				{
// 					int nx = ox + neighbor4[m][0] ;
// 					int ny = oy + neighbor4[m][1] ;
// 					int nz = oz ;
// 					if ( getDataAt( nx, ny, nz ) == 0 )
// 					{
// 						// setDataAt( nx, ny, nz, curwid + 1 ) ;
// 						queue2->prepend( nx, ny, nz ) ;
// 					}
// 				}
// 				
// 				/* Commented out for debugging
// 			
// 				// Find complex nodes in its 3x3 neighborhood
// 				// move them to queue2
// 				for ( i = -1 ; i < 2 ; i ++ )
// 					for ( j = -1 ; j < 2 ; j ++ )
// 						for ( k = -1 ; k < 2 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							// Check simple 
// 							if ( getDataAt( nx, ny, nz ) == curwid && 
// 								// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
// #ifndef NOISE_DIS_HELIX
// 								( isHelixEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
// #else
// 								( isHelixEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
// #endif
// 
// 							{
// 								// Complex, set to next layer
// 								setDataAt( nx, ny, nz, curwid + 1 ) ;
// 								queue2->prepend( nx, ny, nz ) ;
// 								numComplex ++ ;
// 							}
// 						}
// 				*/
// 
// 				// Update scores for nodes in its 5x5 neighborhood
// 				// insert them back into priority queue
// 						
// 				for ( i = -2 ; i < 3 ;i ++ )
// 					for ( j = -2 ; j < 3 ; j ++ )
// 						for ( k = -2 ; k < 3 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							if ( getDataAt( nx, ny, nz ) == curwid )
// 							{
// 								// Compute score
// 								score = getNumPotComplex2( nx, ny, nz ) ;
// 								//score = getNumNeighbor6( nx, ny, nz ) ;
// 								
// 								if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
// 								{
// 									// printf("Update\n") ;
// 									scrvol->setDataAt( nx, ny, nz, score ) ;
// 									// Push to queue
// 									gp = new gridPoint ;
// 									gp->x = nx ;
// 									gp->y = ny ;
// 									gp->z = nz ;
// 									// queue->add( gp, -score ) ;
// 									queue->add( gp, score ) ;
// 								}
// 							}
// 						}
// 						
// 
// 			}
// 
// 			#ifdef VERBOSE
// 			printf("%d complex, %d simple\n", numComplex, numSimple) ;
// 			#endif
// 			
// 			if ( numSimple == 0 )
// 			{
// 					break ;
// 			}
// 		}
// 
// 		// Finally, clean up
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to 0/1...\n") ;
// 		#endif
// 		threshold( 0, 0, 1 ) ;		
// 		delete scrvol;
// 		delete queue;
// 		delete queue2;
// 		delete queue3;
// 		delete queue4;
// 	}
// 
// void Volume::erodeHelix(int disthr)
// 	{
// 		int i, j, k ;
// 		// First, threshold the volume
// 		//printf("Thresholding the volume to -1/0...\n") ;
// 		threshold( 0.1f, -1, 0 ) ;
// 		
// 		/* Debug: remove faces */
// 		//Volume* facevol = markFaceEdge() ;
// 		/* End debugging */
// 
// 		// Next, initialize the linked queue
// 		//printf("Initializing queue...\n") ;
// 		Volume * fvol = new Volume( sizex, sizey, sizez ) ;
// 		GridQueue2* queue2 = new GridQueue2( ) ;
// 		GridQueue2* queue3 = new GridQueue2( ) ;
// 		GridQueue2** queues = new GridQueue2* [ 10000 ] ;
// 
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) >= 0 )
// 					{
// 						if ( ! hasCompleteHelix( i, j, k ) )
// 						// if ( ! hasCompleteHelix( i, j, k, facevol ) )
// 						{
// 							queue2->prepend( i, j, k ) ;
// 							fvol->setDataAt( i, j, k, -1 ) ;
// 						}
// 					}
// 				}
// 		//printf("Total %d nodes\n", queue2->getNumElements() ) ;
// 
// 		// Start erosion
// 		gridQueueEle* ele ;
// 		int dis = -1 ;
// 		while ( queue2->getNumElements() > 0 )
// 		{
// 			// First, set distance
// 			dis -- ;
// 			queues[ - dis ] = new GridQueue2( ) ;
// 			//printf("Distance transform to %d...", dis) ;
// 			queue2->reset() ;
// 			while( ( ele = queue2->getNext() ) != NULL )
// 			{
// 				setDataAt( ele->x, ele->y, ele->z, dis ) ;
// 				queues[ -dis ]->prepend( ele->x, ele->y, ele->z ) ;
// 			}
// 			//printf("%d nodes\n", queues[-dis]->getNumElements()) ;
// 
// 			// Next, find next layer
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				for ( int m = 0 ; m < 6 ; m ++ )
// 				{
// 					int nx = ele->x + neighbor6[m][0] ;
// 					int ny = ele->y + neighbor6[m][1] ;
// 					int nz = ele->z + neighbor6[m][2] ;
// 					if ( getDataAt( nx, ny, nz ) == 0 )
// 					{
// 						fvol->setDataAt( nx, ny, nz, dis ) ;
// 
// 						if ( ! hasCompleteHelix( nx, ny, nz ) )
// 						// if ( ! hasCompleteHelix( nx, ny, nz, facevol ) )
// 						{
// 							setDataAt( nx, ny, nz, 1 ) ;
// 							queue3->prepend( nx, ny, nz ) ;
// 						}
// 					}
// 				}
// 				
// 				ele = queue2->remove() ;
// 			}
// 
// 			// Next, swap queues
// 			GridQueue2* temp = queue2 ;
// 			queue2 = queue3 ;
// 			queue3 = temp ;
// 		}
// 
// 		// Deal with closed rings
// 		dis -- ;
// 		queues[ - dis ] = new GridQueue2( ) ;
// 		#ifdef VERBOSE
// 		printf("Detecting closed rings %d...", dis) ;
// 		#endif
// 		int ftot = 0 ;
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) == 0 )
// 					{
// 						setDataAt( i, j, k, dis ) ;
// 						queues[ -dis ]->prepend( i, j, k ) ;
// /*
// 						int fval = (int) fvol->getDataAt( i, j, k ) ;
// 						if ( fval == 0)
// 						{
// 							// queues[ -dis ]->prepend( i, j, k ) ;
// 						}
// 						else
// 						{
// 							setDataAt( i, j, k, fval - 1 ) ;
// 							// queues[ -fval + 1 ]->prepend( i, j, k ) ;
// 						}
// */
// 						ftot ++ ;
// 					}
// 				}
// 		#ifdef VERBOSE
// 		printf("%d nodes\n", ftot) ;
// 		#endif
// 
// 
// 		// return ;
// 
// 		/* Find local minimum: to help determine erosion level
// 		int cts[ 64 ] ;
// 		for ( i = 0 ; i <= - dis ; i ++ )
// 		{
// 			cts[ i ] = 0 ;
// 		}
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					double val = getDataAt( i, j, k ) ;
// 					if ( val < -1 )
// 					{
// 						int min = 1 ;
// 						for ( int m = 0 ; m < 6 ; m ++ )
// 						{
// 							int nx = i + neighbor6[m][0] ;
// 							int ny = j + neighbor6[m][1] ;
// 							int nz = k + neighbor6[m][2] ;
// 							if ( getDataAt( nx, ny, nz ) < val )
// 							{
// 								min = 0 ;
// 								break ;
// 							}
// 						}
// 
// 						if ( min )
// 						{
// 							cts[ (int) - val ] ++ ;
// 						}
// 					}
// 				}
// 
// 		for ( i = 2 ; i <= - dis ; i ++ )
// 		{
// 			printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
// 		}
// 		*/
// 
// 		// Dilate back
// 		// Starting from nodes with distance - 2 - disthr
// 
// 		if ( disthr + 2 > - dis )
// 		{
// 			disthr = - dis - 2 ;
// 		}
// 		int d; 
// 		
// 		for ( d = - dis ; d > disthr + 1 ; d -- )
// 		{
// 			queues[ d ]->reset() ;
// 			while ( (ele = queues[ d ]->getNext() ) != NULL )
// 			{
// 				setDataAt( ele->x, ele->y, ele->z, d ) ;
// 			}
// 		}
// 		
// 
// 		for ( d = disthr + 1 ; d >= 2 ; d -- )
// 		{
// 			//delete queue3 ;
// 			//queue3 = new GridQueue2( ) ;
// 			queues[ d ]->reset() ;
// 			ele = queues[ d ]->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				int dilatable = 0 ;
// 				for ( int m = 0 ; m < 6 ; m ++ )
// 						{
// 							int nx = ele->x + neighbor6[m][0] ;
// 							int ny = ele->y + neighbor6[m][1] ;
// 							int nz = ele->z + neighbor6[m][2] ;
// 							if ( getDataAt( nx, ny, nz ) == d + 1 )
// 							{
// 								dilatable = 1 ;
// 								break ;
// 							}
// 						}
// 
// 				
// 				if ( ! dilatable )
// 				{
// 					/*
// 					setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
// 					*/
// 					
// 					setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
// 					if ( d > 2 )
// 					{
// 						queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
// 					}
// 					ele = queues[ d ]->remove() ;
// 				}
// 				else
// 				{
// 					setDataAt( ele->x, ele->y, ele->z, d ) ;
// 					ele = queues[ d ]->getNext() ;
// 				}
// 				
// 			}
// 
// 			/* Debug: extract minimal set */
// 			while ( 1 )
// 			{
// 				int num = 0 ;
// 				queues[ d ]->reset() ;
// 				ele = queues[ d ]->getNext() ;
// 				while ( ele != NULL )
// 				{
// 					int critical = 0 ;
// 					setDataAt( ele->x, ele->y, ele->z, -1 ) ;
// 					
// 					for ( i = -1 ; i < 2 ; i ++ )
// 					{
// 						for ( j = -1 ; j < 2 ; j ++ )
// 						{
// 							for ( k = -1 ; k < 2 ; k ++ )
// 							{
// 								if ( i != 0 && j != 0 && k != 0 )
// 								{
// 									continue ;
// 								}
// 								int nx = ele->x + i ;
// 								int ny = ele->y + j ;
// 								int nz = ele->z + k ;
// 								if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteHelix( nx,ny,nz ) ) //, facevol ) )
// 								{
// 									critical = 1 ;
// 									break ;
// 								}
// 							}
// 							if ( critical )
// 							{
// 								break ;
// 							}
// 						}
// 						if ( critical )
// 						{
// 							break ;
// 						}
// 					}
// 					
// 					if ( critical )
// 					{
// 						setDataAt( ele->x, ele->y, ele->z, d ) ;
// 						ele = queues[ d ]->getNext() ;
// 					}
// 					else
// 					{
// 						setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
// 						if ( d > 2 )
// 						{
// 							queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
// 						}
// 						ele = queues[ d ]->remove() ;
// 						num ++ ;
// 					}
// 					
// 				}
// 
// 				#ifdef VERBOSE
// 				printf("Non-minimal: %d\n", num) ;
// 				#endif
// 
// 				if ( num == 0 )
// 				{
// 					break ;
// 				}
// 			}
// 
// 
// 			/* End of debugging */
// 	
// 			/*
// 			queue3->reset() ;
// 			ele = queue3->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
// 				ele = queue3->remove() ;
// 			}
// 			*/
// 		}
// 
// 		// Finally, threshold the volume
// 		#ifdef VERBOSE
// 		//printf("Thresholding the volume to 0/1...\n") ;
// 		#endif
// 		//threshold( -1, 1, 0, 0 ) ;
// 		threshold( 0, 0, 1 ) ;
// 		delete fvol ;
// 		delete queue2;
// 		delete queue3;
// 		for ( d = -dis ; d >= 2 ; d -- ) {
// 			delete queues[d];
// 		}
// 		delete [] queues;
// 
// 	}
// 
// int Volume::erodeSheet(int disthr)
// 	{
// 		int i, j, k ;
// 		// First, threshold the volume
// 		//printf("Thresholding the volume to -1/0...\n") ;
// 		threshold( 0.1f, -1, 0 ) ;
// 
// 		/* Debug: remove cells */
// 		Volume* facevol = markCellFace() ;
// 		/* End debugging */
// 		
// 		// Next, initialize the linked queue
// 		//printf("Initializing queue...\n") ;
// 		Volume * fvol = new Volume( sizex, sizey, sizez ) ;
// 		GridQueue2* queue2 = new GridQueue2( ) ;
// 		GridQueue2* queue3 = new GridQueue2( ) ;
// 		GridQueue2** queues = new GridQueue2* [ 10000 ] ;
// 
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) >= 0 )
// 					{
// 						if ( ! hasCompleteSheet( i, j, k ) )
// 						//if ( ! hasCompleteSheet( i, j, k, facevol ) )
// 						{
// 							queue2->prepend( i, j, k ) ;
// 							fvol->setDataAt( i, j, k, -1 ) ;
// 						}
// 					}
// 				}
// 		#ifdef VERBOSE
// 		printf("Total %d nodes\n", queue2->getNumElements() ) ;
// 		#endif
// 
// 		// Start erosion
// 		gridQueueEle* ele ;
// 		int dis = -1 ;
// 		while ( queue2->getNumElements() > 0 )
// 		{
// 			// First, set distance
// 			dis -- ;
// 			queues[ - dis ] = new GridQueue2( ) ;
// 			//printf("Distance transform to %d...", dis) ;
// 			queue2->reset() ;
// 			while( ( ele = queue2->getNext() ) != NULL )
// 			{
// 				setDataAt( ele->x, ele->y, ele->z, dis ) ;
// 				queues[ -dis ]->prepend( ele->x, ele->y, ele->z ) ;
// 			}
// 			//printf("%d nodes\n", queues[-dis]->getNumElements()) ;
// 
// 			// Next, find next layer
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				// for ( int m = 0 ; m < 6 ; m ++ )
// 				for ( int mx = -1 ; mx < 2 ; mx ++ )
// 					for ( int my = -1 ; my < 2 ; my ++ )
// 						for ( int mz = -1 ; mz < 2 ; mz ++ )
// 						{
// 							if ( mx != 0 && my != 0 && mz != 0 )
// 							{
// 								continue ;
// 							}
// 							//int nx = ele->x + neighbor6[m][0] ;
// 							//int ny = ele->y + neighbor6[m][1] ;
// 							//int nz = ele->z + neighbor6[m][2] ;
// 							int nx = ele->x + mx ;
// 							int ny = ele->y + my ;
// 							int nz = ele->z + mz ;
// 
// 							if ( getDataAt( nx, ny, nz ) == 0 )
// 							{
// 								fvol->setDataAt( nx, ny, nz, dis ) ;
// 
// 								if  ( ! hasCompleteSheet( nx, ny, nz ) )
// 								// if  ( ! hasCompleteSheet( nx, ny, nz, facevol ) )
// 								{
// 									setDataAt( nx, ny, nz, 1 ) ;
// 									queue3->prepend( nx, ny, nz ) ;
// 								}
// 							}
// 						}
// 				
// 				ele = queue2->remove() ;
// 			}
// 
// 			// Next, swap queues
// 			GridQueue2* temp = queue2 ;
// 			queue2 = queue3 ;
// 			queue3 = temp ;
// 		}
// 
// 		/* Deal with closed rings */
// 		
// 		dis -- ;
// 		queues[ - dis ] = new GridQueue2( ) ;
// 		#ifdef VERBOSE
// 		printf("Detecting closed rings %d...", dis) ;
// 		#endif
// 		int ftot = 0 ;
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) == 0 )
// 					{
// 						/*
// 						int fval = (int) fvol->getDataAt( i, j, k ) ;
// 						if ( fval == 0)
// 						{
// 							setDataAt( i, j, k, dis - 2 ) ;
// 							// queues[ -dis ]->prepend( i, j, k ) ;
// 						}
// 						else
// 						{
// 							setDataAt( i, j, k, fval - 1 ) ;
// 							queues[ -fval + 1 ]->prepend( i, j, k ) ;
// 						}
// 						*/
// 						setDataAt( i, j, k, dis ) ;
// 						queues[ -dis ]->prepend( i, j, k ) ;
// 
// 						ftot ++ ;
// 					}
// 				}
// 		#ifdef VERBOSE
// 		printf("%d nodes\n", ftot) ;
// 		#endif
// 		
// 
// 		/* Find local minimum: to help determine erosion level 
// 		int cts[ 64 ] ;
// 		for ( i = 0 ; i <= - dis ; i ++ )
// 		{
// 			cts[ i ] = 0 ;
// 		}
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					double val = getDataAt( i, j, k ) ;
// 					if ( val < -1 )
// 					{
// 						int min = 1 ;
// 						for ( int m = 0 ; m < 6 ; m ++ )
// 						{
// 							int nx = i + neighbor6[m][0] ;
// 							int ny = j + neighbor6[m][1] ;
// 							int nz = k + neighbor6[m][2] ;
// 							if ( getDataAt( nx, ny, nz ) < val )
// 							{
// 								min = 0 ;
// 								break ;
// 							}
// 						}
// 
// 						if ( min )
// 						{
// 							cts[ (int) - val ] ++ ;
// 						}
// 					}
// 				}
// 
// 		for ( i = 2 ; i <= - dis ; i ++ )
// 		{
// 			printf("Local minima: %d with %d nodes.\n", -i, cts[ i ] ) ;
// 		}
// 		*/
// 
// 		// return ;
// 
// 		// Dilate back
// 		// Starting from nodes with distance - 2 - disthr
// 		int d ;
// 		if ( disthr + 2 > - dis )
// 		{
// 			disthr = - dis - 2 ;
// 
// 		}
// 		for ( d = - dis ; d > disthr + 1 ; d -- )
// 		{
// 			queues[ d ]->reset() ;
// 			while ( (ele = queues[ d ]->getNext() ) != NULL )
// 			{
// 				setDataAt( ele->x, ele->y, ele->z, d ) ;
// 			}
// 		}
// 
// 		for (d = disthr + 1 ; d >= 2 ; d -- )
// 		{
// 			
// 			//delete queue3 ;
// 			//queue3 = new GridQueue2( ) ;
// 			queues[ d ]->reset() ;
// 			ele = queues[ d ]->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				int dilatable = 0 ;
// 				// for ( int m = 0 ; m < 6 ; m ++ )
// 				/*
// 				for ( int mx = -1 ; mx < 2 ; mx ++ )
// 				for ( int my = -1 ; my < 2 ; my ++ )
// 				for ( int mz = -1 ; mz < 2 ; mz ++ )
// 				{
// 				if ( mx == 0 || my == 0 || mz == 0 )
// 				{
// 				int nx = ele->x + mx ; // neighbor6[m][0] ;
// 				int ny = ele->y + my ; // neighbor6[m][1] ;
// 				int nz = ele->z + mz ; // neighbor6[m][2] ;
// 				if ( getDataAt( nx, ny, nz ) == - d - 1 )
// 				{
// 				dilatable = 1 ;
// 				break ;
// 				}
// 				}
// 				}
// 				*/
// 				for ( i = 0 ; i < 12 ; i ++ )
// 				{	
// 					int flag = 1, flag2 = 0 ;
// 					for ( j = 0 ; j < 4 ; j ++ )
// 					{
// 						int nx = ele->x + sheetNeighbor[i][j][0] ;
// 						int ny = ele->y + sheetNeighbor[i][j][1] ;
// 						int nz = ele->z + sheetNeighbor[i][j][2] ;
// 						
// 						double val = getDataAt( nx, ny, nz ) ;
// 						
// 						if ( val > - d && val < 0)
// 						{
// 							flag = 0 ;
// 							break ;
// 						}
// 						if ( val == d + 1 )
// 						{
// 							flag2 ++ ;
// 						}
// 					}
// 					
// 					if ( flag && flag2 )
// 					{
// 						dilatable = 1 ;
// 						break ;
// 					}
// 				}
// 				
// 				if ( ! dilatable )
// 				{
// 					// setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
// 					// queue3->prepend( ele->x, ele->y, ele->z ) ;
// 					
// 					setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
// 					if ( d > 2 )
// 					{
// 						queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
// 					}
// 					ele = queues[ d ]->remove() ;
// 				}
// 				else
// 				{
// 					setDataAt( ele->x, ele->y, ele->z, d ) ;
// 					ele = queues[ d ]->getNext() ;
// 				}
// 			}
// 
// 			/* Debug: extract minimal set */
// 			while ( 1 )
// 			{
// 				int num = 0 ;
// 				queues[ d ]->reset() ;
// 				ele = queues[ d ]->getNext() ;
// 				while ( ele != NULL )
// 				{
// 					int critical = 0 ;
// 					setDataAt( ele->x, ele->y, ele->z, -1 ) ;
// 					
// 					for ( i = -1 ; i < 2 ; i ++ )
// 					{
// 						for ( j = -1 ; j < 2 ; j ++ )
// 						{
// 							for ( k = -1 ; k < 2 ; k ++ )
// 							{
// 								if ( i != 0 && j != 0 && k != 0 )
// 								{
// 									continue ;
// 								}
// 								int nx = ele->x + i ;
// 								int ny = ele->y + j ;
// 								int nz = ele->z + k ;
// 								// if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteSheet( nx,ny,nz, facevol ) )
// 								if ( getDataAt(nx,ny,nz) == d + 1 && !hasCompleteSheet( nx,ny,nz ) )
// 								{
// 									critical = 1 ;
// 									break ;
// 								}
// 							}
// 							if ( critical )
// 							{
// 								break ;
// 							}
// 						}
// 						if ( critical )
// 						{
// 							break ;
// 						}
// 					}
// 					
// 					if ( critical )
// 					{
// 						setDataAt( ele->x, ele->y, ele->z, d ) ;
// 						ele = queues[ d ]->getNext() ;
// 					}
// 					else
// 					{
// 						setDataAt( ele->x, ele->y, ele->z, - d + 1 ) ;
// 						if ( d > 2 )
// 						{
// 							queues[ d - 1 ]->prepend( ele->x, ele->y, ele->z ) ;
// 						}
// 						ele = queues[ d ]->remove() ;
// 						num ++ ;
// 					}
// 					
// 				}
// 				#ifdef VERBOSE
// 				printf("Non-minimal: %d\n", num) ;
// 				#endif
// 
// 				if ( num == 0 )
// 				{
// 					break ;
// 				}
// 			}
// 
// 
// 			/* End of debugging */
// 			
// 			/*
// 			queue3->reset() ;
// 			ele = queue3->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				setDataAt( ele->x, ele->y, ele->z, - 1 ) ;
// 				ele = queue3->remove() ;
// 			}
// 			*/
// 		}
// 		
// 
// 		// Finally, threshold the volume
// 		#ifdef VERBOSE
// 		//printf("Thresholding the volume to 0/1...\n") ;
// 		#endif
// 		//threshold( -1, 1, 0, 0 ) ;
// 		threshold( 0, 0, 1 ) ;
// 
// 		delete facevol ;
// 		delete fvol ;
// 		delete queue2;
// 		delete queue3;
// 		for (d = -dis ; d >= 2 ; d -- ) {
// 			delete queues[d];
// 		}
// 		delete [] queues;
// 	
// 		return - dis - 1 ;
// 	}
// 
// void Volume::pad(int padBy, double padValue)
//  {
// 		int newSizeX = sizex + 2*padBy;
// 		int newSizeY = sizey + 2*padBy;
// 		int newSizeZ = sizez + 2*padBy;
// 
// 		float * newData = new float[newSizeX * newSizeY * newSizeZ];
// 		double value;
// 
// 
// 		for(int x = 0; x < newSizeX; x++) {
// 			for(int y = 0; y < newSizeY; y++) {
// 				for(int z = 0; z < newSizeZ; z++) {
// 					if ((x < padBy) || (y < padBy) || (z < padBy) || (x >= padBy + sizex) || (y >= padBy + sizey) || (z >= padBy + sizez)) {
// 						value = padValue;
// 					} else {
// 						value = getDataAt(x-padBy, y-padBy, z-padBy);
// 					}
// 
// 					newData[x * newSizeY * newSizeZ + y * newSizeZ + z] = (float)value;
// 				}
// 			}
// 		}
// 		delete [] data;
// 		data = newData;
// 		sizex = newSizeX;
// 		sizey = newSizeY;
// 		sizez = newSizeZ;
// }
//
// void Volume::skeleton(float thr, Volume* svol, Volume* hvol)
// 	{
// 		int i, j, k ;
// 		// First, threshold the volume
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to -1/0...\n") ;
// 		#endif
// 		threshold( thr, -1, 0 ) ;
// 
// 		// Next, apply convergent erosion 
// 		// by preserving: complex nodes, curve end-points, and sheet points
// 
// 		// Next, initialize the linked queue
// 		#ifdef VERBOSE
// 		printf("Initializing queue...\n") ;
// 		#endif
// 		GridQueue2* queue2 = new GridQueue2( ) ;
// 		GridQueue2* queue3 = new GridQueue2( ) ;
// 		PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );
// 
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) >= 0 )
// 					{
// 						if ( svol->getDataAt(i,j,k) > 0 || hvol->getDataAt(i,j,k) > 0 )
// 						{
// 							setDataAt( i, j, k, MAX_ERODE ) ;
// 						}
// 						else
// 						{
// 							for ( int m = 0 ; m < 6 ; m ++ )
// 							{
// 								if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
// 								{
// 									setDataAt( i, j, k, 1 ) ;
// 									queue2->prepend( i, j, k ) ;
// 									break ;
// 								}
// 							}
// 						}
// 					}
// 				}
// 		int wid = MAX_ERODE ;
// 		#ifdef VERBOSE
// 		printf("Total %d nodes\n", queue2->getNumElements() ) ;
// 
// 
// 		// Perform erosion 
// 		printf("Start erosion to %d...\n", wid) ;
// 		#endif
// 		gridQueueEle* ele ;
// 		gridPoint* gp ;
// 		int ox, oy, oz ;
// 		int score ;
// 		Volume* scrvol = new Volume( this->sizex , this->sizey, this->sizez ) ;
// 		for ( i = 0 ; i < sizex * sizey * sizez ; i ++ )
// 		{
// 			scrvol->setDataAt( i, -1 ) ;
// 		}
// 
// 
// 		for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
// 		{
// 			// At the start of each iteration, 
// 			// queue2 holds all the nodes for this layer
// 			// queue3 and queue are empty
// 
// 			int numComplex = 0, numSimple = 0 ;
// 			#ifdef VERBOSE
// 			printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
// 			#endif
// 			
// 
// 			// Next,
// 			// Compute score for each node left in queue2
// 			// move them into priority queue
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Compute score
// 				score = getNumPotComplex2( ox, oy, oz ) ;
// 				scrvol->setDataAt( ox, oy, oz, score ) ;
// 
// 				// Push to queue
// 				gp = new gridPoint ;
// 				gp->x = ox ;
// 				gp->y = oy ;
// 				gp->z = oz ;
// 				// queue->add( gp, -score ) ;
// 				queue->add( gp, score ) ;
// 								
// 				ele = queue2->remove() ;
// 			}
// 
// 			// Rename queue3 to be queue2, 
// 			// Clear queue3
// 			// From now on, queue2 holds nodes of next level
// 			delete queue2 ;
// 			queue2 = queue3 ;
// 			queue3 = new GridQueue2( ) ;
// 
// 			// Next, start priority queue iteration
// 			while ( ! queue->isEmpty() )
// 			{
// 				// Retrieve the node with the highest score
// 				queue->remove( gp, score ) ;
// 				ox = gp->x ;
// 				oy = gp->y ;
// 				oz = gp->z ;
// 				delete gp ;
// //				score = -score ;
// 
// 				// Ignore the node 
// 				// if it has been processed before
// 				// or it has an updated score
// 				if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
// 				{
// 					continue ;
// 				}
// 
// 				/* Added for debugging */
// 				// Check simple 
// 				if ( ! isSimple( ox, oy, oz ) ) 
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue2->prepend( ox, oy, oz ) ;
// 					numComplex ++ ;
// 				}
// 				else
// 				{
// 					setDataAt( ox, oy, oz, -1 ) ;
// 					numSimple ++ ;
// 				}
// 				/* Adding ends */
// 
// 				// Move its neighboring unvisited node to queue2
// 				for ( int m = 0 ; m < 6 ; m ++ )
// 				{
// 					int nx = ox + neighbor6[m][0] ;
// 					int ny = oy + neighbor6[m][1] ;
// 					int nz = oz + neighbor6[m][2] ;
// 					if ( getDataAt( nx, ny, nz ) == 0 )
// 					{
// 						setDataAt( nx, ny, nz, curwid + 1 ) ;
// 						queue2->prepend( nx, ny, nz ) ;
// 					}
// 				}
// 				
// 				// Update scores for nodes in its 5x5 neighborhood
// 				// insert them back into priority queue
// 						
// 				for ( i = -2 ; i < 3 ;i ++ )
// 					for ( j = -2 ; j < 3 ; j ++ )
// 						for ( k = -2 ; k < 3 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							if ( getDataAt( nx, ny, nz ) == curwid )
// 							{
// 								// Compute score
// 								score = getNumPotComplex2( nx, ny, nz ) ;
// 								
// 								if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
// 								{
// 									// printf("Update\n") ;
// 									scrvol->setDataAt( nx, ny, nz, score ) ;
// 									// Push to queue
// 									gp = new gridPoint ;
// 									gp->x = nx ;
// 									gp->y = ny ;
// 									gp->z = nz ;
// 									// queue->add( gp, -score ) ;
// 									queue->add( gp, score ) ;
// 								}
// 							}
// 						}
// 						
// 
// 			}
// 
// 			#ifdef VERBOSE
// 			printf("%d complex, %d simple\n", numComplex, numSimple) ;
// 			#endif
// 			
// 			if ( numSimple == 0 )
// 			{
// 				delete queue2 ;
// 				break ;
// 			}
// 		}
// 
// 		// Finally, clean up
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to 0/1...\n") ;
// 		#endif
// 		threshold( 0, 0, 1 ) ;	
// 		delete scrvol;
// 		delete queue;
// 		delete queue3;
// 	}
// 
// void Volume::surfaceSkeletonPres(float thr, Volume* svol)
// 	{
// 		int i, j, k ;
// 		// First, threshold the volume
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to -MAX_ERODE/0...\n") ;
// 		#endif
// 		threshold( thr, -MAX_ERODE, 0 ) ;
// 
// 		// Next, initialize the linked queue
// 		#ifdef VERBOSE
// 		printf("Initializing queue...\n") ;
// 		#endif
// 		GridQueue2* queue2 = new GridQueue2( ) ;
// 		GridQueue2* queue3 = new GridQueue2( ) ;
// 		GridQueue2* queue4 = new GridQueue2( ) ;
// 
// 		PriorityQueue <gridPoint,int> * queue = new PriorityQueue <gridPoint,int> ( MAX_QUEUELEN );
// 
// 		for ( i = 0 ; i < sizex ; i ++ )
// 			for ( j = 0 ; j < sizey ; j ++ )
// 				for ( k = 0 ; k < sizez ; k ++ )
// 				{
// 					if ( getDataAt( i, j, k ) >= 0 ) {
// 						if(preserve->getDataAt(i, j, k) > 0) {
// 							setDataAt(i, j, k, MAX_ERODE);
// 						} else {
// 							for ( int m = 0 ; m < 6 ; m ++ )
// 							{
// 								if ( getDataAt( i + neighbor6[m][0], j + neighbor6[m][1], k + neighbor6[m][2] ) < 0 )
// 								{
// 									// setDataAt( i, j, k, 1 ) ;
// 									queue2->prepend( i, j, k ) ;
// 									break ;
// 								}
// 							}
// 						}
// 					}
// 				}
// 		int wid = MAX_ERODE ;
// 		#ifdef VERBOSE
// 		printf("Total %d nodes\n", queue2->getNumElements() ) ;
// 		printf("Start erosion to %d...\n", wid) ;
// 		#endif
// 
// 
// 		// Perform erosion 
// 		gridQueueEle* ele ;
// 		gridPoint* gp ;
// 		int ox, oy, oz ;
// 		int score ;
// 		Volume* scrvol = new Volume( this->sizex , this->sizey, this->sizez ) ;
// 		for ( i = 0 ; i < sizex * sizey * sizez ; i ++ )
// 		{
// 			scrvol->setDataAt( i, -1 ) ;
// 		}
// 
// #ifdef  NOISE_DIS_SHEET
// 		Volume* noisevol = new Volume( sizex, sizey, sizez ) ;
// #endif
// 
// 		for ( int curwid = 1 ; curwid <= wid ; curwid ++ )
// 		{
// 			// At the start of each iteration, 
// 			// queue2 and queue4 holds all the nodes for this layer
// 			// queue3 and queue are empty
// 
// 			int numComplex = 0, numSimple = 0 ;
// 			#ifdef VERBOSE
// 			printf("Processing %d nodes in layer %d\n", queue2->getNumElements(), curwid) ;
// 			#endif			
// 
// 			/*
// 			We first need to assign curwid + 1 to every node in this layer
// 			*/
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				if ( getDataAt(ox,oy,oz) == curwid )
// 				{
// 					ele = queue2->remove() ;
// 				}
// 				else
// 				{
// 					setDataAt(ox,oy,oz, curwid) ;
// 					ele = queue2->getNext() ;
// 				}
// 			}
// 			queue4->reset() ;
// 			ele = queue4->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 
// 				queue2->prepend(ox,oy,oz) ;
// 				ele = queue4->remove() ;
// 			}
// 
// 			// Now queue2 holds all the nodes for this layer
// 
// #ifdef NOISE_DIS_SHEET
// 			/* Extra step: classify nodes in queue2 into noise and non-noise nodes */
// 			queue2->reset() ;
// 			
// 			// First run
// 			int flag = 0 ;
// 			while ( ( ele = queue2->getNext() ) != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				if ( NOISE_DIS_SHEET <= 1 )
// 				{
// 					noisevol->setDataAt( ox, oy, oz, 0 ) ;
// 				}
// 				else
// 				{
// 					flag = 0 ;
// 					for ( int m = 0 ; m < 6 ; m ++ )
// 					{
// 						int nx = ox + neighbor6[m][0] ;
// 						int ny = oy + neighbor6[m][1] ;
// 						int nz = oz + neighbor6[m][2] ;
// 						if ( getDataAt( nx, ny, nz ) == 0 )
// 						{
// 							noisevol->setDataAt( ox, oy, oz, 1 ) ;
// 							flag = 1 ;
// 							break ;
// 						}
// 					}
// 					if ( ! flag )
// 					{
// 						noisevol->setDataAt( ox, oy, oz, 0 ) ;
// 					}
// 				}
// 			}
// 
// 			for ( int cur = 1 ; cur < NOISE_DIS_SHEET ; cur ++ )
// 			{
// 				queue2->reset() ;
// 				int count = 0 ;
// 
// 				while ( ( ele = queue2->getNext() ) != NULL )
// 				{
// 					ox = ele->x ;
// 					oy = ele->y ;
// 					oz = ele->z ;
// 					
// 					if ( noisevol->getDataAt( ox, oy, oz ) == 1 )
// 					{
// 						continue ;
// 					}
// 
// 					flag = 0 ;
// 					for ( int m = 0 ; m < 6 ; m ++ )
// 					{
// 						int nx = ox + neighbor6[m][0] ;
// 						int ny = oy + neighbor6[m][1] ;
// 						int nz = oz + neighbor6[m][2] ;
// 						if ( getDataAt( nx, ny, nz ) > 0 && noisevol->getDataAt( nx, ny, nz ) == 1 )
// 						{
// 							noisevol->setDataAt( ox, oy, oz, 1 ) ;
// 							count ++ ;
// 							break ;
// 						}
// 					}
// 				}
// 
// 				if ( count == 0 )
// 				{
// 					break ;
// 				}
// 			}
// 
// 
// #endif
// 
// 			/* Commented for debugging
// 
// 			// First, 
// 			// check for complex nodes in queue2 
// 			// move them from queue2 to queue3
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Check simple 
// #ifndef NOISE_DIS_SHEET
// 				if ( isSheetEnd( ox, oy, oz ) || ! isSimple( ox, oy, oz ) ) 
// #else
// 				if ( isSheetEnd( ox, oy, oz, noisevol ) || ! isSimple( ox, oy, oz ) ) 
// #endif
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue3->prepend( ox, oy, oz ) ;
// 					ele = queue2->remove() ;
// 					
// 					numComplex ++ ;
// 				}
// 				else
// 				{
// 					ele = queue2->getNext() ;
// 				}
// 			}
// 			*/
// 
// 
// 			// Next,
// 			// Compute score for each node left in queue2
// 			// move them into priority queue
// 			queue2->reset() ;
// 			ele = queue2->getNext() ;
// 			while ( ele != NULL )
// 			{
// 				ox = ele->x ;
// 				oy = ele->y ;
// 				oz = ele->z ;
// 				
// 				// Compute score
// 				score = getNumPotComplex( ox, oy, oz ) ;
// 				scrvol->setDataAt( ox, oy, oz, score ) ;
// 
// 				// Push to queue
// 				gp = new gridPoint ;
// 				gp->x = ox ;
// 				gp->y = oy ;
// 				gp->z = oz ;
// 				// queue->add( gp, -score ) ;
// 				queue->add( gp, score ) ;
// 				
// 				ele = queue2->remove() ;				
// 			}
// 
// 			// Rename queue3 to be queue2, 
// 			// Clear queue3
// 			// From now on, queue2 holds nodes of next level
// 			delete queue2 ;
// 			queue2 = queue3 ;
// 			queue3 = new GridQueue2( ) ;
// 
// 
// 			// Next, start priority queue iteration
// 			while ( ! queue->isEmpty() )
// 			{
// 				// Retrieve the node with the highest score
// 				queue->remove( gp, score ) ;
// 				ox = gp->x ;
// 				oy = gp->y ;
// 				oz = gp->z ;
// 				delete gp ;
// 				// printf("%d\n", score);
// //				score = -score ;
// 
// 				// Ignore the node 
// 				// if it has been processed before
// 				// or it has an updated score
// 				if ( getDataAt( ox, oy, oz ) != curwid || (int) scrvol->getDataAt( ox, oy, oz ) != score )
// 				{
// 					continue ;
// 				}
// 
// 				/* Commented for debugging
// 
// 				// Remove this simple node
// 				setDataAt( ox, oy, oz, -1 ) ;
// 				numSimple ++ ;
// 				// printf("Highest score: %d\n", score) ;
// 				*/
// 
// 				/* Added for debugging */
// 				// Check simple 
// #ifndef NOISE_DIS_SHEET
// 				// if ( hasFeatureFace( ox, oy, oz ) )
// 				if ( (! isSimple( ox, oy, oz )) || isSheetEnd( ox, oy, oz ) )
// 				// if ( hasIsolatedFace(ox,oy,oz)  && (! isNoiseSheetEnd(ox,oy,oz))) 
// #else
// 				// if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) ) 
// 				if ( ! isSimple( ox, oy, oz ) || isSheetEnd( ox, oy, oz, noisevol ) || isHelixEnd( ox, oy, oz, noisevol )) 
// 				// if ( isBertrandEndPoint( ox, oy, oz ) ) 
// #endif
// 				{
// 					// Complex, set to next layer
// 					setDataAt( ox, oy, oz, curwid + 1 ) ;
// 					queue4->prepend( ox, oy, oz ) ;
// 					numComplex ++ ;
// 
// 				}
// 				else
// 				{
// 					setDataAt( ox, oy, oz, -1 ) ;
// 					numSimple ++ ;
// 
// 				}
// 				/* Adding ends */
// 
// 				// Move its neighboring unvisited node to queue2
// 				for ( int m = 0 ; m < 6 ; m ++ )
// 				{
// 					int nx = ox + neighbor6[m][0] ;
// 					int ny = oy + neighbor6[m][1] ;
// 					int nz = oz + neighbor6[m][2] ;
// 					if ( getDataAt( nx, ny, nz ) == 0 )
// 					{
// 						// setDataAt( nx, ny, nz, curwid + 1 ) ;
// 						queue2->prepend( nx, ny, nz ) ;
// 					}
// 				}
// 				
// 				/* Commented for debugging
// 				
// 				// Find complex nodes in its 3x3 neighborhood
// 				// move them to queue2
// 				for ( i = -1 ; i < 2 ; i ++ )
// 					for ( j = -1 ; j < 2 ; j ++ )
// 						for ( k = -1 ; k < 2 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							// Check simple 
// 							if ( getDataAt( nx, ny, nz ) == curwid && 
// 								// ( isSheetEnd( ox, oy, oz ) || ! isSimple( nx, ny, nz )) )
// #ifndef NOISE_DIS_SHEET
// 								( isSheetEnd( nx, ny, nz ) || ! isSimple( nx, ny, nz ) ) )
// #else
// 								( isSheetEnd( nx, ny, nz, noisevol ) || ! isSimple( nx, ny, nz ) ) )
// #endif
// 
// 							{
// 								// Complex, set to next layer
// 								setDataAt( nx, ny, nz, curwid + 1 ) ;
// 								queue2->prepend( nx, ny, nz ) ;
// 								numComplex ++ ;
// 							}
// 						}
// 				*/
// 
// 				// Update scores for nodes in its 5x5 neighborhood
// 				// insert them back into priority queue
// 					
// 				for ( i = -2 ; i < 3 ;i ++ )
// 					for ( j = -2 ; j < 3 ; j ++ )
// 						for ( k = -2 ; k < 3 ; k ++ )
// 						{
// 							int nx = ox + i ;
// 							int ny = oy + j ;
// 							int nz = oz + k ;
// 
// 							if ( getDataAt( nx, ny, nz ) == curwid )
// 							{
// 								// Compute score
// 								score = getNumPotComplex( nx, ny, nz ) ;
// 								
// 								if ( score != (int) scrvol->getDataAt( nx, ny, nz ) )
// 								{
// 									// printf("Update\n") ;
// 									scrvol->setDataAt( nx, ny, nz, score ) ;
// 									// Push to queue
// 									gp = new gridPoint ;
// 									gp->x = nx ;
// 									gp->y = ny ;
// 									gp->z = nz ;
// 									// queue->add( gp, -score ) ;
// 									queue->add( gp, score ) ;
// 								}
// 							}
// 						}
// 						
// 
// 			}
// 
// 			#ifdef VERBOSE
// 			printf("%d complex, %d simple\n", numComplex, numSimple) ;
// 			#endif			
// 
// 			if ( numSimple == 0 )
// 			{
// 					break ;
// 			}
// 		}
// 
// 		// Finally, clean up
// 		#ifdef VERBOSE
// 		printf("Thresholding the volume to 0/1...\n") ;
// 		#endif
// 		threshold( 0, 0, 1 ) ;
// 
// 		delete scrvol;
// 		delete queue;
// 		delete queue2;
// 		delete queue3;
// 		delete queue4;
// 
// 	}
