// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju (taoju@cse.wustl.edu)
// Description:   A Priority queue implementation

#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <cstdlib>
#include <cstdio>

namespace EMAN {

/**
 * Template class for a priority queue. The smallest element is at the front
 */
template < class ValueT, class KeyT >
class PriorityQueue
{
public:
	/// Number of elements in queue
	int queueLength ;

	/// Maximum number of elements int he queue
	int maxLength ;

	/// Queue of elements
	ValueT ** valueQueue ;

	/// Queue of keys
	KeyT * keyQueue ;

public:
	
	/** 
	 * Constructor
	 */
	PriorityQueue ( int max )
	{
		this->maxLength = max ;
		this->queueLength = 0 ;
		this->valueQueue = new ValueT* [ max ] ;
		this->keyQueue = new KeyT [ max ] ;

	};

	/**
	 * Destructor
	 */
	~PriorityQueue ()
	{
		delete [] keyQueue ;
		for ( int i = 0 ; i < queueLength ; i ++ )
		{
			delete valueQueue [ i ] ;
		}
		delete [] valueQueue ;
	};

	/**
	 * Get current length
	 */
	int getLength ( )
	{
		return this->queueLength ;
	};

	/**
	 * Test whether empty 
	 */
	bool isEmpty ( )
	{
		return ( this->queueLength == 0 ) ;
	};

	/**
	 * Test whether full 
	 */
	bool isFull ( )
	{
		return ( this->queueLength == this->maxLength  ) ;
	};

	/**
	 * Add an element
	 */
	void add ( ValueT * v, KeyT k )
	{
		if ( this->isFull() )
		{
			printf("PRIORITY QUEUE FILLED UP !!! \n");
			return ;
		}

		int ind = queueLength ;
		int tind ;
		queueLength ++ ;

		while ( ind > 0 )
		{
			tind = ( ind + 1 ) / 2 - 1 ;
			if ( k < keyQueue[tind] )
			{
				keyQueue[ ind ] = keyQueue [ tind ] ;
				valueQueue [ ind ] = valueQueue [ tind ] ;
				ind = tind ;
			}
			else
			{
				break;
			}
		}

		valueQueue[ ind ] = v ;
		keyQueue [ ind ] = k ;
	};

	/**
	 * Remove an element
	 */
	void remove ( ValueT *& v, KeyT & k )
	{
		if ( this->isEmpty() )
		{
			v = NULL ;
			k = 0 ;
			return ;
		}

		v = valueQueue[0] ;
		k = keyQueue[0] ;
		queueLength -- ;

		if ( queueLength == 0 )
		{
			valueQueue[0] = NULL ;
			return ;
		}

		ValueT * vv = valueQueue [ queueLength ] ;
		KeyT kk = keyQueue [ queueLength ], lowk ;
		int ind = 0, tind, ind2, ind3 ;
		while ( 1 )
		{
			ind2 = 2 * ( ind + 1 ) - 1 ;
			ind3 = ind2 + 1 ;
			tind = ind ; 
			lowk = kk ;

			if ( ind2 >= queueLength )
			{
				break ;
			}
			else 
			{
				if ( keyQueue [ ind2 ] < lowk )
				{
					tind = ind2 ;
					lowk = keyQueue [ ind2 ] ;
				}

				if ( ind3 < queueLength )
				{
					if ( keyQueue [ ind3 ] < lowk )
					{
						tind = ind3 ;
					}
				}

				if ( ind != tind )
				{
					valueQueue [ ind ] = valueQueue [ tind ] ;
					keyQueue [ ind ] = keyQueue [ tind ] ;
					ind = tind ;
				}
				else
				{
					break ;
				}
			}
		}

		valueQueue [ ind ] = vv ;
		keyQueue [ ind ] = kk ;
		valueQueue [ queueLength ] = NULL ;
		keyQueue [ queueLength ] = NULL ;
	};

};


}

#endif
