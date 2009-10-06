#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <cstdlib>
#include <cstdio>

namespace EMAN {

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

	PriorityQueue ( int max )
	{
		this->maxLength = max ;
		this->queueLength = 0 ;
		this->valueQueue = new ValueT* [ max ] ;
		this->keyQueue = new KeyT [ max ] ;

	}

	int get_que_length()
	{
		return queueLength;
	}

	int get_max_length()
	{
		return maxLength;
	}
};


}

#endif
