#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <cstdlib>
#include <cstdio>

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
	
	PriorityQueue ( int max );
	int get_que_length();
	int get_max_length();
};

#endif