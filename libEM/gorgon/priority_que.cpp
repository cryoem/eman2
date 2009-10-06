#include <cstdlib>
#include <cstdio>
#include "priority_que.h"

using namespace std;

PriorityQue::PriorityQueue ( int max )
{
	this->maxLength = max ;
	this->queueLength = 0 ;
	this->valueQueue = new ValueT* [ max ] ;
	this->keyQueue = new KeyT [ max ] ;

};

int PriorityQueue::get_que_length()
{
	return queLength;
}

int PriorityQueue::get_max_length()
{
	return maxLength;
}