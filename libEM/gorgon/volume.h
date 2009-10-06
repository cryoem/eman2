#ifndef VOLUME_H
#define VOLUME_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "priority_que.h"
#include <vector>

class Volume
{
private:
	PriorityQueue* que;
public:
	Volume(int que_max);
};

#endif