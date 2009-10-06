#ifndef VOLUME_H
#define VOLUME_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "priority_que.h"
#include <vector>

namespace EMAN {

class Volume
{
private:
	PriorityQueue<int, int>* que;
public:
	Volume(int que_max);
};

}
#endif
