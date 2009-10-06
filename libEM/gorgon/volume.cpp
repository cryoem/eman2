#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include "priority_queue.h"
#include <vector>
#include "volume.h"

using namespace std;
using namespace EMAN;

Volume::Volume(int que_max)
{
	que = new PriorityQueue<int, int>(que_max);
}
