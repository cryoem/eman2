#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "priority_que.h"
#include <vector>
#include "volume.h"

using namespace std;
using namespace EMAN;

Volume::Volume(int que_max)
{
	que = new PriorityQueue<int, int>(que_max);
}
