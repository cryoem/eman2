#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "priority_que.h"
#include <vector>
#include "volume.h"

using namespace std;

Volume::Volume(int que_max)
{
	que = new PriorityQue(que_max);
}
