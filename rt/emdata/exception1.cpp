#include "emdata.h"
#include "log.h"

using namespace EMAN;

int main()
{
	try {
		EMData e1;
		e1.set_size(10,10,10);
		e1.dot(0, false);
	}
	catch (Exception & e) {
		LOGERR("%s", e.what());
	}
	
	return 0;
}

