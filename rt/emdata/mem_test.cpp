#include <iostream>
#include "emdata.h"

using namespace std;
using namespace EMAN;

int main()
{
	cout << "Hello, memory leak valgrind test!" << endl;
	
	EMData *volume = new EMData(); // initial volume
	volume->read_image("vol001.tfc");
	if (volume) delete volume;

//	ENTERFUNC;
	
//	EXITFUNC;
	
	return 0;
}
