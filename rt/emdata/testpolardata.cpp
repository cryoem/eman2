/**
 * $Id$
 */
#include "emdata.h"
#include "polardata.h"

#include <iostream>

using namespace std;
using namespace EMAN;

void test_polar_data() {
	cout << "Enter test_polar_data() function..." << endl;
	
	PolarData * pd = new PolarData();
	pd->print_polar();
	
	UnevenMatrix * um = pd;
	
	delete um;
	
	cout << "Leave test_polar_data() function..." << endl;
}

int main()
{
	cout << "--- Starting to test PolarData ---" << endl;
	
	try {
		test_polar_data();
	}
	catch (E2Exception & e)	{
		cout << e.what();
	}
	catch (...) {
		cout << "Unknown exception ???" << endl;
	}
	
	return 0;
}
