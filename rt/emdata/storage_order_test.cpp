#include "emdata.h"

#include <iostream>

using namespace EMAN;
using namespace std;

int main()
{
	std::cout << "EMAN2 storage order test ..." << endl;
	
	int nx=3, ny=2;
	
	EMData * e = new EMData();
	e->set_size(nx,ny);
	e->set_value_at(0,0,0.0);
	e->set_value_at(1,0,1.0);
	e->set_value_at(2,0,2.0);
	e->set_value_at(0,1,3.0);
	e->set_value_at(1,1,4.0);
	e->set_value_at(2,1,5.0);
	
	cout << "print out array in a for loop: " << endl;
	for(int j=0; j<ny; j++) {
		for(int i=0; i<nx; i++) {
			cout << e->get_value_at(i, j) << " ";
		}
		cout << endl;
	}
	
	cout << endl << "Print the storage order in memory: " << endl;
	float * rdata = e->get_data();
	for(int i=0; i<nx*ny; i++) {
		cout << rdata[i] << "\t";
	}
	cout << endl;
	
	return 0;
}
