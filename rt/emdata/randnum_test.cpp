#include "randnum.h"
#include <iostream>

using namespace std;
using namespace EMAN;

int main()
{
	cout << "Start random number generator..." << endl;
	
	Randnum r = Randnum(gsl_rng_rand);
	for(int i=0; i<100; ++i) {
		cout << r.get_irand(0,99) << "\t" << r.get_frand() <<  "\t" << r.get_gauss_rand(1.0, 1.0) << endl;
	}
	
	return 0;
}
