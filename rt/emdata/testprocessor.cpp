#include "emdata.h"
using namespace EMAN;

/*
 *@author Wen Jiang
 *@date 2006-7-18
 */
 
int main()
{
	EMData *d;
	int size=512;
	d = new EMData();
	d->set_size(size,size);
	
	// create a random plane
	srand(time(0));	// get_frand() generates the same sequence of numbers everytime
	float xc=Util::get_frand(0,size), yc=Util::get_frand(0,size), zc=Util::get_frand(0,size);
	float norm[3] = {Util::get_frand(0,1),Util::get_frand(0,1),Util::get_frand(0,1)};
	for(int i=0; i<3; i++) norm[i]/=Util::hypot3(norm[0], norm[1], norm[2]); 
	
	for(int j=0; j<size; j++){
		for(int i=0; i<size; i++){
			d->set_value_at(i,j,0,-((i-xc)*norm[0]+(j-yc)*norm[1])/norm[2]+zc);
		}
	}
	float mean = d->get_attr("mean");
	
	float epsilon = 1e-5; // percent

	EMData* dmask=d->copy();
	dmask->process_inplace("eman1.mask.sharp", Dict("inner_radius", size/3-5, "outer_radius", size/3+5, "value", 0));
	
	d->process_inplace("filter.gradientPlaneRemover", Dict("mask", dmask));
	
	int err = 0;
	float max_residual = fabs(d->get_attr("maximum"));
	if ( max_residual > epsilon * mean) {
		err = 1;
		printf("FAILED: max residual=%g percent > error limit=%g\n", max_residual / mean, epsilon);
	}
	else {
		err = 0;
		printf("SUCCESS: max residual=%g percent < error limit=%g\n", max_residual / mean, epsilon);
	}

	return err;
}
	

	
		
		
