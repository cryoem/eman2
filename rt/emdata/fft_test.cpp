#include <iostream>
#include "emdata.h"

using namespace std;
using namespace EMAN;

int main()
{
	cout << "Hello fft_test..." << endl;
	/*
	EMData * img = new EMData();
	img->set_size(5,1,1);
	img->set_value_at(0,0,0,1.0f);
	img->set_value_at(1,0,0,2.0f);
	img->set_value_at(2,0,0,3.0f);
	img->set_value_at(3,0,0,4.0f);
	img->set_value_at(4,0,0,5.0f);
	
	EMData * img2 = img->do_fft();
	*/
	for (int i=0; i<1; ++i) {
		EMData * w = new EMData();
		w->set_size(1024,1024);
		w->process_inplace("testimage.noise.gauss", Dict("mean", 0.0f, "sigma", 1.0f, "seed", 8888));
		
		w->do_fft_inplace();
		w->do_ift_inplace();
		w->postift_depad_corner_inplace();
		
		delete w;
		w = 0;
	}
	
	return 0;
}
/*
void Util::test1(EMData *PROJ)
{
  //int nx = PROJ->get_xsize();
  //int ny = PROJ->get_ysize();
  //int nz = PROJ->get_zsize();
  //int n = nx*ny*nz;
  for(int i = 0;i<20;i++)
  {
    EMData* W =  PROJ->pad_fft();
    W->do_fft_inplace();// W = W->do_fft_inplace();
    W->do_ift_inplace();// W = W->do_ift_inplace();
    W->postift_depad_corner_inplace();
    //float *src = W->get_data();
    //float *dest = PROJ->get_data();
    //memcpy(dest,src,sizeof(float)*n);
    // W->set_size(1,1);
    delete W;W = 0;
   }  
}
*/
