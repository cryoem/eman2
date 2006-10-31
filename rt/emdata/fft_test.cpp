/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

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
� //int nx = PROJ->get_xsize();
� //int ny = PROJ->get_ysize();
� //int nz = PROJ->get_zsize();
� //int n = nx*ny*nz;
� for(int i = 0;i<20;i++)
� {
� � EMData* W = �PROJ->pad_fft();
� � W->do_fft_inplace();// W = W->do_fft_inplace();
� � W->do_ift_inplace();// W = W->do_ift_inplace();
� � W->postift_depad_corner_inplace();
� � //float *src = W->get_data();
� � //float *dest = PROJ->get_data();
� � //memcpy(dest,src,sizeof(float)*n);
� � // W->set_size(1,1);
� � delete W;W = 0;
� �} �
}
*/
