/*
 * Author: John Flanagan, 27/10/2010 (jfflanag@bcm.edu)
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
#include "emdata.h"
#include "testutil.h"


using namespace EMAN;

int main(int   argc, char *argv[])
{
#ifdef EMAN2_USING_CUDA
	EMData *a = new EMData();
	EMData *b = new EMData();
	//time_t start,end;
	
	a->read_image(argv[1]);
	b->read_image(argv[2]);
	
	//EMData *sa = a->process("xform.scale", Dict("scale",2.90909090909090909090909));
	//EMData *sb = b->process("xform.scale", Dict("scale",2.90909090909090909090909));

	//sa->copy_to_cudaro();
	//sb->copy_to_cudaro();
	//sb->copy_to_cuda();

	//EMData *sc = sa->align("rt.3d.sphere", sb, Dict("delta",10,"dphi",10), "ccc.tomo");
	
	//time (&start);
	//for(int i = 0; i < 100000; i++)
	//{
	//	EMData* fft = a->do_fft();
	//	delete fft;
		//a->elementaccessed();
	//}
	//time (&end);
	//float diff = difftime (end,start);
	//printf ("EMAN used %f seconds\n", diff);
	
	a->copy_to_cuda();
	b->copy_to_cuda();
	//int size = (a->get_xsize() + 2)*a->get_ysize()*sizeof(float);
	//EMData::usemempool(size);
	
	//time (&start);
	for(int i = 0; i < 25000; i++)
	{
		//cout << "i " << i << endl;

		EMData* ra = a->align("rotational", b, Dict(), "dot");
		//delete ra;
		a->elementaccessed();
		
	}
	//time (&end);
	//diff = difftime (end,start);
	//printf ("CUDA used %f seconds\n", diff);

	//cout << argv[1] << endl;
	//EMData *b = new EMData();]
        //Transform t(Dict("type","2d","alpha",45));
        //EMData *b = a->process("xform",Dict("transform",(Transform*)&t));
	
	//EMData *c = a->calc_ccf(b);
	//if( b )
	//{
	//	delete b;
	//	b = 0;
	//}

	

	//EMData *fft = a->do_fft();
	//EMData *c = fft->copy(false);

	
	//c->write_image("aatest2.mrc");

	//if( a )
	//{ 
	//	delete a;
	//	a = 0;
	//}
#endif
	return 0;
}
