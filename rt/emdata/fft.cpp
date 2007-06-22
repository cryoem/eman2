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

#include "emdata.h"
#include "util.h"
#include <cstdio>

using namespace EMAN;

void create_lattice_image(const char* filename)
{
    int size = 256;
    int nslices = 4;
    int width = 5;
    
    EMData e;
    e.set_size(size, size, 1);

    float * data = e.get_data();
	
    for (int i = 2; i < nslices; i++) {
		int start = size/nslices*i - width;
		int end = size/nslices*i + width;
	
		for (int j = start; j < end; j++) {
			for (int k = 0; k < size/2; k++) {
				data[j*size+k]  = 1;
			}
		}
    }
    
    for (int l = 0; l < size; l++) {
		for (int i = 1; i < nslices-1; i++) {
			int start = size/nslices*i - width;
			int end = size/nslices*i + width;
	    
			for (int j = start; j < end; j++) {
				data[l*size+j] = 1;
			}
		}
    }
    
    e.update();
    e.write_image(filename, 0, EMUtil::IMAGE_MRC);
    
    return;
}
    


int main(int argc, char* argv[])
{
    int err = 0;
    Util::set_log_level(argc, argv);

    const char* filename = "lattice.mrc";
    create_lattice_image(filename);

	try {
		EMData e1;
		e1.read_image(filename);

		//EMData* fft1 = e1.do_fft();
		e1.do_fft_inplace();
		EMData *fft1 = &e1;
		fft1->write_image("lattice_fft.mrc", 0, EMUtil::IMAGE_MRC);

		EMData* fft1_copy = fft1->copy();
		
		EMData* ift1 = fft1->do_ift();	    
		ift1->write_image("lattice_fft_ift.mrc", 0, EMUtil::IMAGE_MRC);

		EMData* ift2 = fft1_copy->do_ift();	    
		ift2->write_image("lattice_fft_ift2.mrc", 0, EMUtil::IMAGE_MRC);

	}
	catch(...) {
		err = 1;
	}
    
    return err;
}
