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
#include <cstdio>

using namespace EMAN;

int main()
{
    int err = 0;
    EMData e1;

    char filename[1024];
    sprintf(filename, "%s/images/search.dm3", getenv("HOME"));

	try {
		e1.read_image(filename); 
 
		EMData e2;
		e2.set_size(e1.get_xsize(), e1.get_ysize(), e1.get_zsize());
		float* data = e2.get_data();
		int size = e1.get_xsize() * e1.get_ysize() * e1.get_zsize();
    
		for (int i = 0; i < size; i++) {
			data[i] = 0;
		}

		e2.update();
		e2.write_image("search.mrc", 0, EMUtil::IMAGE_MRC);
		Dict d = e2.get_attr_dict();
		EMUtil::dump_dict(d);
    }
	catch(...) {
		err = 1;
		throw;
	}
    return err;
}
