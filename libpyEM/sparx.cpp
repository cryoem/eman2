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

#include "sparx.h"

using namespace EMAN;

void em2flex(EMData& image, Dict & header_dict, scitbx::af::shared<float> & flex_array)
{
	header_dict = image.get_attr_dict();
	size_t size = image.get_xsize() * image.get_ysize() * image.get_zsize();
	
	flex_array.reserve(size);
	float *data = image.get_data();
	
	for(size_t i = 0; i < size; i++) {
		flex_array.push_back(data[i]);
	}
}

void flex2em(const scitbx::af::shared<float> & flex_array, Dict & header_dict, EMData & image)
{
	int nx = header_dict["nx"];
	int ny = header_dict["ny"];
	int nz = header_dict["nz"];
	image.set_size(nx, ny, nz);

	size_t size = (size_t)nx * ny * nz;
	float * data = image.get_data();
	
	for(size_t i = 0; i < size; ++i) {
		data[i] = flex_array[i];
	}

	image.set_attr_dict(header_dict);
}
