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
#include <sstream>
using std::stringstream;

#include <iostream>
using std::cout;
using std::endl;

using namespace EMAN;

int main(int argc, char *argv[])
{
	int SIZE = 96;
	int NTT = 500;

	EMData pat;
	pat.set_size(SIZE, SIZE, 1);
	float *d = pat.get_data();

	for (int i = -SIZE / 2; i < SIZE / 2; i++) {
		for (int j = -SIZE / 2; j < SIZE / 2; j++) {
			d[(i + SIZE / 2) + (j + SIZE / 2) * SIZE] =
				-3.0f * exp(-Util::square((fabs((float)i) + fabs((float)j))) / 10.0f) +
				exp(-Util::square((fabs((float)i) + fabs((float)j) / 2.0f)) / 100.0f) *
				(abs(i) < 2 ? 2.0f : 1.0f);
		}
	}
	pat.update();
	pat.process_inplace("normalize.circlemean");
	pat.process_inplace("mask.sharp", Dict("outer_radius", pat.get_xsize()/2));

	EMData *data;
	data = pat.copy();
	float *dd = data->get_data();

	for (int j = 0; j < SIZE * SIZE; j++) {
		dd[j] += Util::get_gauss_rand(0, 1.0);
	}
	data->update();
	data->process_inplace("normalize.circlemean");
	data->process_inplace("mask.sharp", Dict("outer_radius", data->get_xsize()/2));

	EMData *tmp = 0;

	tmp = data->align("rotate_translate_flip", data, Dict());
	if(tmp)	{delete tmp;tmp = 0;}

	if(data) {delete data; data=0;}

    return 0;
}


