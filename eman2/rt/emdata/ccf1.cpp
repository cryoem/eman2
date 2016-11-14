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

#include "EMData.h"

const char* get_test_image(const char* img)
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/%s", getenv("HOME"), img);
		done = true;
	}
	return filename;
}

int main()
{
	EMData *a = new EMData();
	const char * img1 = get_test_image("3d86_1.mrc");
	const char * img2 = get_test_image("3d86_2.mrc");
	
	a->readImage(img1, -1);

	EMData *b = new EMData();
	b->readImage(img2, -1);

	EMData *c = a->calcCCF(b, 0);
	
	printf("nx = %d, ny = %d, nz = %d\n", a->xSize(), a->ySize(), a->zSize());

	printf("nx = %d, ny = %d, nz = %d\n", c->xSize(), c->ySize(), c->zSize());
	
	c->writeImage("aaatest1.mrc", -1);

	if( a )
	{
		delete a;
		a = 0;
	}

	if( b )
	{
		delete b;
		b = 0;
	}
	
	return 0;
}
