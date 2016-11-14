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
#include "Euler.h"

const char* get_test_image()
{
	static char filename[256];
	static bool done = false;
	if (!done) {
		sprintf(filename, "%s/images/3d99.hed", getenv("HOME"));
		done = true;
	}
	return filename;
}

void rotate(EMData * image, float phi)
{
	char outfile[128];
	sprintf(outfile, "test_%d.hed", (int)phi);
	float f = (float)M_PI / 180;
	const char* imagefile = get_test_image();
	image->readImage(imagefile, -1);
	image->setRAlign(0,0,phi*f);
	image->rotateAndTranslate();
	image->writeIMAGIC3D(outfile);
}

int main()
{
	int err = 0;
	
	EMData *image = new EMData();
#if 0
	for (int i = 0; i <= 180; i += 45) {
		for (int j = 0; j < 360; j += 60) {
			for (int z = 0; z < 360; z += 45) {
				rotate(image, i, j, z);
			}
		}
	}
#endif

	rotate(image, 45);

	if( image )
	{
		delete image;
		image = 0;
	}
	
	return  err;
}
