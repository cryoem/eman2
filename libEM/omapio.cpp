/**
 * $Id$
 */

/*
 * Author: Grant Tang, 06/07/2011 (gtang@bcm.edu)
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

#include "omapio.h"
#include "portable_fileio.h"

using namespace EMAN;

OmapIO::OmapIO(const string & omapname, IOMode rw) :
		filename(omapname), rw_mode(rw), omapfile(0),
		is_big_endian(false), initialized(false), is_new_file(false)
{
	memset(&omaph, 0, sizeof(OmapHeader));
	is_big_endian = ByteOrder::is_host_big_endian();
}

OmapIO::~OmapIO()
{
	if (omapfile) {
		fclose(omapfile);
		omapfile = 0;
	}
}

void OmapIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	initialized = true;
	omapfile = sfopen(filename, rw_mode, &is_new_file);

	char record[512];	//record size 512 bytes

	if (!is_new_file) {
		if (fread(record, 512, 1, omapfile) != 1) {
			throw ImageReadException(filename, "OMAP header");
		}

		for(int i=0; i<512; i++) {
			if(!isprint(record[i])) {
				portable_fseek(omapfile, 0, SEEK_SET);
				break;
			}

			if(record[i] == '\0') break;
		}

		if (fread(&omaph, sizeof(OmapHeader), 1, omapfile) != 1) {
			throw ImageReadException(filename, "OMAP header");
		}

		if (!is_valid(&omaph)) {
			throw ImageReadException(filename, "invalid OMAP");
		}

		if(!ByteOrder::is_host_big_endian()) {	//omap first record is always  big endian
			ByteOrder::swap_bytes((short *) &omaph, 256);	//each record 512 bytes, 256 short intergers
		}
	}

	EXITFUNC;
}

int OmapIO::read_header(EMAN::Dict& dict, int, EMAN::Region const*, bool)
{
	ENTERFUNC;
	init();

	dict["OMAP.xstart"] = (int)omaph.xstart;
	dict["OMAP.ystart"] = (int)omaph.ystart;
	dict["OMAP.zstart"] = (int)omaph.zstart;
	dict["nx"] = (int)omaph.nx;
	dict["ny"] = (int)omaph.ny;
	dict["nz"] = (int)omaph.nz;
	dict["apix_x"] = (int)omaph.apix_x;
	dict["apix_y"] = (int)omaph.apix_y;
	dict["apix_z"] = (int)omaph.apix_z;
	dict["OMAP.cellA"] = (int)(omaph.header10/omaph.scale);
	dict["OMAP.cellB"] = (int)(omaph.header11/omaph.scale);
	dict["OMAP.cellC"] = (int)(omaph.header12/omaph.scale);
	dict["alpha"] = (int)(omaph.alpha/omaph.scale);
	dict["beta"] = (int)(omaph.beta/omaph.scale);
	dict["gamma"] = (int)(omaph.gamma/omaph.scale);
	dict["OMAP.iprod"] = (int)omaph.iprod;
	dict["OMAP.iplus"] = (int)omaph.iplus;
	dict["OMAP.scale"] = (int)omaph.scale?omaph.scale:100;
	dict["OMAP.scale2"] = (int)omaph.scale2?omaph.scale2:100;

	float prod = (float)(omaph.iprod/(int)dict["OMAP.scale2"]);
	float plus = (float)omaph.iplus;

	dict["OMAP.min"] = (omaph.imin - plus)/prod;
	dict["OMAP.imax"] = (omaph.imax - plus)/prod;
	dict["OMAP.sigma"] = (omaph.isigma - plus)/prod;
	dict["OMAP.mean"] = (omaph.imean - plus)/prod;

	if((float)dict["OMAP.sigma"] < 0.001f || (float)dict["OMAP.sigma"] > 10.0f) {
		std::cout << "Warning : Suspect value of sigma : " << (float)dict["OMAP.sigma"] << std::endl;
		dict["OMAP.sigma"] = 0;		//flag bad sigma
	}

	EXITFUNC;
	return 0;
}

int OmapIO::read_data(float *rdata, int, EMAN::Region const*, bool)
{
	ENTERFUNC;
//	std::cout << "file pointer location = " << ftell(omapfile) << std::endl;

	// cubes in each direction
	int inx = omaph.nx/8;
	int iny = omaph.ny/8;
	int inz = omaph.nz/8;

	//include partial cube
	if(omaph.nx%8 > 0) ++inx;
	if(omaph.ny%8 > 0) ++iny;
	if(omaph.nz%8 > 0) ++inz;

	// How much element of last cube to read ?
	int xtraX = omaph.nx%8;
	int xtraY = omaph.ny%8;
	int xtraZ = omaph.nz%8;

//	std::cout << "Total records = " << inx*iny*inz << std::endl;

	float prod = (float)omaph.iprod/omaph.scale2?omaph.scale2:100;
	float plus = omaph.iplus;
	float pixel = 0.0f;
	unsigned char record[512];
	for (int k=0; k < inz; k++) {
		for (int j=0;  j < iny; j++) {
			for (int i=0; i < inx; i++) {
				if (fread(record, 512, 1, omapfile) != 1) {
					throw ImageReadException(filename, "OMAP data");
				}

				//swap bytes for little endian
				bool byteswap = false;
				if(!ByteOrder::is_host_big_endian()) {
					byteswap = true;
				}
				if(byteswap) {
					for (int ii=0; ii < 511; ii+=2) {
						char tempchar = record[ii];
						record[ii] = record[ii+1];
						record[ii+1] = tempchar;
					}
				}

				int cubieSizeX = 8;
				int cubieSizeY = 8;
				int cubieSizeZ = 8;

				//for partial cube
				if (xtraX > 0) if (i == inx-1) cubieSizeX = xtraX;
				if (xtraY > 0) if (j == iny-1) cubieSizeY = xtraY;
				if (xtraZ > 0) if (k == inz-1) cubieSizeZ = xtraZ;

				for (int n=0; n < cubieSizeZ; n++) {
					for (int m=0;  m < cubieSizeY; m++) {
						for (int l=0; l < cubieSizeX; l++) {
							unsigned char sboxLMN = record[8*8*n+8*m+l];

							pixel = ((float)sboxLMN - plus)/prod;

							int pt3 = k*8 + n;
							int pt2 = j*8 + m;
							int pt1 = i*8 + l;

							rdata[pt3*omaph.nx*omaph.ny + pt2*omaph.nx + pt1] = (pixel-plus)/prod;
							if(omaph.isigma>0) rdata[pt3*omaph.nx*omaph.ny + pt2*omaph.nx + pt1] /= (float)omaph.isigma;
						}
					}
				}

			}
		}
	}

	EXITFUNC;
	return 0;
}

int OmapIO::write_header(EMAN::Dict const&, int, EMAN::Region const*, EMAN::EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	throw ImageWriteException("N/A", "No writing for Omap images");

	EXITFUNC;
	return 0;
}

int OmapIO::write_data(float*, int, EMAN::Region const*, EMAN::EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	EXITFUNC;
	return 0;
}

bool OmapIO::is_valid(const void *first_block, off_t file_size)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const short *data = static_cast < const short *>(first_block);
	short xstart = data[0];
	short ystart = data[1];
	short zstart = data[2];
	short nx = data[3];
	short ny = data[4];
	short nz = data[5];
	short const_value = data[18];

	if(!ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&xstart);
		ByteOrder::swap_bytes(&ystart);
		ByteOrder::swap_bytes(&zstart);
		ByteOrder::swap_bytes(&nx);
		ByteOrder::swap_bytes(&ny);
		ByteOrder::swap_bytes(&nz);
		ByteOrder::swap_bytes(&const_value);
	}

	if(const_value != 100) return false;
	if(nx<=0 || ny<=0 || nz<=0 || nx>10000 || ny>10000 || nz>10000) return false;

	EXITFUNC;
	return true;
}

bool OmapIO::is_image_big_endian()
{
	return true;
}

bool OmapIO::is_complex_mode()
{
	return false;
}

void OmapIO::flush()
{
	fflush(omapfile);
}
