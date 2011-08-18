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

	dict["OMAP.xstart"] = omaph.xstart;
	dict["OMAP.ystart"] = omaph.ystart;
	dict["OMAP.zstart"] = omaph.zstart;
	dict["nx"] = omaph.nx;
	dict["ny"] = omaph.ny;
	dict["nz"] = omaph.nz;
	dict["apix_x"] = omaph.apix_x;
	dict["apix_y"] = omaph.apix_y;
	dict["apix_z"] = omaph.apix_z;
	dict["OMAP.header10"] = omaph.header10;
	dict["OMAP.header11"] = omaph.header11;
	dict["OMAP.header12"] = omaph.header12;
	dict["alpha"] = omaph.alpha;
	dict["beta"] = omaph.beta;
	dict["gamma"] = omaph.gamma;
	dict["OMAP.header16"] = omaph.header16;
	dict["OMAP.header17"] = omaph.header17;
	dict["OMAP.scale"] = omaph.scale;
	dict["OMAP.header19"] = omaph.header19;

	EXITFUNC;
	return 0;
}

int OmapIO::read_data(float *rdata, int, EMAN::Region const*, bool)
{
	ENTERFUNC;
	std::cout << "file pointer location = " << ftell(omapfile) << std::endl;

	unsigned char * cdata = (unsigned char *) rdata;
	size_t size = (size_t)omaph.nx*omaph.ny*omaph.nz;
	if (fread(cdata, size, 1, omapfile) != 1) {
		throw ImageReadException(filename, "OMAP data");
	}

	float density_factor = (float)omaph.header16/omaph.header19 + omaph.header17;

	std::cout << "density_factor = " << density_factor << std::endl;

	for (size_t i = 0; i < size; ++i) {
		size_t j = size - 1 - i;
		rdata[j] = static_cast < float >(cdata[j]) * density_factor;
	}

	EXITFUNC;
	return 0;
}

int OmapIO::write_header(EMAN::Dict const&, int, EMAN::Region const*, EMAN::EMUtil::EMDataType, bool)
{
	ENTERFUNC;

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

//	std::cout << "const_value = " << const_value
//			<< ", xstart = " << xstart
//			<< ", ystart = " << ystart
//			<< ", zstart = " << zstart
//			<< ", nx = " << nx
//			<< ", ny = " << ny
//			<< ", nz = " << nz
//			<< std::endl;

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
