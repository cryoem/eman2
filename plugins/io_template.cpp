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

#include "io_template.h"
#include "portable_fileio.h"

using namespace EMAN;

XYZIO::XYZIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), xyz_file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
}

XYZIO::~XYZIO()
{
	if (xyz_file) {
		fclose(xyz_file);
		xyz_file = 0;
	}
}

void XYZIO::init()
{
	if (initialized) {
		return ;
	}

	ENTERFUNC;
	
	initialized = true;
	bool is_new_file = false;
	xyz_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {

	}

	EXITFUNC;
}

bool XYZIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;
	if (!first_block) {
		result = false;
	}

	// check image format validality here
	
	EXITFUNC;
	return result;
}

int XYZIO::read_header(Dict & , int image_index, const Region * , bool )
{
	ENTERFUNC;
	check_read_access(image_index);

	// read header info here
	
	EXITFUNC;
	return 0;
}

int XYZIO::write_header(const Dict & , int image_index, const Region *,
						EMUtil::EMDataType , bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index);
	// write header info here
	EXITFUNC;
	return 0;
}

int XYZIO::read_data(float *data, int image_index, const Region * , bool )
{
	ENTERFUNC;
	check_read_access(image_index, data);

	// read image data here

	EXITFUNC;
	return 0;
}

int XYZIO::write_data(float *data, int image_index, const Region *,
					  EMUtil::EMDataType , bool)
{
	ENTERFUNC;
	check_write_access(rw_mode, image_index, 0, data);

	// write image data here
	
	EXITFUNC;
	return 0;
}

void XYZIO::flush()
{
	if (xyz_file) {
		fflush(xyz_file);
	}
}

bool XYZIO::is_complex_mode()
{
	return false;
}

bool XYZIO::is_image_big_endian()
{
	return is_big_endian;
}

int XYZIO::get_nimg()
{
	init();

	return 1;
}
