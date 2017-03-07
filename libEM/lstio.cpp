/**
 * $Id$
 */

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

#ifndef WIN32
#include <sys/param.h>
#include <unistd.h>
#else
#include <direct.h>
#include <windows.h>
#define M_PI 3.14159265358979323846f
#define MAXPATHLEN (MAX_PATH*4)
#endif

#include <cstdio>
#include <cstring>
#include "lstio.h"
#include "util.h"


using namespace EMAN;

const char *LstIO::MAGIC = "#LST";

LstIO::LstIO(const string & file, IOMode rw)
:	filename(file), rw_mode(rw), lst_file(0)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	initialized = false;
	nimg = 0;
	imageio = 0;
	ref_filename = "";
	last_lst_index = -1;
	last_ref_index = -1;
}

LstIO::~LstIO()
{
	if (lst_file) {
		fclose(lst_file);
		lst_file = 0;
	}
	ref_filename = "";
	if(imageio) {
		delete imageio;
		imageio = 0;
	}
}

void LstIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return ;
	}

	initialized = true;

	bool is_new_file = false;
	lst_file = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {

		char buf[MAXPATHLEN];

		if (!fgets(buf, MAXPATHLEN, lst_file)) {
			throw ImageReadException(filename, "first block");
		}

		if (!is_valid(&buf)) {
			throw ImageReadException(filename, "invalid LST file");
		}

		for (nimg = 0; fgets(buf, MAXPATHLEN, lst_file) != 0; nimg++) {
			if (buf[0] == '#') {
				nimg--;
			}
		}
		rewind(lst_file);
	}
	EXITFUNC;
}

bool LstIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;

	if (!first_block) {
		result = false;
	}
	else {
		result = Util::check_file_by_magic(first_block, MAGIC);
	}

	EXITFUNC;
	return result;
}

int LstIO::calc_ref_image_index(int image_index)
{
	if (image_index == last_lst_index) {
		return last_ref_index;
	}
	else {
		char buf[MAXPATHLEN];
		int step = image_index - last_lst_index;

		if (step < 0) {
			rewind(lst_file);
			step = image_index + 1;
		}

		for (int i = 0; i < step; i++) {
			if (!fgets(buf, MAXPATHLEN, lst_file)) {
				LOGERR("reach EOF in file '%s' before reading %dth image",
					   filename.c_str(), image_index);
				return 1;
			}
			if (buf[0] == '#') {
				i--;
			}
		}
		int ref_image_index = 0;
		char ref_image_path[MAXPATHLEN];
		char unused[256];
		sscanf(buf, " %d %s %[ .,0-9-]", &ref_image_index, ref_image_path, unused);

		char fullpath[MAXPATHLEN];

		char sep = '/';
#ifdef WIN32
		sep = '\\';
#endif
// 		if (ref_image_path[0] == sep) {
			strcpy(fullpath, ref_image_path);
// 		}
// 		else {
// 			if (strrchr(filename.c_str(), sep)) {
// 				strcpy(fullpath, filename.c_str());
// 			}
// 			else {
// #ifndef WIN32
// 				getcwd(fullpath, MAXPATHLEN);
// #else
// 				//GetCurrentDirectory(MAXPATHLEN, fullpath);
// #endif
// 			}
// 
// 			char *p_basename = strrchr(fullpath, sep);
// 			if (p_basename) {
// 				//p_basename++;
// 				//*p_basename = '\0';
// 				char ssep[2];
// 				ssep[0] = sep;
// 				ssep[1] = '\0';
// 				strcat(fullpath, ssep);
// 				strcat(fullpath, ref_image_path);
// 			}
// 		}

		ref_filename = string(fullpath);
		imageio = EMUtil::get_imageio(ref_filename, rw_mode);

		last_ref_index = ref_image_index;
	}

	last_lst_index = image_index;

	return last_ref_index;
}


int LstIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;
	init();
	check_read_access(image_index);
	int ref_image_index = calc_ref_image_index(image_index);
	int err = imageio->read_header(dict, ref_image_index, area, is_3d);
	dict["source_path"] = ref_filename;
	EXITFUNC;
	return err;
}

int LstIO::write_header(const Dict &, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	init();
	fprintf(lst_file, "%s\n", MAGIC);
	EXITFUNC;
	return 0;
}

int LstIO::read_data(float *data, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;
	check_read_access(image_index, data);
	int ref_image_index = calc_ref_image_index(image_index);
	int err = imageio->read_data(data, ref_image_index, area, is_3d);
	EXITFUNC;
	return err;
}

int LstIO::write_data(float *data, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	fprintf(lst_file, "%s\n", (char*)data);
	EXITFUNC;
	return 0;
}

void LstIO::flush()
{
	fflush(lst_file);
}

bool LstIO::is_complex_mode()
{
	return false;
}

bool LstIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int LstIO::get_nimg()
{
	init();
	return nimg;
}
