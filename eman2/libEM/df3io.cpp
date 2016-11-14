/**
 * $Id$
 */

/*
 * Author: Grant Tang, 08/06/2010 (gtang@bcm.edu)
 * Copyright (c) 2000-2010 Baylor College of Medicine
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

#include <algorithm>
#include <climits>

#include "df3io.h"
#include "portable_fileio.h"

using namespace EMAN;

Df3IO::Df3IO(const string & df3_filename, IOMode rw)
:	filename(df3_filename), rw_mode(rw), df3file(0),
 	initialized(false), is_new_file(false), rendermin(0.0), rendermax(0.0)
{
}

Df3IO::~Df3IO()
{
	if (df3file) {
		fclose(df3file);
		df3file = 0;
	}
}

void Df3IO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	initialized = true;
	df3file = sfopen(filename, rw_mode, &is_new_file);

	EXITFUNC;
}

// header is 3 MSB unsigned shorts (nx,ny,nz) in big-endian
int Df3IO::read_header(Dict & dict, int, const Region *, bool )
{
	ENTERFUNC;
	init();
	size_t nr;

	if (!is_new_file) {
		if (fread(&nx, sizeof(unsigned short), 1, df3file) != 1) {
			throw ImageReadException(filename, "DF3 header");
		}

		nr = fread(&ny, sizeof(unsigned short), 1, df3file); nr++;
		nr = fread(&nz, sizeof(unsigned short), 1, df3file); nr++;

		if(!ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&nx);
			ByteOrder::swap_bytes(&ny);
			ByteOrder::swap_bytes(&nz);
		}

		dict["nx"] = (int)nx;
		dict["ny"] = (int)ny;
		dict["nz"] = (int)nz;
	}

	EXITFUNC;
	return 0;
}

int Df3IO::write_header(const Dict & dict, int, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	init();

	nx = (unsigned short)((int)dict["nx"]);
	ny = (unsigned short)((int)dict["ny"]);
	nz = (unsigned short)((int)dict["nz"]);

	portable_fseek(df3file, 0, SEEK_SET);

	unsigned short df3header[3];
	df3header[0] = nx;
	df3header[1] = ny;
	df3header[2] = nz;
	ByteOrder::become_big_endian(df3header, 3);

	EMUtil::getRenderLimits(dict, rendermin, rendermax);

	if(fwrite(df3header, sizeof(unsigned short), 3, df3file) != 3) {
		throw ImageWriteException(filename, "DF3 header");
	}

	EXITFUNC;
	return 0;
}

int Df3IO::read_data(float *rdata, int, const Region *, bool)
{
	ENTERFUNC;

	size_t image_size = (size_t)nx*ny*nz;
	size_t nr;

	// obtain file size:
	portable_fseek (df3file , 0 , SEEK_END);
	size_t fsize = ftell (df3file);
	rewind (df3file);

	unsigned int * uidata = 0;
	unsigned short * usdata = 0;
	unsigned char * ucdata = 0;

	portable_fseek(df3file, sizeof(unsigned short)*3, SEEK_SET);	//skip header
	switch(fsize/image_size) {
	case sizeof(unsigned int):
		uidata = new unsigned int[image_size];
		nr = fread(uidata, sizeof(unsigned int), image_size, df3file); nr++;
		become_host_endian < unsigned int >(uidata, image_size);
		std::copy(uidata, uidata+image_size, rdata);
		if(uidata) {delete [] uidata; uidata=0;}
		break;
	case sizeof(unsigned short):
		usdata = new unsigned short[image_size];
		nr = fread(usdata, sizeof(unsigned short), image_size, df3file); nr++;
		become_host_endian < unsigned short >(usdata, image_size);
		std::copy(usdata, usdata+image_size, rdata);
		if(usdata) {delete [] usdata; usdata=0;}
		break;
	case sizeof(unsigned char):
		ucdata = new unsigned char[image_size];
		nr = fread(ucdata, sizeof(unsigned char), image_size, df3file); nr++;
		std::copy(ucdata, ucdata+image_size, rdata);
		if(ucdata) {delete [] ucdata; ucdata=0;}
		break;
	default:
		throw ImageReadException(filename, "DF3 does not support this kind of data type.");
	}

	EXITFUNC;
	return 0;
}

int Df3IO::write_data(float *data, int, const Region*,
					  EMUtil::EMDataType dt, bool)
{
	ENTERFUNC;

	size_t img_size = (size_t)nx*ny*nz;
	unsigned int * uidata = 0;
	unsigned short * usdata = 0;
	unsigned char * ucdata = 0;

	EMUtil::getRenderMinMax(data, nx, ny, rendermin, rendermax, nz);

	switch(dt) {
	case EMUtil::EM_UINT:
		uidata = new unsigned int[img_size];
		for (size_t i = 0; i < img_size; ++i) {
			if(data[i] <= rendermin) {
				uidata[i] = 0;
			}
			else if(data[i] >= rendermax) {
				uidata[i] = UINT_MAX;
			}
			else {
				uidata[i]=(unsigned int)((data[i]-rendermin)/(rendermax-rendermin)*UINT_MAX);
			}
		}
		ByteOrder::become_big_endian(uidata, img_size);
		if(fwrite(uidata, sizeof(unsigned int), img_size, df3file) != img_size) {
			throw ImageWriteException(filename, "DF3 unsigned int data");
		}
		if(uidata) {delete [] uidata; uidata=0;}
		break;
	case EMUtil::EM_USHORT:
		usdata = new unsigned short[img_size];
		for (size_t i = 0; i < img_size; ++i) {
			if(data[i] <= rendermin) {
				usdata[i] = 0;
			}
			else if(data[i] >= rendermax) {
				usdata[i] = USHRT_MAX;
			}
			else {
				usdata[i]=(unsigned short)((data[i]-rendermin)/(rendermax-rendermin)*USHRT_MAX);
			}
		}
		ByteOrder::become_big_endian(usdata, img_size);
		if(fwrite(usdata, sizeof(unsigned short), img_size, df3file) != img_size) {
			throw ImageWriteException(filename, "DF3 unsigned short data");
		}
		if(usdata) {delete [] usdata; usdata=0;}
		break;
	case EMUtil::EM_UCHAR:
		ucdata = new unsigned char[img_size];
		for (size_t i = 0; i < img_size; ++i) {
			if(data[i] <= rendermin) {
				ucdata[i] = 0;
			}
			else if(data[i] >= rendermax){
				ucdata[i] = UCHAR_MAX;
			}
			else {
				ucdata[i]=(unsigned char)((data[i]-rendermin)/(rendermax-rendermin)*UCHAR_MAX);
			}
		}
		if(fwrite(ucdata, sizeof(unsigned char), img_size, df3file) != img_size) {
			throw ImageWriteException(filename, "DF3 unsigned char data");
		}
		if(ucdata) {delete [] ucdata; ucdata=0;}
		break;
	default:
		throw ImageWriteException(filename,"DF3 does not support this data format");
	}

	EXITFUNC;
	return 0;
}

void Df3IO::flush()
{
	fflush(df3file);
}

bool Df3IO::is_image_big_endian()
{
	return true;	//DF3 file is always big endian
}

bool Df3IO::is_complex_mode()
{
	return false;	//DF3 only support real image
}

bool Df3IO::is_valid(const void *first_block, off_t)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const unsigned short *data = static_cast < const unsigned short *>(first_block);
	unsigned short nx = data[0];
	unsigned short ny = data[1];
	unsigned short nz = data[2];

	if(ByteOrder::is_data_big_endian(&nx) && ByteOrder::is_data_big_endian(&ny) && ByteOrder::is_data_big_endian(&nz)) {
		return true;
	}

	EXITFUNC;
	return false;
}
