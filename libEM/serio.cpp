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

#include "serio.h"
#include "util.h"
#include "portable_fileio.h"

#include <cassert>

using namespace EMAN;

static const short SER_BYTE_ORDER		= 0x4949;
static const short SER_SERIES_ID 		= 0x0197;
//static const short SER_SERIES_VERSION 	= 0x0210;		// can also be 0x220 now
static const int   ValidNumberElementsOffset = 18;	//the offset to ValidNumberElements, which is the number of images in the file

SerIO::SerIO(const string & file, IOMode rw) :
		filename(file), rw_mode(rw), serfile(0), initialized(false),
		is_new_file(false), data_offset_array(0),tag_offset_array(0),
		nimg(0), nx(0), ny(0), nz(0), datatypeid(0), datamode(0)
{
}

SerIO::~SerIO()
{
	if (serfile) {
		fclose(serfile);
		serfile = 0;
	}

	if (data_offset_array) {
		delete [] data_offset_array;
		data_offset_array = 0;
	}

	if (tag_offset_array) {
		delete [] tag_offset_array;
		tag_offset_array = 0;
	}
}

void SerIO::init()
{
	ENTERFUNC;

	if (initialized) {
		return;
	}

	initialized = true;
	serfile = sfopen(filename, rw_mode, &is_new_file);

	if (!is_new_file) {
		if (fread(&serh, sizeof(SerHeader), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}

		if (!is_valid(&serh)) {
			throw ImageReadException(filename, "invalid SER");
		}
	}

	EXITFUNC;
}

bool SerIO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const short *data = static_cast < const short *>(first_block);
	short ByteOrder = data[0];
	short SeriesID = data[1];

	if(ByteOrder != 0x4949 || SeriesID != 0x0197) {
		return false;
	}

	EXITFUNC;
	return true;
}

int SerIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	ENTERFUNC;
	init();
	rewind(serfile);

	short hitem1[3];
	if (fread(hitem1, sizeof(short), 3, serfile) != 3) {
		throw ImageReadException(filename, "SER header");
	}

	dict["SER.ByteOrder"] 		= hitem1[0];
	dict["SER.SeriesID"]		= hitem1[1];
	dict["SER.SeriesVersion"] 	= hitem1[2];

	if((hitem1[0]!=SER_BYTE_ORDER) || (hitem1[1]!=SER_SERIES_ID) || (hitem1[2]!=0x210 && hitem1[2]!=0x220) ) {
		throw ImageReadException(filename, "SER header");
	}
	int hitem2[6];
	uint64_t offset_array;
	if (hitem1[2]==0x210) {
		if (fread(hitem2, sizeof(int), 6, serfile) != 6) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.OffsetArrayOffset"]	= hitem2[4];
		offset_array=(uint64_t)hitem2[4];
	}
	else {
		if (fread(hitem2, sizeof(int), 4, serfile) != 4) {
			throw ImageReadException(filename, "SER header");
		}
		fread(&offset_array, 8,1,serfile);
		fread(&hitem2[5],4,1,serfile);
	}
		
	
	
	dict["SER.DataTypeID"]		= hitem2[0];
	dict["SER.TagTypeID"]		= hitem2[1];
	dict["SER.TotalNumberElements"] = hitem2[2];
	dict["SER.ValidNumberElements"] = hitem2[3];
	dict["SER.NumberDimensions"]	= hitem2[5];

	nimg = (int)dict["SER.ValidNumberElements"];

	if(image_index >= (int)dict["SER.ValidNumberElements"]) {
		throw ImageReadException(filename, "Image index out of bound");
	}

	for(int idx=0; idx<(int)dict["SER.NumberDimensions"]; idx++) {
		read_dim_arr(dict, idx);
	}

// 	long pos = ftell(serfile);
// 	assert(pos == (int)dict["SER.OffsetArrayOffset"]);
	portable_fseek(serfile,(off_t)offset_array,SEEK_SET);

	int tot = (int)dict["SER.TotalNumberElements"];

	data_offset_array = new uint64_t[tot];
	tag_offset_array = new uint64_t[tot];
	if (hitem1[2]==0x210) {
		int *off=new int[tot*2];
		if (fread(off, 4, tot*2, serfile) != (unsigned int)tot*2) {
			throw ImageReadException(filename, "SER header");
		}
		for (int i=0; i<tot; i++) {
			data_offset_array[i]=(uint64_t)off[i];
			tag_offset_array[i]=(uint64_t)off[i+tot];
		}
		delete[] off;
	}
	else {
		if (fread(data_offset_array, sizeof(uint64_t), tot, serfile) != (unsigned int)tot) {
			throw ImageReadException(filename, "SER header");
		}
		if (fread(tag_offset_array, sizeof(uint64_t), tot, serfile) != (unsigned int)tot) {
			throw ImageReadException(filename, "SER header");
		}
	}

	this->datatypeid = (int)dict["SER.DataTypeID"];

	off_t dataoffset = data_offset_array[image_index];
	portable_fseek(serfile, dataoffset, SEEK_SET);

	//To read the attribute in data element(not the actual data)
	read_data_element(dict);

	off_t tagoffset = tag_offset_array[image_index];
	portable_fseek(serfile, tagoffset, SEEK_SET);

	//To read the data tag appended after data values
	read_data_tag(dict);

	EXITFUNC;
	return 0;
}

int SerIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType filestoragetype, bool use_host_endian)
{
	ENTERFUNC;

	throw ImageWriteException(filename, "SER Writing not supported");

	EXITFUNC;
	return 0;
}

int SerIO::read_data(float *rdata, int image_index, const Region *, bool )
{
	ENTERFUNC;

	if(!data_offset_array) {
		throw ImageReadException(filename, "SER header, empty DataOffsetarray");
	}

	size_t size = (size_t)nx * ny * nz;
	off_t data_offset = data_offset_array[image_index];

	size_t i;	//loop index
	unsigned char * puchar = 0;
	char * pchar = 0;
	int * pint = 0;
	unsigned int * puint = 0;
	short * pshort = 0;
	unsigned short * pushort = 0;
	double * pdouble = 0;
	switch(this->datatypeid) {
	case oneD:
		portable_fseek(serfile, data_offset+26, SEEK_SET);	//offset 26 to actual data values
		break;
	case twoD:
		portable_fseek(serfile, data_offset+50, SEEK_SET);	//offset 50 to actual data values
		break;
	default:
		throw ImageReadException(filename, "SER header, wrong DataTypeID");
	}

	switch(this->datamode) {
	case SER_UCHAR:
		puchar = new unsigned char[size];
		if (fread(puchar, sizeof(unsigned char), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(puchar[i]);
		}
		delete [] puchar;
		break;
	case SER_USHORT:
		pushort = new unsigned short[size];
		if (fread(pushort, sizeof(unsigned short), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(pushort[i]);
		}
		delete [] pushort;
		break;
	case SER_UINT:
		puint = new unsigned int[size];
		if (fread(puint, sizeof(unsigned int), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(puint[i]);
		}
		delete [] puint;
		break;
	case SER_CHAR:
		pchar = new char[size];
		if (fread(pchar, sizeof(unsigned char), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(pchar[i]);
		}
		delete [] pchar;
		break;
	case SER_SHORT:
		pshort = new short[size];
		if (fread(pshort, sizeof(short), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(pshort[i]);
		}
		delete [] pshort;
		break;
	case SER_INT:
		pint = new int[size];
		if (fread(pint, sizeof(int), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(pint[i]);
		}
		delete [] pint;
		break;
	case SER_FLOAT:
		if (fread(rdata, sizeof(float), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		break;
	case SER_DOUBLE:
		pdouble = new double[size];
		if (fread(pdouble, sizeof(double), size, serfile) != size) {
			throw ImageReadException(filename, "SER data");
		}
		for (i = 0; i<size; ++i) {
			rdata[i] = static_cast<float>(pdouble[i]);
		}
		delete [] pdouble;
		break;
	case SER_COMPLEX8:
	case SER_COMPLEX16:
		throw ImageReadException(filename, "complex data not supported.");
		break;
	default:
		throw ImageReadException(filename, "Unknown data value type");
	}

	EXITFUNC;
	return 0;
}

int SerIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	throw ImageWriteException(filename, "SER Writing not supported");

	EXITFUNC;
	return 0;
}

bool SerIO::is_complex_mode()
{
	return false;	//not support complex ser for now
}

void SerIO::flush()
{
	fflush(serfile);
}

bool SerIO::is_image_big_endian()
{
	return false;	//ser image is always little endian
}

int SerIO::get_nimg() {
	init();

	portable_fseek(serfile, ValidNumberElementsOffset, SEEK_SET);
	int nimg;
	if (fread(&nimg, sizeof(int), 1, serfile) != 1) {
		throw ImageReadException(filename, "SER header");
	}

	return nimg;
}

void SerIO::read_dim_arr(Dict & dict, int idx)
{
	int dimsize;
	if (fread(&dimsize, sizeof(int), 1, serfile) != 1) {
		throw ImageReadException(filename, "SER header");
	}

	string sidx = Util::int2str(idx);
	dict["SER.DimensionSize"+sidx]	= dimsize;

	double hitem3[2];
	if (fread(hitem3, sizeof(double), 2, serfile) != 2) {
		throw ImageReadException(filename, "SER header");
	}
	dict["SER.CalibrationOffset"+sidx]	= hitem3[0];
	dict["SER.CalibrationDelta"+sidx]	= hitem3[1];

	int celement;
	if (fread(&celement, sizeof(int), 1, serfile) != 1) {
		throw ImageReadException(filename, "SER header");
	}
	dict["SER.CalibrationElement"+sidx]	= celement;

	int desclen;
	if (fread(&desclen, sizeof(int), 1, serfile) != 1) {
		throw ImageReadException(filename, "SER header");
	}
	dict["SER.DescriptionLength"+sidx]	= desclen;

	if(desclen != 0) {
		char * descr = new char[desclen+1];
		//char descr[desclen+1];
		if (fread(descr, sizeof(char), desclen, serfile) != (unsigned int)desclen) {
			throw ImageReadException(filename, "SER header");
		}
		descr[desclen] = '\0';
		string sdescr(descr);
		dict["SER.Description"+sidx] = sdescr;
		delete [] descr;
	}

	int unitslen;
	if (fread(&unitslen, sizeof(int), 1, serfile) != 1) {
		throw ImageReadException(filename, "SER header");
	}
	dict["SER.UnitsLength"+sidx] = unitslen;

	if(unitslen != 0) {
		char * units = new char[unitslen+1];
		//char units[unitslen+1];
		if (fread(units, sizeof(int), unitslen, serfile) != (unsigned int)unitslen) {
			throw ImageReadException(filename, "SER header");
		}
		units[unitslen] = '\0';
		string sunits(units);
		dict["SER.Units"+sidx] = sunits;
		delete [] units;
	}

}

void SerIO::read_data_element(Dict & dict)
{
	if(this->datatypeid == oneD) {	//1D image
		double hitem4[2];
		if (fread(hitem4, sizeof(double), 2, serfile) != 2) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.CalibrationOffset"] 	= hitem4[0];
		dict["SER.CalibrationDelta"]	= hitem4[1];

		int cali;
		if (fread(&cali, sizeof(int), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.CalibrationElement"] = cali;

		short datatype;
		if (fread(&datatype, sizeof(short), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.DataType"] = datatype;

		int arrlen;
		if (fread(&arrlen, sizeof(int), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["nx"] = arrlen;
		dict["ny"] = 1;
		dict["nz"] = 1;

		nx = arrlen;
		ny = 1;
		nz = 1;
	}
	else if(this->datatypeid == twoD) {	//2D image
		double hitem4[2];
		if (fread(hitem4, sizeof(double), 2, serfile) != 2) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.CalibrationOffsetX"] 	= hitem4[0];
		dict["SER.CalibrationDeltaX"]	= hitem4[1];

		int calix;
		if (fread(&calix, sizeof(int), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.CalibrationElementX"] = calix;

		double hitem5[2];
		if (fread(hitem5, sizeof(double), 2, serfile) != 2) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.CalibrationOffsetX"] 	= hitem5[0];
		dict["SER.CalibrationDeltaX"]	= hitem5[1];

		int caliy;
		if (fread(&caliy, sizeof(int), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.CalibrationElementY"] = caliy;

		short datatype;
		if (fread(&datatype, sizeof(short), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.DataType"] = datatype;
		this->datamode = datatype;

		int arrsize[2];
		if (fread(&arrsize, sizeof(int), 2, serfile) != 2) {
			throw ImageReadException(filename, "SER header");
		}
		dict["nx"] = arrsize[0];
		dict["ny"] = arrsize[1];
		dict["nz"] = 1;

		nx = arrsize[0];
		ny = arrsize[1];
		nz = 1;
	}
}

void SerIO::read_data_tag(Dict & dict)
{
	int tag_type = (int)dict["SER.TagTypeID"];
	if( tag_type == timeOnly ) {
		short tagtype;
		if (fread(&tagtype, sizeof(short), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		assert((int)tagtype == tag_type);

		int sertime;
		if (fread(&sertime, sizeof(int), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.Time"] = sertime;
	}
	else if( tag_type == posTime ) {
		short tagtype;
		if (fread(&tagtype, sizeof(short), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		assert((int)tagtype == tag_type);

		int sertime;
		if (fread(&sertime, sizeof(int), 1, serfile) != 1) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.Time"] = sertime;

		double pos[2];
		if (fread(&pos, sizeof(double), 2, serfile) != 2) {
			throw ImageReadException(filename, "SER header");
		}
		dict["SER.PosionX"] = pos[0];
		dict["SER.PosionY"] = pos[1];
	}
	else {
		throw ImageReadException(filename, "SER header, wrong TagTypeID");
	}
}
