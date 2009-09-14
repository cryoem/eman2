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

#include "spiderio.h"
#include "geometry.h"
#include "portable_fileio.h"
#include "emassert.h"
#include "util.h"
#include "transform.h"
#include <iostream>
#include <ctime>
#include <algorithm>

using namespace EMAN;

SpiderIO::SpiderIO(const string & spider_filename, IOMode rw)
:	filename(spider_filename), rw_mode(rw),
	spider_file(0), first_h(0), cur_h(0),
	is_big_endian(ByteOrder::is_host_big_endian()),
	initialized(false)
{
	is_new_file = !Util::is_file_exist(filename);
}

SpiderIO::~SpiderIO()
{
	if (spider_file) {
		fclose(spider_file);
		spider_file = 0;
	}

	if (first_h) {
		free(first_h);
		first_h = 0;
	}

	if (cur_h) {
		free(cur_h);
		cur_h = 0;
	}
}

void SpiderIO::init()
{
	if (initialized) {
		return;
	}
	ENTERFUNC;
	spider_file = sfopen(filename, rw_mode, &is_new_file);
	initialized = true;

	if (!is_new_file) {
		first_h = static_cast < SpiderHeader * >(calloc(1, sizeof(SpiderHeader)));

		if (fread(first_h, sizeof(SpiderHeader), 1, spider_file) != 1) {
			throw ImageReadException(filename, "SPIDER header");
		}

		if (!is_valid_spider(first_h)) {
			throw ImageReadException(filename, "invalid SPIDER");
		}

		float nslice = first_h->nslice;

		is_big_endian = ByteOrder::is_float_big_endian(nslice);
		become_host_endian((float *) first_h, NUM_FLOATS_IN_HEADER);

		if (first_h->istack == SINGLE_IMAGE_HEADER && rw_mode == WRITE_ONLY) {
			fclose(spider_file);
			spider_file = 0;

			spider_file = fopen(filename.c_str(), "wb");
		}
	}

	EXITFUNC;
}

bool SpiderIO::is_valid_spider(const void *first_block)
{
	return SpiderIO::is_valid(first_block);
}

bool SpiderIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;

	if (first_block) {
		const float *data = static_cast < const float *>(first_block);
		float nslice = data[0];
		float nrow = data[1];
		float iform = data[4];
		float nsam = data[11];
		float labrec = data[12];	//NO. of records in file header
		float labbyt = data[21];	//total NO. of bytes in header
		float lenbyt = data[22];	//record length in bytes
		float istack = data[23];

		bool big_endian = ByteOrder::is_float_big_endian(nslice);
		if (big_endian != ByteOrder::is_host_big_endian()) {
			ByteOrder::swap_bytes(&nslice);
			ByteOrder::swap_bytes(&nrow);
			ByteOrder::swap_bytes(&iform);
			ByteOrder::swap_bytes(&nsam);
			ByteOrder::swap_bytes(&labrec);
			ByteOrder::swap_bytes(&labbyt);
			ByteOrder::swap_bytes(&lenbyt);
			ByteOrder::swap_bytes(&istack);
		}

		if( int(nslice) != nslice || int(nrow) != nrow
				|| int(iform) != iform || int(nsam) != nsam
				|| int(labrec) != labrec || int(labbyt) != labbyt
				|| int(lenbyt) != lenbyt ) {
			result =  false;
		}
		else {
			//now we expect this header to be an overall header for SPIDER
			if( int(istack) > 0 ) {
				result = true; //istack>0 for overall header, istack<0 for indexed stack of image
			}
		}

		int ilabrec = static_cast<int>(labrec);
		int ilabbyt = static_cast<int>(labbyt);
		int ilenbyt = static_cast<int>(lenbyt);
		if( ilabbyt != ilabrec * ilenbyt ) {
			result = false;
		}
	}

	EXITFUNC;
	return result;
}

/** If image_index = -1, read the overall spider header.
 */
int SpiderIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

    bool is_read_overall_header = false;
    if (image_index == -1) {
        is_read_overall_header = true;
        image_index = 0;
    }

	check_read_access(image_index);

	if (!first_h) {
		throw ImageReadException(filename, "empty spider header");
	}

	check_region(area, FloatSize(first_h->nsam, first_h->nrow,first_h->nslice), is_new_file,false);

	float size = first_h->nsam * first_h->nrow * first_h->nslice;
	int overall_headlen = 0;

    if (!is_read_overall_header) {
        if (first_h->istack > 0) {	//stack image
            overall_headlen = (int) first_h->headlen;
        }
        else if(first_h->istack == SINGLE_IMAGE_HEADER) {	//single 2D/3D image
            if(image_index != 0) {
                char desc[1024];
                sprintf(desc, "For a single image, index must be 0. Your image index = %d.", image_index);
                throw ImageReadException(filename, desc);
            }
        }
        else {	//complex spider image not supported
        	throw ImageFormatException("complex spider image not supported.");
        }
    }

	size_t single_image_size = (size_t) (first_h->headlen + size * sizeof(float));
	size_t offset = overall_headlen + single_image_size * image_index;

	SpiderHeader *cur_image_hed;
	if (offset == 0) {
		cur_image_hed = first_h;
	}
	else {
		cur_image_hed = static_cast < SpiderHeader * >(calloc(1, sizeof(SpiderHeader)));
		portable_fseek(spider_file, offset, SEEK_SET);

		if (fread(cur_image_hed, sizeof(SpiderHeader), 1, spider_file) != 1) {
			char desc[1024];
			sprintf(desc, "read spider header with image_index = %d failed", image_index);
			throw ImageReadException(filename, desc);
		}

		become_host_endian((float *) cur_image_hed, NUM_FLOATS_IN_HEADER);

		if (cur_image_hed->nsam != first_h->nsam || cur_image_hed->nrow != first_h->nrow
			|| cur_image_hed->nslice != first_h->nslice) {
			char desc[1024];
			sprintf(desc, "%dth image size %dx%dx%d != overall size %dx%dx%d",
					image_index, (int)cur_image_hed->nsam, (int)cur_image_hed->nrow,
					(int)cur_image_hed->nslice,
					(int)first_h->nsam, (int)first_h->nrow, (int)first_h->nslice);
			throw ImageReadException(filename, desc);
		}
	}


	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, (int) cur_image_hed->nsam, &xlen, (int) cur_image_hed->nrow,
							&ylen, (int) cur_image_hed->nslice, &zlen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = zlen;

	dict["datatype"] = EMUtil::EM_FLOAT;

	if(cur_image_hed->mmvalid == 1) {
		dict["minimum"] = cur_image_hed->min;
		dict["maximum"] = cur_image_hed->max;
		dict["mean"] = cur_image_hed->mean;
		dict["sigma"] = cur_image_hed->sigma;
	}

	dict["SPIDER.nslice"] = (int) cur_image_hed->nslice;
	dict["SPIDER.type"] = (int) cur_image_hed->type;

	dict["SPIDER.irec"] = cur_image_hed->irec;

	dict["SPIDER.angvalid"] = (int)cur_image_hed->angvalid;
	if((int)dict["SPIDER.angvalid"] != 0) {
		dict["SPIDER.phi"] = cur_image_hed->phi;
		dict["SPIDER.theta"] = cur_image_hed->theta;
		dict["SPIDER.gamma"] = cur_image_hed->gamma;
	}

	dict["SPIDER.headrec"] = (int) cur_image_hed->headrec;
	dict["SPIDER.headlen"] = (int) cur_image_hed->headlen;
	dict["SPIDER.reclen"] = (int) cur_image_hed->reclen;

	dict["SPIDER.dx"] = cur_image_hed->dx;
	dict["SPIDER.dy"] = cur_image_hed->dy;
	dict["SPIDER.dz"] = cur_image_hed->dz;

	dict["SPIDER.istack"] = (int) cur_image_hed->istack;
	if((int)dict["SPIDER.istack"] > 0) {	//maxim only for overall header
		dict["SPIDER.maxim"] = (int)cur_image_hed->maxim;
	}
	dict["SPIDER.imgnum"] = (int)cur_image_hed->imgnum;

	dict["SPIDER.Kangle"] = (int)cur_image_hed->Kangle;
	if((int)dict["SPIDER.Kangle"] == 1) {
		dict["SPIDER.phi1"] = cur_image_hed->phi1;
		dict["SPIDER.theta1"] = cur_image_hed->theta1;
		dict["SPIDER.psi1"] = cur_image_hed->psi1;
	}
	else if((int)dict["SPIDER.Kangle"] == 2) {
		dict["SPIDER.phi1"] = cur_image_hed->phi1;
		dict["SPIDER.theta1"] = cur_image_hed->theta1;
		dict["SPIDER.psi1"] = cur_image_hed->psi1;
		dict["SPIDER.phi2"] = cur_image_hed->phi2;
		dict["SPIDER.theta2"] = cur_image_hed->theta2;
		dict["SPIDER.psi2"] = cur_image_hed->psi2;
	}

	dict["SPIDER.date"] = string(cur_image_hed->date).substr(0, 11);
	dict["SPIDER.time"] = string(cur_image_hed->time).substr(0, 8);

	dict["SPIDER.title"] = string(cur_image_hed->title);

	if(cur_image_hed->scale>0) {
		dict["SPIDER.scale"] = cur_image_hed->scale;
		Dict dic;
		dic.put("type", "spider");
		dic.put("phi", cur_image_hed->phi);
		dic.put("theta", cur_image_hed->theta);
		dic.put("psi", cur_image_hed->gamma);
		dic.put("tx", cur_image_hed->dx);
		dic.put("ty", cur_image_hed->dy);
		dic.put("tz", cur_image_hed->dz);
		dic.put("scale", cur_image_hed->scale);
		Transform * trans = new Transform(dic);
		if(zlen<=1) {
			dict["xform.projection"] = trans;
		}
		else {
			dict["xform.align3d"] = trans;
		}
		if(trans) {delete trans; trans=0;}
	}


/**	No complex image support for SPIDER in EMAN2
 * 	// complex check is in is_complex_mode, but even/odd check needs to be here
	int type = dict["SPIDER.type"];
	if ( type == IMAGE_2D_FFT_ODD
			|| type == IMAGE_3D_FFT_ODD )
		dict["is_fftodd"] = 1;
	else dict["is_fftodd"] = 0;
	// complex spider files are always ri, so we just set it unconditionally
	dict["is_complex_ri"] = 1;
*/
	if (offset != 0) {
		if( cur_image_hed )
		{
			free(cur_image_hed);
			cur_image_hed = 0;
		}
	}
	EXITFUNC;
	return 0;
}


int SpiderIO::write_header(const Dict & dict, int image_index, const Region* area,
						   EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	if(is_new_file) {	//for a new file write overall header first
		write_single_header(dict, area, image_index, 0, first_h, OVERALL_STACK_HEADER, 1, 0, use_host_endian);
	}
	else {
		swap_header(first_h);
	}

	if(!initialized) {
		init();
	}

	if(!((int)dict["nx"]==first_h->nsam && (int)dict["ny"]==first_h->nrow &&
			(int)dict["nz"]==first_h->nslice)) {
		char desc[1024];
		sprintf(desc, "%dth image size %dx%dx%d != overall size %dx%dx%d",
				image_index, (int)dict["nx"], (int)dict["ny"], (int)dict["nz"],
				(int)first_h->nsam, (int)first_h->nrow, (int)first_h->nslice);
		throw ImageReadException(filename, desc);
	}

	if (!cur_h) {
		cur_h = (SpiderHeader *) calloc(1, static_cast<size_t>(first_h->headlen));
	}

	int MAXIM;
	float size = first_h->nsam * first_h->nrow * first_h->nslice;
	size_t single_image_size = (int) (first_h->headlen + size * sizeof(float));
	size_t offset;
	if(image_index == -1) {	//append image
		offset = (int) first_h->headlen + single_image_size * (int)first_h->maxim;
		MAXIM = (int)first_h->maxim + 1;
	}
	else {
		offset = (int) first_h->headlen + single_image_size * image_index;
		MAXIM = image_index>=(int)first_h->maxim ? image_index+1 : (int)first_h->maxim;
	}

	//update overall header
	if(MAXIM > (int)first_h->maxim) {
		portable_fseek(spider_file, 0, SEEK_SET);
		write_single_header(dict, area, image_index, 0, first_h, OVERALL_STACK_HEADER, MAXIM, 0, use_host_endian);
	}

	portable_fseek(spider_file, offset, SEEK_SET);
	write_single_header(dict, area, image_index, offset, cur_h, SINGLE_IMAGE_HEADER, 0, image_index+1, use_host_endian);

	EXITFUNC;
	return 0;
}

int SpiderIO::write_single_header(const Dict & dict, const Region *area, int image_index, size_t offset,
					SpiderHeader *& hp, int ISTACK, int MAXIM, int IMGNUM, bool use_host_endian)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index);

	if (area) {
		check_region(area, FloatSize(hp->nsam, hp->nrow, hp->nslice), is_new_file);
		EXITFUNC;
		return 0;
	}

	int nx = dict["nx"];
	int ny = dict["ny"];
	int nz = dict["nz"];

	size_t header_size = sizeof(SpiderHeader);
	size_t record_size = nx * sizeof(float);
    size_t num_records = (header_size % record_size) == 0 ? header_size / record_size : header_size / record_size + 1;
	size_t header_length = num_records * record_size;

	if (!hp) {
		hp = static_cast < SpiderHeader * >(calloc(1, header_size));
	}

	hp->angvalid = 0;
	hp->scale = 1.0;
	hp->istack = (float)ISTACK;
	hp->nslice = (float)nz;
	hp->nsam = (float)nx;
	hp->nrow = (float)ny;

	hp->max = dict["maximum"];
	hp->min = dict["minimum"];
	hp->mean = dict["mean"];
	hp->sigma = dict["sigma"];
	hp->mmvalid = 1;

	if(nz<=1 && dict.has_key("xform.projection")) {
		hp->angvalid = 1;
		Transform * t = dict["xform.projection"];
		Dict d = t->get_params("spider");
		hp->phi = d["phi"];
		hp->theta = d["theta"];
		hp->gamma = d["psi"];
		hp->dx = d["tx"];
		hp->dy = d["ty"];
		hp->dz = d["tz"];
		hp->scale = d["scale"];
		if(t) {delete t; t=0;}
	}
	else if(nz>1 && dict.has_key("xform.align3d")) {
		hp->angvalid = 1;
		Transform * t = dict["xform.align3d"];
		Dict d = t->get_params("spider");
		hp->phi = d["phi"];
		hp->theta = d["theta"];
		hp->gamma = d["psi"];
		hp->dx = d["tx"];
		hp->dy = d["ty"];
		hp->dz = d["tz"];
		hp->scale = d["scale"];
		if(t) {delete t; t=0;}
	}

	if(nz == 1) {
		hp->type = IMAGE_2D;
	}
	else {
		hp->type = IMAGE_3D;
	}

// complex image file not supported in EMAN2

	hp->reclen = (float)record_size;
	hp->headrec = (float)num_records;
	hp->headlen = (float)header_length;

	if(ISTACK == OVERALL_STACK_HEADER) {
		hp->maxim = (float)MAXIM;
	}
	hp->irec = (float)(num_records + ny*nz);

	hp->imgnum = (float)IMGNUM;

	time_t tod;
 	time(&tod);
 	struct tm * ttt = localtime(&tod);
 	char ctime[9];
 	char cdate[12];
 	strftime(ctime, 9, "%H:%M:%S", ttt);
 	std::copy(&ctime[0], &ctime[8], hp->time);
 	strftime(cdate, 12, "%d-%b-%Y", ttt);
	std::copy(&cdate[0], &cdate[11], hp->date);

	if(dict.has_key("SPIDER.title")) {
		string title = static_cast<string>(dict["SPIDER.title"]);
		std::copy(&(title[0]), &(title[title.length()]), hp->title);
	}

	portable_fseek(spider_file, offset, SEEK_SET);

	if(use_host_endian) {
		if (fwrite(hp, header_size, 1, spider_file) != 1) {
			throw ImageWriteException(filename, "write spider header failed");
		}
	}
	else {	//swap byte order
		SpiderHeader * hp2 = new SpiderHeader(*hp);
		ByteOrder::swap_bytes((float *) hp2, NUM_FLOATS_IN_HEADER);
		if (fwrite(hp2, header_size, 1, spider_file) != 1) {
			throw ImageWriteException(filename, "write spider header failed");
		}
		if(hp2) {delete hp2; hp2=0;}
	}

	size_t pad_size = header_length - header_size;
	char *pad = static_cast < char *>(calloc(pad_size, 1));
	fwrite(pad, pad_size, 1, spider_file);
	if( pad )
	{
		free(pad);
		pad = 0;
	}

	EXITFUNC;
	return 0;
}


int SpiderIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index, data);

	check_region(area, FloatSize((int) first_h->nsam, (int) first_h->nrow,
								 (int) first_h->nslice), is_new_file,false);

	int overall_headlen = 0;
	if (first_h->istack >0) {	//for over all header length
		overall_headlen = (int) first_h->headlen;
	}
	else {
		if(image_index != 0) {
			char desc[256];
			sprintf(desc, "For single image, index must be 0. Your image index = %d.", image_index);
			throw ImageReadException(filename, desc);
		}
	}

	size_t size = static_cast < size_t > (first_h->nsam * first_h->nrow * first_h->nslice);
	size_t single_image_size = static_cast < size_t > (first_h->headlen + size * sizeof(float));
	off_t offset = overall_headlen + single_image_size * image_index;
	portable_fseek(spider_file, offset, SEEK_SET);

	portable_fseek(spider_file, (int) first_h->headlen, SEEK_CUR);

#if 1
	EMUtil::process_region_io(data, spider_file, READ_ONLY, 0, sizeof(float),
							  (int) first_h->nsam, (int) first_h->nrow,
							  (int) first_h->nslice, area);
#endif
#if 0
	unsigned int nz = static_cast < unsigned int >(first_h->nslice);
	int sec_size = static_cast < int >(first_h->nsam * first_h->nrow * sizeof(float));

	if (fread(data, sec_size, nz, spider_file) != nz) {
		LOGERR("Incomplete SPIDER data read");
		return 1;
	}
#endif

	int xlen = 0, ylen = 0, zlen = 0;
	EMUtil::get_region_dims(area, (int) first_h->nsam, &xlen, (int) first_h->nrow, &ylen,
							(int) first_h->nslice, &zlen);

	int data_size = xlen * ylen * zlen;
	become_host_endian(data, data_size);
	EXITFUNC;
	return 0;
}

int SpiderIO::write_data(float *data, int image_index, const Region* area,
						 EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	if(!cur_h) {
		throw ImageWriteException(filename, "Please write header before write data");
	}

	if(first_h->istack == SINGLE_IMAGE_HEADER) {
		throw ImageWriteException(filename, "Cannot mix single spider and stack spider");
	}

	size_t offset;
	float size = first_h->nsam * first_h->nrow * first_h->nslice;
	size_t single_image_size = (int) (first_h->headlen + size * sizeof(float));
	offset = (int) first_h->headlen + single_image_size * image_index + (int) first_h->headlen;

	swap_data(data, (size_t)size);

	write_single_data(data, area, cur_h, offset, image_index, (int)first_h->maxim+1, use_host_endian);

	EXITFUNC;
	return 0;
}


void SpiderIO::flush()
{
	fflush(spider_file);
}

int SpiderIO::write_single_data(float *data, const Region * area, SpiderHeader *& hp,
								size_t offset, int img_index, int max_nimg, bool use_host_endian)
{
	ENTERFUNC;

	check_write_access(rw_mode, img_index, max_nimg, data);

	if (area) {
		check_region(area, FloatSize(hp->nsam, hp->nrow, hp->nslice), is_new_file);
	}

	if (!hp) {
		throw ImageWriteException(filename, "NULL image header");
	}

	portable_fseek(spider_file, offset, SEEK_SET);

	int size = (int)(hp->nsam * hp->nrow * hp->nslice);
	if(!use_host_endian) {
		ByteOrder::swap_bytes(data, size);
	}

	// remove the stuff in #if 0 ... #endif if the following works.
	EMUtil::process_region_io(data, spider_file, WRITE_ONLY,0, sizeof(float),
							  (int) hp->nsam, (int) hp->nrow,
							  (int) hp->nslice, area);

	EXITFUNC;
	return 0;
}

void SpiderIO::swap_data(float *data, size_t size)
{
	if (data && need_swap()) {
		ByteOrder::swap_bytes(data, size);
	}
}

void SpiderIO::swap_header(SpiderHeader * header)
{
	if (header && need_swap()) {
		ByteOrder::swap_bytes((float *) header, NUM_FLOATS_IN_HEADER);
	}
}

bool SpiderIO::is_complex_mode()
{
	int type = static_cast<int>(first_h->type);
	if (type == IMAGE_2D_FFT_ODD
		|| type == IMAGE_2D_FFT_EVEN
		|| type == IMAGE_3D_FFT_ODD
		|| type == IMAGE_3D_FFT_EVEN
	   ) return true;
	return false;
}

bool SpiderIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int SpiderIO::get_nimg()
{
	init();
	if (!first_h) {
		Assert(is_new_file);
		return 0;
	}
	else if (first_h->istack > 0) {						//image stack
		return static_cast < int >(first_h->maxim);
	}
	else if (first_h->istack == SINGLE_IMAGE_HEADER) {	//single 2D/3D image
		return 1;
	}
	else {												//complex image
		throw ImageFormatException("complex spider image not supported.");
	}
}

bool SpiderIO::need_swap() const
{
	if (!is_new_file && (is_big_endian != ByteOrder::is_host_big_endian())) {
		return true;
	}
	return false;
}
/* vim: set ts=4 noet: */
