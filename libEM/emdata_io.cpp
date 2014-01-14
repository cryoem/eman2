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

#include <iomanip>

#include "emdata.h"
#include "all_imageio.h"
#include "ctf.h"

#include <iostream>
using std::cout;
using std::endl;

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

using namespace EMAN;

void EMData::read_image(const string & filename, int img_index, bool nodata,
						const Region * region, bool is_3d)
{
	ENTERFUNC;

	ImageIO *imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY);

	if (!imageio) {
		throw ImageFormatException("cannot create an image io");
	}
	else {
		int err = imageio->read_header(attr_dict, img_index, region, is_3d);
		if (err) {
			throw ImageReadException(filename, "imageio read header failed");
		}
		else {
			LstIO * myLstIO = dynamic_cast<LstIO *>(imageio);
			if(!myLstIO)	attr_dict["source_path"] = filename;	//"source_path" is set to full path of reference image for LstIO, so skip this statement
			attr_dict["source_n"] = img_index;
			if (imageio->is_complex_mode()) {
				set_complex(true);
				set_fftpad(true);
			}
			if (attr_dict.has_key("is_fftodd") && (int)attr_dict["is_fftodd"] == 1) {
				set_fftodd(true);
			}
			if ((int) attr_dict["is_complex_ri"] == 1) {
				set_ri(true);
			}
			save_byteorder_to_dict(imageio);

			nx = attr_dict["nx"];
			ny = attr_dict["ny"];
			nz = attr_dict["nz"];
			attr_dict.erase("nx");
			attr_dict.erase("ny");
			attr_dict.erase("nz");

			if (!nodata) {

				if (region) {
					nx = (int)region->get_width();
					if (nx <= 0) nx = 1;
					ny = (int)region->get_height();
					if (ny <= 0) ny = 1;
					nz = (int)region->get_depth();
					if (nz <= 0) nz = 1;
					set_size(nx,ny,nz);
					to_zero(); // This could be avoided in favor of setting only the regions that were not read to to zero... but tedious
				} // else the dimensions of the file being read match those of this
				else {
					set_size(nx, ny, nz);
				}

				// If GPU features are enabled there is  danger that rdata will
				// not be allocated, but set_size takes care of this, so this
				// should be safe.
				int err = imageio->read_data(get_data(), img_index, region, is_3d);
				if (err) {
					throw ImageReadException(filename, "imageio read data failed");
				}
				else {
					update();
				}
			}
			else {
				if (rdata!=0) EMUtil::em_free(rdata);
				rdata=0;
			}

		}
	}
    
	EMUtil::close_imageio(filename, imageio);
	imageio = 0;
	EXITFUNC;
}

void EMData::read_binedimage(const string & filename, int img_index, int binfactor, bool fast, bool is_3d)
{
	ENTERFUNC;

	ImageIO *imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY);

	if (!imageio) {
		throw ImageFormatException("cannot create an image io");
	}
	else {
		int err = imageio->read_header(attr_dict, img_index, 0, is_3d);
		if (err) {
			throw ImageReadException(filename, "imageio read header failed");
		}
		else {
			attr_dict["source_path"] = filename;
			attr_dict["source_n"] = img_index;
			if (imageio->is_complex_mode()) {
				set_complex(true);
				set_fftpad(true);
			}
			if (attr_dict.has_key("is_fftodd") && (int)attr_dict["is_fftodd"] == 1) {
				set_fftodd(true);
			}
			if ((int) attr_dict["is_complex_ri"] == 1) {
				set_ri(true);
			}
			save_byteorder_to_dict(imageio);

			int ori_nx = nx = attr_dict["nx"];
			int ori_ny = ny = attr_dict["ny"];
			int ori_nz = nz = attr_dict["nz"];
			attr_dict.erase("nx");
			attr_dict.erase("ny");
			attr_dict.erase("nz");

			// At this point nx, ny and nz are all reduced by binfactor
			set_size(nx/binfactor, ny/binfactor, nz/binfactor);

			//here is where we read in the binned data
			EMData* tempdata = new EMData();
			size_t sizeofslice = nx*ny*sizeof(float);

			//zbin factor use 1 to speed binning(but don't benfit by averaging in Z)
			int zbin = binfactor;
			if(fast) zbin = 1;
			//verbose
			float percent = 0.1f;
			for(int k = 0; k < ori_nz; k+=binfactor){
				if(k > ori_nz*percent){
					cout << float(k)/float(ori_nz) << "% Done!" << endl;
					percent+=0.1f;
				}
				// read in a slice region
				const Region* binregion = new Region(0,0,k,ori_nx,ori_ny,zbin);
				tempdata->read_image(filename, 0, false, binregion);
				// shrink the slice
				if (binfactor > 1) tempdata->process_inplace("math.meanshrink",Dict("n",binfactor));
				size_t offset = nx*ny*k/binfactor;
				//add slice to total
				EMUtil::em_memcpy(get_data()+offset,tempdata->get_data(),sizeofslice);
				delete binregion;
			}

			delete tempdata;
			update();
		}
	}

    EMUtil::close_imageio(filename, imageio);
	imageio = 0;
	EXITFUNC;
}

#include <sys/stat.h>

void EMData::write_image(const string & filename, int img_index,
						 EMUtil::ImageType imgtype,
						 bool header_only, const Region * region,
						 EMUtil::EMDataType filestoragetype,
						 bool use_host_endian)
{
	ENTERFUNC;

	struct stat fileinfo;
	if ( region && stat(filename.c_str(),&fileinfo) != 0 ) throw UnexpectedBehaviorException("To write an image using a region the file must already exist and be the correct dimensions");

	if (is_complex() && is_shuffled())
		fft_shuffle();

	if (imgtype == EMUtil::IMAGE_UNKNOWN) {
		const char *ext = strrchr(filename.c_str(), '.');
		if (ext) {
			ext++;
			imgtype = EMUtil::get_image_ext_type(ext);
		}
	}
	ImageIO::IOMode rwmode = ImageIO::READ_WRITE;

	//set "nx", "ny", "nz" and "changecount" in attr_dict, since they are taken out of attribute dictionary

	attr_dict["nx"] = nx;
	attr_dict["ny"] = ny;
	attr_dict["nz"] = nz;
	attr_dict["changecount"] = changecount;

    // Check if this is a write only format.
    if (Util::is_file_exist(filename)) {
        if (!header_only && region == 0) {
            ImageIO * tmp_imageio = EMUtil::get_imageio(filename, ImageIO::READ_ONLY, imgtype);
            if (tmp_imageio->is_single_image_format()) {
                rwmode = ImageIO::WRITE_ONLY;
            }
            EMUtil::close_imageio(filename, tmp_imageio);
            tmp_imageio = 0;
        }
    }

	LOGVAR("getimageio %d",rwmode);
	ImageIO *imageio = EMUtil::get_imageio(filename, rwmode, imgtype);
	if (!imageio) {
		throw ImageFormatException("cannot create an image io");
	}
	else {
		update_stat();
		/* Let each image format decide how to deal with negative image_index*/
//		if (img_index < 0) {
//			img_index = imageio->get_nimg();
//		}
		LOGVAR("header write %d",img_index);

		switch(filestoragetype) {
		case EMUtil::EM_UINT:
			attr_dict["datatype"] = (int)EMUtil::EM_UINT;
			break;
		case EMUtil::EM_USHORT:
			attr_dict["datatype"] = (int)EMUtil::EM_USHORT;
			break;
		case EMUtil::EM_SHORT:
			attr_dict["datatype"] = (int)EMUtil::EM_SHORT;
			break;
		case EMUtil::EM_CHAR:
			attr_dict["datatype"] = (int)EMUtil::EM_CHAR;
			break;
		case EMUtil::EM_UCHAR:
			attr_dict["datatype"] = (int)EMUtil::EM_UCHAR;
			break;
		default:
			attr_dict["datatype"] = (int)EMUtil::EM_FLOAT;;	//default float
		}

		int err = imageio->write_header(attr_dict, img_index, region, filestoragetype,
										use_host_endian);
		if (err) {
			throw ImageWriteException(filename, "imageio write header failed");
		}
		else {
			if (!header_only) {
				if (imgtype == EMUtil::IMAGE_LST) {
					const char *reffile = attr_dict["LST.reffile"];
					if (strcmp(reffile, "") == 0) {
						reffile = path.c_str();
					}
					int refn = attr_dict["LST.refn"];
					if (refn < 0) {
						refn = pathnum;
					}

					const char *comment = attr_dict["LST.comment"];
					char *lstdata = new char[1024];
					sprintf(lstdata, "%d\t%s", refn, reffile);
					if(strcmp(comment, "") != 0) {
						sprintf(lstdata+strlen(lstdata), "\t%s\n", comment);
					}
					else {
						strcat(lstdata, "\n");
					}
					err = imageio->write_data((float*)lstdata, img_index,
											  region, filestoragetype, use_host_endian);
					if( lstdata )
					{
						delete [] lstdata;
						lstdata = 0;
					}
				}
				if (imgtype == EMUtil::IMAGE_LSTFAST) {
					const char *reffile = attr_dict["LST.reffile"];
					if (strcmp(reffile, "") == 0) {
						reffile = path.c_str();
					}
					int refn = attr_dict["LST.refn"];
					if (refn < 0) {
						refn = pathnum;
					}

					const char *comment = attr_dict["LST.comment"];
					char *lstdata = new char[1024];
					sprintf(lstdata, "%d\t%s", refn, reffile);
					if(strcmp(comment, "") != 0) {
						sprintf(lstdata+strlen(lstdata), "\t%s\n", comment);
					}
					else {
						strcat(lstdata, "\n");
					}
					err = imageio->write_data((float*)lstdata, img_index,
											  region, filestoragetype, use_host_endian);
					if( lstdata )
					{
						delete [] lstdata;
						lstdata = 0;
					}
				}
				else {
					err = imageio->write_data(get_data(), img_index, region, filestoragetype,
											  use_host_endian);
				}
				if (err) {
					imageio->flush();
					throw ImageWriteException(filename, "imageio write data failed");
				}
			}
		}
	}
	//PNG image already do cleaning in write_data function.
	if (!(imgtype == EMUtil::IMAGE_PNG)) {
		imageio->flush();
	}

    EMUtil::close_imageio(filename, imageio);
	imageio = 0;
	EXITFUNC;
}


void EMData::append_image(const string & filename,
						  EMUtil::ImageType imgtype, bool header_only)
{
	ENTERFUNC;
	write_image(filename, -1, imgtype, header_only, 0);
	EXITFUNC;
}


void EMData::write_lst(const string & filename, const string & reffile,
					   int refn, const string & comment)
{
	ENTERFUNC;
	attr_dict["LST.reffile"] = reffile;
	attr_dict["LST.refn"] = refn;
	attr_dict["LST.comment"] = comment;
	write_image(filename, -1, EMUtil::IMAGE_LST, false);
	EXITFUNC;
}


void EMData::print_image(const string str, ostream& out) {
	out << "Printing EMData object: " << str << std::endl;
	int nx = get_xsize();
	int ny = get_ysize();
	int nz = get_zsize();
	for (int iz = 0; iz < nz; iz++) {
		out << "(z = " << iz << " slice)" << std::endl;
		for (int ix = 0; ix < nx; ix++) {
			for (int iy = 0; iy < ny; iy++) {
				out << std::setiosflags(std::ios::fixed)
					<< std::setiosflags(std::ios_base::scientific)
					<< std::setw(12)
					 << std::setprecision(5) << (*this)(ix,iy,iz) << "  ";
				if (((iy+1) % 6) == 0) {
					out << std::endl << "   ";
				}
			}
			out << std::endl;
		}
	}
}

vector < shared_ptr<EMData> > EMData::read_images(const string & filename, vector < int >img_indices,
									   bool header_only)
{
	ENTERFUNC;

	int total_img = EMUtil::get_image_count(filename);
	size_t num_img = img_indices.size();

	for (size_t i = 0; i < num_img; i++) {
		if (img_indices[i] < 0 || img_indices[i] >= total_img) {
			throw OutofRangeException(0, total_img, img_indices[i], "image index");
		}
	}

	size_t n = (num_img == 0 ? total_img : num_img);

	vector< shared_ptr<EMData> > v;
	for (size_t j = 0; j < n; j++) {
		shared_ptr<EMData> d(new EMData());
		size_t k = (num_img == 0 ? j : img_indices[j]);
		try {
			d->read_image(filename, (int)k, header_only);
		}
		catch(E2Exception &e) {
			throw(e);
		}
		if ( d != 0 )
		{
			v.push_back(d);
		}
		else
			throw ImageReadException(filename, "imageio read data failed");
	}

	EXITFUNC;
	return v;
}


vector < shared_ptr<EMData> >EMData::read_images_ext(const string & filename, int img_index_start,
										   int img_index_end, bool header_only,
										   const string & ext)
{
	ENTERFUNC;

	if (img_index_end < img_index_start) {
		throw InvalidValueException(img_index_end, "image index end < image index start");
	}
	string new_filename = filename;
	new_filename = new_filename.insert(new_filename.rfind("."), ext);
	int num_img = EMUtil::get_image_count(new_filename);

	if (img_index_start < 0 || img_index_start >= num_img) {
		throw OutofRangeException(0, num_img-1, img_index_start, "image index start");
	}

	if (img_index_end >= num_img) {
		img_index_end = num_img - 1;
	}

	vector < shared_ptr<EMData> >v;

	for (int i = img_index_start; i < img_index_end; i++) {
		shared_ptr<EMData> d(new EMData());
		try {
			d->read_image(new_filename, i, header_only);
		}
		catch(E2Exception &e) {
			throw(e);
		}
		v.push_back(d);
	}
	EXITFUNC;
	return v;
}

