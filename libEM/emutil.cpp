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

#ifdef WIN32
#include <windows.h>
	#define MAXPATHLEN (MAX_PATH*4)
#else
#include <sys/param.h>
#endif	// WIN32

#include "all_imageio.h"
#include "portable_fileio.h"
#include "emcache.h"
#include "emdata.h"
#include "ctf.h"
#include "emassert.h"
#include "exception.h"
#include "hdf_filecache.h"

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

//#ifdef EMAN2_USING_CUDA_MALLOC
//#include "cuda/cuda_util.h"
//#endif

using namespace EMAN;

static const int ATTR_NAME_LEN = 128;

EMUtil::ImageType EMUtil::get_image_ext_type(const string & file_ext)
{
	ENTERFUNC;

	static bool initialized = false;
	static map < string, ImageType > imagetypes;

	if (!initialized) {
		imagetypes["rec"] = IMAGE_MRC;
		imagetypes["mrc"] = IMAGE_MRC;
		imagetypes["MRC"] = IMAGE_MRC;
		imagetypes["ali"] = IMAGE_MRC;
		imagetypes["st"] = IMAGE_MRC;		// IMOD stack file

		imagetypes["mrcs"] = IMAGE_MRC;
		imagetypes["MRCS"] = IMAGE_MRC;

		imagetypes["raw"] = IMAGE_MRC;
		imagetypes["RAW"] = IMAGE_MRC;

		imagetypes["tnf"] = IMAGE_MRC;
		imagetypes["TNF"] = IMAGE_MRC;

		imagetypes["ccp4"] = IMAGE_MRC;
		imagetypes["map"] = IMAGE_MRC;

		imagetypes["dm3"] = IMAGE_DM3;
		imagetypes["DM3"] = IMAGE_DM3;

		imagetypes["dm4"] = IMAGE_DM4;
		imagetypes["DM4"] = IMAGE_DM4;

		imagetypes["dat"] = IMAGE_SPIDER;
		imagetypes["spi"] = IMAGE_SPIDER;
		imagetypes["SPI"] = IMAGE_SPIDER;
		imagetypes["dat"] = IMAGE_SPIDER;
		imagetypes["DAT"] = IMAGE_SPIDER;

		imagetypes["spider"] = IMAGE_SPIDER;
		imagetypes["SPIDER"] = IMAGE_SPIDER;

		imagetypes["spidersingle"] = IMAGE_SINGLE_SPIDER;
		imagetypes["SPIDERSINGLE"] = IMAGE_SINGLE_SPIDER;

		imagetypes["singlespider"] = IMAGE_SINGLE_SPIDER;
		imagetypes["SINGLESPIDER"] = IMAGE_SINGLE_SPIDER;

		imagetypes["img"] = IMAGE_IMAGIC;
		imagetypes["IMG"] = IMAGE_IMAGIC;

		imagetypes["hed"] = IMAGE_IMAGIC;
		imagetypes["HED"] = IMAGE_IMAGIC;

		imagetypes["imagic"] = IMAGE_IMAGIC;
		imagetypes["IMAGIC"] = IMAGE_IMAGIC;

		imagetypes["pgm"] = IMAGE_PGM;
		imagetypes["PGM"] = IMAGE_PGM;

		imagetypes["lst"] = IMAGE_LST;
		imagetypes["LST"] = IMAGE_LST;

		imagetypes["lsx"] = IMAGE_LSTFAST;	// but .lst or another extension would also be ok
		imagetypes["LSX"] = IMAGE_LSTFAST;

		imagetypes["pif"] = IMAGE_PIF;
		imagetypes["PIF"] = IMAGE_PIF;

		imagetypes["png"] = IMAGE_PNG;
		imagetypes["PNG"] = IMAGE_PNG;

		imagetypes["h5"] = IMAGE_HDF;
		imagetypes["H5"] = IMAGE_HDF;

		imagetypes["hd5"] = IMAGE_HDF;
		imagetypes["HD5"] = IMAGE_HDF;

		imagetypes["hdf"] = IMAGE_HDF;
		imagetypes["HDF"] = IMAGE_HDF;

		imagetypes["tif"] = IMAGE_TIFF;
		imagetypes["TIF"] = IMAGE_TIFF;

		imagetypes["tiff"] = IMAGE_TIFF;
		imagetypes["TIFF"] = IMAGE_TIFF;

		imagetypes["fts"] = IMAGE_FITS;
		imagetypes["FTS"] = IMAGE_FITS;

		imagetypes["vtk"] = IMAGE_VTK;
		imagetypes["VTK"] = IMAGE_VTK;

		imagetypes["hdr"] = IMAGE_SAL;
		imagetypes["HDR"] = IMAGE_SAL;

		imagetypes["sal"] = IMAGE_SAL;
		imagetypes["SAL"] = IMAGE_SAL;

		imagetypes["map"] = IMAGE_ICOS;
		imagetypes["MAP"] = IMAGE_ICOS;

		imagetypes["icos"] = IMAGE_ICOS;
		imagetypes["ICOS"] = IMAGE_ICOS;

		imagetypes["am"] = IMAGE_AMIRA;
		imagetypes["AM"] = IMAGE_AMIRA;

		imagetypes["amira"] = IMAGE_AMIRA;
		imagetypes["AMIRA"] = IMAGE_AMIRA;

		imagetypes["emim"] = IMAGE_EMIM;
		imagetypes["EMIM"] = IMAGE_EMIM;

		imagetypes["xplor"] = IMAGE_XPLOR;
		imagetypes["XPLOR"] = IMAGE_XPLOR;

		imagetypes["em"] = IMAGE_EM;
		imagetypes["EM"] = IMAGE_EM;

		imagetypes["dm2"] = IMAGE_GATAN2;
		imagetypes["DM2"] = IMAGE_GATAN2;

		imagetypes["v4l"] = IMAGE_V4L;
		imagetypes["V4L"] = IMAGE_V4L;

		imagetypes["jpg"] = IMAGE_JPEG;
		imagetypes["JPG"] = IMAGE_JPEG;
		imagetypes["jpeg"] = IMAGE_JPEG;
		imagetypes["JPEG"] = IMAGE_JPEG;

		imagetypes["df3"] = IMAGE_DF3;
		imagetypes["DF3"] = IMAGE_DF3;

		imagetypes["Omap"] = IMAGE_OMAP;
		imagetypes["omap"] = IMAGE_OMAP;
		imagetypes["OMAP"] = IMAGE_OMAP;
		imagetypes["BRIX"] = IMAGE_OMAP;
		imagetypes["brix"] = IMAGE_OMAP;
		imagetypes["DSN6"] = IMAGE_OMAP;

		imagetypes["situs"] = IMAGE_SITUS;
		imagetypes["SITUS"] = IMAGE_SITUS;

		imagetypes["ser"] = IMAGE_SER;
		imagetypes["SER"] = IMAGE_SER;

		initialized = true;
	}

	ImageType result = IMAGE_UNKNOWN;

	if (imagetypes.find(file_ext) != imagetypes.end()) {
		result = imagetypes[file_ext];
	}

	EXITFUNC;

	return result;
}

bool EMUtil::is_valid_filename(const string & filename)
{
	ImageType type = get_image_ext_type(Util::get_filename_ext(filename));

	return (type != IMAGE_UNKNOWN);
}

EMUtil::ImageType EMUtil::fast_get_image_type(const string & filename,
											  const void *first_block,
											  off_t file_size)
{
	ENTERFUNC;

	Assert(filename != "");
	Assert(first_block != 0);
	Assert(file_size > 0);

#ifdef ENABLE_V4L2
	if (filename.compare(0,5,"/dev/")==0) return IMAGE_V4L;
#endif

	string ext = Util::get_filename_ext(filename);

	if (ext == "") {
		return IMAGE_UNKNOWN;
	}

	ImageType image_type = get_image_ext_type(ext);

	switch (image_type) {
	case IMAGE_MRC:
		if (MrcIO::is_valid(first_block, file_size)) {
			return IMAGE_MRC;
		}
		break;
	case IMAGE_DM3:
		if (DM3IO::is_valid(first_block)) {
			return IMAGE_DM3;
		}
		break;
	case IMAGE_DM4:
		if (DM4IO::is_valid(first_block)) {
			return IMAGE_DM4;
		}
		break;
#ifdef EM_HDF5
	case IMAGE_HDF:
		if (HdfIO2::is_valid(first_block)) {
			return IMAGE_HDF;
		}
		break;
#endif
	case IMAGE_LST:
		if (LstIO::is_valid(first_block)) {
			return IMAGE_LST;
		}
		break;
	case IMAGE_LSTFAST:
		if (LstFastIO::is_valid(first_block)) {
			return IMAGE_LSTFAST;
		}
		break;
#ifdef EM_TIFF
	case IMAGE_TIFF:
		if (TiffIO::is_valid(first_block)) {
			return IMAGE_TIFF;
		}
		break;
#endif
	case IMAGE_SPIDER:
		if (SpiderIO::is_valid(first_block)) {
			return IMAGE_SPIDER;
		}
		break;
	case IMAGE_SINGLE_SPIDER:
		if (SingleSpiderIO::is_valid(first_block)) {
			return IMAGE_SINGLE_SPIDER;
		}
		break;
	case IMAGE_PIF:
		if (PifIO::is_valid(first_block)) {
			return IMAGE_PIF;
		}
		break;
#ifdef EM_PNG
	case IMAGE_PNG:
		if (PngIO::is_valid(first_block)) {
			return IMAGE_PNG;
		}
		break;
#endif
	case IMAGE_VTK:
		if (VtkIO::is_valid(first_block)) {
			return IMAGE_VTK;
		}
		break;
	case IMAGE_PGM:
		if (PgmIO::is_valid(first_block)) {
			return IMAGE_PGM;
		}
		break;
	case IMAGE_ICOS:
		if (IcosIO::is_valid(first_block)) {
			return IMAGE_ICOS;
		}
		break;
	case IMAGE_SAL:
		if (SalIO::is_valid(first_block)) {
			return IMAGE_SAL;
		}
		break;
	case IMAGE_AMIRA:
		if (AmiraIO::is_valid(first_block)) {
			return IMAGE_AMIRA;
		}
		break;
	case IMAGE_XPLOR:
		if (XplorIO::is_valid(first_block)) {
			return IMAGE_XPLOR;
		}
		break;
	case IMAGE_GATAN2:
		if (Gatan2IO::is_valid(first_block)) {
			return IMAGE_GATAN2;
		}
		break;
	case IMAGE_EM:
		if (EmIO::is_valid(first_block, file_size)) {
			return IMAGE_EM;
		}
		break;
	case IMAGE_DF3:
		if (EmIO::is_valid(first_block, file_size)) {
			return IMAGE_DF3;
		}
		break;
	case IMAGE_OMAP:
		if (OmapIO::is_valid(first_block, file_size)) {
			return IMAGE_OMAP;
		}
		break;
	case IMAGE_SITUS:
		if (SitusIO::is_valid(first_block)) {
			return IMAGE_SITUS;
		}
		break;
	case IMAGE_SER:
		if (SerIO::is_valid(first_block)) {
			return IMAGE_SER;
		}
		break;
	case IMAGE_IMAGIC:
		if (ImagicIO::is_valid(first_block)) {
			return IMAGE_IMAGIC;
		}
		break;
	default:
		return IMAGE_UNKNOWN;
	}

	EXITFUNC;

	return IMAGE_UNKNOWN;
}

EMUtil::ImageType EMUtil::get_image_type(const string & in_filename)
{
	ENTERFUNC;

	Assert(in_filename != "");

#ifdef ENABLE_V4L2
	if (in_filename.compare(0,5,"/dev/")==0) return IMAGE_V4L;
#endif

	string filename = in_filename;

	string old_ext = Util::get_filename_ext(filename);

	if (old_ext == ImagicIO::IMG_EXT) {
		filename = Util::change_filename_ext(filename, ImagicIO::HED_EXT);
	}
	else if (old_ext == "hdf") {
		return IMAGE_HDF;
	}

	FILE *in = fopen(filename.c_str(), "rb");

	if (! in) {
		throw FileAccessException(filename);
	}

	char first_block[1024];
	size_t n = fread(first_block, sizeof(char), sizeof(first_block), in);
	portable_fseek(in, 0, SEEK_END);
	off_t file_size = portable_ftell(in);

	if (n == 0) {
//		This produces annoying console messages		
//		LOGERR("file '%s' is an empty file", filename.c_str());
		fclose(in);
		return IMAGE_UNKNOWN;
	}

	fclose(in);

	ImageType image_type = fast_get_image_type(filename, first_block, file_size);

	if (image_type != IMAGE_UNKNOWN) {
		return image_type;
	}

	if (SpiderIO::is_valid(first_block)) {
		image_type = IMAGE_SPIDER;
	}
	else if (SingleSpiderIO::is_valid(first_block)) {
		image_type = IMAGE_SINGLE_SPIDER;
	}
	else if (MrcIO::is_valid(first_block, file_size)) {
		image_type = IMAGE_MRC;
	}
	else if (DM3IO::is_valid(first_block)) {
		image_type = IMAGE_DM3;
	}
	else if (DM4IO::is_valid(first_block)) {
		image_type = IMAGE_DM4;
	}
#ifdef EM_HDF5
	else if (HdfIO2::is_valid(first_block)) {
		image_type = IMAGE_HDF;
	}
#endif
	else if (LstIO::is_valid(first_block)) {
		image_type = IMAGE_LST;
	}
	else if (LstFastIO::is_valid(first_block)) {
		image_type = IMAGE_LSTFAST;
	}
#ifdef EM_TIFF
	else if (TiffIO::is_valid(first_block)) {
		image_type = IMAGE_TIFF;
	}
#endif
	else if (PifIO::is_valid(first_block)) {
		image_type = IMAGE_PIF;
	}
#ifdef EM_PNG
	else if (PngIO::is_valid(first_block)) {
		image_type = IMAGE_PNG;
	}
#endif
	else if (VtkIO::is_valid(first_block)) {
		image_type = IMAGE_VTK;
	}
	else if (PgmIO::is_valid(first_block)) {
		image_type = IMAGE_PGM;
	}
	else if (IcosIO::is_valid(first_block)) {
		image_type = IMAGE_ICOS;
	}
	else if (SalIO::is_valid(first_block)) {
		image_type = IMAGE_SAL;
	}
	else if (AmiraIO::is_valid(first_block)) {
		image_type = IMAGE_AMIRA;
	}
	else if (XplorIO::is_valid(first_block)) {
		image_type = IMAGE_XPLOR;
	}
	else if (Gatan2IO::is_valid(first_block)) {
		image_type = IMAGE_GATAN2;
	}
	else if (FitsIO::is_valid(first_block)) {
		image_type = IMAGE_FITS;
	}
	else if (EmIO::is_valid(first_block, file_size)) {
		image_type = IMAGE_EM;
	}
	else if(OmapIO::is_valid(first_block, file_size)) {
		image_type = IMAGE_OMAP;
	}
	else if(SitusIO::is_valid(first_block)) {
		image_type = IMAGE_SITUS;
	}
	else if(SerIO::is_valid(first_block)) {
		image_type = IMAGE_SER;
	}
	else if (ImagicIO::is_valid(first_block)) {
		image_type = IMAGE_IMAGIC;
	}
	else if (Df3IO::is_valid(first_block)) {
		image_type = IMAGE_DF3;
	}
	else {
		// LOGERR("I don't know this image's type: '%s'", filename.c_str());

		throw ImageFormatException("invalid image type");
	}

	EXITFUNC;

	return image_type;
}

int EMUtil::get_image_count(const string & filename)
{
	ENTERFUNC;

	Assert(filename != "");

	int nimg = 0;

	ImageIO *imageio = get_imageio(filename, ImageIO::READ_ONLY);

	if (imageio) {
		nimg = imageio->get_nimg();
	}

   EMUtil::close_imageio(filename, imageio);

	imageio = NULL;

	EXITFUNC;

	return nimg;
}

ImageIO *EMUtil::get_imageio(const string & filename, int rw,
							 ImageType image_type)
{
	ENTERFUNC;

   // printf("EMUtil::get_imageio\n");

	Assert(filename != "");
	Assert(rw == ImageIO::READ_ONLY ||
		   rw == ImageIO::READ_WRITE ||
		   rw == ImageIO::WRITE_ONLY);

	ImageIO *imageio = 0;
   int persist = 0;

   #ifdef IMAGEIO_CACHE
   imageio = GlobalCache::instance()->get_imageio(filename, rw);
   if (imageio) {
    	return imageio;
   }
   #endif

	ImageIO::IOMode rw_mode = static_cast < ImageIO::IOMode > (rw);

	if (image_type == IMAGE_UNKNOWN) {
		image_type = get_image_type(filename);
	}

	if (image_type == IMAGE_UNKNOWN) {
		if(rw == ImageIO::WRITE_ONLY || rw == ImageIO::READ_WRITE) {
			throw ImageFormatException("writing to this image format not supported.");
		}
	}

	switch (image_type) {
#ifdef ENABLE_V4L2
	case IMAGE_V4L:
		imageio = new V4L2IO(filename, rw_mode);
		break;
#endif
	case IMAGE_MRC:
		imageio = new MrcIO(filename, rw_mode);
		break;
	case IMAGE_IMAGIC:
		imageio = new ImagicIO2(filename, rw_mode);
		if (rw_mode==ImageIO::READ_ONLY && ((ImagicIO2 *)imageio)->init_test()==-1 ) {
			delete imageio;
			imageio = new ImagicIO(filename, rw_mode);
		}
		break;
	case IMAGE_DM3:
		imageio = new DM3IO(filename, rw_mode);
		break;
	case IMAGE_DM4:
		imageio = new DM4IO(filename, rw_mode);
		break;
#ifdef EM_TIFF
	case IMAGE_TIFF:
		imageio = new TiffIO(filename, rw_mode);
		break;
#endif
#ifdef EM_HDF5
	case IMAGE_HDF:
        persist = 30;
        if (rw_mode != ImageIO::READ_ONLY) {
            persist = 3;
        }
		imageio = new HdfIO2(filename, rw_mode);
		if (((HdfIO2 *)imageio)->init_test()==-1) {
			delete imageio;
			imageio = new HdfIO(filename, rw_mode);
		}
		break;
#endif	//EM_HDF5
	case IMAGE_LST:
		imageio = new LstIO(filename, rw_mode);
		break;
	case IMAGE_LSTFAST:
		imageio = new LstFastIO(filename, rw_mode);
		break;
	case IMAGE_PIF:
		imageio = new PifIO(filename, rw_mode);
		break;
	case IMAGE_VTK:
		imageio = new VtkIO(filename, rw_mode);
		break;
	case IMAGE_SPIDER:
		imageio = new SpiderIO(filename, rw_mode);
		break;
	case IMAGE_SINGLE_SPIDER:
		imageio = new SingleSpiderIO(filename, rw_mode);
		break;
	case IMAGE_PGM:
		imageio = new PgmIO(filename, rw_mode);
		break;
#ifdef EM_JPEG
	case IMAGE_JPEG:
		imageio = new JpegIO(filename,rw_mode);
		break;
#endif
	case IMAGE_ICOS:
		imageio = new IcosIO(filename, rw_mode);
		break;
#ifdef EM_PNG
	case IMAGE_PNG:
		imageio = new PngIO(filename, rw_mode);
		break;
#endif
	case IMAGE_SAL:
		imageio = new SalIO(filename, rw_mode);
		break;
	case IMAGE_AMIRA:
		imageio = new AmiraIO(filename, rw_mode);
		break;
	case IMAGE_GATAN2:
		imageio = new Gatan2IO(filename, rw_mode);
		break;
	case IMAGE_EM:
		imageio = new EmIO(filename, rw_mode);
		break;
	case IMAGE_XPLOR:
		imageio = new XplorIO(filename, rw_mode);
		break;
	case IMAGE_FITS:
		imageio = new FitsIO(filename, rw_mode);
		break;
	case IMAGE_DF3:
		imageio = new Df3IO(filename, rw_mode);
		break;
	case IMAGE_OMAP:
		imageio = new OmapIO(filename, rw_mode);
		break;
	case IMAGE_SITUS:
		imageio = new SitusIO(filename, rw_mode);
		break;
	case IMAGE_SER:
		imageio = new SerIO(filename, rw_mode);
		break;
	default:
		break;
	}

   #ifdef IMAGEIO_CACHE
	if (persist > 0) {
		GlobalCache::instance()->add_imageio(filename, rw, persist, imageio);
	}
   #endif

	EXITFUNC;

	return imageio;
}

void EMUtil::close_imageio(const string & filename, const ImageIO * io)
{
    //printf("EMUtil::close_imageio\n");
    #ifdef IMAGEIO_CACHE
    if (GlobalCache::instance()->contains(filename)) {
        GlobalCache::instance()->close_imageio(filename);
    } else {
        delete io;
    }
    #else
    delete io;
    #endif
}

const char *EMUtil::get_imagetype_name(ImageType t)
{
	switch (t) {
	case IMAGE_V4L:
		return "V4L2";
		break;
	case IMAGE_MRC:
		return "MRC";
		break;
	case IMAGE_SPIDER:
		return "SPIDER";
		break;
	case IMAGE_SINGLE_SPIDER:
		return "Single-SPIDER";
		break;
	case IMAGE_IMAGIC:
		return "IMAGIC";
		break;
	case IMAGE_PGM:
		return "PGM";
		break;
	case IMAGE_LST:
		return "LST";
		break;
	case IMAGE_LSTFAST:
		return "Fast LST";
		break;
	case IMAGE_PIF:
		return "PIF";
		break;
	case IMAGE_PNG:
		return "PNG";
		break;
	case IMAGE_HDF:
		return "HDF5";
		break;
	case IMAGE_DM3:
		return "GatanDM3";
		break;
	case IMAGE_DM4:
		return "GatanDM4";
		break;
	case IMAGE_TIFF:
		return "TIFF";
		break;
	case IMAGE_VTK:
		return "VTK";
		break;
	case IMAGE_SAL:
		return "HDR";
		break;
	case IMAGE_ICOS:
		return "ICOS_MAP";
		break;
	case IMAGE_EMIM:
		return "EMIM";
		break;
	case IMAGE_GATAN2:
		return "GatanDM2";
		break;
	case IMAGE_JPEG:
		return "JPEG";
		break;
	case IMAGE_AMIRA:
		return "AmiraMesh";
		break;
	case IMAGE_XPLOR:
		return "XPLOR";
		break;
	case IMAGE_EM:
		return "EM";
		break;
	case IMAGE_FITS:
		return "FITS";
		break;
	case IMAGE_DF3:
		return "DF3";
		break;
	case IMAGE_OMAP:
		return "OMAP";
		break;
	case IMAGE_SITUS:
		return "SITUS";
		break;
	case IMAGE_SER:
		return "SER";
		break;
	case IMAGE_UNKNOWN:
		return "unknown";
	}

	return "unknown";
}

const char *EMUtil::get_datatype_string(EMDataType type)
{
	switch (type) {
	case EM_CHAR:
		return "CHAR";
	case EM_UCHAR:
		return "UNSIGNED CHAR";
	case EM_SHORT:
		return "SHORT";
	case EM_USHORT:
		return "UNSIGNED SHORT";
	case EM_INT:
		return "INT";
	case EM_UINT:
		return "UNSIGNED INT";
	case EM_FLOAT:
		return "FLOAT";
	case EM_DOUBLE:
		return "DOUBLE";
	case EM_SHORT_COMPLEX:
		return "SHORT_COMPLEX";
	case EM_USHORT_COMPLEX:
		return "USHORT_COMPLEX";
	case EM_FLOAT_COMPLEX:
		return "FLOAT_COMPLEX";
	case EM_UNKNOWN:
		return "UNKNOWN";
	}

	return "UNKNOWN";
}

void EMUtil::get_region_dims(const Region * area, int nx, int *area_x,
							 int ny, int *area_y, int nz, int *area_z)
{
	Assert(area_x);
	Assert(area_y);

	if (! area) {
		*area_x = nx;
		*area_y = ny;

		if (area_z) {
			*area_z = nz;
		}
	}
	else {
		Vec3i size = area->get_size();
		*area_x = size[0];
		*area_y = size[1];

		if (area_z) {
			if (area->get_ndim() > 2 && nz > 1) {
				*area_z = size[2];
			}
			else {
				*area_z = 1;
			}
		}

	}
}

void EMUtil::get_region_origins(const Region * area, int *p_x0, int *p_y0, int *p_z0,
								int nz, int image_index)
{
	Assert(p_x0);
	Assert(p_y0);

	if (area) {
		*p_x0 = static_cast < int >(area->origin[0]);
		*p_y0 = static_cast < int >(area->origin[1]);

		if (p_z0 && nz > 1 && area->get_ndim() > 2) {
			*p_z0 = static_cast < int >(area->origin[2]);
		}
	}
	else {
		*p_x0 = 0;
		*p_y0 = 0;

		if (p_z0) {
			*p_z0 = nz > 1 ? 0 : image_index;
		}
	}
}

double EMUtil::mode_size_product(size_t factor, size_t mode_size)
{
	const size_t mode_size_half = 11111111;

	double product;

	if (mode_size == mode_size_half) {
		product = factor * 0.5;
	}
	else {
		product = factor * mode_size;
	}

	return product;
}

void EMUtil::process_region_io(void *vdata, FILE * file,
							   int rw_mode, int image_index,
							   size_t mode_size, int nx, int ny, int nz,
							   const Region * area, bool need_flip,
							   ImageType imgtype, int pre_row, int post_row)
{
	Assert(vdata != 0);
	Assert(file != 0);
	Assert(rw_mode == ImageIO::READ_ONLY ||
		   rw_mode == ImageIO::READ_WRITE ||
		   rw_mode == ImageIO::WRITE_ONLY);

	if (mode_size == 0) throw UnexpectedBehaviorException("The mode size was 0?");

	const size_t mode_size_half = 11111111;

	unsigned char * cdata = (unsigned char *)vdata;

	int dx0 = 0; // data x0
	int dy0 = 0; // data y0
	int dz0 = 0; // data z0

	int fx0 = 0; // file x0
	int fy0 = 0; // file y0
	int fz0 = nz > 1 ? 0 : image_index; // file z0

	int xlen = 0;
	int ylen = 0;
	int zlen = 0;

	get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

	bool debug = (getenv("DEBUG_IO") != NULL);

	if (debug) {
		printf ("-------------- process_region_io --------------\n");
		printf ("rw, indx, modsiz, nx, ny, nz = %d %d %ld %d %d %d\n",
					rw_mode, image_index, mode_size, nx, ny, nz);
		printf ("flip, imtyp, prerow, postrow = %d %d %d %d\n",
					(int) need_flip, (int) imgtype, pre_row, post_row);
		printf ("xlen, ylen, zlen = %d %d %d\n", xlen, ylen, zlen);

		if (area != NULL) {
			printf ("x0, y0, z0, mx, my, mz, ndim = %g %g %g %g %g %g %d\n",
						area->x_origin(),
						area->y_origin(),
						area->z_origin(),
						area->get_width(),
						area->get_height(),
						area->get_depth(),
						area->get_ndim());
		}
	}

	if (area) { // Accommodate for all boundary overlaps of the region
		Vec3i origin = area->get_origin();

		fx0 = origin[0]; dx0 = origin[0];
		fy0 = origin[1]; dy0 = origin[1];

		if (nz > 1 && area->get_ndim() > 2) {
			fz0 = origin[2]; dz0 = origin[2];
		}

		if (need_flip) {
			Vec3i size = area->get_size();
			fy0 = ny-(origin[1]+size[1]);
		}

		if (fx0 < 0) {
			dx0 *= -1;
			xlen = xlen + fx0; // because there are less reads
			fx0 = 0;
		} else {
			dx0 = 0;
			//fx0 *= -1;
		}

		if (fy0 < 0) {
			dy0 *= -1;
			ylen = ylen + fy0; // because there are less reads
			fy0 = 0;
		} else {
			if (need_flip){
				dy0*=-1;
			}
			else dy0 = 0;
			//fy0 *= -1;
		}

		if (fz0 < 0) {
			dz0 *= -1;
			zlen = zlen + fz0; // because there are less reads
			fz0 = 0;
		} else {
			dz0 = 0;
			//fz0 *= -1;
		}

		if ((fx0 + xlen)> nx) xlen = nx-fx0;
		if ((fy0 + ylen)> ny) ylen = ny-fy0;
		if ((fz0 + zlen)> nz && nz > 1) zlen = nz-fz0;

		// This is fine - the region was entirely outside the image

		if ( xlen <= 0 || ylen <= 0 || zlen <= 0 ) return;

		if (mode_size == mode_size_half) {
			// Have an area with 8 bit packed MRC format
			// (2 4-bit valus per 8 bit byte),
			// with an effective mode size of half a byte,
			// where most x-pixel parameters need to be even
			// for everything to work well.

			bool error = false;

			if (fx0 % 2 != 0  &&  false) { // Don't check right now.
				cout << "First region x-pixel location "
						  "must be even for packed 8 bit MRC." << endl;
				error = true;
			}

			if ((nx - fx0 - xlen) % 2 != 0  &&  false) { // Don't check right now.
				cout << "No. of x-pixel locations after region "
						  "must be even for packed 8 bit MRC." << endl;
				error = true;
			}

			if (xlen % 2 != 0) {
				cout << "No. of region x-pixels "
						  "must be even for packed 8 bit MRC." << endl;
				error = true;
			}

			if (error) {
				return;
			}
		}
	}

	if (xlen <= 0) {
		cout << "Xlen was too small " << xlen << endl;

		return;
	}

	Vec3i size;

	if (area != 0) size = area->get_size();
	else size = Vec3d(nx,ny,nz);

	//size_t area_sec_size = ylen    * mode_size_product(xlen,    mode_size);
	size_t memory_sec_size = size[1] * mode_size_product(size[0], mode_size);
	size_t img_row_size    = mode_size_product(nx,      mode_size) +
									 pre_row + post_row;
	size_t area_row_size   = mode_size_product(xlen,    mode_size);
	size_t memory_row_size = mode_size_product(size[0], mode_size);

	if ( area_row_size <= 0 ) {
		cout << "Xlen was too small " << xlen << " mode_size " << mode_size << endl;

		return;
	}

	size_t x_pre_gap  = mode_size_product(fx0, mode_size);
	size_t x_post_gap = mode_size_product(nx - fx0 - xlen, mode_size);

	size_t y_pre_gap  = fy0 * img_row_size;
	size_t y_post_gap = (ny - fy0 - ylen) * img_row_size;

	if (x_pre_gap  < 0) x_pre_gap  = 0;
	if (x_post_gap < 0) x_post_gap = 0;

	size_t extra = img_row_size - (x_pre_gap + area_row_size + x_post_gap);

	if (extra > 0) x_post_gap += extra;

	portable_fseek(file, img_row_size * ny * fz0, SEEK_CUR);

	float nxlendata[1];
	int floatsize = (int) sizeof(float);
	nxlendata[0] = (float)(nx * floatsize);

	if (debug) {
		printf ("fz0, dz0, xlen, ylen, zlen = %d %d %d %d %d\n",
					fz0, dz0, xlen, ylen, zlen);
		printf ("x_pre_gap, x_post_gap, y_pre_gap, y_post_gap = %d %d %d %d\n",
					x_pre_gap, x_post_gap, y_pre_gap, y_post_gap);
		printf ("mem_sec_siz, img_row_siz, are_row_siz, mem_row_siz = %ld %ld %ld %ld\n",
					memory_sec_size, img_row_size, area_row_size, memory_row_size);
		printf ("-----------------------------------------------\n");
	}

	for (int k = dz0; k < (dz0+zlen); k++) {
		// k is image/slice number, starting from 0

		if (y_pre_gap > 0) {
			portable_fseek(file, y_pre_gap, SEEK_CUR);
		}

		//long k2 = k * area_sec_size;
		long k2 = k*memory_sec_size;

		for (int j = dy0; j < (dy0+ylen); j++) {
			if (pre_row > 0) {
				if (imgtype == IMAGE_ICOS && rw_mode != ImageIO::READ_ONLY && !area) {
					fwrite(nxlendata, floatsize, 1, file);
				}
				else {
					portable_fseek(file, pre_row, SEEK_CUR);
				}
			}

			if (x_pre_gap > 0) {
				portable_fseek(file, x_pre_gap, SEEK_CUR);
			}

			int jj = j;

			if (need_flip) {
				jj = (dy0+ylen) - 1 - j;

				// region considerations add complications
				// in the flipping scenario (imagic format)

				if (dy0 > 0) {
					jj += dy0;
				}
			}

			if (rw_mode == ImageIO::READ_ONLY) {
				if (fread(&cdata[k2 + jj * memory_row_size +
						  (size_t) mode_size_product(dx0, mode_size)],
						  area_row_size, 1, file) != 1) {

//					cout << jj << " " << k2 << " " << memory_row_size
//						  << " " << dx0 << " "
//						  << mode_size << " " << area_row_size << " "
//						  << cdata << ftell(file) << "done" << endl;

					cout << "Reached premature end-of-file reading "
						  << "region from image/slice number " << k
						  << " of file with " << ftell(file)
						  << " bytes." << endl;

					throw ImageReadException("Unknownfilename",
													 "incomplete data read");
				}
			}
			else {
				if (fwrite(&cdata[k2 + jj * memory_row_size +
							(size_t) mode_size_product(dx0, mode_size)],
						   area_row_size, 1, file) != 1) {
					throw ImageWriteException("", "incomplete data write");
				}
			}

			if (x_post_gap > 0) {
				portable_fseek(file, x_post_gap, SEEK_CUR);
			}

			if (post_row > 0) {
				if (imgtype == IMAGE_ICOS && rw_mode != ImageIO::READ_ONLY && !area) {
					fwrite(nxlendata, floatsize, 1, file);
				}
				else {
					portable_fseek(file, post_row, SEEK_CUR);
				}
			}
		}

		if (y_post_gap > 0) {
			portable_fseek(file, y_post_gap, SEEK_CUR);
		}
	}
}

void EMUtil::dump_dict(const Dict & dict)
{
	vector < string > keys = dict.keys();
	vector < EMObject > values = dict.values();

	for (unsigned int i = 0; i < keys.size(); i++) {
		EMObject obj = values[i];

		if (! obj.is_null()) {
			string val = obj.to_str();

			if (keys[i] == "datatype") {
				val = get_datatype_string((EMDataType) (int) obj);
			}

			fprintf(stdout, "%25s\t%s\n", keys[i].c_str(), val.c_str());
		}
	}
}

bool EMUtil::is_same_size(const EMData * const em1, const EMData * const em2)
{
	if (em1->get_xsize() == em2->get_xsize() &&
		em1->get_ysize() == em2->get_ysize() &&
		em1->get_zsize() == em2->get_zsize()) {
		return true;
	}

	return false;
}

bool EMUtil::is_complex_type(EMDataType datatype)
{
	if (datatype == EM_SHORT_COMPLEX ||
		datatype == EM_USHORT_COMPLEX ||
		datatype == EM_FLOAT_COMPLEX) {

		return true;
	}

	return false;
}

EMData *EMUtil::vertical_acf(const EMData * image, int maxdy)
{
	if (!image) {
		throw NullPointerException("NULL Image");
	}

	EMData *ret = new EMData();

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	if (maxdy <= 1) {
		maxdy = ny / 8;
	}

	ret->set_size(nx, maxdy, 1);

	float *data = image->get_data();
	float *ret_data = ret->get_data();

	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < maxdy; y++) {
			float dot = 0;

			for (int yy = maxdy; yy < ny - maxdy; yy++) {
				dot += data[x + (yy + y) * nx] * data[x + (yy - y) * nx];
			}

			ret_data[x + y * nx] = dot;
		}
	}

	ret->update();

	return ret;
}

EMData *EMUtil::make_image_median(const vector < EMData * >&image_list)
{
	if (image_list.size() == 0) {
		return 0;
	}

	EMData *image0 = image_list[0];

	int image0_nx = image0->get_xsize();
	int image0_ny = image0->get_ysize();
	int image0_nz = image0->get_zsize();
	size_t size = (size_t)image0_nx * image0_ny * image0_nz;

	EMData *result = new EMData();

	result->set_size(image0_nx, image0_ny, image0_nz);

	float *dest = result->get_data();
	int nitems = static_cast < int >(image_list.size());
	float *srt = new float[nitems];
	float **src = new float *[nitems];

	for (int i = 0; i < nitems; i++) {
		src[i] = image_list[i]->get_data();
	}

	for (size_t i = 0; i < size; ++i) {
		for (int j = 0; j < nitems; j++) {
			srt[j] = src[j][i];
		}

		for (int j = 0; j < nitems; j++) {
			for (int k = j + 1; k < nitems; k++) {
				if (srt[j] < srt[k]) {
					float v = srt[j];
					srt[j] = srt[k];
					srt[k] = v;
				}
			}
		}

		int l = nitems / 2;

		if (nitems < 3) {
			dest[i] = srt[l];
		}
		else {
			dest[i] = (srt[l] + srt[l + 1] + srt[l - 1]) / 3.0f;
		}
	}

	if (srt) {
		delete [] srt;
		srt = NULL;
	}

	if (src) {
		delete [] src;
		src = NULL;
	}

	result->update();

	return result;
}

bool EMUtil::is_same_ctf(const EMData * image1, const EMData * image2)
{
	if (!image1) {
		throw NullPointerException("image1 is NULL");
	}

	if (!image2) {
		throw NullPointerException("image2 is NULL");
	}

	Ctf *ctf1 = image1->get_ctf();
	Ctf *ctf2 = image2->get_ctf();

	if ((!ctf1 && !ctf2) && (image1->has_ctff() == false && image2->has_ctff() == false)) {
		return true;
	}

	if (ctf1 && ctf2) {
		bool result = ctf1->equal(ctf2);
		delete ctf1;
		ctf1 = NULL;
		delete ctf2;
		ctf2 = NULL;

		return result;
	}

	return false;
}

static int imgscore_cmp(const void *imgscore1, const void *imgscore2)
{
	Assert(imgscore1 != 0);
	Assert(imgscore2 != 0);

	float c = ((ImageScore *)imgscore1)->score - ((ImageScore *)imgscore2)->score;

	if (c<0) {
		return 1;
	}
	else if (c>0) {
		return -1;
	}

	return 0;
}

ImageSort::ImageSort(int nn)
{
	Assert(nn > 0);
	n = nn;
	image_scores = new ImageScore[n];
}

ImageSort::~ImageSort()
{
	if (image_scores) {
		delete [] image_scores;
		image_scores = NULL;
	}
}

void ImageSort::sort()
{
	qsort(image_scores, n, sizeof(ImageScore), imgscore_cmp);
}

void ImageSort::set(int i, float score)
{
	Assert(i >= 0);

	image_scores[i] = ImageScore(i, score);
}

int ImageSort::get_index(int i) const
{
	Assert(i >= 0);

	return image_scores[i].index;
}


float ImageSort::get_score(int i) const
{
	Assert(i >= 0);

	return image_scores[i].score;
}


int ImageSort::size() const
{
	return n;
}

void EMUtil::process_ascii_region_io(float *data, FILE * file, int rw_mode,
									 int , size_t mode_size, int nx, int ny, int nz,
									 const Region * area, bool has_index_line,
									 int nitems_per_line, const char *outformat)
{
	Assert(data != 0);
	Assert(file != 0);
	Assert(rw_mode == ImageIO::READ_ONLY ||
		   rw_mode == ImageIO::READ_WRITE ||
		   rw_mode == ImageIO::WRITE_ONLY);

	int xlen = 0, ylen = 0, zlen = 0;
	get_region_dims(area, nx, &xlen, ny, &ylen, nz, &zlen);

	int x0 = 0;
	int y0 = 0;
	int z0 = 0;

	if (area) {
		x0 = (int)area->origin[0];
		y0 = (int)area->origin[1];
		z0 = (int)area->origin[2];
	}

	int nlines_per_sec = (nx *ny) / nitems_per_line;
	int nitems_last_line = (nx * ny) % nitems_per_line;

	if (nitems_last_line != 0) {
		nlines_per_sec++;
	}

	if (has_index_line) {
		nlines_per_sec++;
	}

	if (z0 > 0) {
		jump_lines(file, z0 * nlines_per_sec);
	}

	int nlines_pre_sec = (y0 * nx + x0) / nitems_per_line;
	int gap_nitems = nx - xlen;
	int ti = 0;
	int rlines = 0;

	for (int k = 0; k < zlen; k++) {
		EMUtil::jump_lines(file, nlines_pre_sec+1);

		int head_nitems = (y0 * nx + x0) % nitems_per_line;
		int tail_nitems = 0;
		bool is_head_read = false;

		for (int j = 0; j < ylen; j++) {
			if (head_nitems > 0 && !is_head_read) {
				EMUtil::process_numbers_io(file, rw_mode, nitems_per_line, mode_size,
										   nitems_per_line-head_nitems,
										   nitems_per_line-1, data, &ti, outformat);
				rlines++;
			}

			EMUtil::process_lines_io(file, rw_mode, nitems_per_line,
									 mode_size, (xlen - head_nitems),
									 data, &ti, outformat);

			rlines += ((xlen - head_nitems) / nitems_per_line);

			tail_nitems = (xlen - head_nitems) % nitems_per_line;

			if ((gap_nitems + tail_nitems) > 0) {
				head_nitems = nitems_per_line -
					(gap_nitems + tail_nitems) % nitems_per_line;
			}
			else {
				head_nitems = 0;
			}

			is_head_read = false;

			if (tail_nitems > 0) {
				if ((gap_nitems < (nitems_per_line-tail_nitems)) &&
					(j != (ylen-1))) {
					EMUtil::exclude_numbers_io(file, rw_mode, nitems_per_line,
											   mode_size, tail_nitems,
											   tail_nitems+gap_nitems-1, data, &ti, outformat);
					is_head_read = true;
					rlines++;
				}
				else {
					EMUtil::process_numbers_io(file, rw_mode, nitems_per_line, mode_size,
											   0, tail_nitems-1, data, &ti, outformat);
					rlines++;
				}
			}

			if (gap_nitems > (nitems_per_line-tail_nitems)) {
				int gap_nlines = (gap_nitems - (nitems_per_line-tail_nitems)) /
					nitems_per_line;

				if (gap_nlines > 0 && j != (ylen-1)) {
					EMUtil::jump_lines(file, gap_nlines);
				}
			}
		}

		int ytail_nitems = (ny-ylen-y0) * nx + (nx-xlen-x0) - (nitems_per_line-tail_nitems);
		EMUtil::jump_lines_by_items(file, ytail_nitems, nitems_per_line);
	}
}

void EMUtil::jump_lines_by_items(FILE * file, int nitems, int nitems_per_line)
{
	Assert(file);
	Assert(nitems_per_line > 0);

	if (nitems <= 0) {
		return;
	}

	int nlines = nitems / nitems_per_line;

	if ((nitems % nitems_per_line) != 0) {
		nlines++;
	}

	if (nlines > 0) {
		jump_lines(file, nlines);
	}
}

void EMUtil::jump_lines(FILE * file, int nlines)
{
	Assert(file);

	if (nlines > 0) {
		char line[MAXPATHLEN];

		for (int l = 0; l < nlines; l++) {
			if (!fgets(line, sizeof(line), file)) {
				Assert("read xplor file failed");
			}
		}
	}
}

void EMUtil::process_numbers_io(FILE * file, int rw_mode,
								int nitems_per_line, size_t mode_size, int start,
								int end, float *data, int *p_i, const char * outformat)
{
	Assert(file);
	Assert(start >= 0);
	Assert(start <= end);
	Assert(end <= nitems_per_line);
	Assert(data);
	Assert(p_i);
	Assert(outformat);

	char line[MAXPATHLEN];

	if (rw_mode == ImageIO::READ_ONLY) {
		if (!fgets(line, sizeof(line), file)) {
			Assert("read xplor file failed");
		}

		int nitems_in_line = (int) (strlen(line) / mode_size);
		Assert(end <= nitems_in_line);
		vector<float> d(nitems_in_line);
		char * pline = line;

		for (int i = 0; i < nitems_in_line; i++) {
			sscanf(pline, "%f", &d[i]);
			pline += (int)mode_size;
		}

		for (int i = start; i <= end; i++) {
			data[*p_i] = d[i];
			(*p_i)++;
		}
	}
	else {
		portable_fseek(file, mode_size * start, SEEK_CUR);

		for (int i = start; i <= end; i++) {
			fprintf(file, outformat, data[*p_i]);
			(*p_i)++;
		}

		portable_fseek(file, mode_size * (nitems_per_line - end-1)+1, SEEK_CUR);
	}
}

void EMUtil::exclude_numbers_io(FILE * file, int rw_mode,
								int nitems_per_line, size_t mode_size, int start,
								int end, float * data, int *p_i, const char * outformat)
{
	Assert(file);
	Assert(mode_size > 0);
	Assert(start >= 0);
	Assert(end <= nitems_per_line);
	Assert(data);
	Assert(p_i);
	Assert(outformat);

	char line[MAXPATHLEN];

	if (rw_mode == ImageIO::READ_ONLY) {

		if (!fgets(line, sizeof(line), file)) {
			Assert("read xplor file failed");
		}

		int nitems_in_line =  (int) (strlen(line) / mode_size);
		Assert(end <= nitems_in_line);

		vector<float> d(nitems_in_line);
		char *pline = line;

		for (int i = 0; i < nitems_in_line; i++) {
			sscanf(pline, "%f", &d[i]);
			pline = pline + (int)mode_size;
		}


		for (int i = 0; i < start; i++) {
			data[*p_i] = d[i];
			(*p_i)++;
		}

		for (int i = end+1; i < nitems_in_line; i++) {
			data[*p_i] = d[i];
			(*p_i)++;
		}
	}
	else {
		for (int i = 0; i < start; i++) {
			fprintf(file, outformat, data[*p_i]);
			(*p_i)++;
		}

		portable_fseek(file, (end-start+1) * mode_size, SEEK_CUR);

		for (int i = end+1; i < nitems_per_line; i++) {
			fprintf(file, outformat, data[*p_i]);
			(*p_i)++;
		}
		portable_fseek(file, 1, SEEK_CUR);
	}
}

void EMUtil::process_lines_io(FILE * file, int rw_mode,
							  int nitems_per_line, size_t mode_size,
							  int nitems, float *data, int *p_i,
							  const char * outformat)
{
	Assert(file);
	Assert(data);
	Assert(p_i);

	if (nitems > 0) {
		int nlines = nitems / nitems_per_line;

		for (int i = 0; i < nlines; i++) {
			EMUtil::process_numbers_io(file, rw_mode, nitems_per_line, mode_size, 0,
									   nitems_per_line-1, data, p_i, outformat);
		}
	}
}

vector<string> EMUtil::get_euler_names(const string & euler_type)
{
    vector<string> v;
    string b = "euler_";

    if (euler_type == "EMAN") {
        v.push_back(b + "alt");
        v.push_back(b + "az");
        v.push_back(b + "phi");
    }
    else if (euler_type == "MRC") {
        v.push_back(b + "theta");
        v.push_back(b + "phi");
        v.push_back(b + "omega");
    }
    else if (euler_type == "IMAGIC") {
        v.push_back(b + "alpha");
        v.push_back(b + "beta");
        v.push_back(b + "gamma");
    }
    else if (euler_type == "SPIDER") {
        v.push_back(b + "phi");
        v.push_back(b + "theta");
        v.push_back(b + "gamma");
    }
    else if (euler_type == "SPIN" ||
             euler_type == "SGIROT") {
        v.push_back(b + "q");
        v.push_back(b + "n1");
        v.push_back(b + "n2");
        v.push_back(b + "n3");
    }

    else if (euler_type == "QUATERNION") {
        v.push_back(b + "e0");
        v.push_back(b + "e1");
        v.push_back(b + "e2");
        v.push_back(b + "e3");
    }

    return v;
}

vector<EMObject> EMUtil::get_all_attributes(const string & file_name, const string & attr_name)
{
	vector<EMObject> v;

	Assert(file_name != "");
	Assert(attr_name != "");

	vector< shared_ptr<EMData> > vpImg = EMData::read_images(file_name, vector<int>(), true);
	vector< shared_ptr<EMData> >::const_iterator iter;

	for (iter = vpImg.begin(); iter!=vpImg.end(); ++iter) {
		v.push_back((*iter)->get_attr_default(attr_name));
	}

	return v;
}

void EMUtil::getRenderLimits(const Dict & dict, float & rendermin, float & rendermax)
{
	const char    min_or_max_flag     = 'm';
	const char    std_devs_flag       = 's';

	const float   flag_base           =  1.0e30;
	const float   use_data_min_or_max = -flag_base;

	string        svalue;
	const char *  str;
	char          flag;
	float         num_std_devs;

	bool          debug = (getenv("DEBUG_RENDER_LIMITS") != NULL);

	rendermin = 0.0;
	rendermax = 0.0;

	if (debug) {
		printf ("into RenderLimits, rmin = %g, rmax = %g\n", rendermin, rendermax);
	}

	if (dict.has_key("render_min")) {
		svalue = static_cast<string>(dict["render_min"]);
		str    = svalue.c_str();
		flag   = str[0];

		if (debug) printf ("render_min = %s\n", str);

		if (flag == min_or_max_flag) {
			rendermin = use_data_min_or_max;
		}
		else if (flag == std_devs_flag) {
			num_std_devs = 4.0;
			sscanf (str+1, "%g", & num_std_devs);
			rendermin = flag_base * num_std_devs;
		}
		else {
			rendermin = (float) dict["render_min"];
		}
	}

	if (dict.has_key("render_max")) {
		svalue = static_cast<string>(dict["render_max"]);
		str    = svalue.c_str();
		flag   = str[0];

		if (debug) printf ("render_max = %s\n", str);

		if (flag == min_or_max_flag) {
			rendermax = use_data_min_or_max;
		}
		else if (flag == std_devs_flag) {
			num_std_devs = 4.0;
			sscanf (str+1, "%g", & num_std_devs);
			rendermax = flag_base * num_std_devs;
		}
		else {
			rendermax = (float) dict["render_max"];
		}
	}

	if (debug) {
		printf ("out of RenderLimits, rmin = %g, rmax = %g\n", rendermin, rendermax);
	}
}

void EMUtil::getRenderMinMax(float * data, const int nx, const int ny,
				float & rendermin, float & rendermax, const int nz)
{
	const float   flag_base           =  1.0e30;
	const float   use_data_min_or_max = -flag_base;
	const float   use_num_std_devs    =  flag_base / 100.0;

	float         num_std_devs;

	bool          debug = (getenv("DEBUG_RENDER_LIMITS") != NULL);

	if (debug) {
		printf ("into RenderMinMax, rmin = %g, rmax = %g\n", rendermin, rendermax);
	}

	if (rendermax <= rendermin ||
		Util::is_nan(rendermin) || Util::is_nan(rendermax) ||
		fabs(rendermin) > use_num_std_devs ||
		fabs(rendermax) > use_num_std_devs) {

		double m = 0.0f, s = 0.0f;

		size_t size = (size_t)nx*ny*nz;
		float min = data[0], max = data[0];

		for (size_t i = 0; i < size; ++i) {
			m += data[i];
			s += data[i] * data[i];

			min = data[i] < min ? data[i] : min;
			max = data[i] > max ? data[i] : max;
//			if (!Util::goodf(&data[i])) printf("NAN in image at pixel %ld\n",i);
		}

		m /= (float)(size);
		s = sqrt(s/(float)(size)-m*m);

		if (debug) printf ("min, mean, max, s.d. = %g %g %g %g\n", min, m, max, s);

		if (s <= 0 || Util::is_nan(s)) s = 1.0; // this means all data values are the same

		if (rendermin == use_data_min_or_max) {
			rendermin = min;
		}
		else if (rendermin > use_num_std_devs) {
			num_std_devs = rendermin / flag_base;
			rendermin = m - s * num_std_devs;
		}
		else {
			rendermin = m - s * 5.0f;
		}

		if (rendermax == use_data_min_or_max) {
			rendermax = max;
		}
		else if (rendermax > use_num_std_devs) {
			num_std_devs = rendermax / flag_base;
			rendermax = m + s * num_std_devs;
		}
		else {
			rendermax = m + s * 5.0f;
		}

		if (rendermin <= min) rendermin = min;
		if (rendermax >= max) rendermax = max;

//	printf("rendermm %f %f %f %f\n",rendermin,rendermax,m,s);
	}

	if (debug) {
		printf ("out of RenderMinMax, rmin = %g, rmax = %g\n", rendermin, rendermax);
	}
}

#ifdef EM_HDF5
EMObject EMUtil::read_hdf_attribute(const string & filename, const string & key, int image_index)
{
	ImageType image_type = get_image_type(filename);

	if (image_type != IMAGE_HDF) {
		throw ImageFormatException("This function only applies to HDF5 file.");
	}

	HdfIO2* imageio = new HdfIO2(filename, ImageIO::READ_ONLY);
	imageio->init();

	// Each image is in a group for later expansion. Open the group

	hid_t file = imageio->get_fileid();

	char ipath[50];
	sprintf(ipath, "/MDF/images/%d", image_index);

	hid_t igrp = H5Gopen(file, ipath);

	if (igrp < 0) { // group not existed
		throw _NotExistingObjectException(string(ipath));
	}

	string s("EMAN.");
	s += key;
	hid_t attr = H5Aopen_name(igrp, s.c_str());
	EMObject emobj = imageio->read_attr(attr);

	H5Aclose(attr);
	H5Gclose(igrp);
	delete imageio;

	return emobj;
}

int EMUtil::write_hdf_attribute(const string & filename, const string & key, EMObject value, int image_index)
{
	ImageType image_type = get_image_type(filename);

	if (image_type != IMAGE_HDF) {
		throw ImageFormatException("This function only applies to HDF5 file.");
	}

	HdfIO2* imageio = new HdfIO2(filename, ImageIO::WRITE_ONLY);
	imageio->init();

	// Each image is in a group for later expansion. Open the group

	hid_t file = imageio->get_fileid();

	char ipath[50];
	sprintf(ipath, "/MDF/images/%d", image_index);

	hid_t igrp = H5Gopen(file, ipath);

	if (igrp<0) { // group not existed
		throw _NotExistingObjectException(string(ipath));
	}

	string s("EMAN.");
	s += key;
	int ret = imageio->write_attr(igrp, s.c_str(), value);

	H5Gclose(igrp);
	delete imageio;

	return ret;
}

int EMUtil::delete_hdf_attribute(const string & filename, const string & key, int image_index)
{
	ImageType image_type = get_image_type(filename);

	if (image_type != IMAGE_HDF) {
		throw ImageFormatException("This function only applies to HDF5 file.");
	}

	HdfIO2* imageio = new HdfIO2(filename, ImageIO::READ_WRITE);
	imageio->init();

	// Each image is in a group for later expansion. Open the group

	hid_t file = imageio->get_fileid();

	char ipath[50];
	sprintf(ipath, "/MDF/images/%d", image_index);

	hid_t igrp = H5Gopen(file,ipath);

	if (igrp < 0) { // group not existed
		throw _NotExistingObjectException(string(ipath));
	}

	string s("EMAN.");
	s += key;
	herr_t ret = H5Adelete(igrp, s.c_str());

	H5Gclose(igrp);
	delete imageio;

	if (ret >= 0) return 0;
	else return -1;
}
#endif	//EM_HDF5
