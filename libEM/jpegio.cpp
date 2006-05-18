/**
 * $Id$
 */
#ifdef EM_PNG

#include <stdlib.h>
#include <stdio.h>
#include "jpegio.h"
#include "geometry.h"
#include "util.h"


using namespace EMAN;

JpegIO::JpegIO(const string & file, IOMode rw):	filename(file), rw_mode(rw), jpeg_file(0), initialized(false),rendermin(0),rendermax(0)
{}

JpegIO::~JpegIO()
{
	if (jpeg_file) {
		fclose(jpeg_file);
		jpeg_file = 0;
	}
	
}

void JpegIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}

	initialized = true;

	bool is_new_file = false;
	jpeg_file = sfopen(filename, rw_mode, &is_new_file, true);

	if (!is_new_file) {
		throw ImageReadException(filename, "JPEG reading not supported");
	}
	else {
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_compress(&cinfo);
		jpeg_stdio_dest(&cinfo, jpeg_file);
	}

	EXITFUNC;
}

bool JpegIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	return false;
}

int JpegIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	throw ImageReadException(filename, "JPEG reading not supported");

	EXITFUNC;
	return 0;
}

int JpegIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	check_write_access(rw_mode, image_index);

	if ((int) dict["nz"] != 1) {
		LOGERR("Only support 2D JPEG file write");
		return 1;
	}

	rendermin=(float)dict["render_min"];	// float value representing black in the output
	rendermax=(float)dict["render_max"];	// float value representign white in the output

	cinfo.image_width = (int) dict["nx"];      /* image width and height, in pixels */
	cinfo.image_height = (int) dict["ny"];
	cinfo.input_components = 1;     /* # of color components per pixel */
	cinfo.in_color_space = JCS_GRAYSCALE; /* colorspace of input image */
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo,75,true);

	EXITFUNC;
	return 0;
}

int JpegIO::read_data(float *data, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	EXITFUNC;
	return 0;
}

int JpegIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	if (image_index>0) throw ImageReadException("N/A", "JPEG files are single-image only");
	if (area->size[0]!=cinfo.image_width || area->size[1]!=cinfo.image_height)
		 throw ImageReadException("N/A", "No region writing for JPEG images");
	int nx=(int)area->size[0],ny=(int)area->size[1];

	// If we didn't get any parameters in 'render_min' or 'render_max', we need to find some good ones
	if (rendermax<=rendermin) {
		float m=0,s=0;
		
		for (int i=0; i<nx*ny; i++) { m+=data[i]; s+=data[i]*data[i]; }
		m/=(float)(nx*ny);
		s=sqrt(m/(float)(nx*ny)-m*m);
		if (s<=0) s=1.0;	// this means all data values are the same
		rendermin=m-s*3.0;
		rendermax=m+s*3.0;
	}

	unsigned char *cdata=(unsigned char *)malloc(nx+1);

	// convert and write the data 1 scanline at a time
	JSAMPROW rp[1];
	rp[0]=cdata;
	jpeg_start_compress(&cinfo, TRUE);
	for (int i=0; i<ny; i++) {
		for (int j=0; j<nx; j++) {
			if (data[i*nx+j]<=rendermin) cdata[j]=0;
			else if (data[i*nx+j]>=rendermax) cdata[j]=255;
			else cdata[j]=(int)((data[i*nx+j]-rendermin)/(rendermax-rendermin)*256.0);
		}
		jpeg_write_scanlines(&cinfo, rp, 1);
	}	

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	free(cdata);
	EXITFUNC;
	return 0;
}

void JpegIO::flush()
{

}

bool JpegIO::is_complex_mode()
{
	return false;
}

bool JpegIO::is_image_big_endian()
{
	return true;
}


#endif	//EM_PNG
