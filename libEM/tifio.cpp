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

#ifdef EM_TIFF

#include "tifio.h"
#include "util.h"
#include "geometry.h"
#include <tiffio.h>
#include <climits>

using namespace EMAN;

namespace {
	/* x% weighting -> fraction of full color */
	#define PCT(x)	(((x)*255+127)/100)
	const int RED = PCT(30);		/* 30% */
	const int GREEN = PCT(59);		/* 59% */
	const int BLUE = PCT(11);		/* 11% */
}

TiffIO::TiffIO(string tiff_filename, IOMode rw)
:	filename(tiff_filename), rw_mode(rw), tiff_file(0),
	bitspersample(0), photometric(0), initialized(false),
	rendermin(0.0f), rendermax(0.0f)
{
	is_big_endian = ByteOrder::is_host_big_endian();
}


TiffIO::~TiffIO()
{
	if (tiff_file) {
		TIFFClose(tiff_file);
		tiff_file = 0;
	}
}

void TiffIO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}
	initialized = true;

	bool is_new_file = false;
	FILE *tmp_in = sfopen(filename, rw_mode, &is_new_file, true);
	if (!tmp_in) {
		throw ImageReadException(filename, "open TIFF");
	}

	if( !is_new_file ) {
		char buf[64];
		if (fread(buf, sizeof(buf), 1, tmp_in) != 1) {
			throw ImageReadException(filename, "first block");
		}

		if (!is_valid(&buf)) {
			throw ImageReadException(filename, "invalid TIFF");
		}

		if (buf[0] == TIFF_BIG_ENDIAN) {
			is_big_endian = true;
		}
		else {
			is_big_endian = false;
		}
	}

	fclose(tmp_in);
	tmp_in = 0;

	TIFFSetWarningHandler(0);

	if( rw_mode == ImageIO::READ_ONLY ) {
		tiff_file = TIFFOpen(filename.c_str(), "r");

		if (!tiff_file) {
			throw ImageReadException(filename, "open TIFF");
		}

		TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &bitspersample);

		if (bitspersample != CHAR_BIT &&
				bitspersample != (CHAR_BIT * sizeof(short)) &&
				bitspersample != (CHAR_BIT * sizeof(float)) ) {
			char desc[256];
			sprintf(desc, "invalid %d bits. only %d-bit and %d-bit TIFF are supported",
					bitspersample, CHAR_BIT, (int)(CHAR_BIT * sizeof(short)));
			throw ImageReadException(filename, desc);
		}
	}
	else {
		tiff_file = TIFFOpen(filename.c_str(), "w");
		if (!tiff_file) {
			throw ImageReadException(filename, "open TIFF");
		}
	}

	EXITFUNC;
}

bool TiffIO::is_valid(const void *first_block)
{
	ENTERFUNC;
	bool result = false;

	if (!first_block) {
		result = false;
	}
	else {
		const char *data = static_cast < const char *>(first_block);

		if ((data[0] == data[1]) && (data[0] == TIFF_LITTLE_ENDIAN || data[1] == TIFF_BIG_ENDIAN)) {
			result = true;
		}
	}
	EXITFUNC;
	return result;
}

int TiffIO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	init();
	
	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}

	if(image_index != 0) {
		throw ImageReadException(filename, "no stack allowed for MRC image. For take 2D slice out of 3D image, read the 3D image first, then use get_clip().");
	}

	int nx = 0;
	int ny = 0;
	TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &ny);

	check_region(area, IntSize(nx, ny));

	float min = 0;
	float max = 0;
	TIFFDataType data_type = TIFF_NOTYPE;
	float resolution_x = 0;
	float resolution_y = 0;

	TIFFGetField(tiff_file, TIFFTAG_MINSAMPLEVALUE, &min);
	TIFFGetField(tiff_file, TIFFTAG_MAXSAMPLEVALUE, &max);

	TIFFGetField(tiff_file, TIFFTAG_PHOTOMETRIC, &photometric);

	TIFFGetField(tiff_file, TIFFTAG_SAMPLEFORMAT, &data_type);
	TIFFGetField(tiff_file, TIFFTAG_XRESOLUTION, &resolution_x);
	TIFFGetField(tiff_file, TIFFTAG_YRESOLUTION, &resolution_y);

	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;

	dict["minimum"] = min;
	dict["maximum"] = max;

	if (bitspersample == CHAR_BIT) {
		dict["datatype"] = EMUtil::EM_UCHAR;
	}
	else if (bitspersample == sizeof(unsigned short) * CHAR_BIT) {
		dict["datatype"] = EMUtil::EM_USHORT;
	}
	else if (bitspersample == sizeof(float) * CHAR_BIT) {
		dict["datatype"] = EMUtil::EM_FLOAT;
	}

	dict["bitspersample"] = bitspersample;
	dict["resolution_x"] = resolution_x;
	dict["resolution_y"] = resolution_y;
	EXITFUNC;
	return 0;
}

int TiffIO::read_data(float *rdata, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, rdata);

	int nx = 0;
	int ny = 0;
	TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &ny);

	int err = 0;

	/* for grey scale image, use TIFFReadEncodedStrip() and TIFFReadEncodedTile()
	 * because the reading of strip image is twice time faster than TIFFReadRGBAImage -Grant*/
	if(photometric == PHOTOMETRIC_MINISWHITE || photometric == PHOTOMETRIC_MINISBLACK) {
		unsigned char *cdata;	//buffer for strip or tile

		if(TIFFIsTiled(tiff_file)) {
			tsize_t tileSize = TIFFTileSize(tiff_file);
			tsize_t tileMax = TIFFNumberOfTiles(tiff_file);
			tsize_t tileCount;

			uint32 tileWidth, tileLength;
			TIFFGetField(tiff_file, TIFFTAG_TILEWIDTH, &tileWidth);
			TIFFGetField(tiff_file, TIFFTAG_TILELENGTH, &tileLength);

			if((cdata=(unsigned char*)_TIFFmalloc(tileSize))==NULL){
				fprintf(stderr,"Error: Could not allocate enough memory\n");
				return(-1);
			}

			int tilePerLine = nx/tileWidth + 1;
			int NX, NY;	//(NX, NY) is the cordinates of tile
			int xpos, ypos; //(xpos, ypos) is the actual coordinates of pixel (j,i) in image

			for(tileCount=0; tileCount<tileMax; tileCount++) {
				if(TIFFReadEncodedTile(tiff_file, tileCount, cdata, tileSize) == -1) {
					fprintf(stderr,"Error reading tiled image\n");return(-1);
				}
				else {
					NX = tileCount%tilePerLine;
					NY = tileCount/tilePerLine;
					uint32 i, j;
					for(i=0; i<tileLength; i++) {
						for(j=0; j<tileWidth; j++) {
							xpos = NX*tileWidth + j;
							ypos = NY*tileLength + i;

							if(bitspersample == CHAR_BIT) {
								if(xpos<nx && ypos<ny) {	//discard those pixel in tile which is out of actual image's boundary
									photometric == PHOTOMETRIC_MINISWHITE ?
										rdata[nx*(ny-1)-(ypos*nx)+xpos] = -(float) ((unsigned char*)cdata)[i*tileWidth+j] :
										rdata[nx*(ny-1)-(ypos*nx)+xpos] = (float) ((unsigned char*)cdata)[i*tileWidth+j];
								}
							}
							else if(bitspersample == sizeof(unsigned short) * CHAR_BIT) {
								if(xpos<nx && ypos<ny) {	//discard those pixel in tile which is out of actual image's boundary
									photometric == PHOTOMETRIC_MINISWHITE ?
										rdata[nx*(ny-1)-(ypos*nx)+xpos] = -(float) ((unsigned short*)cdata)[i*tileWidth+j] :
										rdata[nx*(ny-1)-(ypos*nx)+xpos] = (float) ((unsigned short*)cdata)[i*tileWidth+j];
								}
							}
							else if (bitspersample == sizeof(float) * CHAR_BIT) {
								photometric == PHOTOMETRIC_MINISWHITE ?
										rdata[nx*(ny-1)-(ypos*nx)+xpos] = -(float) ((float*)cdata)[i*tileWidth+j] :
										rdata[nx*(ny-1)-(ypos*nx)+xpos] = (float) ((float*)cdata)[i*tileWidth+j];
							}
							else {
								fprintf(stderr,"BAILING OUT:Allow only 8- or 16-bits image\n");
								return(-1);
							}
						}
					}
				}
			}

		}
		else {
			check_region(area, IntSize(nx, ny));
			int xlen = 0, ylen = 0, x0 = 0, y0 = 0;
			EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);
			EMUtil::get_region_origins(area, &x0, &y0);

			int strip_size = TIFFStripSize(tiff_file);
			uint32 num_strips = TIFFNumberOfStrips(tiff_file);

			if((cdata = static_cast < unsigned char *>(_TIFFmalloc(strip_size)))==NULL) {
				fprintf(stderr,"Error: Could not allocate enough memory\n");
				return(-1);
			}

			int k = 0;
			int num_read = 0;
			int mode_size = bitspersample / CHAR_BIT;
			int total_rows = 0;

			for (uint32 i = 0; i < num_strips; i++) {
				if ((num_read = TIFFReadEncodedStrip(tiff_file, i, cdata, strip_size)) == -1) {
					LOGERR("reading stripped TiFF image '%s' failed", filename.c_str());
					err = 1;
					break;
				}

				int nitems = num_read / mode_size;
				int nrows = nitems / nx;
				total_rows += nrows;

				int y_start = 0;
				int y_end = nrows;

				if (area) {
					if (total_rows >= y0 && total_rows < y0 + nrows) {
						y_start = nrows - (total_rows - y0);
					}
					else if (total_rows >= (y0 + ylen) && total_rows < (y0 + ylen + nrows)) {
						y_end = y0 + ylen - total_rows + nrows;
					}
					else if (total_rows >= (y0 + ylen + nrows)) {
						break;
					}
				}

				for (int l = y_start; l < y_end; l++) {
					for (int j = x0; j < x0 + xlen; j++) {
						if (bitspersample == CHAR_BIT) {
							photometric == PHOTOMETRIC_MINISWHITE ?
								rdata[k] = -(float) ((unsigned char*)cdata)[l * nx + j] :
								rdata[k] = (float) ((unsigned char*)cdata)[l * nx + j];
						}
						else if (bitspersample == sizeof(unsigned short) * CHAR_BIT) {
							photometric == PHOTOMETRIC_MINISWHITE ?
								rdata[k] = -(float)((unsigned short*)cdata)[l * nx + j] :
								rdata[k] = (float)((unsigned short*)cdata)[l * nx + j];
						}
						else if (bitspersample == sizeof(float) * CHAR_BIT) {
							photometric == PHOTOMETRIC_MINISWHITE ?
								rdata[k] = -(float)((float*)cdata)[l * nx + j] :
								rdata[k] = (float)((float*)cdata)[l * nx + j];
						}
						k++;
					}
				}
			}
			Util::flip_image(rdata, xlen, ylen);
		}
		_TIFFfree(cdata);
	}
	else {	//process color image, convert to greyscale
		size_t npixels = nx * ny;
		uint32 * raster = (uint32*) _TIFFmalloc(npixels * sizeof(uint32));
		if(raster != NULL) {
			if(TIFFReadRGBAImage(tiff_file, nx, ny, raster, 0)) {
				int abgr = 0;	//raw ABGR pixel value
				int red=0, green=0, blue=0;
				for (int i=0; i<nx; ++i) {
					for (int j=0; j<ny; ++j) {
						abgr 	= raster[j+ny*i];
						red 	= TIFFGetR(abgr);
						green 	= TIFFGetG(abgr);
						blue 	= TIFFGetB(abgr);
						rdata[j+ny*i] = static_cast<float>(red*RED+green*GREEN+blue*BLUE);
					}
				}
				_TIFFfree(raster);
			}
		}
	}

	EXITFUNC;
	return err;
}


int TiffIO::write_header(const Dict & dict, int image_index, const Region*, EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	//support single image TIFF only, index must be zero
	if(image_index != 0) {
		throw ImageWriteException(filename, "MRC file does not support stack.");
	}
	check_write_access(rw_mode, image_index);

	nx = (unsigned int) (int) dict["nx"];
	ny = (unsigned int) (int) dict["ny"];
	nz = (unsigned int) (int)dict["nz"];
	if (nz != 1) {
		LOGERR("Only support 2D TIFF file write");
		return 1;
	}

	EMUtil::EMDataType datatype = (EMUtil::EMDataType) (int) dict["datatype"];
	if (datatype == EMUtil::EM_UCHAR) {
		bitspersample = CHAR_BIT;
	}
	else if(datatype == EMUtil::EM_USHORT) {
		bitspersample = CHAR_BIT * sizeof(short);
	}
	else if(datatype == EMUtil::EM_FLOAT) {
		bitspersample = CHAR_BIT * sizeof(float);
	}
	else {
		LOGWARN("Don't support data type '%s' in TIFF. Convert to '%s'.",
				EMUtil::get_datatype_string(datatype),
				EMUtil::get_datatype_string(EMUtil::EM_USHORT));
		bitspersample = CHAR_BIT * sizeof(short);
	}
	TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(tiff_file, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, nx);
	TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, ny);
	TIFFSetField(tiff_file, TIFFTAG_ROWSPERSTRIP, ny);
	TIFFSetField(tiff_file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tiff_file, TIFFTAG_SOFTWARE, "EMAN2" );
	//TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, NO_COMPRESSION);
	//TIFFSetField(tiff_file, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);

	EXITFUNC;
	return 0;
}

int TiffIO::write_data(float * data, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;

	// If we didn't get any parameters in 'render_min' or 'render_max', we need to find some good ones
	getRenderMinMax(data, nx, ny, rendermin, rendermax);

	if(bitspersample == CHAR_BIT) {
		unsigned char *cdata = new unsigned char[nx*ny];

		int src_idx, dst_idx;
		for (unsigned int i = 0; i < ny; ++i) {
			for (unsigned int j = 0; j < nx; ++j) {
				src_idx = i*nx+j;
				dst_idx = nx*(ny-1) - (i*nx) +j;
				if(data[src_idx] < rendermin){
					cdata[dst_idx] = 0;
				}
				else if(data[src_idx] > rendermax) {
					cdata[dst_idx] = UCHAR_MAX;
				}
				else {
					cdata[dst_idx] = (unsigned char)((data[src_idx] - rendermin) / (rendermax - rendermin) * UCHAR_MAX);
				}
			}
		}

		if(TIFFWriteEncodedStrip(tiff_file, 0, cdata, nx*ny) == -1)
			{ printf("Fail to write tiff file.\n"); return -1; }

		if( cdata )	{ delete[] cdata; cdata = 0; }
	}
	else if(bitspersample == CHAR_BIT*sizeof(short)) {
		unsigned short *sdata = new unsigned short[nx*ny];

		int src_idx, dst_idx;
		for (unsigned int i = 0; i < ny; ++i) {
			for (unsigned int j = 0; j < nx; ++j) {
				src_idx = i*nx+j;
				dst_idx = nx*(ny-1) - (i*nx) +j;
				if(data[src_idx] < rendermin){
					sdata[dst_idx] = 0;
				}
				else if(data[src_idx] > rendermax) {
					sdata[dst_idx] = USHRT_MAX;
				}
				else {
					sdata[dst_idx] = (unsigned short)((data[src_idx] - rendermin) / (rendermax - rendermin) * USHRT_MAX);
				}
			}
		}

		if(TIFFWriteEncodedStrip(tiff_file, 0, sdata, nx*ny*sizeof(short)) == -1)
			{ printf("Fail to write tiff file.\n"); return -1; }

		if( sdata )	{ delete[] sdata; sdata = 0; }
	}
	else if(bitspersample == CHAR_BIT*sizeof(float)) {
		float *fdata = new float[nx*ny];

		int src_idx, dst_idx;
		for (unsigned int i = 0; i < ny; ++i) {
			for (unsigned int j = 0; j < nx; ++j) {
				src_idx = i*nx+j;
				dst_idx = nx*(ny-1) - (i*nx) +j;
				fdata[dst_idx] = data[src_idx];
			}
		}

		if(TIFFWriteEncodedStrip(tiff_file, 0, fdata, nx*ny*sizeof(float)) == -1)
					{ printf("Fail to write tiff file.\n"); return -1; }

		if( fdata )	{ delete[] fdata; fdata = 0; }
	}
	else {
		LOGWARN("TIFF in EMAN2 only support data type 8 bit, 16 bit or 32 bit.");
	}

	EXITFUNC;
	return 0;
}

void TiffIO::flush()
{
	TIFFFlush(tiff_file);
}

bool TiffIO::is_complex_mode()
{
	return false;
}

bool TiffIO::is_image_big_endian()
{
	init();
	return is_big_endian;
}


#endif	//EM_TIFF
