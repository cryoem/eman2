/*
 * Author: tunay.durmaz@bcm.edu, 08/20/2020
 * Copyright (c) 2020- Baylor College of Medicine
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

#include "eerio.h"

#include <tiffio.h>

using namespace EMAN;


auto Decoder::operator()(unsigned int count, unsigned int sub_pix) const {
	return std::make_pair(x(count, sub_pix), y(count, sub_pix));
}

typedef vector<pair<int, int>> COORDS;

auto decode_eer_data(EerWord *data, Decoder &decoder) {
	EerStream is((data));
	EerRle    rle;
	EerSubPix sub_pix;

	is >> rle >> sub_pix;
	int count = rle;

	COORDS coords;

	while (count < decoder.camera_size * decoder.camera_size) {
		coords.push_back(decoder(count, sub_pix));

		is >> rle >> sub_pix;

		count += rle+1;
	}

	return coords;
}

void TIFFOutputWarning(const char* module, const char* fmt, va_list ap)
{
	if(string(fmt).rfind("Unknown field with tag %d (0x%x) encountered") != string::npos)
		return;

	cout << module << ": ";
	vprintf(fmt, ap);
	cout << endl;
}

string read_acquisition_metadata(TIFF *tiff) {
	char *metadata_c = nullptr;
	uint32_t count = 0;

	TIFFSetDirectory(tiff, 0);
	TIFFGetField(tiff, 65001, &count, &metadata_c);

	return string(metadata_c, count);
}

auto read_compression(TIFF *tiff) {
	uint16_t compression = 0;

	TIFFGetField(tiff, TIFFTAG_COMPRESSION, &compression);

	return compression;
}

auto read_raw_data(TIFF *tiff) {
	auto num_strips = TIFFNumberOfStrips(tiff);
	vector<unsigned int> strip_sizes(num_strips);

	for(size_t i=0; i<num_strips; ++i)
		strip_sizes[i] = TIFFRawStripSize(tiff, i);

	std::vector<unsigned char> data;

	for(size_t i=0; i<num_strips; ++i) {
		auto prev_size = data.size();
		data.resize(prev_size + strip_sizes[i]);
		TIFFReadRawStrip(tiff, i, data.data() + prev_size, strip_sizes[i]);
	}

	return data;
}


EerIO::EerIO(const string & fname, IOMode rw, Decoder &dec)
:	ImageIO(fname, rw), decoder(dec)
{
	TIFFSetWarningHandler(TIFFOutputWarning);

	tiff_file = TIFFOpen(filename.c_str(), "r");

	acquisition_metadata = read_acquisition_metadata(tiff_file);

	for( ; TIFFReadDirectory(tiff_file); )
		num_frames = TIFFCurrentDirectory(tiff_file) + 1;

}

EerIO::~EerIO()
{
	TIFFClose(tiff_file);
}

void EerIO::init()
{
	ENTERFUNC;

	EXITFUNC;
}

int EerIO::get_nimg()
{
	return num_frames;
}

bool EerIO::is_image_big_endian()
{
	return is_big_endian;
}


int EerIO::read_header(Dict & dict, int image_index, const Region * area, bool is_3d)
{
	TIFFSetDirectory(tiff_file, image_index);

	int nx = 0;
	int ny = 0;

	TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &ny);

	dict["nx"] = decoder.num_pix();
	dict["ny"] = decoder.num_pix();
	dict["nz"] = 1;

	return 0;
}


int EerIO::write_header(const Dict & dict, int image_index, const Region* area,
						EMUtil::EMDataType filestoragetype, bool use_host_endian)
{
	ENTERFUNC;

	EXITFUNC;

	return 0;
}

int EerIO::read_data(float *rdata, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	TIFFSetDirectory(tiff_file, image_index);
	
	auto data = read_raw_data(tiff_file);

	auto coords = decode_eer_data((EerWord *) data.data(), decoder);
	for(auto &c : coords)
		rdata[c.first + c.second * decoder.num_pix()] += 1;

	EXITFUNC;

	return 0;
}

int EerIO::write_data(float *data, int image_index, const Region* area,
					  EMUtil::EMDataType, bool use_host_endian)
{
	ENTERFUNC;

	EXITFUNC;
	return 0;
}


bool EerIO::is_complex_mode()
{
	return false;
}

void EerIO::flush()
{
}

bool EerIO::is_single_image_format() const
{
	return false;
}
