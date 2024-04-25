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

#include <algorithm>
#include <tiffio.h>
#include <boost/property_tree/xml_parser.hpp>


using boost::property_tree::ptree;
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
#ifdef DEBUG
	cout << module << ": ";
	vprintf(fmt, ap);
	cout << endl;
#endif  //DEBUG
}

string read_acquisition_metadata(TIFF *tiff) {
	char *metadata_c = nullptr;
	uint32_t count = 0;

	TIFFSetDirectory(tiff, 0);
	TIFFGetField(tiff, 65001, &count, &metadata_c);

	return string(metadata_c, count);
}

string to_snake_case(const string &s) {
	auto ret(s);
	int sh = 0;
	for(int i=0; i<s.size(); i++) {
		if(isupper(s[i])) {
			if(!isupper(s[i-1])) {
				ret.insert(i + sh, "_");
				sh++;
			}
			ret[i + sh] = ::tolower(s[i]);
		}
	}
	
	return ret;
}

Dict parse_acquisition_data(string metadata) {
	Dict dict;
	std::istringstream ins(metadata);
	ptree pt;

	read_xml(ins, pt);

	for(auto &v : pt.get_child("metadata")) {
		auto xml_item_name = v.second.get_child("<xmlattr>.name").data();
		auto xml_item_value = v.second.data();

		auto key = "EER." + to_snake_case(xml_item_name);

		dict[key] = xml_item_value;
	}
	
	return dict;
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

	auto acquisition_metadata = read_acquisition_metadata(tiff_file);
	acquisition_data_dict = parse_acquisition_data(acquisition_metadata);

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

	dict["EER.compression"] = read_compression(tiff_file);

	for(auto &d : acquisition_data_dict)
		dict[d.first] = d.second;

	if(auto it = dict.find("EER.sensor_pixel_size.width"); it != dict.end()) {
		if ((float)it->second == 0.0)
			dict["apix_x"] = 1.0f;
		else
			dict["apix_x"] = (float)it->second * (float)1.0e10;
	}

	if(auto it = dict.find("EER.sensor_pixel_size.height"); it != dict.end()) {
		if ((float)it->second == 0.0)
			dict["apix_y"] = 1.0f;
		else
			dict["apix_y"] = (float)it->second * (float)1.0e10;
	}

	dict["apix_z"] = dict["apix_x"];

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

	std::fill(rdata, rdata + decoder.camera_size * decoder.camera_size, 0);

	auto coords = decode_eer_data((EerWord *) data.data(), decoder);
	for(const auto &c : coords)
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
