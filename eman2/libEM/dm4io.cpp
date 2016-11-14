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

#include "dm4io.h"
#include "portable_fileio.h"
#include "geometry.h"

using namespace EMAN;
using namespace EMAN::GatanDM4;

const char *TagTable::IMAGE_WIDTH_TAG = "Dimensions #0";
const char *TagTable::IMAGE_HEIGHT_TAG = "Dimensions #1";
const char *TagTable::IMAGE_NIMG_TAG = "Dimensions #2";
const char *TagTable::IMAGE_DATATYPE_TAG = "DataType";
const char *TagTable::IMAGE_THUMB_INDEX_TAG = "ImageIndex";


TagTable::TagTable()
	:	img_index(0), is_big_endian(true), img_counted(1), num_images_found(0)
{
}

TagTable::~TagTable()
{
	for (unsigned int i = 0; i < data_list.size(); i++) {
		if (data_list[i]) {
			delete[]data_list[i];
			data_list[i] = 0;
		}
	}
}

void TagTable::add(const string & name, const string & value)
{
	const char *value_str = value.c_str();

	if (name == IMAGE_WIDTH_TAG) {
		x_list.push_back(atoi(value_str));
	}
	else if (name == IMAGE_HEIGHT_TAG) {
		y_list.push_back(atoi(value_str));
	}
	else if (name == IMAGE_DATATYPE_TAG) {
		datatype_list.push_back(atoi(value_str));
	}
	else if (name == IMAGE_THUMB_INDEX_TAG) {
		set_thumb_index(atoi(value_str));
	}
	else if(name == IMAGE_NIMG_TAG){
		img_counted=atoi(value_str);
	}
	else {
		tags[name] = value;
	}
}

void TagTable::add_data(char *data)
{
	if (!data) {
		throw NullPointerException("DM4 data is NULL");
	}
	else {
		data_list.push_back(data);
	}
}

string TagTable::get_string(const string & name)
{
	return tags[name];
}

int TagTable::get_int(const string & name)
{
	return atoi(tags[name].c_str());
}

float TagTable::get_float(const string & name)
{
	return static_cast < float >(atof(tags[name].c_str()));
}

double TagTable::get_double(const string & name)
{
	return atof(tags[name].c_str());
}

void TagTable::dump() const
{
	map < string, string >::const_iterator p;

	for (p = tags.begin(); p != tags.end(); p++) {
//		LOGDEBUG("  %s: %s", (*p).first.c_str(), (*p).second.c_str());
		printf("  %s: %s\n", (*p).first.c_str(), (*p).second.c_str());
	}
}

int TagTable::get_xsize() const
{
	return x_list[img_index];
}

int TagTable::get_image_counted() const
{
	return img_counted;
}

int TagTable::get_ysize() const
{
	return y_list[img_index];
}

int TagTable::get_datatype() const
{
	return datatype_list[img_index];
}

char *TagTable::get_data() const
{
	return data_list[img_index];
}

void TagTable::set_thumb_index(int i)
{
	if (i != 0 && i != 1) {
		throw OutofRangeException(0, 1, i, "image index");
	}
	else {
		if (i == 0) {
			img_index = 1;
		}
		else {
			img_index = 0;
		}
	}
}

int TagTable::get_num_images_found() const
{
	return num_images_found;
}

void TagTable::set_num_images_found(int num_found)
{
	num_images_found = num_found;
}

/////////////////////////////////////////////////////////////////////////////////////////////

TagData::TagData(FILE * data_file, TagTable * table, const string & tagname)
	:	in(data_file), tagtable(table), name(tagname), tag_type(UNKNOWN)
{
}

TagData::~TagData()
{
}

string TagData::read_native(bool is_value_stored)
{
	size_t sz = typesize();
	size_t nr;
	char val_str[32];

	if (tag_type == SHORT) {
		short val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == USHORT) {
		unsigned short val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == INT) {
		int val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == CHAR || tag_type == OCTET) {
		char val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == BOOLEAN) {
		bool val = false;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == UINT) {
		unsigned int val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%u", (int) val);
	}
	else if (tag_type == FLOAT) {
		float val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%f", val);
	}
	else if (tag_type == DOUBLE) {
		double val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%10e", val);
	}
	else if (tag_type == OCTEU) {
		long long val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%lld", val);
	}
	else if (tag_type == OCTEV) {
		unsigned long long val = 0;
		nr = fread(&val, sz, 1, in); nr++;
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%lld", val);
	}
	else {
		LOGERR("invalid tag type: '%lld'", tag_type);
		exit(1);
	}

	if (is_value_stored) {
		tagtable->add(name, val_str);
	}

	LOGVAR("value = '%s'", val_str);

	return string(val_str);
}

vector < int >TagData::read_array_types()
{
	LOGVAR("TagData::read_array_types()");

	long long array_type = 0;
	size_t nr;
	nr = fread(&array_type, sizeof(array_type), 1, in); nr++;

	ByteOrder::become_big_endian(&array_type);

	LOGVAR("array data type = '%s'", GatanDM4::to_str((Type) array_type));

	vector < int >item_types;

	if (array_type == STRUCT) {
		item_types = read_struct_types();
	}
	else if (array_type == ARRAY) {
		item_types = read_array_types();
		LOGERR("DM4: don't know how to handle this array type");
	}
	else {
		item_types.push_back(array_type);
	}

	return item_types;
}

// string tag data are stored in unicode
string TagData::read_string(int size)
{
	if (size <= 0) {
		return string("");
	}

	unsigned short *buf = new unsigned short[size];
	char *str = new char[size + 1];

	size_t nr;
	nr = fread(buf, size * sizeof(unsigned short), 1, in); nr++;
	tagtable->become_host_endian < unsigned short >(buf, size);

	for (int i = 0; i < size; i++) {
		str[i] = static_cast < char >(buf[i]);
	}

	str[size] = '\0';
	string str1 = string(str);

	if (str) {
		delete [] str;
		str = NULL;
	}

	if (buf) {
		delete [] buf;
		buf = NULL;
	}

	return str1;
}

int TagData::read_array_data(vector < int >item_types, bool nodata, int image_index, int num_images)
{
	ENTERFUNC;
	if (item_types.size() == 0) {
		LOGERR("DM4 item types cannot be empty");
		return 1;
	}

	int err = 0;
	long long array_size = 0;

	size_t nr;
	nr = fread(&array_size, sizeof(array_size), 1, in); nr++;
	ByteOrder::become_big_endian(&array_size);

	LOGVAR("array size = %lld\n", array_size);

	size_t item_size = 0;
	for (size_t i = 0; i < item_types.size(); i++) {
		item_size += typesize(item_types[i]);
	}

	LOGVAR("%s array item size = %lld\n", name.c_str(), item_size);

	size_t buf_size = item_size * array_size;

	if (item_types.size() == 1 && item_types[0] == USHORT && nodata) {
		string val = read_string(array_size);
		tagtable->add(name, val);
		LOGVAR("value: %s", val.c_str());
	}
	else if (!nodata && name == "Data") {
		int num_found = tagtable->get_num_images_found();
		num_found++;
		tagtable->set_num_images_found(num_found);

		char * data;

		if (image_index < 0  ||  buf_size % num_images != 0  ||  num_found == 1) {
			data = new char[buf_size];
			nr = fread(data, item_size, array_size, in); nr++;
		}
		else {
			size_t image_size = buf_size / num_images;

			data = new char[image_size];
			portable_fseek(in, image_index * image_size, SEEK_CUR);
			nr = fread(data, image_size, 1, in); nr++;
			portable_fseek(in, (num_images - image_index - 1) * image_size, SEEK_CUR);
			array_size = array_size / num_images;
		}

		if (item_size == sizeof(short)) {
			tagtable->become_host_endian((short *) data, array_size);
		}
		else if (item_size == sizeof(int)) {
			tagtable->become_host_endian((int *) data, array_size);
		}
		else if (item_size == sizeof(double)) {
			tagtable->become_host_endian((double *) data, array_size);
		}
		else {
			LOGERR("cannot handle this type of DM4 image data");
			return 1;
		}

		tagtable->add_data(data);
	}
	else {
		portable_fseek(in, buf_size, SEEK_CUR);
	}
	EXITFUNC;
	return err;
}

vector < int >TagData::read_struct_types()
{
	LOGVAR("TagData::read_struct_types()");

	long long namelength = 0;
	long long nfields = 0;

	size_t nr;
	nr = fread(&namelength, sizeof(namelength), 1, in); nr++;
	ByteOrder::become_big_endian(&namelength);

	nr = fread(&nfields, sizeof(nfields), 1, in); nr++;
	ByteOrder::become_big_endian(&nfields);

	LOGVAR("namelength = %lld\n", namelength);
	LOGVAR("num fields = %lld\n", nfields);

	vector < int >field_types;

	for (unsigned int i = 0; i < nfields; i++) {
		nr = fread(&namelength, sizeof(namelength), 1, in); nr++;
		ByteOrder::become_big_endian(&namelength);

		long long field_type = 0;
		nr = fread(&field_type, sizeof(field_type), 1, in); nr++;
		ByteOrder::become_big_endian(&field_type);
		
		LOGVAR("%dth namelength = %lld, type = '%s'",
			   i, namelength, GatanDM4::to_str((Type) field_type));
		field_types.push_back(field_type);
	}

	return field_types;
}

int TagData::read_any(bool nodata, int image_index, int num_images)
{
	int err = 0;

	size_t nr;
	nr = fread(&tag_type, sizeof(tag_type), 1, in); nr++;
	
	ByteOrder::become_big_endian(&tag_type);
	LOGVAR("TagData::read_any tag type = '%s'\n", GatanDM4::to_str((Type) tag_type));


	if (tag_type == ARRAY) {
		vector < int >item_types = read_array_types();
		err = read_array_data(item_types, nodata, image_index, num_images);
	}
	else if (tag_type == STRUCT) {
		vector < int >field_types = read_struct_types();

		for (unsigned int i = 0; i < field_types.size(); i++) {
			
			tag_type = static_cast < Type > (field_types[i]);
			
			string val = read_native(false);
			
			char int_str[32];
			sprintf(int_str, " #%d", i);
			string fieldname = name + string(int_str);
			tagtable->add(fieldname, val);
		}
	}
	else if (tag_type == STRING) {

		int str_sz = 0;
		size_t nr;
		nr = fread(&str_sz, sizeof(str_sz), 1, in); nr++;
		ByteOrder::become_big_endian(&str_sz);

		char *val = new char[str_sz + 1];
		nr = fread(val, str_sz, 1, in); nr++;
		val[str_sz] = '\0';
		string val_str = string(val);

		if (val) {
			delete [] val;
			val = NULL;
		}

		tagtable->add(name, val_str);
	}
	else {
		read_native(true);
	}

	return err;
}

int TagData::read_tag_data(bool nodata, int image_index, int num_images)
{
	LOGVAR("TagData::read_tag_data()");
	int err = 0;

	const char *DATA_TYPE_MARK = "%%%%";
	const size_t mark_sz = strlen(DATA_TYPE_MARK);
	char *mark = new char[mark_sz + 1];

	long long interval;
	
	size_t nr;
	nr = fread(&interval, sizeof(interval), 1, in); nr++;

	ByteOrder::become_big_endian(&interval);

	nr = fread(mark, mark_sz, 1, in); nr++;
	mark[mark_sz] = '\0';

	if (strcmp(mark, DATA_TYPE_MARK) != 0) {
		LOGERR("data type label has been changed from '%s' to '%s'",
			   DATA_TYPE_MARK, mark);
		return 1;
	}

	if (mark) {
		delete [] mark;
		mark = NULL;
	}

	long long encoded_types_size = 0;
	nr = fread(&encoded_types_size, sizeof(long long), 1, in); nr++;
	ByteOrder::become_big_endian(&encoded_types_size);

	LOGVAR("encoded types size = %lld\n", encoded_types_size);

	err = read_any(nodata, image_index, num_images);

	return err;
}

size_t TagData::typesize() const
{
	return typesize((int) tag_type);
}

size_t TagData::typesize(int t) const
{
	size_t size = 0;
	Type type = static_cast < Type > (t);

	
	switch (type) {
	case SHORT:
		size = sizeof(short);
		break;
	case USHORT:
		size = sizeof(unsigned short);
		break;
	case INT:
		size = sizeof(int);
		break;
	case UINT:
		size = sizeof(unsigned int);
		break;
	case FLOAT:
		size = sizeof(float);
		break;
	case DOUBLE:
		size = sizeof(double);
		break;
	case BOOLEAN:
		size = sizeof(bool);
		break;
	case CHAR:
	case OCTET:
		size = sizeof(char);
		break;
	case OCTEU:
		size = sizeof(double);
		break;
	case OCTEV:
		size = sizeof(double);
		break;
	default:
		LOGERR("no such type: '%d'\n", type);
		break;
	}

	return size;
}

/////////////////////////////////////////////////////////////////////////////////////////////

TagEntry::TagEntry(FILE * data_file, TagTable * table, TagGroup * parent)
	:	in(data_file), tagtable(table), parent_group(parent), name("")
{
}

TagEntry::~TagEntry()
{
}

int TagEntry::read_tag_entry(bool nodata, int image_index, int num_images)
{
	LOGVAR("TagEntry::read_tag_entry()");
	int err = 0;
	long long pos = 0;
	char tagtype = 0;
	char *tmp_name = 0;

	pos = ftell(in);
	size_t nr;
	nr = fread(&tagtype, sizeof(char), 1, in); nr++;

	if (tagtype != GROUP_TAG && tagtype != DATA_TAG) {
		portable_fseek(in, sizeof(char) * 7, SEEK_CUR);
		nr = fread(&tagtype, sizeof(char), 1, in); nr++;
	}

	if (tagtype != GROUP_TAG && tagtype != DATA_TAG) {
		tagtype = fgetc(in);
		if (tagtype == EOF)
		{
			return 1;
		}
		else{
			LOGERR("TagEntry::read_tag_entry() invalid tag type: %d @ position %lld", tagtype, pos);
			return 1;
		}
	}

	short name_len = 0;
	nr = fread(&name_len, sizeof(short), 1, in); nr++;

	ByteOrder::become_big_endian(&name_len);

	if (name_len != 0) {
		tmp_name = new char[name_len + 1];
		nr = fread(tmp_name, name_len, 1, in); nr++;
		tmp_name[name_len] = '\0';
	}
	else {
		string parent_name = parent_group->get_name();
		name_len = static_cast < short >(parent_name.size() + 4);
		tmp_name = new char[name_len + 1];
		sprintf(tmp_name, "%s #%d", parent_name.c_str(), parent_group->get_entry_id());
	}

	name = string(tmp_name);

	if (tmp_name) {
		delete [] tmp_name;
		tmp_name = NULL;
	}

	LOGVAR("\ntag name: '%s', len: %d, type: '%s'",
		   name.c_str(), name_len, GatanDM4::to_str((EntryType) tagtype));

	if (tagtype == DATA_TAG) {
		TagData tag_data(in, tagtable, name);
		err = tag_data.read_tag_data(nodata, image_index, num_images);
	}
	else if (tagtype == GROUP_TAG) {
		long long tot_size = 0;	//size of DataType record + size of data
		nr = fread(&tot_size, sizeof(long long), 1, in); nr++;
		ByteOrder::become_big_endian(&tot_size);

		TagGroup group(in, tagtable, name);
		err = group.read_tag_group(nodata, image_index, num_images);
	}

/*
	long long tot_size = 0;	//size of DataType record + size of data
	nr = fread(&tot_size, sizeof(long long), 1, in); nr++;
*/
	return err;
}

/////////////////////////////////////////////////////////////////////////////////////////////

TagGroup::TagGroup(FILE * data_file, TagTable * table, const string & groupname)
	:	in(data_file), tagtable(table), name(groupname), entry_id(0)
{
}

TagGroup::~TagGroup()
{
}

int TagGroup::read_tag_group(bool nodata, int image_index, int num_images)
{
	LOGVAR("TagGroup::read_tag_group()");
	char is_sorted, is_open;

	long long ntags = 0;
	
//	portable_fseek(in, sizeof(char) * 2, SEEK_CUR);
	size_t nr;
	nr = fread(&is_sorted, sizeof(is_sorted), 1, in); nr++;
	nr = fread(&is_open,   sizeof(is_open),   1, in); nr++;

	nr = fread(&ntags, sizeof(ntags), 1, in); nr++;
	
	ByteOrder::become_big_endian(&ntags);

	LOGVAR("DM4: ntags = %d\n", ntags);
	
	int err = 0;
	// char flagend;
	for (int i = 0; i < ntags; i++) {
		/*
		portable_fseek(in, sizeof(char) * 9, SEEK_CUR);
		nr = fread(&flagend, sizeof(char), 1, in); nr++;
		
		if (flagend ==EOF){
			break;
		}
		else{
			portable_fseek(in, -sizeof(char) * 10, SEEK_CUR);;
		}
		
		*/
		
		TagEntry tag_entry(in, tagtable, this);
		err = tag_entry.read_tag_entry(nodata, image_index, num_images);

		if (err) {
			break;
		}
	}

	return err;
}

string TagGroup::get_name() const
{
	return name;
}

int TagGroup::get_entry_id()
{
	int id = entry_id;
	entry_id++;
	return id;
}

/////////////////////////////////////////////////////////////////////////////////////////////

DM4IO::DM4IO(const string & dm4_filename, IOMode rw)
	:	filename(dm4_filename), rw_mode(rw), dm4file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	tagtable = new TagTable();
}

DM4IO::~DM4IO()
{
	if (dm4file) {
		fclose(dm4file);
		dm4file = 0;
	}
	if (tagtable) {
		delete tagtable;
		tagtable = 0;
	}
}

void DM4IO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}
	initialized = true;

	if (rw_mode != READ_ONLY) {
		throw ImageReadException(filename, "only support DM4 read-only");
	}

	dm4file = sfopen(filename, READ_ONLY);

	int stream_version = 0;
	if (fread(&stream_version, sizeof(stream_version), 1, dm4file) != 1) {
		throw ImageReadException(filename, "read stream version of DM4 file");
	}

	long long recsize;
	if (fread(&recsize, sizeof(recsize), 1, dm4file) != 1) {
		throw ImageReadException(filename, "read size of TagGroup recoed of DM4 file");
	}

	int endianness;
	if (fread(&endianness, sizeof(endianness), 1, dm4file) != 1) {
		throw ImageReadException(filename, "read endianness indicator of DM4 file");
	}

	if (ByteOrder::is_data_big_endian(&stream_version) != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&stream_version);
		ByteOrder::swap_bytes(&recsize);
		ByteOrder::swap_bytes(&endianness);
	}

	if (endianness == 0) {
		is_big_endian = true;
	}
	else {
		is_big_endian = false;
	}

	tagtable->set_endian(is_big_endian);
	
	LOGDEBUG("dm4 ver = %d, image size = %d, is_big_endian = %d",
			stream_version, recsize, (int) is_big_endian);

	EXITFUNC;
}

bool DM4IO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const int *data = static_cast < const int *>(first_block);

	int img_ver = data[0];
	int byte_order = data[3];
	if (ByteOrder::is_data_big_endian(&img_ver) != ByteOrder::is_host_big_endian()) {
		ByteOrder::swap_bytes(&img_ver);
		ByteOrder::swap_bytes(&byte_order);
	}

	if (img_ver != 4) {
		return false;
	}

	if (byte_order != 0 && byte_order != 1) {
		return false;
	}

	return true;
}

bool DM4IO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int DM4IO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	int err = 0;

	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}
	image_index = 0;
	check_read_access(image_index);

	portable_fseek(dm4file, NUM_ID_INT * sizeof(int), SEEK_SET);
	TagGroup root_group(dm4file, tagtable, "");
	root_group.read_tag_group(true, 0, 1);

	int nx = tagtable->get_xsize();
	int ny = tagtable->get_ysize();

	check_region(area, IntSize(nx, ny));
	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;

	dict["nimg"] = tagtable->get_image_counted();
	dict["DM4.acq_date"] = tagtable->get_string("Acquisition Date");
	dict["DM4.acq_time"] = tagtable->get_string("Acquisition Time");
	dict["DM4.source"] = tagtable->get_string("Source");
	dict["DM4.exposure_number"] = tagtable->get_int("Exposure Number");
	dict["DM4.exposure_time"] = tagtable->get_double("Exposure (s)");
	dict["DM4.zoom"] = tagtable->get_double("Zoom");
	dict["DM4.antiblooming"] = tagtable->get_int("Antiblooming");
	dict["DM4.indicated_mag"] = tagtable->get_double("Indicated Magnification");
	dict["DM4.actual_mag"] = tagtable->get_double("Actual Magnification");
	dict["DM4.pixel_size"] = tagtable->get_double("Pixel Size (um) #0");
	dict["DM4.name"] = tagtable->get_string("Name");

	dict["DM4.voltage"] = tagtable->get_double("Voltage")/1000.0;
	dict["microscope_voltage"]=(float)dict["DM4.voltage"];
	dict["DM4.cs"] = tagtable->get_double("Cs(mm)");
	dict["microscope_cs"]=(float)dict["DM4.cs"];

	dict["DM4.frame_type"] = tagtable->get_string("Processing");
	dict["DM4.camera_x"] = tagtable->get_int("Active Size (pixels) #0");
	dict["DM4.camera_y"] = tagtable->get_int("Active Size (pixels) #1");
	dict["DM4.binning_x"] = tagtable->get_int("Binning #0");
	dict["DM4.binning_y"] = tagtable->get_int("Binning #1");
	dict["datatype"] = to_em_datatype(tagtable->get_datatype());

	if ((float)dict["DM4.actual_mag"] >0.0) {
		float apix=10000.0*(float)dict["DM4.pixel_size"]/(float)dict["DM4.actual_mag"];
		dict["apix_x"]=apix;
		dict["apix_y"]=apix;
		dict["apix_z"]=apix;
	}

	EXITFUNC;
	return err;
}

int DM4IO::read_data(float *rdata, int image_index, const Region * area, bool)
{
	ENTERFUNC;

	check_read_access(image_index, rdata);

	portable_fseek(dm4file, NUM_ID_INT * sizeof(int), SEEK_SET);
	TagGroup root_group(dm4file, tagtable, "");
	root_group.read_tag_group(true, 0, 1);

	int nx = tagtable->get_xsize();
	int ny = tagtable->get_ysize();
	int num_images = tagtable->get_image_counted();

	portable_fseek(dm4file, NUM_ID_INT * sizeof(int), SEEK_SET);
	root_group.read_tag_group(false, image_index, num_images);

	check_region(area, IntSize(nx, ny));

	int xlen = 0, ylen = 0, x0 = 0, y0 = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);
	EMUtil::get_region_origins(area, &x0, &y0);

	char *data = tagtable->get_data();
	int data_type = tagtable->get_datatype();

	long long offset = 0;
	long long k = 0;
	long long off = 0;
	int xlast = x0 + xlen;
	int ylast = y0 + ylen;

//	bool flip_vertical = (getenv ("FLIP_VERTICAL") != NULL);
	bool flip_vertical = true;

	for (int iy = y0; iy < ylast; iy++) {
		if (flip_vertical) {
			off = offset + (y0 + ylast - 1 - iy) * nx;
		}
		else {
			off = offset + iy * nx;
		}

		switch (data_type) {
		case GatanDM4::DataType::SIGNED_INT8_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((char *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::UNSIGNED_INT8_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((unsigned char *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::SIGNED_INT16_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((short *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::UNSIGNED_INT16_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((unsigned short *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::SIGNED_INT32_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((int *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::UNSIGNED_INT32_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((unsigned int *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::REAL4_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((float *) data)[off + ix];
				k++;
			}
			break;
		case GatanDM4::DataType::REAL8_DATA:
			for (int ix = x0; ix < xlast; ix++) {
				rdata[k] = (float) ((double *) data)[off + ix];
				k++;
			}
			break;				
		default:
			string desc = string("unsupported DM4 data type") +
				string(GatanDM4::to_str((GatanDM4::DataType::GatanDataType) data_type));
			throw ImageReadException(filename, desc);
			k += xlen;
		}
	}

// #include "debug_read_data.h"

	EXITFUNC;
	return 0;

	EXITFUNC;
}

bool DM4IO::is_complex_mode()
{
	return false;
}

int DM4IO::write_header(const Dict &, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	LOGWARN("DM4 write is not supported.");
	EXITFUNC;
	return 1;
}

int DM4IO::write_data(float *, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	LOGWARN("DM4 write is not supported.");
	EXITFUNC;
	return 1;
}

int DM4IO::get_nimg()
{
	init();

	TagGroup root_group(dm4file, tagtable, "");
	root_group.read_tag_group(true, 0, 1);

	return tagtable->get_image_counted();
}

void DM4IO::flush()
{
}

int GatanDM4::to_em_datatype(int gatan_datatype)
{
	DataType::GatanDataType type = static_cast < DataType::GatanDataType > (gatan_datatype);
	int t = 0;

	switch (type) {
	case GatanDM4::DataType::SIGNED_INT8_DATA:
		t = EMUtil::EM_CHAR;
		break;
	case GatanDM4::DataType::UNSIGNED_INT8_DATA:
		t = EMUtil::EM_UCHAR;
		break;
	case GatanDM4::DataType::SIGNED_INT16_DATA:
		t = EMUtil::EM_SHORT;
		break;
	case GatanDM4::DataType::UNSIGNED_INT16_DATA:
		t = EMUtil::EM_USHORT;
		break;
	case GatanDM4::DataType::SIGNED_INT32_DATA:
		t = EMUtil::EM_INT;
		break;
	case GatanDM4::DataType::UNSIGNED_INT32_DATA:
		t = EMUtil::EM_UINT;
		break;
	default:
		t = EMUtil::EM_UNKNOWN;
	}

	return t;
}


const char *GatanDM4::to_str(TagData::Type type)
{
	const char *str = "unknown";

	switch (type) {
	case TagData::SHORT:
		str = "short";
		break;
	case TagData::INT:
		str = "int";
		break;
	case TagData::USHORT:
		str = "unsigned short";
		break;
	case TagData::UINT:
		str = "unsigned int";
		break;
	case TagData::FLOAT:
		str = "float";
		break;
	case TagData::DOUBLE:
		str = "double";
		break;
	case TagData::BOOLEAN:
		str = "boolean";
		break;
	case TagData::CHAR:
		str = "char";
		break;
	case TagData::OCTET:
		str = "octet";
		break;
	case TagData::STRUCT:
		str = "struct";
		break;
	case TagData::STRING:
		str = "string";
		break;
	case TagData::ARRAY:
		str = "array";
		break;
	default:
		str = "unknown";
	}
	return str;
}

const char *GatanDM4::to_str(TagEntry::EntryType type)
{
	const char *str = "unknown";

	switch (type) {
	case TagEntry::GROUP_TAG:
		str = "Group";
		break;
	case TagEntry::DATA_TAG:
		str = "Data";
		break;
	default:
		str = "unknown";
	}
	return str;
}


const char *GatanDM4::to_str(GatanDM4::DataType::GatanDataType type)
{
	switch (type) {
	case GatanDM4::DataType::SIGNED_INT16_DATA:
		return "SIGNED_INT16_DATA";
	case GatanDM4::DataType::REAL4_DATA:
		return "REAL4_DATA";
	case GatanDM4::DataType::COMPLEX8_DATA:
		return "COMPLEX8_DATA";
	case GatanDM4::DataType::OBSELETE_DATA:
		return "OBSELETE_DATA";
	case GatanDM4::DataType::PACKED_DATA:
		return "PACKED_DATA";
	case GatanDM4::DataType::UNSIGNED_INT8_DATA:
		return "UNSIGNED_INT8_DATA";
	case GatanDM4::DataType::SIGNED_INT32_DATA:
		return "SIGNED_INT32_DATA";
	case GatanDM4::DataType::RGB_DATA:
		return "RGB_DATA";
	case GatanDM4::DataType::SIGNED_INT8_DATA:
		return "SIGNED_INT8_DATA";
	case GatanDM4::DataType::UNSIGNED_INT16_DATA:
		return "UNSIGNED_INT16_DATA";
	case GatanDM4::DataType::UNSIGNED_INT32_DATA:
		return "UNSIGNED_INT32_DATA";
	case GatanDM4::DataType::REAL8_DATA:
		return "REAL8_DATA";
	case GatanDM4::DataType::COMPLEX16_DATA:
		return "COMPLEX16_DATA";
	case GatanDM4::DataType::BINARY_DATA:
		return "BINARY_DATA";
	case GatanDM4::DataType::RGB_UINT8_0_DATA:
		return "RGB_UINT8_0_DATA";
	case GatanDM4::DataType::RGB_UINT8_1_DATA:
		return "RGB_UINT8_1_DATA";
	case GatanDM4::DataType::RGB_UINT16_DATA:
		return "RGB_UINT16_DATA";
	case GatanDM4::DataType::RGB_FLOAT32_DATA:
		return "RGB_FLOAT32_DATA";
	case GatanDM4::DataType::RGB_FLOAT64_DATA:
		return "RGB_FLOAT64_DATA";
	case GatanDM4::DataType::RGBA_UINT8_0_DATA:
		return "RGBA_UINT8_0_DATA";
	case GatanDM4::DataType::RGBA_UINT8_1_DATA:
		return "RGBA_UINT8_1_DATA";
	case GatanDM4::DataType::RGBA_UINT8_2_DATA:
		return "RGBA_UINT8_2_DATA";
	case GatanDM4::DataType::RGBA_UINT8_3_DATA:
		return "RGBA_UINT8_3_DATA";
	case GatanDM4::DataType::RGBA_UINT16_DATA:
		return "RGBA_UINT16_DATA";
	case GatanDM4::DataType::RGBA_FLOAT32_DATA:
		return "RGBA_FLOAT32_DATA";
	case GatanDM4::DataType::RGBA_FLOAT64_DATA:
		return "RGBA_FLOAT64_DATA";
	case GatanDM4::DataType::POINT2_SINT16_0_DATA:
		return "POINT2_SINT16_0_DATA";
	case GatanDM4::DataType::POINT2_SINT16_1_DATA:
		return "POINT2_SINT16_1_DATA";
	case GatanDM4::DataType::POINT2_SINT32_0_DATA:
		return "POINT2_SINT32_0_DATA";
	case GatanDM4::DataType::POINT2_FLOAT32_0_DATA:
		return "POINT2_FLOAT32_0_DATA";
	case GatanDM4::DataType::RECT_SINT16_1_DATA:
		return "RECT_SINT16_1_DATA";
	case GatanDM4::DataType::RECT_SINT32_1_DATA:
		return "RECT_SINT32_1_DATA";
	case GatanDM4::DataType::RECT_FLOAT32_1_DATA:
		return "RECT_FLOAT32_1_DATA";
	case GatanDM4::DataType::RECT_FLOAT32_0_DATA:
		return "RECT_FLOAT32_0_DATA";
	case GatanDM4::DataType::SIGNED_INT64_DATA:
		return "SIGNED_INT64_DATA";
	case GatanDM4::DataType::UNSIGNED_INT64_DATA:
		return "UNSIGNED_INT64_DATA";
	default:
		break;
	}
	return "Unknown Type";
}
