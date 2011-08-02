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

#include <cstring>

#include "dm3io.h"
#include "portable_fileio.h"
#include "geometry.h"

using namespace EMAN;
using namespace EMAN::Gatan;

const char *TagTable::IMAGE_WIDTH_TAG = "Dimensions #0";
const char *TagTable::IMAGE_HEIGHT_TAG = "Dimensions #1";
const char *TagTable::IMAGE_DATATYPE_TAG = "DataType";
const char *TagTable::IMAGE_THUMB_INDEX_TAG = "ImageIndex";

TagTable::TagTable()
	:	img_index(0), is_big_endian(true)
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
	else {
		tags[name] = value;
	}
}

void TagTable::add_data(char *data)
{
	if (!data) {
		throw NullPointerException("DM3 data is NULL");
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
	char val_str[32];

	if (tag_type == SHORT) {
		short val = 0;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == USHORT) {
		unsigned short val = 0;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == INT) {
		int val = 0;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == CHAR || tag_type == OCTET) {
		char val = 0;
		fread(&val, sz, 1, in);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == BOOLEAN) {
		bool val = false;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", val);
	}
	else if (tag_type == UINT) {
		unsigned int val = 0;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%d", (int) val);
	}
	else if (tag_type == FLOAT) {
		float val = 0;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%f", val);
	}
	else if (tag_type == DOUBLE) {
		double val = 0;
		fread(&val, sz, 1, in);
		tagtable->become_host_endian(&val);
		sprintf(val_str, "%10e", val);
	}
	else {
		LOGERR("invalid tag type: '%d'", tag_type);
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

	int array_type = 0;
	fread(&array_type, sizeof(array_type), 1, in);

	ByteOrder::become_big_endian(&array_type);

	LOGVAR("array data type = '%s'", Gatan::to_str((Type) array_type));

	vector < int >item_types;

	if (array_type == STRUCT) {
		item_types = read_struct_types();
	}
	else if (array_type == ARRAY) {
		item_types = read_array_types();
		LOGERR("DM3: don't know how to handle this array type");
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

	fread(buf, size * sizeof(unsigned short), 1, in);
	tagtable->become_host_endian < unsigned short >(buf, size);

	for (int i = 0; i < size; i++) {
		str[i] = static_cast < char >(buf[i]);
	}

	str[size] = '\0';
	string str1 = string(str);

	if( str )
	{
		delete[]str;
		str = 0;
	}
	if( buf )
	{
		delete[]buf;
		buf = 0;
	}

	return str1;
}


int TagData::read_array_data(vector < int >item_types, bool nodata)
{
	ENTERFUNC;
	if (item_types.size() == 0) {
		LOGERR("DM3 item types cannot be empty");
		return 1;
	}

	int err = 0;
	int array_size = 0;

	fread(&array_size, sizeof(array_size), 1, in);
	ByteOrder::become_big_endian(&array_size);

	LOGVAR("array size = %d\n", array_size);

	size_t item_size = 0;
	for (size_t i = 0; i < item_types.size(); i++) {
		item_size += typesize(item_types[i]);
	}

	LOGVAR("%s array item size = %d\n", name.c_str(), item_size);

	size_t buf_size = item_size * array_size;

	if (item_types.size() == 1 && item_types[0] == USHORT) {
		string val = read_string(array_size);
		tagtable->add(name, val);
		LOGVAR("value: %s", val.c_str());
	}
	else if (!nodata && name == "Data") {
		char *data = new char[buf_size];
		fread(data, buf_size, 1, in);

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
			LOGERR("cannot handle this type of DM3 image data");
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

	unsigned int namelength = 0;
	unsigned int nfields = 0;

	fread(&namelength, sizeof(namelength), 1, in);
	ByteOrder::become_big_endian(&namelength);

	fread(&nfields, sizeof(nfields), 1, in);
	ByteOrder::become_big_endian(&nfields);

	LOGVAR("namelength = %d\n", namelength);
	LOGVAR("num fields = %d\n", nfields);

	vector < int >field_types;

	for (unsigned int i = 0; i < nfields; i++) {
		fread(&namelength, sizeof(namelength), 1, in);
		ByteOrder::become_big_endian(&namelength);

		int field_type = 0;
		fread(&field_type, sizeof(field_type), 1, in);
		ByteOrder::become_big_endian(&field_type);

		LOGVAR("%dth namelength = %d, type = '%s'",
			   i, namelength, Gatan::to_str((Type) field_type));
		field_types.push_back(field_type);
	}

	return field_types;
}



int TagData::read_any(bool nodata)
{
	int err = 0;

	fread(&tag_type, sizeof(tag_type), 1, in);
	ByteOrder::become_big_endian(&tag_type);
	LOGVAR("tag type = '%s'\n", Gatan::to_str((Type) tag_type));

	if (tag_type == ARRAY) {
		vector < int >item_types = read_array_types();
		err = read_array_data(item_types, nodata);
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
		fread(&str_sz, sizeof(str_sz), 1, in);
		ByteOrder::become_big_endian(&str_sz);

		char *val = new char[str_sz + 1];
		fread(val, str_sz, 1, in);
		val[str_sz] = '\0';
		string val_str = string(val);
		if( val )
		{
			delete[]val;
			val = 0;
		}

		tagtable->add(name, val_str);
	}
	else {
		read_native(true);
	}

	return err;
}

int TagData::read(bool nodata)
{
	LOGVAR("TagData::read()");
	int err = 0;

	const char *DATA_TYPE_MARK = "%%%%";
	const size_t mark_sz = strlen(DATA_TYPE_MARK);
	char *mark = new char[mark_sz + 1];
	fread(mark, mark_sz, 1, in);
	mark[mark_sz] = '\0';

	if (strcmp(mark, DATA_TYPE_MARK) != 0) {
		LOGERR("data type label has been changed from '%s' to '%s'",
			   DATA_TYPE_MARK, mark);
		return 1;
	}

	if( mark )
	{
		delete[]mark;
		mark = 0;
	}

	int encoded_types_size = 0;
	fread(&encoded_types_size, sizeof(int), 1, in);
	ByteOrder::become_big_endian(&encoded_types_size);

	LOGVAR("encoded types size = %d\n", encoded_types_size);

	err = read_any(nodata);

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
	default:
		LOGERR("no such type: '%d'\n", type);
		break;
	}

	return size;
}

/////////////////////////////////////////////

TagGroup::TagGroup(FILE * data_file, TagTable * table, const string & groupname)
	:	in(data_file), tagtable(table), name(groupname), entry_id(0)
{
}

TagGroup::~TagGroup()
{
}

int TagGroup::read(bool nodata)
{
	LOGVAR("TagGroup::read()");

	int ntags = 0;
	portable_fseek(in, sizeof(char) * 2, SEEK_CUR);
	fread(&ntags, sizeof(ntags), 1, in);
	ByteOrder::become_big_endian(&ntags);
	LOGVAR("DM3: ntags = %d\n", ntags);

	int err = 0;

	for (int i = 0; i < ntags; i++) {
		TagEntry tag_entry(in, tagtable, this);
		err = tag_entry.read(nodata);

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



/////////////////////////////////////////////

TagEntry::TagEntry(FILE * data_file, TagTable * table, TagGroup * parent)
	:	in(data_file), tagtable(table), parent_group(parent), name("")
{
}

TagEntry::~TagEntry()
{
}


int TagEntry::read(bool nodata)
{
	LOGVAR("TagEntry::read()");
	int err = 0;
	char tag_type = 0;
	char *tmp_name = 0;

	fread(&tag_type, sizeof(char), 1, in);

	if (tag_type != GROUP_TAG && tag_type != DATA_TAG) {
		LOGERR("TagEntry::read() invalid tag type: %d", tag_type);
		return 1;
	}

	short name_len = 0;
	fread(&name_len, sizeof(short), 1, in);
	ByteOrder::become_big_endian(&name_len);

	if (name_len != 0) {
		tmp_name = new char[name_len + 1];
		fread(tmp_name, name_len, 1, in);
		tmp_name[name_len] = '\0';
	}
	else {
		string parent_name = parent_group->get_name();
		name_len = static_cast < short >(parent_name.size() + 4);
		tmp_name = new char[name_len + 1];
		sprintf(tmp_name, "%s #%d", parent_name.c_str(), parent_group->get_entry_id());
	}

	name = string(tmp_name);
	if( tmp_name )
	{
		delete[]tmp_name;
		tmp_name = 0;
	}

	LOGVAR("\ntag name: '%s', len: %d, type: '%s'",
		   name.c_str(), name_len, Gatan::to_str((EntryType) tag_type));

	if (tag_type == DATA_TAG) {
		TagData tag_data(in, tagtable, name);
		err = tag_data.read(nodata);
	}
	else if (tag_type == GROUP_TAG) {
		TagGroup group(in, tagtable, name);
		err = group.read(nodata);
	}

	return err;
}


////////////////////////////////////////////

DM3IO::DM3IO(const string & dm3_filename, IOMode rw)
	:	filename(dm3_filename), rw_mode(rw), dm3file(0), initialized(false)
{
	is_big_endian = ByteOrder::is_host_big_endian();
	tagtable = new TagTable();
}

DM3IO::~DM3IO()
{
	if (dm3file) {
		fclose(dm3file);
		dm3file = 0;
	}
	if (tagtable) {
		delete tagtable;
		tagtable = 0;
	}
}

void DM3IO::init()
{
	ENTERFUNC;
	if (initialized) {
		return;
	}
	initialized = true;

	if (rw_mode != READ_ONLY) {
		throw ImageReadException(filename, "only support DM3 read-only");
	}

	dm3file = sfopen(filename, READ_ONLY);

	int buf[NUM_ID_INT];
	if (fread(buf, sizeof(buf), 1, dm3file) != 1) {
		throw ImageReadException(filename, "read first block of DM3 file");
	}

	if (!is_valid(&buf)) {
		throw ImageReadException(filename, "invalid DM3 file");
	}

	int byte_order = buf[2];

	if (byte_order == 0) {
		is_big_endian = true;
	}
	else {
		is_big_endian = false;
	}

	tagtable->set_endian(is_big_endian);
	ByteOrder::become_big_endian(buf, 3);

	LOGDEBUG("dm3 ver = %d, image size = %d, is_big_endian = %d",
			 buf[0], buf[1], (int) is_big_endian);

	EXITFUNC;
}



bool DM3IO::is_valid(const void *first_block)
{
	ENTERFUNC;

	if (!first_block) {
		return false;
	}

	const int *data = static_cast < const int *>(first_block);

	int img_ver = data[0];
	int img_size = data[1];
	int byte_order = data[2];

	ByteOrder::become_big_endian(&img_ver);

	if (img_ver != 3) {
		return false;
	}

	ByteOrder::become_big_endian(&img_size);
	ByteOrder::become_big_endian(&byte_order);

	if (byte_order != 0 && byte_order != 1) {
		return false;
	}

	return true;
}

bool DM3IO::is_image_big_endian()
{
	init();
	return is_big_endian;
}

int DM3IO::read_header(Dict & dict, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	int err = 0;

	//single image format, index can only be zero
	if(image_index == -1) {
		image_index = 0;
	}
	image_index = 0;
	check_read_access(image_index);

	portable_fseek(dm3file, NUM_ID_INT * sizeof(int), SEEK_SET);
	TagGroup root_group(dm3file, tagtable, "");
	root_group.read(true);

	int nx = tagtable->get_xsize();
	int ny = tagtable->get_ysize();

	check_region(area, IntSize(nx, ny));
	int xlen = 0, ylen = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);

	dict["nx"] = xlen;
	dict["ny"] = ylen;
	dict["nz"] = 1;

	dict["DM3.exposure_number"] = tagtable->get_int("Exposure Number");
	dict["DM3.exposure_time"] = tagtable->get_double("Exposure (s)");
	dict["DM3.zoom"] = tagtable->get_double("Zoom");
	dict["DM3.antiblooming"] = tagtable->get_int("Antiblooming");
	dict["DM3.magnification"] = tagtable->get_double("Indicated Magnification");

	dict["DM3.voltage"] = tagtable->get_double("Voltage")/1000.0;
	dict["DM3.cs"] = tagtable->get_double("Cs(mm)");

	dict["DM3.frame_type"] = tagtable->get_string("Processing");
	dict["DM3.camera_x"] = tagtable->get_int("Active Size (pixels) #0");
	dict["DM3.camera_y"] = tagtable->get_int("Active Size (pixels) #1");
	dict["DM3.binning_x"] = tagtable->get_int("Binning #0");
	dict["DM3.binning_y"] = tagtable->get_int("Binning #1");
	dict["datatype"] = to_em_datatype(tagtable->get_datatype());

//tagtable->dump();

	EXITFUNC;
	return err;
}

int DM3IO::read_data(float *rdata, int image_index, const Region * area, bool)
{
	ENTERFUNC;
	//single image format, index can only be zero
	image_index = 0;
	check_read_access(image_index, rdata);

	portable_fseek(dm3file, NUM_ID_INT * sizeof(int), SEEK_SET);

	TagGroup root_group(dm3file, tagtable, "");
	root_group.read(false);

	int nx = tagtable->get_xsize();
	int ny = tagtable->get_ysize();

	check_region(area, IntSize(nx, ny));

	int xlen = 0, ylen = 0, x0 = 0, y0 = 0;
	EMUtil::get_region_dims(area, nx, &xlen, ny, &ylen);
	EMUtil::get_region_origins(area, &x0, &y0);

	char *data = tagtable->get_data();
	int data_type = tagtable->get_datatype();

	int k = 0;
	for (int i = y0; i < y0 + ylen; i++) {
		for (int j = x0; j < x0 + xlen; j++) {
			switch (data_type) {
			case Gatan::DataType::SIGNED_INT8_DATA:
				rdata[k] = (float) ((char *) data)[i * nx + j];
				break;
			case Gatan::DataType::UNSIGNED_INT8_DATA:
				rdata[k] = (float) ((unsigned char *) data)[i * nx + j];
				break;
			case Gatan::DataType::SIGNED_INT16_DATA:
				rdata[k] = (float) ((short *) data)[i * nx + j];
				break;
			case Gatan::DataType::UNSIGNED_INT16_DATA:
				rdata[k] = (float) ((unsigned short *) data)[i * nx + j];
				break;
			case Gatan::DataType::SIGNED_INT32_DATA:
				rdata[k] = (float) ((int *) data)[i * nx + j];
				break;
			case Gatan::DataType::UNSIGNED_INT32_DATA:
				rdata[k] = (float) ((unsigned int *) data)[i * nx + j];
				break;
			default:
				string desc = string("unsupported DM3 data type") +
					string(Gatan::to_str((Gatan::DataType::GatanDataType) data_type));
				throw ImageReadException(filename, desc);
			}
			k++;
		}
	}
	EXITFUNC;
	return 0;
}

bool DM3IO::is_complex_mode()
{
	return false;
}

int DM3IO::write_header(const Dict &, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	LOGWARN("DM3 write is not supported.");
	EXITFUNC;
	return 1;
}

int DM3IO::write_data(float *, int, const Region* , EMUtil::EMDataType, bool)
{
	ENTERFUNC;
	LOGWARN("DM3 write is not supported.");
	EXITFUNC;
	return 1;
}

void DM3IO::flush()
{
}

int Gatan::to_em_datatype(int gatan_datatype)
{
	DataType::GatanDataType type = static_cast < DataType::GatanDataType > (gatan_datatype);
	int t = 0;

	switch (type) {
	case Gatan::DataType::SIGNED_INT8_DATA:
		t = EMUtil::EM_CHAR;
		break;
	case Gatan::DataType::UNSIGNED_INT8_DATA:
		t = EMUtil::EM_UCHAR;
		break;
	case Gatan::DataType::SIGNED_INT16_DATA:
		t = EMUtil::EM_SHORT;
		break;
	case Gatan::DataType::UNSIGNED_INT16_DATA:
		t = EMUtil::EM_USHORT;
		break;
	case Gatan::DataType::SIGNED_INT32_DATA:
		t = EMUtil::EM_INT;
		break;
	case Gatan::DataType::UNSIGNED_INT32_DATA:
		t = EMUtil::EM_UINT;
		break;
	default:
		t = EMUtil::EM_UNKNOWN;
	}

	return t;
}


const char *Gatan::to_str(TagData::Type type)
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

const char *Gatan::to_str(TagEntry::EntryType type)
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


const char *Gatan::to_str(Gatan::DataType::GatanDataType type)
{
	switch (type) {
	case Gatan::DataType::SIGNED_INT16_DATA:
		return "SIGNED_INT16_DATA";
	case Gatan::DataType::REAL4_DATA:
		return "REAL4_DATA";
	case Gatan::DataType::COMPLEX8_DATA:
		return "COMPLEX8_DATA";
	case Gatan::DataType::OBSELETE_DATA:
		return "OBSELETE_DATA";
	case Gatan::DataType::PACKED_DATA:
		return "PACKED_DATA";
	case Gatan::DataType::UNSIGNED_INT8_DATA:
		return "UNSIGNED_INT8_DATA";
	case Gatan::DataType::SIGNED_INT32_DATA:
		return "SIGNED_INT32_DATA";
	case Gatan::DataType::RGB_DATA:
		return "RGB_DATA";
	case Gatan::DataType::SIGNED_INT8_DATA:
		return "SIGNED_INT8_DATA";
	case Gatan::DataType::UNSIGNED_INT16_DATA:
		return "UNSIGNED_INT16_DATA";
	case Gatan::DataType::UNSIGNED_INT32_DATA:
		return "UNSIGNED_INT32_DATA";
	case Gatan::DataType::REAL8_DATA:
		return "REAL8_DATA";
	case Gatan::DataType::COMPLEX16_DATA:
		return "COMPLEX16_DATA";
	case Gatan::DataType::BINARY_DATA:
		return "BINARY_DATA";
	case Gatan::DataType::RGB_UINT8_0_DATA:
		return "RGB_UINT8_0_DATA";
	case Gatan::DataType::RGB_UINT8_1_DATA:
		return "RGB_UINT8_1_DATA";
	case Gatan::DataType::RGB_UINT16_DATA:
		return "RGB_UINT16_DATA";
	case Gatan::DataType::RGB_FLOAT32_DATA:
		return "RGB_FLOAT32_DATA";
	case Gatan::DataType::RGB_FLOAT64_DATA:
		return "RGB_FLOAT64_DATA";
	case Gatan::DataType::RGBA_UINT8_0_DATA:
		return "RGBA_UINT8_0_DATA";
	case Gatan::DataType::RGBA_UINT8_1_DATA:
		return "RGBA_UINT8_1_DATA";
	case Gatan::DataType::RGBA_UINT8_2_DATA:
		return "RGBA_UINT8_2_DATA";
	case Gatan::DataType::RGBA_UINT8_3_DATA:
		return "RGBA_UINT8_3_DATA";
	case Gatan::DataType::RGBA_UINT16_DATA:
		return "RGBA_UINT16_DATA";
	case Gatan::DataType::RGBA_FLOAT32_DATA:
		return "RGBA_FLOAT32_DATA";
	case Gatan::DataType::RGBA_FLOAT64_DATA:
		return "RGBA_FLOAT64_DATA";
	case Gatan::DataType::POINT2_SINT16_0_DATA:
		return "POINT2_SINT16_0_DATA";
	case Gatan::DataType::POINT2_SINT16_1_DATA:
		return "POINT2_SINT16_1_DATA";
	case Gatan::DataType::POINT2_SINT32_0_DATA:
		return "POINT2_SINT32_0_DATA";
	case Gatan::DataType::POINT2_FLOAT32_0_DATA:
		return "POINT2_FLOAT32_0_DATA";
	case Gatan::DataType::RECT_SINT16_1_DATA:
		return "RECT_SINT16_1_DATA";
	case Gatan::DataType::RECT_SINT32_1_DATA:
		return "RECT_SINT32_1_DATA";
	case Gatan::DataType::RECT_FLOAT32_1_DATA:
		return "RECT_FLOAT32_1_DATA";
	case Gatan::DataType::RECT_FLOAT32_0_DATA:
		return "RECT_FLOAT32_0_DATA";
	case Gatan::DataType::SIGNED_INT64_DATA:
		return "SIGNED_INT64_DATA";
	case Gatan::DataType::UNSIGNED_INT64_DATA:
		return "UNSIGNED_INT64_DATA";
	default:
		break;
	}
	return "Unknown Type";
}

