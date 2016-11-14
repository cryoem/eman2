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
 
#ifndef __dm3io_h__
#define __dm3io_h__

#include "imageio.h"

using std::vector;
using std::map;

namespace EMAN
{
	namespace Gatan
	{
		class TagTable
		{
		  public:
			TagTable();
			~TagTable();

			void add(const string & name, const string & value);
			void add_data(char *data);

			string get_string(const string & name);
			int get_int(const string & name);
			float get_float(const string & name);
			double get_double(const string & name);

			int get_xsize() const;
			int get_ysize() const;
			int get_datatype() const;
			char *get_data() const;

			void dump() const;

			template < class T > void become_host_endian(T * data, int n = 1) {
				if (is_big_endian != ByteOrder::is_host_big_endian()) {
					ByteOrder::swap_bytes(data, n);
				}
			}

			void set_endian(bool big_endian)
			{
				is_big_endian = big_endian;
			}

		  private:
			static const char *IMAGE_WIDTH_TAG;
			static const char *IMAGE_HEIGHT_TAG;
			static const char *IMAGE_DATATYPE_TAG;
			static const char *IMAGE_THUMB_INDEX_TAG;
			void set_thumb_index(int i);

		  private:
			int img_index;
			bool is_big_endian;
			std::map < string, string > tags;
			vector < int >x_list;
			vector < int >y_list;
			vector < int >datatype_list;
			vector < char *>data_list;
		};

		class TagData
		{
		  public:
			enum Type
			{
				UNKNOWN = 0,
				SHORT = 2,
				INT = 3,
				USHORT = 4,
				UINT = 5,
				FLOAT = 6,
				DOUBLE = 7,
				BOOLEAN = 8,
				CHAR = 9,
				OCTET = 10,
				STRUCT = 15,
				STRING = 18,
				ARRAY = 20
			};

			TagData(FILE * data_file, TagTable * tagtable, const string & tagname);
			~TagData();

			int read(bool nodata = false);

		  private:
			size_t typesize() const;
			size_t typesize(int type) const;
			int read_any(bool nodata = false);

			vector < int >read_array_types();
			int read_array_data(vector < int >item_types, bool nodata = false);
			vector < int >read_struct_types();
			string read_native(bool is_value_stored);
			string read_string(int size);

		  private:
			FILE * in;
			TagTable *tagtable;
			string name;
			Type tag_type;
		};


		class TagGroup
		{
		  public:
			TagGroup(FILE * data_file, TagTable * tagtable, const string & groupname);
			~TagGroup();

			int read(bool nodata = false);
			string get_name() const;
			int get_entry_id();

		  private:
			FILE * in;
			TagTable *tagtable;
			string name;
			int entry_id;
		};


		class TagEntry
		{
		  public:
			enum EntryType
			{
				GROUP_TAG = 20,
				DATA_TAG = 21
			};

			TagEntry(FILE * data_file, TagTable * tagtable, TagGroup * parent_group);
			~TagEntry();

			int read(bool nodata = false);

		  private:
			FILE * in;
			TagTable *tagtable;
			TagGroup *parent_group;
			string name;
		};


		class DataType
		{
		  public:
			enum GatanDataType
			{
				NULL_DATA,
				SIGNED_INT16_DATA,
				REAL4_DATA,
				COMPLEX8_DATA,
				OBSELETE_DATA,
				PACKED_DATA,
				UNSIGNED_INT8_DATA,
				SIGNED_INT32_DATA,
				RGB_DATA,
				SIGNED_INT8_DATA,
				UNSIGNED_INT16_DATA,
				UNSIGNED_INT32_DATA,
				REAL8_DATA,
				COMPLEX16_DATA,
				BINARY_DATA,
				RGB_UINT8_0_DATA,
				RGB_UINT8_1_DATA,
				RGB_UINT16_DATA,
				RGB_FLOAT32_DATA,
				RGB_FLOAT64_DATA,
				RGBA_UINT8_0_DATA,
				RGBA_UINT8_1_DATA,
				RGBA_UINT8_2_DATA,
				RGBA_UINT8_3_DATA,
				RGBA_UINT16_DATA,
				RGBA_FLOAT32_DATA,
				RGBA_FLOAT64_DATA,
				POINT2_SINT16_0_DATA,
				POINT2_SINT16_1_DATA,
				POINT2_SINT32_0_DATA,
				POINT2_FLOAT32_0_DATA,
				RECT_SINT16_1_DATA,
				RECT_SINT32_1_DATA,
				RECT_FLOAT32_1_DATA,
				RECT_FLOAT32_0_DATA,
				SIGNED_INT64_DATA,
				UNSIGNED_INT64_DATA,
				LAST_DATA
			};
		};

		int to_em_datatype(int gatan_datatype);
		const char *to_str(Gatan::TagData::Type type);
		const char *to_str(Gatan::TagEntry::EntryType type);
		const char *to_str(Gatan::DataType::GatanDataType type);
	}

	/** Gatan DM3 file is a hierarchical binary image
	 * format. Everything in the image is a <key, value> pair, where
	 * key may be a container-type key which contains more key/value
	 * pairs. To read its header information, the whole file has to be
	 * parsed. During parsing, we check the keys that we are
	 * interested in and get their values.
	 *
	 * The real binary data itself is also in this key/value
	 * hierarchy.
	 *
	 * 1 Gatan DM3 file contains 1 2D image.
	*/

	class DM3IO:public ImageIO
	{
	  public:
		explicit DM3IO(const string & filename, IOMode rw_mode = READ_ONLY);
		~DM3IO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		enum
		{ NUM_ID_INT = 3 };

	  private:
		string filename;
		IOMode rw_mode;
		FILE *dm3file;
		bool is_big_endian;
		bool initialized;
		Gatan::TagTable * tagtable;
	};

}




#endif	//__dm3io_h__
