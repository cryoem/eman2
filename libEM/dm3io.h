#ifndef __dm3io_h__
#define __dm3io_h__

#include "imageio.h"
#include <map>
#include <vector>
#include <stdio.h>

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

	    void add(string name, string value);
	    void add_data(char *data);

	    string get_string(string name);
	    int get_int(string name);
	    float get_float(string name);
	    double get_double(string name);

	    int get_xsize() const;
	    int get_ysize() const;
	    int get_datatype() const;
	    char *get_data() const;

	    void dump() const;

	    template <class T> void become_platform_endian(T * data, int n = 1) {
		if (is_big_endian != ByteOrder::is_machine_big_endian()) {
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
	    std::map<string, string> tags;
	    vector<int> x_list;
	    vector<int> y_list;
	    vector<int> datatype_list;
	    vector<char *> data_list;
	};

	class TagData
	{
	public:
	    enum Type {
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

	    TagData(FILE * data_file, TagTable * tagtable, string tagname);
	    ~TagData();

	    int read(bool nodata = false);

	private:
	    size_t typesize() const;
	    size_t typesize(int type) const;
	    int read_any(bool nodata = false);

	    vector<int> read_array_types();
	    int read_array_data(vector<int> item_types, bool nodata = false);
	    vector<int> read_struct_types();
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
	    TagGroup(FILE * data_file, TagTable * tagtable, string groupname);
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
	    enum EntryType {
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
	    enum GatanDataType {
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


    class DM3IO : public ImageIO
    {
    public:
	DM3IO(string filename, IOMode rw_mode = READ_ONLY);
	~DM3IO();

	DEFINE_IMAGEIO_FUNC;
	static bool is_valid(const void *first_block);
	
    private:
	enum { NUM_ID_INT = 3 };

    private:
	string filename;
	IOMode rw_mode;
	FILE *dm3file;
	bool is_big_endian;
	bool initialized;
	Gatan::TagTable * tagtable;
    };

}




#endif
