#ifndef __obejct__em__
#define __obejct__em__

#include <string>
#include <map>

using std::string;
using std::map;

namespace EMAN {

    enum MapInfoType {
	NORMAL,
	ICOS2F_FIRST_OCTANT,
	ICOS2F_FULL,
	ICOS2F_HALF,
	ICOS3F_HALF,
	ICOS3F_FULL,
	ICOS5F_HALF,
	ICOS5F_FULL,
	ICOS_UNKNOWN
    };

    class EMData;
    
    class EMObject {
    public:
	EMObject();
	EMObject(int n);
	EMObject(float f);
	EMObject(double d);
	EMObject(string str);
	EMObject(EMData* em);
	
	~EMObject();

	int get_int() const;
	float get_float() const;
	double get_double() const;
	string get_string() const;
	EMData* get_EMData() const;

	bool is_null() const;
	string to_str() const;
	
    private:
	enum ObjectType {
	    INT_OBJECT,
	    FLOAT_OBJECT,
	    DOUBLE_OBJECT,
	    STRING_OBJECT,
	    EMDATA_OBJECT,
	    UNKNOWN_OBJECT
	};

	const char* get_object_type_name(ObjectType t) const;
	
    private:
	union {
	    int n;
	    float f;
	    double d;
	};
	
	EMData* emdata;
	string str;
	ObjectType type;
    };

    
   
    
}

#endif
