/**
 * $Id$
 */

#ifndef eman__testutil_h__
#define eman__testutil_h__

#include "geometry.h"
#include "vec3.h"
#include "emobject.h"
#include "emdata.h"

#include <map>
#include <vector>
#include <string>

using std::map;
using std::vector;
using std::string;

namespace EMAN
{
	/**TestUtil defines function assisting testing of EMAN2.
	 */

	class Dict;
	
	class TestUtil
	{
	public:
		static int get_debug_int(int i);
		static float get_debug_float(int i);
		static string get_debug_string(int i);
		
		static string get_debug_image(const string & imagename);

		static void to_emobject(const Dict & d);
		
		static IntPoint test_IntPoint(const IntPoint & p);
		static FloatPoint test_FloatPoint(const FloatPoint & p);
		static IntSize test_IntSize(const IntSize & p);
		static FloatSize test_FloatSize(const FloatSize & p);
		static Vec3i test_Vec3i(const Vec3i & p);
		static Vec3f test_Vec3f(const Vec3f & p);

		static vector<int> test_vector_int(const vector<int> & v);
		static vector<float> test_vector_float(const vector<float> & v);
		static vector<long> test_vector_long(const vector<long> & v);
		static vector<string> test_vector_string(const vector<string> & v);
		static vector<EMData*> test_vector_emdata(const vector<EMData*> & v);
		static vector<Pixel> test_vector_pixel(const vector<Pixel> & v);

		static map<string, int> test_map_int(const map<string, int>& d);
		static map<string, long> test_map_long(const map<string, long>& d);
		static map<string, float> test_map_float(const map<string, float>& d);
		static map<string, string> test_map_string(const map<string, string>& d);
		static map<string, EMObject> test_map_emobject(const map<string, EMObject>& d);
		static map<string, vector<string> > test_map_vecstring(const map<string,
															   vector<string> >& d);

		static Dict test_dict(const Dict & d);
		
	private:
		static float tf[10];
		static int ti[10];
	};
}

#endif


