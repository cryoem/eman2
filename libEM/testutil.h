/**
 * $Id$
 */

#ifndef eman__testutil_h__
#define eman__testutil_h__

#include "geometry.h"
#include "vec3.h"
#include <map>
#include <iostream>
using namespace std;

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
		static const char* get_debug_string(int i);
		
		static const char* get_debug_image(const char* imagename);

		static void to_emobject(const Dict & d);
		
		static IntPoint test_IntPoint(const IntPoint & p);
		static FloatPoint test_FloatPoint(const FloatPoint & p);
		static IntSize test_IntSize(const IntSize & p);
		static FloatSize test_FloatSize(const FloatSize & p);
		static Vec3i test_Vec3i(const Vec3i & p);
		static Vec3f test_Vec3f(const Vec3f & p);

		template<class T>
		static map<string, T> test_map_dict(const map<string, T>& d)
		{
			for (int i = 0; i < 3; i++) {
				cout << "map[\"" << ts[i] << "\"]=" << d[ts[i]] << ";";
			}
			cout << endl;
		}
		

	private:
		static float tf[10];
		static int ti[10];
		static const char *ts[10];
	};
}

#endif
