/**
 * $Id$
 */

#ifndef eman__testutil_h__
#define eman__testutil_h__

#include "geometry.h"
#include "vec3.h"

namespace EMAN
{
	/**TestUtil defines function assisting testing of EMAN2.
	 */

	class Dict;
	
	class TestUtil
	{
	public:
		static const char* get_debug_image(const char* imagename);
		static void to_emobject(const Dict & d);		

		static void to_IntPoint(const IntPoint & p);
		static IntPoint from_IntPoint();

		static void to_FloatPoint(const FloatPoint & p);
		static  FloatPoint from_FloatPoint();
		
		static void to_IntSize(const IntSize & p);
		static  IntSize from_IntSize();

		static void to_FloatSize(const FloatSize & p);
		static  FloatSize from_FloatSize();

		static void to_Vec3f(const Vec3f & p);
		static  Vec3f from_Vec3f();

		static void to_Vec3i(const Vec3i & p);
		static  Vec3i from_Vec3i();

		
	};
}

#endif
