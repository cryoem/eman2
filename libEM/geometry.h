/**
 * $Id$
 */
#ifndef eman_geometry_h_
#define eman_geometry_h_

#include <string>
using std::string;

namespace EMAN
{
	/** Size is used to describe a 2D or 3D rectangular size.
     */
	class Size
	{
	public:
		Size():x(0), y(0), z(0)
		{
		}
		Size(int xx, int yy):x(x), y(yy), z(0)
		{
		}
		Size(int xx, int yy, int zz):x(xx), y(yy), z(zz)
		{
		}

		int get_ndim() const
		{
			if (z > 1) {
				return 3;
			}
			return 2;
		}

		int x;
		int y;
		int z;
	};

	/** Template Point defines a point in a 2D/3D space. It stores the
     * point's coordinates. The coordinates can be in any native
     * numerical type: int, float, double, etc. 
     */
	template < class T > class Point {
	public:
		Point():x(0), y(0), z(0), ndim(0) {
		}
		Point(T xx, T yy):x(xx), y(yy), z(0), ndim(2) {
		}
		Point(T xx, T yy, T zz):x(xx), y(yy), z(zz), ndim(3) {
		}
		
		int get_ndim() const
		{
			return ndim;
		}

		T x;
		T y;
		T z;
		
	private:
		int ndim;
	};

	class Pixel {
	public:
		Pixel(int xx, int yy, int zz, float vv) : x(xx), y(yy), z(zz), value(vv) { }	
		int x;
		int y;
		int z;
		float value;
	};

	
	bool operator<(const Pixel& p1, const Pixel& p2);
	

	
	/** Region defines a 2D or 3D rectangular region specified by origins and sizes.
     */
	class Region
	{
	public:
		Region()
		{
			origin = Point < float >(0, 0);
			size = Size(0, 0);
		}

		Region(float x, float y, int xsize, int ysize)
		{
			origin = Point < float >(x, y);
			size = Size(xsize, ysize);
		}

		Region(float x, float y, float z, int xsize, int ysize, int zsize)
		{
			origin = Point < float >(x, y, z);
			size = Size(xsize, ysize, zsize);
		}

		Region(const Point < float >&o, const Size & s):origin(o), size(s)
		{
		}

		~Region() {
		}

		/** to check whether a point is inside this region
		 */
		bool inside_region() const;
		bool inside_region(const Point < float >&p) const;
		bool inside_region(float x, float y) const;
		bool inside_region(float x, float y, float z) const;

		int get_ndim() const
		{
			return origin.get_ndim();
		}
		string get_string() const;

		Point < float >origin;
		Size size;
	};
}

#endif
