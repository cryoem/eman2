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


	/** IntPoint defines a integer-coordinate point in a 2D/3D space. */
	class IntPoint {
	public:
		IntPoint():x(0), y(0), z(0), ndim(0) {
		}
		IntPoint(int xx, int yy):x(xx), y(yy), z(0), ndim(2) {
		}
		IntPoint(int xx, int yy, int zz):x(xx), y(yy), z(zz), ndim(3) {
		}
		
		int get_ndim() const
		{
			return ndim;
		}

		int x;
		int y;
		int z;
		
	private:
		int ndim;
	};

	/** FloatPoint defines a float-coordinate point in a 2D/3D space. */
	class FloatPoint {
	public:
		FloatPoint():x(0), y(0), z(0), ndim(0) {
		}
		FloatPoint(float xx, float yy):x(xx), y(yy), z(0), ndim(2) {
		}
		FloatPoint(float xx, float yy, float zz):x(xx), y(yy), z(zz), ndim(3) {
		}
		
		int get_ndim() const
		{
			return ndim;
		}

		float x;
		float y;
		float z;
	private:
		int ndim;
	};

	
	class Pixel {
	public:
		Pixel(int xx, int yy, int zz, float vv) : x(xx), y(yy), z(zz), value(vv) { }

		IntPoint get_point() const
		{
			return IntPoint(x, y, z);
		}

		float get_value() const
		{
			return value;
		}
		
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
			origin = FloatPoint (0, 0);
			size = Size(0, 0);
		}

		Region(float x, float y, int xsize, int ysize)
		{
			origin = FloatPoint (x, y);
			size = Size(xsize, ysize);
		}

		Region(float x, float y, float z, int xsize, int ysize, int zsize)
		{
			origin = FloatPoint (x, y, z);
			size = Size(xsize, ysize, zsize);
		}

		Region(const FloatPoint &o, const Size & s):origin(o), size(s)
		{
		}

		~Region() {
		}

		/** to check whether a point is inside this region
		 */
		bool inside_region() const;
		bool inside_region(const FloatPoint &p) const;
		bool inside_region(float x, float y) const;
		bool inside_region(float x, float y, float z) const;

		int get_ndim() const
		{
			return origin.get_ndim();
		}
		string get_string() const;

		FloatPoint origin;
		Size size;
	};
}

#endif
