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
	Size() : xsize(0), ysize(0), zsize(0)
	{
	}
	Size(int x, int y) : xsize(x), ysize(y), zsize(0)
	{
	}
	Size(int x, int y, int z) : xsize(x), ysize(y), zsize(z)
	{
	}

	int get_ndim() const
	{
	    if (zsize > 1) {
		return 3;
	    }
	    return 2;
	}

	int xsize;
	int ysize;
	int zsize;
    };

    /** Template Point defines a point in a 2D/3D space. It stores the
     * point's coordinates. The coordinates can be in any native
     * numerical type: int, float, double, etc. 
     */
    template <class T> class Point {
    public:
	Point() : x(0), y(0), z(0), ndim(0) {
	}
	Point(T xx, T yy) : x(xx), y(yy), z(0), ndim(2) {
	}
	Point(T xx, T yy, T zz) : x(xx), y(yy), z(zz), ndim(3) {
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


    /** Region defines a 2D or 3D rectangular region specified by origins and sizes.
     */    
    class Region
    {
    public:
	Region()
	{
	    origin = Point<float>(0, 0);
	    size = Size(0, 0);
	}

	Region(float x, float y, int xsize, int ysize)
	{
	    origin = Point<float>(x, y);
	    size = Size(xsize, ysize);
	}

	Region(float x, float y, float z, int xsize, int ysize, int zsize)
	{
	    origin = Point<float>(x, y, z);
	    size = Size(xsize, ysize, zsize);
	}

	Region(const Point<float>&o, const Size & s) : origin(o), size(s)
	{
	}

	~Region() {
	}

	/** to check whether a point is inside this region
	 */
	bool inside_region() const;
	bool inside_region(const Point<float> & p) const;
	bool inside_region(float x, float y) const;
	bool inside_region(float x, float y, float z) const;

	int get_ndim() const
	{
	    return origin.get_ndim();
	}
	string get_string() const;

	Point<float> origin;
	Size size;
    };
}

#endif
