/**
 * $Id$
 */
#ifndef eman_geometry_h_
#define eman_geometry_h_

#include <string>
using std::string;

namespace EMAN
{

	/** IntSize is used to describe a 2D or 3D rectangular size.
     */
	class IntSize
	{
	public:
		IntSize(int xx=0, int yy=0, int zz=0)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
		}

		int get_ndim() const
		{
			if (data[2] > 1) {
				return 3;
			}
			else if (data[1] > 1) {
				return 2;
			}
			return 1;
		}

		int operator[] (int i) const
		{
			return data[i];
		}

		int & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		int data[3];
	};

	/** FloatSize is used to describe a 2D or 3D rectangular size.
     */
	class FloatSize
	{
	public:
		FloatSize(float xx=0, float yy=0, float zz=0)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
		}

		int get_ndim() const
		{
			if (data[2] > 1) {
				return 3;
			}
			else if (data[1] > 1) {
				return 2;
			}
			return 1;
		}
		
		float operator[] (int i) const
		{
			return data[i];
		}

		float & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		float data[3];
	};
	
	/** IntPoint defines a integer-coordinate point in a 2D/3D space. */
	class IntPoint {
	public:
		IntPoint()
		{
			data[0] = 0;
			data[1] = 0;
			data[2] = 0;
			ndim = 0;
		}
		IntPoint(int xx)
		{
			data[0] = xx;
			data[1] = 0;
			data[2] = 0;
			ndim = 1;
		}
		IntPoint(int xx, int yy)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = 0;
			ndim = 2;
		}
		IntPoint(int xx, int yy, int zz)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
			ndim = 3;
		}
		
		int get_ndim() const
		{
			return ndim;
		}
		
		int operator[] (int i) const
		{
			return data[i];
		}

		int & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		int data[3];
		int ndim;
	};

	/** FloatPoint defines a float-coordinate point in a 2D/3D space. */
	class FloatPoint {
	public:
		FloatPoint()
		{
			data[0] = 0;
			data[1] = 0;
			data[2] = 0;
			ndim = 0;
		}
		FloatPoint(float xx)
		{
			data[0] = xx;
			data[1] = 0;
			data[2] = 0;
			ndim = 1;
		}
		FloatPoint(float xx, float yy)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = 0;
			ndim = 2;
		}
		FloatPoint(float xx, float yy, float zz)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
			ndim = 3;
		}
		
		int get_ndim() const
		{
			return ndim;
		}

		float operator[] (int i) const
		{
			return data[i];
		}

		float & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		float data[3];
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
			size = FloatSize(0, 0);
		}

		Region(float x, float y, float xsize, float ysize)
		{
			origin = FloatPoint (x, y);
			size = FloatSize(xsize, ysize);
		}

		Region(float x, float y, float z, float xsize, float ysize, float zsize)
		{
			origin = FloatPoint(x, y, z);
			size = FloatSize(xsize, ysize, zsize);
		}

		Region(const FloatPoint &o, const FloatSize & s):origin(o), size(s)
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
		FloatSize size;
	};
}

#endif
