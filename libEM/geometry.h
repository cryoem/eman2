/**
 * $Id$
 */
#ifndef eman_geometry_h_
#define eman_geometry_h_

#include <string>
using std::string;

namespace EMAN
{

	/** IntSize is used to describe a 1D, 2D or 3D rectangular size in integers.
     */
	class IntSize
	{
	public:
		
		/** Construct an IntSize object.
		 * @param xx The x direction size. Default is 0.
		 * @param yy The y direction size. Default is 0.
		 * @param zz The z direction size. Default is 0.
		 */
		IntSize(int xx=0, int yy=0, int zz=0)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
		}

		/** Get its dimension, 1D, 2D, or 3D.
		 * @return The dimension.
		 */
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

		/** Get the ith direction's size. Used as a rvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's size. 
		 */
		int operator[] (int i) const
		{
			return data[i];
		}

		/** Get the ith direction's size. Used as a lvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's size. 
		 */
		int & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		int data[3];
	};

	/** FloatSize is used to describe a 1D, 2D or 3D rectangular size
		in floating numbers.
     */
	class FloatSize
	{
	public:
		
		/** Construct a FloatSize object.
		 * @param xx The x direction size. Default is 0.
		 * @param yy The y direction size. Default is 0.
		 * @param zz The z direction size. Default is 0.
		 */
		FloatSize(float xx=0, float yy=0, float zz=0)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
		}

		/** Construct a FloatSize object.
		 * @param xx The x direction size. Default is 0.
		 * @param yy The y direction size. Default is 0.
		 * @param zz The z direction size. Default is 0.
		 */
		FloatSize(int xx, int yy=0, int zz=0)
		{
			data[0] = (float)xx;
			data[1] = (float)yy;
			data[2] = (float)zz;
		}

		/** Construct a FloatSize object.
		 * @param xx The x direction size. Default is 0.
		 * @param yy The y direction size. Default is 0.
		 * @param zz The z direction size. Default is 0.
		 */
		FloatSize(double xx, double yy=0, double zz=0)
		{
			data[0] = (float)xx;
			data[1] = (float)yy;
			data[2] = (float)zz;
		}

		/** Get its dimension, 1D, 2D, or 3D.
		 * @return The dimension.
		 */
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
		
		/** Get the ith direction's size. Used as a rvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's size. 
		 */
		float operator[] (int i) const
		{
			return data[i];
		}
		
		/** Get the ith direction's size. Used as a lvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's size. 
		 */
		float & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		float data[3];
	};
	
	/** IntPoint defines an integer-coordinate point in a 1D/2D/3D space.
	 */
	class IntPoint {
	public:
		/** Construct a point at the origin location.
		 */
		IntPoint()
		{
			data[0] = 0;
			data[1] = 0;
			data[2] = 0;
			ndim = 0;
		}
		
		/** Construct a 1D point.
		 * @param xx The x coordinate value.
		 */
		IntPoint(int xx)
		{
			data[0] = xx;
			data[1] = 0;
			data[2] = 0;
			ndim = 1;
		}
		
		/** Construct a 2D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 */
		IntPoint(int xx, int yy)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = 0;
			ndim = 2;
		}
		
		/** Construct a 3D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 * @param zz The z coordinate value.
		 */
		IntPoint(int xx, int yy, int zz)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
			ndim = 3;
		}

		/** Get the dimension of the point, 1D/2D/3D.
		 * @return The dimension of the point.
		 */
		int get_ndim() const
		{
			return ndim;
		}
		
		/** Get the ith direction's coordinate. Used as a rvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's coordinate. 
		 */
		int operator[] (int i) const
		{
			return data[i];
		}

		/** Get the ith direction's coordinate. Used as a lvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's coordinate. 
		 */
		int & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		int data[3];
		int ndim;
	};

	/** FloatPoint defines a float-coordinate point in a 1D/2D/3D space.
	*/
	class FloatPoint {
	public:

		/** Construct a point at the origin location.
		 */
		FloatPoint()
		{
			data[0] = 0;
			data[1] = 0;
			data[2] = 0;
			ndim = 0;
		}
		
		/** Construct a 1D point.
		 * @param xx The x coordinate value.
		 */
		FloatPoint(float xx)
		{
			data[0] = xx;
			data[1] = 0;
			data[2] = 0;
			ndim = 1;
		}

		/** Construct a 2D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 */
		FloatPoint(float xx, float yy)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = 0;
			ndim = 2;
		}
		
		/** Construct a 3D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 * @param zz The z coordinate value.
		 */
		FloatPoint(float xx, float yy, float zz)
		{
			data[0] = xx;
			data[1] = yy;
			data[2] = zz;
			ndim = 3;
		}
		
		/** Construct a 1D point.
		 * @param xx The x coordinate value.
		 */
		FloatPoint(int xx)
		{
			data[0] = (float)xx;
			data[1] = 0;
			data[2] = 0;
			ndim = 1;
		}
		
		/** Construct a 2D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 */
		FloatPoint(int xx, int yy)
		{
			data[0] = (float)xx;
			data[1] = (float)yy;
			data[2] = 0;
			ndim = 2;
		}
		
		/** Construct a 3D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 * @param zz The z coordinate value.
		 */
		FloatPoint(int xx, int yy, int zz)
		{
			data[0] = (float)xx;
			data[1] = (float)yy;
			data[2] = (float)zz;
			ndim = 3;
		}
			
		/** Construct a 1D point.
		 * @param xx The x coordinate value.
		 */
		FloatPoint(double xx)
		{
			data[0] = (float)xx;
			data[1] = 0;
			data[2] = 0;
			ndim = 1;
		}
		
		/** Construct a 2D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 */
		FloatPoint(double xx, double yy)
		{
			data[0] = (float)xx;
			data[1] = (float)yy;
			data[2] = 0;
			ndim = 2;
		}
		
		/** Construct a 3D point.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 * @param zz The z coordinate value.
		 */
		FloatPoint(double xx, double yy, double zz)
		{
			data[0] = (float)xx;
			data[1] = (float)yy;
			data[2] = (float)zz;
			ndim = 3;
		}

		/** Get the dimension of the point, 1D/2D/3D.
		 * @return The dimension of the point.
		 */
		int get_ndim() const
		{
			return ndim;
		}

		/** Get the ith direction's coordinate. Used as a rvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's coordinate. 
		 */
		float operator[] (int i) const
		{
			return data[i];
		}

		/** Get the ith direction's coordinate. Used as a lvalue.
		 * @param i The ith direction, with 0 is x, 1 is y, 2 is z.
		 * @return The ith direction's coordinate. 
		 */
		float & operator[] (int i) 
		{
			return data[i];
		}
		
	private:
		float data[3];
		int ndim;
	};

	/** Pixel describes a 3D pixel's coordinates and its intensity value. 
	 */
	class Pixel {
	public:
		/** Construct a Pixel object given its 3D coordinates and its value.
		 * @param xx The x coordinate value.
		 * @param yy The y coordinate value.
		 * @param zz The z coordinate value.
		 * @param vv The pixel's intensity value.
		 */
		Pixel(int xx, int yy, int zz, float vv) : x(xx), y(yy), z(zz), value(vv) { }

		/** Get the pixel's coordinates as an integer point.
		 * @return An integer point containing the pixel's coordinates.
		*/
		IntPoint get_point() const
		{
			return IntPoint(x, y, z);
		}

		/** Get the pixel's intensity value.
		 * @return The pixel's intensity value.
		 */
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
	bool operator==(const Pixel& p1, const Pixel& p2);
	bool operator!=(const Pixel& p1, const Pixel& p2);

	
	/** Region defines a 2D or 3D rectangular region specified by its
	 * origin coordinates and all edges' sizes. The coordinates and
	 * edge sizes can be integer or floating numbers.
     */
	class Region
	{
	public:
		/** Construct a null region with its origin at coordinate origins
		 * and its sizes to be 0.
		 */
		Region()
		{
			origin = FloatPoint ();
			size = FloatSize();
		}

		/** Construct a 2D integer region.
		 */
		Region(int x, int y, int xsize, int ysize)
		{
			origin = FloatPoint (x, y);
			size = FloatSize(xsize, ysize);
		}

		/** Construct a 3D integer region.
		 */
		Region(int x, int y, int z, int xsize, int ysize, int zsize)
		{
			origin = FloatPoint(x, y, z);
			size = FloatSize(xsize, ysize, zsize);
		}
		
		/** Construct a 2D floating-number region.
		 */
		Region(float x, float y, float xsize, float ysize)
		{
			origin = FloatPoint (x, y);
			size = FloatSize(xsize, ysize);
		}
		
		/** Construct a 3D floating-number region.
		 */
		Region(float x, float y, float z, float xsize, float ysize, float zsize)
		{
			origin = FloatPoint(x, y, z);
			size = FloatSize(xsize, ysize, zsize);
		}

		/** Construct a 2D floating-number region.
		 */
		Region(double x, double y, double xsize, double ysize)
		{
			origin = FloatPoint (x, y);
			size = FloatSize(xsize, ysize);
		}

		/** Construct a 3D floating-number region.
		 */
		Region(double x, double y, double z, double xsize, double ysize, double zsize)
		{
			origin = FloatPoint(x, y, z);
			size = FloatSize(xsize, ysize, zsize);
		}

		/** Construct a region given's orginal point and edge sizes.
		 */
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

		/** Get the region's dimension.
		 * @return The region's dimension.
		 */
		int get_ndim() const
		{
			return origin.get_ndim();
		}

		/** Get the description of this region in a string.
		 * @return the description of this region in a string.
		 */
		string get_string() const;

		FloatPoint origin;  /* region original point. */
		FloatSize size;     /* region edge sizes. */
	};
}

#endif
