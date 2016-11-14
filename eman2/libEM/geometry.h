/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

#ifndef eman_geometry_h_
#define eman_geometry_h_

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <algorithm> // for copy
		 
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
		explicit IntSize(int xx=0, int yy=0, int zz=0)
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
		explicit FloatSize(float xx=0, float yy=0, float zz=0)
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
			else if (data[0] > 1) {
				return 1;
			}
			else {
				return 0;
			}
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
		
		inline operator vector<float>() const {
			vector<float> t(data,data+3);
			return t;
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
		explicit IntPoint(int xx)
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
	
	IntPoint operator -( const IntPoint& p);
// 	{
// 		return IntPoint(-p[0],-p[1],-p[2]);
// 	}

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
		explicit FloatPoint(float xx)
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
		explicit FloatPoint(int xx)
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
		explicit FloatPoint(double xx)
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

		FloatPoint(const FloatPoint & fp)
		{
			data[0] = fp.data[0];
			data[1] = fp.data[1];
			data[2] = fp.data[2];
			ndim = fp.ndim;
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
		
		inline operator vector<float>() const {
			vector<float> t(data,data+3);
			return t;
		}
		
		operator IntPoint () const { return IntPoint((int)data[0],(int)data[1],(int)data[2]); }
		
		inline FloatPoint& operator=(const vector<float>& v) {
			copy(v.begin(),v.end(),data);
			return *this;
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

		Pixel(const Pixel & p) : x(p.x), y(p.y), z(p.z), value(p.value) {}
		
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

		/** Construct a 1D integer region
		 */
		Region(int x, int xsize)
		{
			origin = FloatPoint (x);
			size = FloatSize(xsize);
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
		
		/** Construct a 1D floating-number region.
		 */
		Region(float x, float xsize)
		{
			origin = FloatPoint (x);
			size = FloatSize(xsize);
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

		/** Construct a 1D floating-number region.
		 */
		Region(double x, double xsize)
		{
			origin = FloatPoint (x);
			size = FloatSize(xsize);
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

		Region(const Region & r)
		{
			origin = r.origin;
			size = r.size;
		}

		
		~Region() {
		}

		/** to check whether a point is inside this region
		 */
		bool inside_region() const;
		bool inside_region(const FloatPoint &p) const;
		bool inside_region(float x) const;
		bool inside_region(float x, float y) const;
		bool inside_region(float x, float y, float z) const;


		/** get the width */
		float get_width() const { return size[0]; }
		/** get the height */
		float get_height() const { return size[1]; }
		/** get the depth */
		float get_depth() const { return size[2]; }
		
		/** set the width */
		void set_width(const float& v) { size[0] = v; }
		/** set the height */
		void set_height(const float& v) { size[1] = v; }
		/** set the depth */
		void set_depth(const float& v) { size[2] = v; }
		
		/** get the x element of the origin */
		float x_origin() const { return origin[0]; }
		/** get the y element of the origin */
		float y_origin() const { return origin[1]; }
		/** get the z element of the origin */
		float z_origin() const { return origin[2]; }
		
		/** get the size of each dimension as a vector */
		vector<float> get_size() const { return size; }
		/** get the origin as a vector */
		vector<float> get_origin() const { return origin; }
		/** set the origin using a vector */
		void set_origin(const vector<float>& v) { origin = v; }
		
		
		/** To check whether 'this' region is inside a given box
		 * assuming the box's origins are (0,0,0).
		 * @param box The nD rectangular box.
		 * @return True if 'this' region is inside the box; Otherwise, false.
		 */
		bool is_region_in_box(const FloatSize & box) const;
		
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
