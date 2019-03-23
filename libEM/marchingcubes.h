/*
 * Author: Tao Ju, 5/16/2007 <taoju@cs.wustl.edu>, code ported by Grant Tang
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

#ifndef _MARCHING_CUBES_H_
#define _MARCHING_CUBES_H_

#include <vector>
using std::vector;

#include "vecmath.h"
#include "isosurface.h"

// Marching cubes debug will turn on debug and timing information
#define MARCHING_CUBES_DEBUG 0

#include <ostream>
using std::ostream;

#include <climits>
// for CHAR_BIT

#ifdef __APPLE__
	#include "OpenGL/gl.h"
	#include "OpenGL/glext.h"
#else // WIN32, LINUX
	#include "GL/gl.h"
	#include "GL/glext.h"
#endif	//__APPLE__

namespace EMAN
{
	/** CustomVector has some trivial optimizations of the STL vector.
	* It has get_data() functionality, which gives direct access to
	* the underlying data, which is necessary for handing OpenGL vertex and normal
	* arrays - the native STL vector does not give straightforward access to the data
	* pointers. This class also has a push_back_3 function which does a memcpy
	* of 3 elements - this can be used instead of 3 push_back function calls. Although it wasn't
	* rigorously tested it should save some time by virtue of less function calls. Although the
	* savings may be trivial, I just haven't tested it.
	*
	* This class was not written for general use because- it is designed specifically for
	* use in conjunction with the MarchingCubes class only.
	* @author David Woolford
	*/
	template<typename type>
	class CustomVector
	{
		public:
			/** Constructor
			* @param starting_size the starting size of the underlying data
			*/
			explicit CustomVector(unsigned int starting_size=1024) : data(0), size(0), elements(0)
			{
				resize(starting_size);
			}

			~CustomVector()
			{
				if(data) {free(data); data=0;}
			}

			/** Clear
			* clears all data and resizes
			* @param starting_size the starting size of the underlying data
			*/
			inline void clear(unsigned int starting_size=1024)
			{
				if (data) {free(data); data = 0;}
				size = 0;
				elements = 0;
				resize(starting_size);
			}

			/** Resize
			* Resize the data pointer using realloc
			* In general you want to resize to a larger size...
			* @param n the new size of the underlying data
			*
			*/
			inline void resize(const unsigned int n)
			{
				data = (type*)realloc(data, n*sizeof(type));

				for(unsigned int i = size; i < n; ++i) data[i] = 0;
				size = n;
			}

			/** Push back a value
			* Potentially resizes
			* @param t the value to push back
			*/
			inline void push_back(const type& t)
			{
				if ( elements == size ) resize(2*size);
				data[elements] = t;
				++elements;
			}

			/** Push back 3 values
			 * Potentially resizes
			 * @param p a pointer to 3 contiguous values
			 */
			inline void push_back_3(const type* const p)
			{
				if ( elements+2 >= size ) resize(2*size);
				memcpy( &data[elements], p, 3*sizeof(type));
				elements = elements + 3;
			}

			/** Multiply contiguous groups of 3 by different values
			 */
			inline void mult3(const type& v1, const type& v2, const type& v3)
			{
				for(unsigned int i = 0; (i + 2) < elements; i += 3 ){
					data[i] *= v1;
					data[i+1] *= v2;
					data[i+2] *= v3;
				}
			}

			/** The number of elements, is the same as STL::vector size()
			* Should be called size() but oh well...
			* This is the actual number of stored floating point numbers
			* not the number of 'groups of 3'
			* @return the number of elements stored by this object
			*/
			inline unsigned int elem() { return elements; }

			/** Operator[] - provide convenient set functionality (note NOT get)
			* @param idx the index to access for setting purposes
			* potentially resizes
			*/
			inline type& operator[](const unsigned int idx)
			{
				unsigned int csize = size;
				while (idx >= csize ) {
					csize *= 2;
				}
				if ( csize != size ) resize(csize);
				return data[idx];
			}

			/** get the underlying data as a pointer objects
			* @return the data pointer
			*/
			inline type* get_data() { return data; }

		private:
			/// The data pointer
			type* data;
			/// The size of the associated memory block
			unsigned int size;
			/// The number of stored elements.
			unsigned int elements;
	};

	/** Class to encapsulate an RGB color generator for marching cubes isosurface generator
	For now you can only color by radius, but in the near future you will be able to color by map, etc
	 **/
	class ColorRGBGenerator{
	public:
		ColorRGBGenerator();
		/** Constructor */
		ColorRGBGenerator(EMData* emdata);
		
		/** Generate a color based on pixel coords*/
		float* getRGBColor(int x, int y, int z);
		
		/** set the emdata */
		void set_data(EMData* data);
		
		/** Generate a color map **/
		void generateRadialColorMap();
		
		/** Set min max data */
		void set_cmap_data(EMData* data);
		
		/** Set origin */
		inline void setOrigin(int orix, int oriy, int oriz)
		{
			originx = orix;
			originy = oriy;
			originz = oriz;
			needtorecolor = true;
		}
		/** Set scaling */
		inline void setScale(float i, float o)
		{
			inner = i;
			outer = o;
			needtorecolor = true;
		}
		/** Set RGB mode 0 = none, 1 = color by radius, more to come :) */
		inline void setRGBmode(int mode)
		{
			rgbmode = mode;
			needtorecolor = true;
		}
		
		/** Set the mn max for cmap coloring */
		inline void setMinMax(float min, float max)
		{
			minimum = min;
			maximum = max;
			needtorecolor = true;
		}
		
		/** Return RGB mode */
		inline int getRGBmode()
		{
			return rgbmode;
		}
		
		/** Lets us know if we need to recalclate colors*/
		inline bool getNeedToRecolor()
		{
			return needtorecolor;
		}
		
		inline void setNeedToRecolor(bool recolor)
		{
			needtorecolor = recolor;
		}
	private:
		int rgbmode;
		int originx;
		int originy;
		int originz;
		float inner;
		float outer;
		float minimum;
		float maximum;
	
		bool needtorecolor;	// dirty bit to let the sytem know when we need to recolor
		
		float* colormap; 	// pointer to a colormap
		EMData* em_data;	// pointer to EMdata
		EMData* cmap;		// pointer to colormap data
		
		float rgb[3];		//little array to hold RGB values;
	};
	
	class MarchingCubes : public Isosurface {
		friend class GLUtil;
	public:
		/** Default constructor
		*/
		MarchingCubes();

		/** Most commonly used constructor
		* calls set_data(em)
		* @param em the EMData object to generate triangles and normals for
		*/
		MarchingCubes(EMData * em);
		virtual ~MarchingCubes();

		/** Sets Voxel data for Isosurface implementation
		* Calls calculate_min_max_vals which generates the tree of data
		* @param data the emdata object to be rendered in 3D
		* @exception ImageDimensionException if the image z dimension is 1
		*/
		void set_data(EMData* data);

		/** Set Isosurface value
		* @param value the new isosurface value
		*/
		void set_surface_value(const float value);

		/** Get the current isosurface value being used
		* @return the current isosurface value
		*/
		float get_surface_value() const { return _surf_value; }

		/** Set sampling rates
		* A smaller value means a finer sampling.
		* The finest sampling level is -1
		* Sampling values increment in steps of 1, and a single increment
		* is interpreted as one step up or down the tree stored in minvals
		* and maxvals
		* @param rate the tree level to render
		 */
		void set_sampling(const int rate) { drawing_level = rate; }

		/** Current the current sampling rate
		* Finest sampling is -1.
		*/
		int get_sampling() const { return drawing_level; }

		/** Get the range of feasible sampling rates
		*/
		int get_sampling_range() { return minvals.size()-1; }

		/** Color the vertices
		 */
		void color_vertices();
		
		/** Get the isosurface as dictionary
		* Traverses the tree and marches the cubes
		* @return a dictionary object containing to float pointers (to vertex and normal data), and an int pointer (to face data)
		*/
		Dict get_isosurface();

		void surface_face_z();
		
		/** Functions to control colroing mode
		 * */
		inline void setRGBorigin(int x, int y, int z)
		{
			rgbgenerator.setOrigin(x, y, z);
		}
		
		inline void setRGBscale(float i, float o)
		{
			rgbgenerator.setScale(i, o);
		}
		
		inline void setRGBmode(int mode)
		{
			rgbgenerator.setRGBmode(mode);
		}
		
		/** Return RGB mode */
		inline int getRGBmode()
		{
			return rgbgenerator.getRGBmode();
		}
		
		/** Sets the colormap */
		inline void setCmapData(EMData* data)
		{
			rgbgenerator.set_cmap_data(data);
		}
		
		/** Sets the colormap mix max range */
		inline void setCmapMinMax(float min, float max)
		{
			rgbgenerator.setMinMax(min, max);
		}
		
	private:
		map<int, int> point_map;
		unsigned long _isodl;
		GLuint buffer[4];

		/** Calculate the min and max value trees
		* Stores minimum and maximum cube neighborhood values in a tree structure
		* @exception NullPointerException if _emdata is null... this should not happen but is left for clarity for
		* programmers
		*/
		bool calculate_min_max_vals();
		
	
		/** Clear the minimum and maximum value search trees
		* Frees memory in the minvals and maxvals
		*/
		void clear_min_max_vals();

		/// Vectors for storing the search trees
		vector<EMData*> minvals, maxvals;

		/// The "sampling rate"
		int drawing_level;

		/** The main cube drawing function
		* To start the process of generate triangles call with draw_cube(0,0,0,minvals.size()-1)
		* Once cur_level becomes drawing_level marching_cube is called
		* @param x the current x value, relative to cur_level
		* @param y the current y value, relative to cur_level
		* @param z the current z value, relative to cur_level
		* @param cur_level the current tree traversal level
		*
		*/
		void draw_cube(const int x, const int y, const int z, const int cur_level );


		/** Function for managing cases where a triangles can potentially be rendered
		* Called by draw_cube. Generates vertices, normals, and keeps track of common points
		* @param fX the current x coordinate, relative to cur_level
		* @param fY the current y coordinate, relative to cur_level
		* @param fZ the current z coordinate, relative to cur_level
		* @param cur_level
		*/
		void marching_cube(int fX, int fY, int fZ, const int cur_level);

		/** Calculate and generate the entire set of vertices and normals using current states
		 * Calls draw_cube(0,0,0,minvals.size()-1)
		*/
		void calculate_surface();
		
		/** Find the approximate point of intersection of the surface between two
		 * points with the values fValue1 and fValue2
		 *
		 * @param fValue1
		 * @param fValue2
		 * @param fValueDesired
		 * @return offset
		 */
		float get_offset(float fValue1, float fValue2, float fValueDesired);

		/** Get edge num
		* needs better commenting
		*/
		int get_edge_num(int x, int y, int z, int edge);

		/** Find the gradient of the scalar field at a point. This gradient can
		 * be used as a very accurate vertx normal for lighting calculations.
		 *
		 * THIS FUNCTION IS NO LONGER CALLED - d.woolford
		 * but is retained because it may be useful, perhaps even for saving time
		 * @param normal where to store the normal
		 * @param fX
		 * @param fY
		 * @param fZ
		 */
		void get_normal(Vector3 &normal, int fX, int fY, int fZ);

		///.Custom vectors for storing points, normals and faces
		CustomVector<float> pp;
		CustomVector<float> cc;
		CustomVector<int> vv;
		CustomVector<float> nn;
		CustomVector<unsigned int> ff;
		
		/** COlor by radius generator */
		ColorRGBGenerator rgbgenerator;
		
		bool needtobind;	// A dirty bit to signal when the the MC algorithm or color has chaged and hence a need to update GPU buffers
	};

	/** @author David Woolford
	 * @date March 2009
	 */
//	class CudaMarchingCubes : public Isosurface {
//		public:
//
//			/** Most commonly used constructor
//			* calls set_data(em)
//			* @param em the EMData object to generate triangles and normals for
//			*/
//			CudaMarchingCubes(EMData * em);
//			virtual ~CudaMarchingCubes();
//			/**
//			 * Set Isosurface value
//			 */
//			virtual void set_surface_value(const float value);
//
//			/** Sets Voxel data for Isosurface implementation
//			* Calls calculate_min_max_vals which generates the tree of data
//			* @param data the emdata object to be rendered in 3D
//			* @exception ImageDimensionException if the image z dimension is 1
//			*/
//			void set_data(EMData* data);
//
//			virtual float get_surface_value() const { return _surf_value; }
//
//			/**
//			 * Set Grid Size
//			 */
//			virtual void set_sampling(const int size) {} ;
//
//			virtual int get_sampling() const { return 0; };
//
//			/** Get the number of feasible samplings
//			*
//			 */
//			virtual int get_sampling_range() {return 0; }
//
//			virtual Dict get_isosurface()  { return Dict;}
//
//	}
	
	/** A work in progress by David Woolford
	*/
	class U3DWriter{
	public:
		typedef unsigned int U32;
		typedef long unsigned int U64;
		typedef double F64;
		typedef float F32;
		typedef short int I16;
		typedef short unsigned int U16;
		typedef unsigned char U8;



		U3DWriter();
		~U3DWriter();

		int write(const string& filename);

		ostream& write(ostream&);

	private:
		unsigned int size_of_in_bytes();

		void test_type_sizes();

		ostream& write_header(ostream&);
		ostream& write_clod_mesh_generator_node(ostream& os);


		template<typename type>
		ostream& write(ostream& os, const type& T) {
			os.write( (const char*)(&T), sizeof(type) );
			return os;
		}


// 		U8 reverse_bits_3d( U8 val ) {
// 			U8 ret = 0;
// 			U8 n_bits = sizeof ( val ) * CHAR_BIT;
//
// 			for ( unsigned i = 0; i < n_bits; ++i ) {
// 				ret = ( ret << 1 ) | ( val & 1 );
// 				val >>= 1;
// 			}
//
// 			return ret;
// 		}

// 		U16 reverse_bits_u3d( U16 val ) {
// 			bool is_big_endian = ByteOrder::is_host_big_endian();
//
// 			U16 ret = 0;
//
// 			U8* pu8 = (U8*)&val;
// 			U8* pu8ret = (U8*)&ret;
// 			U8 high_order = 0;
// 			U8 low_order = 0;
// 			if (is_big_endian) {
// 				high_order = pu8[0];
// 				low_order = pu8[1];
// 			}
// 			else {
// 				high_order = pu8[1];
// 				low_order = pu8[0];
// 			}
//
// 			pu8ret[0] = low_order;
// 			pu8ret[1] = higher_order;
//
// 			high_order = reverse_bits_u3d(high_order);
// 			low_order = reverse_bits_u3d(low_order);
//
// 			return ret;
//
// 		}
//
// 		U32 reverse_bits_u3d( U32 val ) {
// 			bool is_big_endian = ByteOrder::is_host_big_endian();
//
// 			U32 ret = 0;
//
// 			U16 * pu16 = (U16*)&val;
// 			U16* pu16ret = (U16*)&ret;
// 			U16 high_order = 0;
// 			U16 low_order = 0;
// 			if (is_big_endian) {
// 				high_order = pu16[0];
// 				low_order = pu16[1];
// 			}
// 			else {
// 				high_order = pu16[1];
// 				low_order = pu16[0];
// 			}
//
// 			high_order = reverse_bits_u3d(high_order);
// 			low_order = reverse_bits_u3d(low_order);
//
// 			pu16ret[0] = low_order;
// 			pu16ret[1] = higher_order;
//
// 			return ret;
//
// 		}
//

		U32 DIFFUSE_COLOR_COUNT;
		U32 SPECULAR_COLOR_COUNT;


		CustomVector<F32> pp;
		CustomVector<F32> nn;
		CustomVector<unsigned int> ff;
	};

	// Template specialization. Have to be careful when dealing with strings
	template<> ostream& U3DWriter::write(ostream& os, const string& );

// 	class U3DBitStreamWriter {
// 	public:
// 		typedef unsigned int U32;
// 		typedef long unsigned int U64;
// 		typedef double F64;
// 		typedef float F32;
// 		typedef short int I16;
// 		typedef short unsigned int U16;
// 		typedef unsigned char U8;
//
//
// 		U3D_BitStreamWriter();
// 		~U3D_BitStreamWriter();
//
// 		// We could do this kind of thing with templates,  but
// 		// we would have to specialize each function anyhow,
// 		// so it makes sense just to do it the C way.
// 		void writeU8( U8 value );
// 		void writeU16( U16 value );
// 		void writeU32( U32 value );
// 		void writeU64( U64 value );
// 		void writeI32( I32 value );
// 		void writeF32( F32 value );
// 		void writeF64( F64 value );
//
// 		void WriteCompressedU32( U32 context, U32 value );
// 		void WriteCompressedU16( U32 context, U32 value );
// 		void WriteCompressedU8( U32 context, U32 value );
//
// 	};



}

#endif
