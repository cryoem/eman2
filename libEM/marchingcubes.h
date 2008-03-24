#ifndef _MARCHING_CUBES_H_
#define _MARCHING_CUBES_H_

#include <vector>

#include "vecmath.h"
#include "isosurface.h"

// Marching cubes debug will turn on debug and timing information
#define MARCHING_CUBES_DEBUG 0

using std::vector;

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
	*/
	template<typename type>
	class CustomVector
	{
		public:
			/** Constructor
			* @param starting_size the starting size of the underlying data
			*/
			explicit CustomVector(int starting_size=1024) : data(0), size(0), elements(0)
			{
				resize(starting_size);
			}
			
			/** Clear
			* clears all data and resizes
			* @param starting_size the starting size of the underlying data 
			*/
			inline void clear(int starting_size=1024)
			{
				if ( data != 0 )
				{
					delete [] data;
					data = 0;
				}
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
			inline void resize(const int n)
			{
				data = (type*)realloc(data, n*sizeof(type));
				
				for(int i = size; i < n; ++i) data[i] = 0;
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
				for(int i = 0; (i + 2) < elements; i += 3 ){
					data[i] *= v1;
					data[i+1] *= v2;
					data[i+2] *= v3;
				}
			}
			
			/** The number of elements, is the same as STL::vector size()
			* Should be called size() but oh well...
			* @return the number of elements stored by this object
			*/
			inline int elem() { return elements; }
			
			/** Operator[] - provide convenient set functionality (note NOT get)
			* @param idx the index to access for setting purposes
			* potentially resizes
			*/
			inline type& operator[](const int idx)
			{
				int csize = size;
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
			int size;
			/// The number of stored elements.
			int elements;
	};
	
	class MarchingCubes : public Isosurface {
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
		
		/** Get the isosurface as dictionary
		* Traverses the tree and marches the cubes
		* @return a dictionary object containing to float pointers (to vertex and normal data), and an int pointer (to face data)
		*/
		Dict get_isosurface();

#ifdef EMAN2_USING_OPENGL
		/** Get an isosurface display list
		* Traverses the tree, marches the cubes, and renders a display list using the associated vertices and normals
		* Uses OpenGL arrays for maximum performance
		* @return an OpenGL display list number
		*/
		unsigned long get_isosurface_dl(unsigned int tex_id = 0);
#endif //EMAN2_USING_OPENGL
		
	private:	
		map<int, int> point_map;
		unsigned long _isodl;
		
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
		* @param x the current y value, relative to cur_level
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
		CustomVector<float> nn;
		CustomVector<int> ff;
	};

}

#endif
