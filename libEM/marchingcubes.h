#ifndef _MARCHING_CUBES_H_
#define _MARCHING_CUBES_H_

#include <vector>

#include "vecmath.h"
#include "isosurface.h"

using std::vector;

namespace EMAN
{

	struct CubeNode {
		int level;
		int size;
		int xsize, ysize, zsize;
		CubeNode* children[8];
		float min, max;
		int num_children;
		int x,y,z;
		bool is_leaf;
	
		~CubeNode() {
			if(!is_leaf) {
				for(int i=0; i<num_children; i++)
					delete children[i];
			}
			//delete &children;
		}
	};
	
	class MarchingCubes : public Isosurface {
	public:
		MarchingCubes();
		MarchingCubes(EMData * em, bool isSmooth = false);
		virtual ~MarchingCubes();
	
		/**
		 * Sets Voxel data for Isosurface implementation
		 */
		void set_data(EMData* data);

		/**
		 * Set Isosurface value
		 */
		void set_surface_value(const float value);
	
		float get_surface_value() const ;
	
		/**
		 * Set Grid Size
		 */
		void set_sample_density(const int size);
	
		float get_sample_density() const ;
		
		Dict get_isosurface(bool smooth);
		
		unsigned long get_isosurface_dl(bool smooth);
		
	private:	
		int _sample;
		CubeNode* _root;
		map<int, int> point_map;
		unsigned long _isodl;
		
		bool calculate_min_max_vals(EMData* em );
		void clear_min_max_vals();
		vector<EMData*> minvals, maxvals;
		EMData* root;
		int levels, drawing_level;
		void draw_cube(const int x, const int y, const int z, const int cur_level );
		void marching_cube(int fX, int fY, int fZ, const EMData* const em);
		
		void calculate_surface(bool smooth);
		
		/** Find the gradient of the scalar field at a point. This gradient can 
		 * be used as a very accurate vertx normal for lighting calculations.
		 * 
		 * @param normal 
		 * @param fX 
		 * @param fY 
		 * @param fZ
		 */
		void get_normal(Vector3 &normal, int fX, int fY, int fZ);
		
		/** Find the approximate point of intersection of the surface between two 
		 * points with the values fValue1 and fValue2
		 * 
		 * @param fValue1 
		 * @param fValue2 
		 * @param fValueDesired 
		 * @return offset
		 */
		float get_offset(float fValue1, float fValue2, float fValueDesired);
		
		int get_edge_num(int x, int y, int z, int edge);

		template<typename type>
		class CustomVector
		{
			public:
			CustomVector(const int starting_size=1024) : data(0), size(0), elements(0)
			{
				resize(starting_size);
			}
			
			inline void clear(const int starting_size=1024)
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
			
			inline void resize(const int n)
			{
// 				cout << "Called resize with arg " << n << endl;
				
				data = (type*)realloc(data, n*sizeof(type));
				
				for(int i = size; i < n; ++i) data[i] = 0;
				size = n;
			}
			
			inline void push_back(const type& t)
			{
				if ( elements == size ) resize(2*size);
				data[elements] = t;
				++elements;
			}
			
			inline int elem() { return elements; }
			
			inline type& operator[](const int idx)
			{
				int csize = size;
				while (idx >= csize ) {
					csize *= 2;
				}
				if ( csize != size ) resize(csize);
				return data[idx];
			}
			
			inline type* get_data() { return data; }
			
			private:
			type* data;
			int size;
			int elements;
		};
		
// 		CustomVector p_map;
		CustomVector<float> pp;
		CustomVector<float> nn;
		CustomVector<int> ff;
	};

}

#endif
