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
		
		Dict get_isosurface(bool smooth) const ;
		
		unsigned long get_isosurface_dl(bool smooth);

		static const int leaf_level;
		
		int get_leaf_level() { return leaf_level; }
		int get_root_level();
	private:	
		int _sample;
		CubeNode* _root;
		map<int, int> point_map;
		unsigned long _isodl;
	
		void calculate_surface(bool smooth);
		
		/** MarchCube performs the Marching Cubes algorithm on a single cube 
		 * and BYPASS cubes without data.
		 * 
		 * @param fX 
		 * @param fY 
		 * @param fZ
		 * @param fxScale
		 * @param fyScale
		 * @param fzScale
		 */
		void marching_cube(int fX, int fY, int fZ, int fxScale, int fyScale, int fzScale );
		
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
		void build_search_tree();
		CubeNode* get_cube_node(int x, int y, int z, int level, int xsize, int ysize, int zsize);
		void draw_cube(CubeNode* node);
		
		/** Compute smooth normals. 
		 * */
		void calculate_smooth_normals();
	};

}

#endif
