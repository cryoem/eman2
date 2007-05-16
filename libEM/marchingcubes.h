#ifndef _MARCHING_CUBES_H_
#define _MARCHING_CUBES_H_

#include <vector>

#include "vecmath.h"
#include "isosurface.h"

using std::vector;

namespace EMAN
{

	struct CubeNode {
		int level, size;
		CubeNode* children[8];
		float min, max;
		int x,y,z;
		bool is_leaf;
	
		~CubeNode() {
			if(!is_leaf) {
				for(int i=0; i<8; i++)
					delete children[i];
			}
			//delete &children;
		}
	};
	
	class MarchingCubes : public Isosurface {
	public:
		MarchingCubes() : _root(0) {};
		MarchingCubes(EMData * em, bool isSmooth = false);
		virtual ~MarchingCubes();
	
		/**
		 * Sets Voxel data for Isosurface implementation
		 */
		void setVolumeData(EMData* data);

		/**
		 * Set Isosurface value
		 */
		void setSurfaceValue(const float value);
	
		float getSurfaceValue() const ;
	
		/**
		 * Set Grid Size
		 */
		void setSampleDensity(const int size);
	
		float getSampleDensity() const ;
		
		Dict get_isosurface(bool smooth) const ;
	
	private:	
		int _sample;
		CubeNode* _root;
		map<int, int> point_map;
	
		void calculateSurface(bool smooth);
		
		/** MarchCube performs the Marching Cubes algorithm on a single cube 
		 * and BYPASS cubes without data.
		 * 
		 * @param fX 
		 * @param fY 
		 * @param fZ
		 * @param fScale
		 */
		void marchingCube(int fX, int fY, int fZ, int fScale);
		
		/** Find the gradient of the scalar field at a point. This gradient can 
		 * be used as a very accurate vertx normal for lighting calculations.
		 * 
		 * @param normal 
		 * @param fX 
		 * @param fY 
		 * @param fZ
		 */
		void getNormal(Vector3 &normal, int fX, int fY, int fZ);
		
		/** Find the approximate point of intersection of the surface between two 
		 * points with the values fValue1 and fValue2
		 * 
		 * @param fValue1 
		 * @param fValue2 
		 * @param fValueDesired 
		 * @return offset
		 */
		float getOffset(float fValue1, float fValue2, float fValueDesired);
		
		int getEdgeNum(int x, int y, int z, int edge);
		void buildSearchTree();
		CubeNode* getCubeNode(int x, int y, int z, int level, int size);
		void drawCube(CubeNode* node);
		
		/** Compute smooth normals. 
		 * */
		void calculate_smooth_normals();
	};

}

#endif
