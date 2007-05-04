#ifndef _MARCHING_CUBES_H_
#define _MARCHING_CUBES_H_

#include <vector>

#include "vecmath.h"
#include "Isosurface.h"

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
		MarchingCubes(EMData * em);
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
	
		/**
		 * Get Isosurface mesh
		 */
//		const Mesh& getMesh() const ;
	
//		void drawMesh(bool smooth) const { _mesh->draw(smooth); }
	
	private:
		MarchingCubes(){};
		
		//Mesh* _mesh;
		float _surf_value;
		int _sample;
		//VoxelData _voxel;
		//EMData * _emdata;
		CubeNode* _root;
		vector<float> *points, *normals, *normalsSm;
		vector<int> *faces;
		map<int, int> point_map;
	
		void calculateSurface();
		void marchingCube(int fX, int fY, int fZ, int fScale);
		void getNormal(Vector3 &normal, int fX, int fY, int fZ);
		float getOffset(float fValue1, float fValue2, float fValueDesired);
		int getEdgeNum(int x, int y, int z, int edge);
		void buildSearchTree();
		CubeNode* getCubeNode(int x, int y, int z, int level, int size);
		void drawCube(CubeNode* node);
		
		void calculate_smooth_normals();
	};

}

#endif
