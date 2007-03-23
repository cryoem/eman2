#ifndef _MARCHING_CUBES_H_
#define _MARCHING_CUBES_H_

#include "Isosurface.h"
#include <vector>
#include <hash_map>

using std::hash_map;
using std::vector;

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
	MarchingCubes();
	~MarchingCubes();

	/**
	 * Sets Voxel data for Isosurface implementation
	 */
	void setVolumeData(const VoxelData& data);

	void loadMRC(std::string file);

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
	const Mesh& getMesh() const ;

	void drawMesh(bool smooth) const { _mesh->draw(smooth); }

private:
	Mesh* _mesh;
	float _surf_value;
	int _sample;
	VoxelData _voxel;
	CubeNode* _root;
	vector<float> *points, *normals;
	vector<int> *faces;
	hash_map<int, int> point_map;

	void calculateSurface();
	void marchingCube(int fX, int fY, int fZ, int fScale);
	void getNormal(Vector3 &normal, int fX, int fY, int fZ);
	float getOffset(float fValue1, float fValue2, float fValueDesired);
	int getEdgeNum(int x, int y, int z, int edge);
	void buildSearchTree();
	CubeNode* getCubeNode(int x, int y, int z, int level, int size);
	void drawCube(CubeNode* node);
};

#endif
