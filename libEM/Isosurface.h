#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_

#include "Mesh.h"
#include "VoxelData.h"

namespace EMAN
{

	class Isosurface {
	public:
		Isosurface(){}
		~Isosurface(){}
		
		/**
		 * Sets Voxel data for Isosurface implementation
		 */
		virtual void setVolumeData(const VoxelData& data) = 0;
	
		/**
		 * Set Isosurface value
		 */
		virtual void setSurfaceValue(const float value) = 0;
	
		virtual float getSurfaceValue() const = 0;
	
		/**
		 * Set Grid Size
		 */
		virtual void setSampleDensity(const int size) = 0;
	
		virtual float getSampleDensity() const = 0;
	
		/**
		 * Get Isosurface mesh
		 */
		virtual const Mesh& getMesh() const = 0;
	
	};

}

#endif
