#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_

#include "emdata.h"

namespace EMAN
{

	class Isosurface {
	public:
		Isosurface() : _emdata(0), points(0), normals(0), normalsSm(0), faces(0), isSmooth(false), _surf_value(1) {}
		virtual ~Isosurface(){
			if(points) {delete points; points=0;}
			if(normals) {delete normals; normals=0;}
			if(normalsSm) {delete normalsSm; normalsSm=0;}
			if(faces) {delete faces; faces=0;}
		}
		
		/**
		 * Sets Voxel data for Isosurface implementation
		 */
		virtual void set_data(EMData* data) {
			_emdata = data;
		}
	
		/**
		 * Set Isosurface value
		 */
		virtual void set_surface_value(const float value) = 0;
	
		virtual float get_surface_value() const = 0;
	
		/**
		 * Set Grid Size
		 */
		virtual void set_sample_density(const int size) = 0;
	
		virtual float get_sample_density() const = 0;
		
		virtual Dict get_isosurface(bool smooth) const =0;
	
	protected:
		EMData * _emdata;
		
		vector<float> *points, *normals, *normalsSm;
		vector<int> *faces;
		
		bool isSmooth;	//boolean to indicate if smooth normals needed
		float _surf_value;	
	};

}

#endif
