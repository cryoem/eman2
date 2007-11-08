#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_

#include "emdata.h"

namespace EMAN
{

	class Isosurface {
	public:
		Isosurface() : _emdata(0), _surf_value(1) {}
		virtual ~Isosurface(){}
		
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
		virtual void set_sampling(const int size) = 0;
	
		virtual int get_sampling() const = 0;
		
		/** Get the number of feasible samplings
		*
		 */
		virtual int get_sampling_range() = 0;
		
		virtual Dict get_isosurface()  = 0;
#ifdef EMAN2_USING_OPENGL
		virtual unsigned long get_isosurface_dl() = 0;
#endif
	
	protected:
		EMData * _emdata;
		
		float _surf_value;	
	};

}

#endif
