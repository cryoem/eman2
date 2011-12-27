/*
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

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include "marchingcubes.h"
#include "emdata.h"

// Using =======================================================================
using namespace boost::python;
//using namespace EMAN;

// Declarations ================================================================
namespace  {

struct EMAN_Isosurface_Wrapper: EMAN::Isosurface, wrapper<EMAN::Isosurface>
{
    void default_set_data(EMAN::EMData* p0) {
        EMAN::Isosurface::set_data(p0);
    }

    void set_data(EMAN::EMData* p0) {
        if(override set_data = this->get_override("set_data")) {
        	EMAN::Isosurface::set_data(p0);
        }
    }

    void set_surface_value(const float p0) {
        this->get_override("set_surface_value")(p0);
    }

    float get_surface_value() {
    	return this->get_override("get_surface_value")();
    }

    void set_sampling(float p0) {
        this->get_override("set_sample_density")(p0);
    }

    void setRGBorigin(int x, int y, int z) {
        this->get_override("set_rgb_scale")(x,y,z);
    }
    
    void set_rgb_scale(float i, float o) {
        this->get_override("set_rgb_scale")(i, o);
    }
    
    void set_rgb_mode(int mode) {
        this->get_override("set_rgb_mode")(mode);
    }
    
    void set_cmap_data(EMAN::EMData* data) {
        this->get_override("set_cmap_data")(data);
    }
    
    void set_cmap_minmax(float min, float max) {
        this->get_override("set_cmap_minmax")(min, max);
    }
    
    EMAN::Dict get_isosurface() {
       return this->get_override("get_isosurface")();
    }

	int get_sampling_range() {
		return this->get_override("get_sampling_range")();
	}

};

}

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyMarchingCubes2)
{
	class_<EMAN_Isosurface_Wrapper, boost::noncopyable>("Isosurface", no_init)
		.def("set_data", &EMAN::Isosurface::set_data, &EMAN_Isosurface_Wrapper::default_set_data)
		.def("set_surface_value", pure_virtual(&EMAN::Isosurface::set_surface_value))
		.def("get_surface_value", pure_virtual(&EMAN::Isosurface::get_surface_value))
		.def("set_sampling", pure_virtual(&EMAN::Isosurface::set_sampling))
		.def("get_sampling", pure_virtual(&EMAN::Isosurface::get_sampling))
		.def("get_isosurface", pure_virtual(&EMAN::Isosurface::get_isosurface))
		.def("get_sampling_range", pure_virtual(&EMAN::Isosurface::get_sampling_range))
		.def("set_rgb_origin", pure_virtual(&EMAN::Isosurface::setRGBorigin))
		.def("set_rgb_scale", pure_virtual(&EMAN::Isosurface::setRGBscale))
		.def("set_rgb_mode",  pure_virtual(&EMAN::Isosurface::setRGBmode))
		.def("set_cmap_data",  pure_virtual(&EMAN::Isosurface::setCmapData))
		.def("set_cmap_minmax",  pure_virtual(&EMAN::Isosurface::setCmapMinMax))
		;

	/* We do not wrap default constructor of MarchingCubes into Python */
	//class_< EMAN::MarchingCubes, bases<EMAN::Isosurface> >("MarchingCubes", init<  >())
	class_< EMAN::MarchingCubes, bases<EMAN::Isosurface> >("MarchingCubes", init< EMAN::EMData *>())
		//.def(init< EMAN::EMData *, optional< bool > >())
		;
}

/*
// Module ======================================================================
BOOST_PYTHON_MODULE(libpyMarchingCubes2)
{
	scope* EMAN_MarchingCubes_scope = new scope(

	);
   class_< EMAN::Isosurface, boost::noncopyable, EMAN::MarchingCubes
    		>("Isosurface", init<  >())
        .def(init< const EMAN::Isosurface& >())
        .def("setSurfaceValue", &EMAN::Isosurface::setSurfaceValue)
        .def("get_isosurface", &EMAN::Isosurface::get_isosurface)
        .def("setVolumeData", &EMAN::Isosurface::setVolumeData)
        .def("getSurfaceValue", &EMAN::Isosurface::getSurfaceValue)
        .def("set_sample_density", &EMAN::Isosurface::set_sample_density)
        .def("getSampleDensity", &EMAN::Isosurface::getSampleDensity)
    ;

    class_< EMAN::MarchingCubes, bases< EMAN::Isosurface >  >("MarchingCubes", init<  >())
        .def(init< const EMAN::MarchingCubes& >())
        .def(init< EMAN::EMData*, bool >())
    ;

}
*/
