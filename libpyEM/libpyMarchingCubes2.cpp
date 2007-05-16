
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include "marchingcubes.h"

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
    
    void set_surface_value(float p0) {
        this->get_override("set_surface_value")(p0);
    }
    
    float get_surface_value() {
    	return this->get_override("get_surface_value")();
    }
    
    void set_sample_density(float p0) {
        this->get_override("set_sample_density")(p0);
    }
    
    float get_sample_density() {
        return this->get_override("get_sample_density")();
    }
    
    EMAN::Dict get_isosurface(bool p0) {
       return this->get_override("get_isosurface")(p0);
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
		.def("set_sample_density", pure_virtual(&EMAN::Isosurface::set_sample_density))
		.def("get_sample_density", pure_virtual(&EMAN::Isosurface::get_sample_density))
		.def("get_isosurface", pure_virtual(&EMAN::Isosurface::get_isosurface))
		;
		
	class_< EMAN::MarchingCubes, bases<EMAN::Isosurface> >("MarchingCubes", init<  >())
		.def(init< EMAN::EMData *, optional< bool > >())
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
        .def("setSampleDensity", &EMAN::Isosurface::setSampleDensity)
        .def("getSampleDensity", &EMAN::Isosurface::getSampleDensity)
    ;

    class_< EMAN::MarchingCubes, bases< EMAN::Isosurface >  >("MarchingCubes", init<  >())
        .def(init< const EMAN::MarchingCubes& >())
        .def(init< EMAN::EMData*, bool >())
    ;

}
*/
