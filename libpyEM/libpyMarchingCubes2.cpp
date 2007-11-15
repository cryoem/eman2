
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
    
    void set_surface_value(const float p0) {
        this->get_override("set_surface_value")(p0);
    }
    
    float get_surface_value() {
    	return this->get_override("get_surface_value")();
    }
    
    void set_sampling(float p0) {
        this->get_override("set_sample_density")(p0);
    }
    
    EMAN::Dict get_isosurface() {
       return this->get_override("get_isosurface")();
    }
#ifdef EMAN2_USING_OPENGL
	unsigned long get_isosurface_dl(unsigned int) {
		return this->get_override("get_isosurface_dl")();
	}
#endif
	
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
#ifdef EMAN2_USING_OPENGL
		.def("get_isosurface_dl", pure_virtual(&EMAN::Isosurface::get_isosurface_dl))
#endif //EMAN2_USING_OPENGL
		.def("get_sampling_range", pure_virtual(&EMAN::Isosurface::get_sampling_range))
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
        .def("setSampleDensity", &EMAN::Isosurface::setSampleDensity)
        .def("getSampleDensity", &EMAN::Isosurface::getSampleDensity)
    ;

    class_< EMAN::MarchingCubes, bases< EMAN::Isosurface >  >("MarchingCubes", init<  >())
        .def(init< const EMAN::MarchingCubes& >())
        .def(init< EMAN::EMData*, bool >())
    ;

}
*/
