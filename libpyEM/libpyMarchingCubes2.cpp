
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include "MarchingCubes.h"
//#include "Isosurface.h"

// Using =======================================================================
using namespace boost::python;
//using namespace EMAN;

// Declarations ================================================================
namespace  {

struct EMAN_Isosurface_Wrapper: EMAN::Isosurface, wrapper<EMAN::Isosurface>
{  
    void default_setVolumeData(EMAN::EMData* p0) {
        EMAN::Isosurface::setVolumeData(p0);
    }
        
	void setVolumeData(EMAN::EMData* p0) {
        if(override setVolumeData = this->get_override("setVolumeData")) {
        	EMAN::Isosurface::setVolumeData(p0);
        }
    }
    
    void setSurfaceValue(float p0) {
        this->get_override("setSurfaceValue")(p0);
    }
    
    float getSurfaceValue() {
    	return this->get_override("getSurfaceValue")();
    }
    
    void setSampleDensity(float p0) {
        this->get_override("setSampleDensity")(p0);
    }
    
    float getSampleDensity() {
        return this->get_override("getSampleDensity")();
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
		.def("setVolumeData", &EMAN::Isosurface::setVolumeData, &EMAN_Isosurface_Wrapper::default_setVolumeData)
		.def("setSurfaceValue", pure_virtual(&EMAN::Isosurface::setSurfaceValue))
		.def("getSurfaceValue", pure_virtual(&EMAN::Isosurface::getSurfaceValue))
		.def("setSampleDensity", pure_virtual(&EMAN::Isosurface::setSampleDensity))
		.def("getSampleDensity", pure_virtual(&EMAN::Isosurface::getSampleDensity))
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
