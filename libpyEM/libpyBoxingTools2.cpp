
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include "boxingtools.h"

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {
	
	// Module ======================================================================
	BOOST_PYTHON_MODULE(libpyBoxingTools2)
	{
		scope* EMAN_BoxingTools_scope = new scope( 
		class_< EMAN::BoxingTools>("BoxingTools", init<  >())
		.def(init< const EMAN::BoxingTools& >())
		.def("get_min_delta_profile", &EMAN::BoxingTools::get_min_delta_profile)
		.def("is_local_maximum", &EMAN::BoxingTools::is_local_maximum)
		.def("auto_correlation_pick",&EMAN::BoxingTools::auto_correlation_pick)
		.def("set_radial_zero",&EMAN::BoxingTools::set_radial_zero)
		.staticmethod("get_min_delta_profile")
		.staticmethod("is_local_maximum")
		.staticmethod("auto_correlation_pick")
		.staticmethod("set_radial_zero")
		);
		
		delete EMAN_BoxingTools_scope;
	}
}
