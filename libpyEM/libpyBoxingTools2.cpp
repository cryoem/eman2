
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
		.def("set_radial_non_zero",&EMAN::BoxingTools::set_radial_non_zero)
		.def("find_radial_max",&EMAN::BoxingTools::find_radial_max)
		.def("classify",&EMAN::BoxingTools::classify)
		.def("get_color",&EMAN::BoxingTools::get_color)
		.def("set_mode",&EMAN::BoxingTools::set_mode)
		.def("set_region",&EMAN::BoxingTools::set_region)
		.staticmethod("get_min_delta_profile")
		.staticmethod("is_local_maximum")
		.staticmethod("auto_correlation_pick")
		.staticmethod("set_radial_non_zero")
		.staticmethod("find_radial_max")
		.staticmethod("classify")
		.staticmethod("get_color")
		.staticmethod("set_mode")
		.staticmethod("set_region")
		);
		
		
		enum_< EMAN::BoxingTools::CmpMode >("CmpMode")
		.value("SWARM_DIFFERENCE", EMAN::BoxingTools::SWARM_DIFFERENCE)
		.value("SWARM_RATIO", EMAN::BoxingTools::SWARM_RATIO)
		;
		
		delete EMAN_BoxingTools_scope;
	}
}
