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
#include "boxingtools.h"

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

	// Module ======================================================================
	BOOST_PYTHON_MODULE(libpyBoxingTools2)
	{
		scope* EMAN_BoxingTools_scope = new scope(
		class_< EMAN::BoxingTools>("BoxingTools",
				"BoxingTools is class for encapsulating common\n"
				"boxing operations that may become expensive if\n"
				"they are implemented in python.",
				init<  >())
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
		.value("SWARM_AVERAGE_RATIO", EMAN::BoxingTools::SWARM_AVERAGE_RATIO)
		;

		delete EMAN_BoxingTools_scope;
	}
}
