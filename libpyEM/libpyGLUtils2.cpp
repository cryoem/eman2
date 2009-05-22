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

#ifdef EMAN2_USING_OPENGL

// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

#include "glutil.h"
#include "emdata.h"
#include "marchingcubes.h"

// Using =======================================================================
using namespace boost::python;

namespace {
	BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_GLUtil_get_isosurface_dl_overloads_1_2, EMAN::GLUtil::get_isosurface_dl, 1, 2)
}

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyGLUtils2)
{
	scope* EMAN_GLUtil_scope = new scope(
		 class_< EMAN::GLUtil>("GLUtil", init<  >())
		.def(init< const EMAN::GLUtil& >())
		.def("gen_glu_mipmaps", &EMAN::GLUtil::gen_glu_mipmaps)
		.def("gen_gl_texture", &EMAN::GLUtil::gen_gl_texture)
		.def("render_amp8_gl_texture", &EMAN::GLUtil::render_amp8_gl_texture)
		.def("nearest_projected_points", &EMAN::GLUtil::nearest_projected_points )
		.def("colored_rectangle", &EMAN::GLUtil::colored_rectangle )
		.def("mx_bbox", &EMAN::GLUtil::mx_bbox )
		.def("render_amp8", &EMAN::GLUtil::render_amp8)
		.def("get_isosurface_dl", &EMAN::GLUtil::get_isosurface_dl, EMAN_GLUtil_get_isosurface_dl_overloads_1_2())
		.staticmethod("gen_glu_mipmaps")
		.staticmethod("gen_gl_texture")
		.staticmethod("render_amp8_gl_texture")
		.staticmethod("nearest_projected_points")
		.staticmethod("colored_rectangle")
		.staticmethod("mx_bbox")
		.staticmethod("render_amp8")
		.staticmethod("get_isosurface_dl")
	);

	delete EMAN_GLUtil_scope;

}

#endif //EMAN2_USING_OPENGL

