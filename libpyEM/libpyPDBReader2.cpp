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
#include <pdbreader.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyPDBReader2)
{
    scope* EMAN_PDBReader_scope = new scope(
    class_< EMAN::PDBReader >("PDBReader", init<  >())
        .def(init< const EMAN::PDBReader& >())
        .def(init< unsigned int >())
        .def("zero", &EMAN::PDBReader::zero)
        .def("copy", &EMAN::PDBReader::copy, return_value_policy< manage_new_object >())
        .def("get_number_points", &EMAN::PDBReader::get_number_points)
        .def("__len__", &EMAN::PDBReader::get_number_points)
        .def("set_number_points", &EMAN::PDBReader::set_number_points)
        .def("get_points", &EMAN::PDBReader::get_points)
        .def("read_from_pdb", &EMAN::PDBReader::read_from_pdb)
        .def("save_to_pdb", &EMAN::PDBReader::save_to_pdb)
        .def("set_points_array", &EMAN::PDBReader::set_points_array)
	.def("right_transform", &EMAN::PDBReader::right_transform)
	.def("makePointArray", &EMAN::PDBReader::makePointArray,return_value_policy< manage_new_object >())
	.def("get_x", &EMAN::PDBReader::get_x)
	.def("get_y", &EMAN::PDBReader::get_y)
	.def("get_z", &EMAN::PDBReader::get_z)
	.def("get_atomName", &EMAN::PDBReader::get_atomName)
	.def("get_resName", &EMAN::PDBReader::get_resName)
	.def("get_resNum", &EMAN::PDBReader::get_resNum)


    );

    delete EMAN_PDBReader_scope;

}

