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
#include <emdata.h>
#include <sparx/fundamentals.h>

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {


}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFundamentals2)
{
    enum_< EMAN::fp_flag >("fp_flag")
        .value("CIRCULANT", EMAN::CIRCULANT)
        .value("CIRCULANT_NORMALIZED", EMAN::CIRCULANT_NORMALIZED)
        .value("PADDED", EMAN::PADDED)
        .value("PADDED_NORMALIZED", EMAN::PADDED_NORMALIZED)
        .value("PADDED_LAG", EMAN::PADDED_LAG)
        .value("PADDED_NORMALIZED_LAG", EMAN::PADDED_NORMALIZED_LAG)
    ;

    enum_< EMAN::kernel_shape >("kernel_shape")
        .value("CROSS", EMAN::CROSS)
        .value("BLOCK", EMAN::BLOCK)
        .value("CIRCULAR", EMAN::CIRCULAR)
    ;

    enum_< EMAN::morph_type >("morph_type")
        .value("BINARY", EMAN::BINARY)
        .value("GRAYLEVEL", EMAN::GRAYLEVEL)
    ;

    def("correlation", &EMAN::correlation, return_value_policy< manage_new_object >());
    def("convolution", &EMAN::convolution, return_value_policy< manage_new_object >());
    def("autocorrelation", &EMAN::autocorrelation, return_value_policy< manage_new_object >());
    def("self_correlation", &EMAN::self_correlation, return_value_policy< manage_new_object >());
    def("periodogram", &EMAN::periodogram, return_value_policy< manage_new_object >());
    def("rsconvolution", &EMAN::rsconvolution, return_value_policy< manage_new_object >());
    def("filt_median_", &EMAN::filt_median_, return_value_policy< manage_new_object >());
    def("filt_dilation_", &EMAN::filt_dilation_, return_value_policy< manage_new_object >());
    def("filt_erosion_", &EMAN::filt_erosion_, return_value_policy< manage_new_object >());
}

