
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <fundamentals.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyFundamentals2)
{
    enum_< EMAN::fp_flag >("fp_flag")
        .value("PADDED_NORMALIZED_LAG", EMAN::PADDED_NORMALIZED_LAG)
        .value("PADDED_NORMALIZED", EMAN::PADDED_NORMALIZED)
        .value("CIRCULANT_NORMALIZED", EMAN::CIRCULANT_NORMALIZED)
        .value("PADDED", EMAN::PADDED)
        .value("CIRCULANT", EMAN::CIRCULANT)
        .value("PADDED_LAG", EMAN::PADDED_LAG)
    ;

    def("correlation", &EMAN::correlation, return_value_policy< manage_new_object >());
    def("convolution", &EMAN::convolution, return_value_policy< manage_new_object >());
    def("autocorrelation", &EMAN::autocorrelation, return_value_policy< manage_new_object >());
    def("self_correlation", &EMAN::self_correlation, return_value_policy< manage_new_object >());
    def("periodogram", &EMAN::periodogram, return_value_policy< manage_new_object >());
    def("norm_pad_ft", &EMAN::norm_pad_ft, return_value_policy< manage_new_object >());
}

