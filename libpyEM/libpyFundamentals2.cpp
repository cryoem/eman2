
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

