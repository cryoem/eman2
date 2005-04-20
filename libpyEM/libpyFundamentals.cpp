#include <boost/python.hpp>
#include <emdata.h>
#include <fundamentals.h>

using namespace boost::python;

BOOST_PYTHON_MODULE(libpyFundamentals) {
	enum_<EMAN::fp_flag>("fp_flag")
		.value("CIRCULANT", EMAN::CIRCULANT)
		.value("CIRCULANT_NORMALIZED", EMAN::CIRCULANT_NORMALIZED)
		.value("PADDED", EMAN::PADDED)
		.value("PADDED_LAG", EMAN::PADDED_LAG)
		.value("PADDED_NORMALIZED_LAG", EMAN::PADDED_NORMALIZED_LAG)
		;
	def("correlation", &EMAN::correlation, return_value_policy< manage_new_object >() );
	def("convolution", &EMAN::convolution, return_value_policy< manage_new_object >() );
	def("autocorrelation", &EMAN::autocorrelation, return_value_policy< manage_new_object >() );
	def("self_correlation", &EMAN::self_correlation, return_value_policy< manage_new_object >() );
	def("periodogram", &EMAN::periodogram, return_value_policy< manage_new_object >() );
	def("norm_pad_ft", &EMAN::norm_pad_ft, return_value_policy< manage_new_object >() );
}
