#ifndef eman__sparx_h__
#define eman__sparx_h__

#include <scitbx/array_family/boost_python/flex_fwd.h> // MUST be the very first include

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace python = boost::python;

namespace EMAN {
	void em2flex(EMData& image, Dict & header_dict, scitbx::af::shared<float> & flex_array);
	void flex2em(const scitbx::af::shared<float> & flex_array, Dict & header_dict, EMData & image);
}


#endif
