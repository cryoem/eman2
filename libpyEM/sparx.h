#ifndef eman__sparx_h__
#define eman__sparx_h__

#include <scitbx/array_family/boost_python/flex_fwd.h> // MUST be the very first include

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace python = boost::python;

namespace EMAN {

	/** return a cctbx flex array storing the raw image data */
	scitbx::af::shared<float> em2flex(const EMData& image);
	

}


#endif
