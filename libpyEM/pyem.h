#ifndef __pyem__h__
#define __pyem__h__

#include <string>
#include <imageio.h>
#include <transform.h>
#include <ctf.h>

#include <boost/python.hpp>

namespace python = boost::python;
using std::string;
#if 0

python::list py_read_images_by_index(string filename, python::list img_indices, bool nodata = false);
python::list py_read_images_by_ext(string filename, int img_index_start, int img_index_end,
				bool nodata = false, string ext = "");
#endif
#endif
