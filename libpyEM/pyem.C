#include "pyem.h"
#include "pylist.h"
#include "emdata.h"

#include <boost/python/detail/api_placeholder.hpp>

using namespace EMAN;


python::list py_read_images_by_index(string filename, python::list img_indices, bool nodata)
{
    vector<int> v = PyList::list2vector<int>(img_indices);
    std::vector<EMData*> images = EMData::read_images_by_index(filename, v, nodata);
    python::list result = PyList::vector2list<EMData*>(images);
    return result;
}

python::list py_read_images_by_ext(string filename, int img_index_start, int img_index_end, bool nodata, string ext)
{
    std::vector<EMData*> images = EMData::read_images_by_ext(filename, img_index_start, img_index_end, nodata, ext);
    return PyList::vector2list<EMData*>(images);
}

