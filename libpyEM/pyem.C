#include "pyem.h"
#include "pylist.h"
#include "emdata.h"

#include <boost/python/detail/api_placeholder.hpp>

using namespace EMAN;

python::list py_read_images_index(string filename, python::list img_indices, int nimg, bool nodata)
{
    int* array = new int[python::len(img_indices)];
    PyList::list2array(img_indices, array);
    
    std::vector<EMData*> images = EMData::read_images_by_index(filename, array, nimg, nodata);
    
    delete [] array;
    array = 0;
    
    return PyList::vector2list(images);
}

python::list py_read_images_ext(string filename, int img_index_start, int img_index_end, bool nodata, string ext)
{
    std::vector<EMData*> images = EMData::read_images_by_ext(filename, img_index_start, img_index_end, nodata, ext);
    return PyList::vector2list(images);
}

