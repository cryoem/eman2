#define NO_IMPORT_ARRAY

#include <Python.h>

#include "pylist.h"
#include "emdata.h"

using namespace EMAN;


static python::numeric::array make_numeric_array(float * data, vector<int> dims)
{
	size_t total = 1;
	vector<int>::iterator iter = dims.begin();
	while(iter != dims.end()){
		total *= *iter;
		++iter;
	}    
#if 1
	python::object obj(python::handle<>(PyArray_FromDimsAndData(dims.size(),&dims[0],
																PyArray_FLOAT, (char*)data)));
#endif
#if 0
	python::object obj(python::handle<>(PyArray_FromDims(dims.size(),&dims[0], PyArray_FLOAT)));
	char *arr_data = ((PyArrayObject*) obj.ptr())->data;
	memcpy(arr_data, data, sizeof(float) * total);
#endif
	return python::extract<python::numeric::array>(obj);
}


python::numeric::array Wrapper::em2numpy(EMData *image)
{
	float * data = image->get_data();
	int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();
	
	vector<int> dims;
	
	if (nz > 1) {
		dims.push_back(nz);
	}
	
	if (ny > 1) {
		dims.push_back(ny);
	}
	
	dims.push_back(nx);
	
	return make_numeric_array(data, dims);
}

void Wrapper::numpy2em(python::numeric::array& array, EMData* image)
{
	int nz = python::len(array);
    int ny = python::len(array[0]);
    int nx = python::len(array[0][0]);
	int nxy = nx * ny;
	
    image->set_size(nx, ny, nz);
    
    float* data = image->get_data();
    
    for (int i = 0; i < nz; i++) {
		int i2 = i * nxy;
		for (int j = 0; j < ny; j++) {
			python::numeric::array array2 =
				python::extract<python::numeric::array>(array[i][j]);
			python::list l = python::list(array2);
			int j2 = i2 + j * nx;
			for (int k = 0; k < nx; k++) {
				data[j2 + k] = python::extract<float>(l[k]);
			}
		}
    }
}
