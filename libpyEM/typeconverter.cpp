#define NO_IMPORT_ARRAY

#include <Python.h>

#include "typeconverter.h"
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

	python::object obj(python::handle<>(PyArray_FromDimsAndData(dims.size(),&dims[0],
																PyArray_FLOAT, (char*)data)));


	return python::extract<python::numeric::array>(obj);
}


python::numeric::array EMNumPy::em2numpy(EMData *image)
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

void EMNumPy::numpy2em(python::numeric::array& array, EMData* image)
{
	if (!PyArray_Check(array.ptr())) {
		PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
		return;
	}
	
	PyArrayObject * array_ptr = (PyArrayObject*) array.ptr();
	int ndim = array_ptr->nd;
	int * dims_ptr = array_ptr->dimensions;
	
	int nx = 1;
	int ny = 1;
	int nz = 1;

	if (ndim <= 0 || ndim > 3) {
		LOGERR("%dD Numeric array to EMData is not supported.", ndim);
		return;
	}

	if (ndim == 1) {
		nx = dims_ptr[0];
	}
	else if (ndim == 2) {
		ny = dims_ptr[0];
		nx = dims_ptr[1];
	}
	else if (ndim == 3) {
		nz = dims_ptr[0];
		ny = dims_ptr[1];
		nx = dims_ptr[2];
	}

	char* array_data = ((PyArrayObject*) array.ptr())->data;
	image->set_shared_data(nx, ny, nz, (float*)array_data);
	
}

#if 0
// need further investigation to make it work
PyObject* EMObject_to_python::convert(EMObject const& emobj)
{

	EMObject::ObjectType t = emobj.get_type();
	PyObject * result = 0;
			
	if (t == EMObject::INT) {
		result = PyInt_FromLong((int)emobj);
	}
	else if (t == EMObject::FLOAT) {
		result = PyFloat_FromDouble((float) emobj);
	}
	else if (t == EMObject::DOUBLE) {
		result = PyFloat_FromDouble((double) emobj);
	}
	else if (t == EMObject::STRING) {
		result = PyString_FromString((const char*) emobj);
	}
	else if (t == EMObject::EMDATA) {
		EMData * img = (EMData*) emobj;
		result = python::object(img).ptr();
	}
	else if (t == EMObject::XYDATA) {
		XYData * xyd = (XYData*) xyd;
		result = python::object(xyd).ptr();
	}
			
	if (result) {
		return python::incref(result);
	}
	return 0;
}


#endif
