#include "pylist.h"
#include "emdata.h"

using namespace EMAN;

python::numeric::array Wrapper::em2numpy(EMData *image)
{
	python::list datalist = python::list();
	
	if (!image) {
		return python::numeric::array(datalist);
	}

	float * data = image->get_data();
    int nx = image->get_xsize();
    int ny = image->get_ysize();
    int nz = image->get_zsize();
	int nxy = nx * ny;
    
    for (int i = 0; i < nz; i++) {
		python::list ll = python::list();
		int i2 = i * nxy;
		
		for (int j = 0; j < ny; j++) {
			python::list l = python::list();
			int j2 = j * nx + i2;
			
			for (int k = 0; k < nx; k++) {
				l.append(data[j2 + k]);
			}
			
			ll.append(l);
		}
		datalist.append(ll);
    }
    
    return python::numeric::array(datalist);
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
			python::numeric::array array2 = python::extract<python::numeric::array>(array[i][j]);
			python::list l = python::list(array2);
			int j2 = i2 + j * nx;
			for (int k = 0; k < nx; k++) {
				data[j2 + k] = python::extract<float>(l[k]);
			}
		}
    }
}
