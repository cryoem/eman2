/*
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

#include <Python.h>
#include <numpy/arrayobject.h>

// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <emutil.h>
#include <sparx/lapackblas.h>
#include <imageio.h>
#include <testutil.h>
#include <xydata.h>
#include <emobject.h>
#include <randnum.h>
#include "ctf.h"
#include "geometry.h"
#include "portable_fileio.h"

// Using =======================================================================
using namespace boost::python;

#if PY_MAJOR_VERSION >= 3
	#define IS_PY3K
#endif

// Declarations ================================================================
namespace  {

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_im_diff_overloads_2_3, EMAN::Util::im_diff, 2, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_TwoDTestFunc_overloads_5_7, EMAN::Util::TwoDTestFunc, 5, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_even_angles_overloads_1_5, EMAN::Util::even_angles, 1, 5)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_decimate_overloads_2_4, EMAN::Util::decimate, 2, 4)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_window_overloads_2_7, EMAN::Util::window, 2, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_pad_overloads_2_8, EMAN::Util::pad, 2, 8)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_tf_overloads_2_7, EMAN::Util::tf, 2, 7)

//BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_ctf_img_overloads_5_12, EMAN::Util::ctf_img, 5, 12)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_histogram_overloads_2_5, EMAN::Util::histogram, 2, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_multiref_polar_ali_helical_overloads_10_11, EMAN::Util::multiref_polar_ali_helical, 10, 11)
BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_multiref_polar_ali_helical_local_overloads_11_13, EMAN::Util::multiref_polar_ali_helical_local, 11, 13)
BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_multiref_polar_ali_helical_90_overloads_10_11, EMAN::Util::multiref_polar_ali_helical_90, 10, 11)
BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_multiref_polar_ali_helical_90_local_overloads_11_13, EMAN::Util::multiref_polar_ali_helical_90_local, 11, 13)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_multiref_polar_ali_helicon_local_overloads_11_13, EMAN::Util::multiref_polar_ali_helicon_local, 11, 13)
BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_multiref_polar_ali_helicon_90_local_overloads_11_13, EMAN::Util::multiref_polar_ali_helicon_90_local, 11, 13)


void read_raw_emdata(EMAN::EMData *ths, const char* path, size_t offset, int rw_mode, int image_index, int mode, const EMAN::Region * area=0) { 
	Py_BEGIN_ALLOW_THREADS
	const char *iomode;
	if (rw_mode==1) iomode="r";	// 1 is READ_ONLY
	else iomode="r+";
	FILE *file = fopen(path,iomode);
	if (file==NULL) { printf("ERROR: cannot open %s\n",path); return; }
	portable_fseek(file,offset,SEEK_SET);
	size_t mode_size;
	if (mode==0) mode_size=4;			// single precision float
	else if (mode==1) mode_size=4;	// unsigned 32 bit int
	else if (mode==2) mode_size=4;	// signed 32 bit int 
	else if (mode==3) mode_size=2;	// unsigned 16 bit int
	else if (mode==4) mode_size=2;	// signed 16 bit int 
		
	EMAN::EMUtil::process_region_io(ths->get_data(),file,rw_mode,image_index,mode_size,ths->get_xsize(),ths->get_ysize(),ths->get_zsize(),area);
	fclose(file);
	
	// Now we convert the data we read into floating point
	size_t n=ths->get_xsize()*ths->get_ysize()*ths->get_zsize();
	float *data=ths->get_data();
	if (mode==0) ;
	else if (mode==1) {
		unsigned int *idata=(unsigned int *)ths->get_data();
		for (size_t i=n; i>0; i--) data[i-1]=(float)idata[i-1];		// not as silly as it looks, size_t is unsigned. Not a perfect solution, but it works
	}
	else if (mode==2) {
		int *idata=(int *)ths->get_data();
		for (size_t i=n; i>0; i--) data[i-1]=(float)idata[i-1];
	}
	else if (mode==3) {
		unsigned short *idata=(unsigned short *)ths->get_data();
		for (size_t i=n; i>0; i--) data[i-1]=(float)idata[i-1];
	}
	else if (mode==4) {
		short *idata=(short *)ths->get_data();
		for (size_t i=n; i>0; i--) data[i-1]=(float)idata[i-1];
	}
	else printf("read_raw_emdata: Unknown mode\n");
	
	ths->update();
	
	Py_END_ALLOW_THREADS
}


struct EMAN_Util_sincBlackman_Wrapper: EMAN::Util::sincBlackman
{
    EMAN_Util_sincBlackman_Wrapper(PyObject* py_self_, const EMAN::Util::sincBlackman& p0):
        EMAN::Util::sincBlackman(p0), py_self(py_self_) {}

    EMAN_Util_sincBlackman_Wrapper(PyObject* py_self_, int p0, float p1):
        EMAN::Util::sincBlackman(p0, p1), py_self(py_self_) {}

    EMAN_Util_sincBlackman_Wrapper(PyObject* py_self_, int p0, float p1, int p2):
        EMAN::Util::sincBlackman(p0, p1, p2), py_self(py_self_) {}

    void build_sBtable() {
        call_method< void >(py_self, "build_sBtable");
    }

    void default_build_sBtable() {
        EMAN::Util::sincBlackman::build_sBtable();
    }


    PyObject* py_self;
};

struct EMAN_Util_KaiserBessel_Wrapper: EMAN::Util::KaiserBessel
{
    EMAN_Util_KaiserBessel_Wrapper(PyObject* py_self_, const EMAN::Util::KaiserBessel& p0):
        EMAN::Util::KaiserBessel(p0), py_self(py_self_) {}

    EMAN_Util_KaiserBessel_Wrapper(PyObject* py_self_, float p0, int p1, float p2, float p3, int p4):
        EMAN::Util::KaiserBessel(p0, p1, p2, p3, p4), py_self(py_self_) {}

    EMAN_Util_KaiserBessel_Wrapper(PyObject* py_self_, float p0, int p1, float p2, float p3, int p4, float p5):
        EMAN::Util::KaiserBessel(p0, p1, p2, p3, p4, p5), py_self(py_self_) {}

    EMAN_Util_KaiserBessel_Wrapper(PyObject* py_self_, float p0, int p1, float p2, float p3, int p4, float p5, int p6):
        EMAN::Util::KaiserBessel(p0, p1, p2, p3, p4, p5, p6), py_self(py_self_) {}

    void build_I0table() {
        call_method< void >(py_self, "build_I0table");
    }

    void default_build_I0table() {
        EMAN::Util::KaiserBessel::build_I0table();
    }

    float sinhwin(float p0) const {
        return call_method< float >(py_self, "sinhwin", p0);
    }

    float default_sinhwin(float p0) const {
        return EMAN::Util::KaiserBessel::sinhwin(p0);
    }

    float i0win(float p0) const {
        return call_method< float >(py_self, "i0win", p0);
    }

    float default_i0win(float p0) const {
        return EMAN::Util::KaiserBessel::i0win(p0);
    }

    PyObject* py_self;
};

struct EMAN_Util_FakeKaiserBessel_Wrapper: EMAN::Util::FakeKaiserBessel
{
    EMAN_Util_FakeKaiserBessel_Wrapper(PyObject* py_self_, const EMAN::Util::FakeKaiserBessel& p0):
        EMAN::Util::FakeKaiserBessel(p0), py_self(py_self_) {}

    EMAN_Util_FakeKaiserBessel_Wrapper(PyObject* py_self_, float p0, int p1, float p2, float p3, int p4):
        EMAN::Util::FakeKaiserBessel(p0, p1, p2, p3, p4), py_self(py_self_) {}

    EMAN_Util_FakeKaiserBessel_Wrapper(PyObject* py_self_, float p0, int p1, float p2, float p3, int p4, float p5):
        EMAN::Util::FakeKaiserBessel(p0, p1, p2, p3, p4, p5), py_self(py_self_) {}

    EMAN_Util_FakeKaiserBessel_Wrapper(PyObject* py_self_, float p0, int p1, float p2, float p3, int p4, float p5, int p6):
        EMAN::Util::FakeKaiserBessel(p0, p1, p2, p3, p4, p5, p6), py_self(py_self_) {}

    float sinhwin(float p0) const {
        return call_method< float >(py_self, "sinhwin", p0);
    }

    float default_sinhwin(float p0) const {
        return EMAN::Util::FakeKaiserBessel::sinhwin(p0);
    }

    float i0win(float p0) const {
        return call_method< float >(py_self, "i0win", p0);
    }

    float default_i0win(float p0) const {
        return EMAN::Util::FakeKaiserBessel::i0win(p0);
    }

    void build_I0table() {
        call_method< void >(py_self, "build_I0table");
    }

    void default_build_I0table() {
        EMAN::Util::FakeKaiserBessel::build_I0table();
    }

    PyObject* py_self;
};


BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMUtil_get_imageio_overloads_2_3, EMAN::EMUtil::get_imageio, 2, 3)

#ifdef EM_HDF5
BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMUtil_read_hdf_attribute_2_3, EMAN::EMUtil::read_hdf_attribute, 2, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMUtil_write_hdf_attribute_3_4, EMAN::EMUtil::write_hdf_attribute, 3, 4)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_EMUtil_delete_hdf_attribute_2_3, EMAN::EMUtil::delete_hdf_attribute, 2, 3)
#endif	//EM_HDF5

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_TestUtil_check_image_overloads_1_2, EMAN::TestUtil::check_image, 1, 2)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_TestUtil_make_image_file_overloads_2_6, EMAN::TestUtil::make_image_file, 2, 6)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_TestUtil_verify_image_file_overloads_2_6, EMAN::TestUtil::verify_image_file, 2, 6)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_TestUtil_make_image_file2_overloads_2_6, EMAN::TestUtil::make_image_file2, 2, 6)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_TestUtil_verify_image_file2_overloads_2_6, EMAN::TestUtil::verify_image_file2, 2, 6)

}// namespace

/*
struct EMAN_Util_Wrapper: EMAN::Util
{
	EMAN_Util_Wrapper(PyObject* py_self_, const EMAN::Util& p):
		EMAN::Util(p), py_self(py_self_) {}

	EMAN_Util_Wrapper(PyObject* py_self_):
		EMAN::Util(), py_self(py_self_) {}

	EMAN::Dict get_stats(const std::vector<double,std::allocator<double> >& data) {
		return call_method< EMAN::Dict>(py_self,"get_stats", data);
	}

	PyObject* py_self;
};*/

using boost::python::numeric::array;

float* get_fptr( array& a )
{
/*
	if (!PyArray_Check(a.ptr())) {
		//PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject for get_fptr");
		//return NULL;
		throw std::runtime_error( "Expected a PyArryaObject for get_fptr" );
	} */

	PyArrayObject* aptr = (PyArrayObject*)a.ptr();
        char datatype = aptr->descr->type;
        if( datatype != 'f' )
        {
		//PyErr_SetString(PyExc_ValueError, "expected a float PyArrayObject for get_fptr");
		//return NULL;
		throw std::runtime_error( "Expected a float array for get_fptr" );

	}

	return (float*)(aptr->data);


}


int* get_iptr( array& a )
{
/*
	if (!PyArray_Check(a.ptr())) {
		PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject for get_iptr");
		return NULL;
	}
*/
	PyArrayObject* aptr = (PyArrayObject*)a.ptr();
        char datatype = aptr->descr->type;
        if( datatype != 'i' && datatype != 'l') // for some reaons on a mac it is 'l'
        {
		//cout << "datatype != 'i's: " << datatype << endl;
        	//PyErr_SetString(PyExc_ValueError, "expected an integer PyArrayObject for get_iptr");
		throw std::runtime_error( "Expected a int array for get_fptr" );
	}

	return (int*)(aptr->data);
}



int pysstevd(const string& jobz, int n, array& diag, array& subdiag, array& qmat, int kstep, array& fwork, int lfwrk, array& iwork, int liwrk )
{
    int info;

    float* d = get_fptr( diag );
    float* e = get_fptr( subdiag );
    float* f = get_fptr( qmat );
    float* g = get_fptr( fwork );
    int*  ih = get_iptr( iwork );

    sstevd_( (char*)jobz.c_str(), &n, d, e, f, &kstep, g, &lfwrk, ih, &liwrk, &info );
    return info;
}

float pysnrm2( int n, array& a, int incx )
{
    float* f = get_fptr( a );
    return snrm2_(&n, f, &incx);
}

int pysgemv( const string& trans, int m, int n, float alpha, array& a, int lda, array& x, int incx, float beta, array& y, int incy )
{
    float* fa = get_fptr( a );
    float* fx = get_fptr( x );
    float* fy = get_fptr( y );
    return sgemv_( trans.c_str(), &m, &n, &alpha, fa, &lda, fx, &incx, &beta, fy, &incy );
}

int pysaxpy( int n, float alpha, array& x, int incx, array& y, int incy )
{
    float* fx = get_fptr( x );
    float* fy = get_fptr( y );
    return saxpy_( &n, &alpha, fx, &incx, fy, &incy );
}

float pysdot( int n, array& x, int incx, array& y, int incy )
{
    float* fx = get_fptr( x );
    float* fy = get_fptr( y );
    assert( fx != NULL && fy != NULL );
    return sdot_( &n, fx, &incx, fy, &incy );
}

void readarray( object& f, array& x, int size)
{
#ifdef IS_PY3K
	extern PyTypeObject PyIOBase_Type;
	if(!PyObject_IsInstance(f.ptr(), (PyObject *)&PyIOBase_Type) )
#else
	if( !PyFile_Check(f.ptr()) )
#endif	//IS_PY3K
    {
        std::cout << "Error: expecting a file object" << std::endl;
        return;
    }

#ifdef IS_PY3K
	int fd = PyObject_AsFileDescriptor( f.ptr() );
	FILE*  fh = fdopen(fd, "r");
#else
    FILE*  fh = PyFile_AsFile( f.ptr() );
#endif	//IS_PY3K
    float* fx = get_fptr( x );

    fread( fx, sizeof(float), size, fh );
}


// k_means_cont_table_ is locate to util_sparx.cpp
int pyk_means_cont_table(array& group1, array& group2, array& stb, long int s1, long int s2, int flag) {
    int* pt_group1 = get_iptr(group1);
    int* pt_group2 = get_iptr(group2);
    int* pt_stb  = get_iptr(stb);
    return EMAN::Util::k_means_cont_table_(pt_group1, pt_group2, pt_stb, s1, s2, flag);
}

// bb_enumerateMPI is locate in util_sparx.cpp
vector<int> pybb_enumerateMPI(array& parts, array& classDims, int nParts, int nClasses, int T, int nguesses,int LARGEST_CLASS,int J, int max_branching, float stmult, int
branchfunc, int LIM) {
    int* pt_parts = get_iptr(parts);
    int* pt_classDims = get_iptr(classDims);
    return EMAN::Util::bb_enumerateMPI_(pt_parts, pt_classDims, nParts, nClasses,T,nguesses,LARGEST_CLASS, J, max_branching, stmult, branchfunc, LIM);
}

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyUtils2)
{
	scope* EMAN_Util_scope = new scope(
		 class_< EMAN::Util>("Util", "Util is a collection of utility functions.", init<  >())
		.def(init< const EMAN::Util& >())
		.def("coveig", &EMAN::Util::coveig, args("n", "covmat", "eigval", "eigvec"))
		.def("coveig_for_py", &EMAN::Util::coveig_for_py, args("ncov", "covmatpy"), "same function than Util::coveig but wrapped to use directly in python code")
		.def("WTM", &EMAN::Util::WTM)
		.def("WTF", &EMAN::Util::WTF)
		.def("CANG", &EMAN::Util::CANG)
		.def("BPCQ", &EMAN::Util::BPCQ)
		.def("infomask", &EMAN::Util::infomask)
		.def("helixshiftali", &EMAN::Util::helixshiftali)
		.def("snakeshiftali", &EMAN::Util::snakeshiftali)
		.def("curhelixshiftali", &EMAN::Util::curhelixshiftali)
		.def("bsplineBase", &EMAN::Util::bsplineBase)
		.def("bsplineBasedu", &EMAN::Util::bsplineBasedu)
		.def("convertTocubicbsplineCoeffs", &EMAN::Util::convertTocubicbsplineCoeffs)
		.def("cluster_pairwise", &EMAN::Util::cluster_pairwise, args("d", "K", "T", "F"))
		.def("cluster_equalsize", &EMAN::Util::cluster_equalsize, args("d"))
		.def("vareas", &EMAN::Util::vareas, args("d"))
		.def("cyclicshift", &EMAN::Util::cyclicshift, args("image", "params"), "Performs inplace integer cyclic shift as specified by the 'dx','dy','dz' parameters on a 3d volume.\nImplements the inplace swapping using reversals as descibed in  also:\nhttp://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Intro/Eg01/\n@author  Phani Ivatury\n@date 18-2006\n@see http://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Intro/Eg01/\n\nA[0] A[1] A[2] A[3] A[4] A[5] A[6] A[7] A[8] A[9]\n\n10   20   30   40   50   60   70   80   90   100\n------------\n  m = 3 (shift left three places)\n\n  Reverse the items from 0..m-1 and m..N-1:\n\n 30   20   10   100  90   80   70   60   50   40\n\n Now reverse the entire sequence:\n\n  40   50   60   70   80   90   100  10   20   30\n\n\n  cycl_shift() in libpy/fundementals.py calls this function\n\n\
Usage:\n EMData *im1 = new EMData();\n im1->set_size(70,80,85);\n im1->to_one();\nDict params; params['dx'] = 10;params['dy'] = 10000;params['dz'] = -10;\nUtils::cyclicshift(im1,params);\nim1.peak_search(1,1)")
		.def("im_diff", &EMAN::Util::im_diff, EMAN_Util_im_diff_overloads_2_3(args("V1", "V2", "mask"), "V1 - \nV2 - \nmask - (default=Null)"))
		.def("TwoDTestFunc", &EMAN::Util::TwoDTestFunc, EMAN_Util_TwoDTestFunc_overloads_5_7(args("Size", "p", "q", "a", "b", "flag", "alphaDeg"), "Creates a Two D Test Pattern\n \nSize - must be odd\np - the x frequency\nq - the y frequency\na - the x falloff\nb - the y falloff\nflag - (default=0)\nalphaDeg - the projection angle(default=0)\n \nreturn The 2D test pattern in real space, fourier space,\nor the projection in real or fourier space\nor the FH of the pattern")[return_value_policy< manage_new_object >()])
		.def("splint", &EMAN::Util::splint, args("xa", "ya", "y2a", "n", "xq", "yq", "m"), "Given the arrays xa(ordered, ya of length n, which tabulate a function\nand given the array y2a which is the output of spline and an unordered array xq,\nthis routine returns a cubic-spline interpolated array yq.\n \nxa - \nya - of x is the tabulated function of length n\ny2a - is returned from spline: second derivs\nn - \nxq - is the x values to be splined: has m points.\nyq - are the splined values\nm -")
		.def("even_angles", &EMAN::Util::even_angles, EMAN_Util_even_angles_overloads_1_5(args("delta", "t1", "t2", "p1", "p2"), "Compute a vector containing quasi-evenly spaced Euler angles.\nThe order of angles in the vector is phi, theta, psi.\n \ndelta - Delta theta (spacing in theta).\nt1 - Starting (min) value of theta in degrees(default=0).\nt2 - Ending (max) value of theta in degrees(default=90).\np1 - Starting (min) value of phi in degrees(default=0)\np2 - Ending (max) value of phi in degrees(default = 359.9)\n \nreturn Vector of angles as a flat list of phi_0, theta_0, psi_0, ..., phi_N, theta_N, psi_N."))
		.def("quadri", &EMAN::Util::quadri, args("x", "y", "nx", "ny", "image"), "Quadratic interpolation (2D).\n \nNote:  This routine starts counting from 1, not 0!\n \nThis routine uses six image points for interpolation:\n \n@see M. Abramowitz & I.E. Stegun, Handbook of Mathematical\nFunctions (Dover, New York, 1964), Sec. 25.2.67.\nhttp://www.math.sfu.ca/~cbm/aands/page_882.htm\n \n@see http://www.cl.cam.ac.uk/users/nad/pubs/quad.pdf\n \n@verbatim\n       f3    fc\n       |\n       | x\nf2-----f0----f1\n       |\n       |\n       f4\n@endverbatim\n \nf0 - f4 are image values near the interpolated point X.\nf0 is the interior mesh point nearest x.\n \nCoords:\nf0 = (x0, y0)\nf1 = (xb, y0)\nf2 = (xa, y0)\nf3 = (x0, yb)\nf4 = (x0, ya)\nfc = (xc, yc)\n \nMesh spacings:\nhxa -- x- mesh spacing to the left of f0\nhxb -- x- mesh spacing to the right of f0\n\
hyb -- y- mesh spacing above f0\nhya -- y- mesh spacing below f0\n \nInterpolant:\n  f = f0 + c1*(x-x0) + c2*(x-x0)*(x-x1)\n		 + c3*(y-y0) + c4*(y-y0)*(y-y1)\n		 + c5*(x-x0)*(y-y0)\n \nx - x-coord value\ny - y-coord value\nnx - \nny - \nimage - Image object (pointer)\n \nreturn Interpolated value.")
		.def("get_pixel_conv_new", &EMAN::Util::get_pixel_conv_new, args("nx", "ny", "nz", "delx", "dely", "delz", "data", "kb"), "Here counting is in C style, so coordinates of the pixel delx should be [0-nx-1]\nCommented by Zhengfan Yang on 04/20/07\nThis function is written to replace get_pixel_conv(), which is too slow to use in practice.\nI made the following changes to get_pixel_conv():1. Use the same data passing scheme as quadri() and move the function from emdata_sparx.cpp to util_sparx.cpp\n2. Reduce usage of i0win_tab (from 98 calls to 14 calls in 2D case, from 1029 calls to 21 calls in 3D case!)\n3. Unfold the 'for' loop\n4. Reduce the usage of multiplications through some bracketing (from 98 times to 57 times in 2D case, from 1029 times to 400 times in 3D case)\nThe shortcoming of this routine is that it only works for window size N=7. In case you want to use other window\nsize, say N=5, you can easily modify it by referring my code.")
		.def("bilinear", &EMAN::Util::bilinear, args("xold", "yold", "nsam", "nrow", "xim"), "")
		.def("Polar2D", &EMAN::Util::Polar2D, return_value_policy< manage_new_object >(), args("image", "numr", "mode"), "")
		.def("Polar2Dm", &EMAN::Util::Polar2Dm, return_value_policy< manage_new_object >(), args("image", "cns2", "cnr2", "numr", "cmode"), "")
		.def("alrl_ms", &EMAN::Util::alrl_ms, args("xim", "nsam", "nrow", "cns2", "cnr2", "numr", "circ", "lcirc", "nring", "mode"), "")
		.def("Polar2Dmi", &EMAN::Util::Polar2Dmi, return_value_policy< manage_new_object >(), args("image", "cns2", "cnr2", "numr", "cmode", "kb"), "")
		.def("fftr_q", &EMAN::Util::fftr_q, args("xcmplx", "nv"), "")
		.def("fftr_d", &EMAN::Util::fftr_d, args("xcmplx", "nv"), "")
		.def("fftc_q", &EMAN::Util::fftc_q, args("br", "bi", "ln", "ks"), "")
		.def("fftc_d", &EMAN::Util::fftc_d, args("br", "bi", "ln", "ks"), "")
		.def("Frngs", &EMAN::Util::Frngs, args("circ", "numr"), "This function conducts the Single Precision Fourier Transform for a set of rings")
		.def("Frngs_inv", &EMAN::Util::Frngs_inv, args("circ", "numr"), "This function conducts the Single Precision Inverse Fourier Transform for a set of rings")
		.def("Applyws", &EMAN::Util::Applyws, args("circ", "numr", "wr"), "This is a copy of Applyws routine from alignment.py")
		.def("Crosrng_e", &EMAN::Util::Crosrng_e, args("circ1", "cir2", "numr", "mirrored"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_e is the original one")
		.def("Crosrng_rand_e", &EMAN::Util::Crosrng_rand_e, args("circ1", "cir2", "numr", "mirrored", "previous_max", "an", "psi_pos"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_e is the original one")
		.def("Crosrng_ew", &EMAN::Util::Crosrng_ew, args("circ1", "cir2", "numr", "w", "mirrored"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_ew is the one that you could apply weights to different rings")
		.def("Crosrng_ms", &EMAN::Util::Crosrng_ms, args("circ1", "circ2", "numr"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_ms assumes the user already applied weights to circ1, it also returns both\nstraight and mirrored positions simultaneously.")
		.def("Crosrng_ms_delta", &EMAN::Util::Crosrng_ms_delta, args("circ1", "circ2", "numr", "delta_start", "delta"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_ms assumes the user already applied weights to circ1, it also returns both\nstraight and mirrored positions simultaneously.")
		.def("Crosrng_sm_psi", &EMAN::Util::Crosrng_sm_psi, args("circ1", "circ2", "numr", "psi", "flag", "psi_max"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\n \nchecks both straight & mirrored position\ninput - fourier transforms of rings!!\ncirc1 already multiplied by weights!")
		.def("Crosrng_ns", &EMAN::Util::Crosrng_ns, args("circ1", "circ2", "numr"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates")
		.def("Crosrng_msg", &EMAN::Util::Crosrng_msg, return_value_policy< manage_new_object >(), args("circ1", "circ2", "numr"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_msg differs from the previous ones in that it returns the cross-correlation\nfunction entirely instead of the peak value and position, thus makes it\npossible to use the gridding method to determine the peak position\n \nchecks both straight & mirrored positions\ninput - fourier transforms of rings!!\ncirc1 already multiplied by weights!\nreturns EM object with 1D ccf")
		.def("Crosrng_msg_s", &EMAN::Util::Crosrng_msg_s, return_value_policy< manage_new_object >(), args("circ1", "circ2", "numr"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_msg_s is same as Crosrng_msg except that it only checks straight position\n \nThis program is half of the Crosrng_msg. It only checks straight position.\ninput - fourier transforms of rings!!\ncirc1 already multiplied by weights!\nreturns EM object with 1D ccf")
		.def("Crosrng_msg_m", &EMAN::Util::Crosrng_msg_m, return_value_policy< manage_new_object >(), args("circ1", "circ2", "numr"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates\nCrosrng_msg_m is same as Crosrng_msg except that it only checks mirrored position\n \nThis program is half of the Crosrng_msg. It only checks mirrored position.\ninput - fourier transforms of rings!!\ncirc1 already multiplied by weights!\nreturns EM object with 1D ccf")
		.def("Crosrng_msg_vec_p",&EMAN::Util::Crosrng_msg_vec_p, args("circ1", "circ2", "numr"), "A little notes about different Crosrng:\n \nBasically, they all do cross-correlation function to two images in polar coordinates")
		.def("update_fav", &EMAN::Util::update_fav, args("ave", "dat", "tot", "mirror", "numr"), "")
		.def("sub_fav", &EMAN::Util::sub_fav, args("ave", "dat", "tot", "mirror", "numr"), "")
		.def("ener", &EMAN::Util::ener, args("ave", "numr"))
		.def("ener_tot", &EMAN::Util::ener_tot, args("data", "numr", "tot"), "")
		.def("min_dist_real", &EMAN::Util::min_dist_real, args("image", "data"), "k-means helper")
		.def("min_dist_four", &EMAN::Util::min_dist_four, args("image", "data"), "k-means helper")
		.def("cml_weights", &EMAN::Util::cml_weights, args("cml"), "new code common-lines\nhelper function for the weights calculation by Voronoi to Cml")
		.def("cml_init_rot", &EMAN::Util::cml_init_rot, args("Ori"), "new code common-lines\n2009-03-25 15:35:05 JB. This function prepare rotation matrix for common-lines")
		.def("cml_update_rot", &EMAN::Util::cml_update_rot, args("Rot", "iprj", "nph", "th", "nps"), "new code common-lines")
		.def("cml_line_insino", &EMAN::Util::cml_line_insino, args("Rot", "i_prj", "n_prj"), "new code common-lines\n2009-03-25 15:35:53 JB. This function calculates common-lines between sinogram")
		.def("cml_line_insino_all", &EMAN::Util::cml_line_insino_all, args("Rot", "seq", "n_prj", "n_lines"), "new code common-lines\n2009-03-30 15:35:07 JB. This function calculates all common-lines between sinogram")
		.def("cml_line_in3d", &EMAN::Util::cml_line_in3d, args("Ori", "seq", "nprj", "nlines"), "new code common-lines\n2009-03-26 10:46:14 JB. This function calculate all common-lines in space for Voronoi")
		.def("cml_spin_psi", &EMAN::Util::cml_spin_psi, args("data", "com", "weights", "iprj", "iw", "n_psi", "d_psi", "n_prj"), "new code common-lines\n2009-03-26 11:37:53 JB. This function spin all angle psi and evaluate the partial discrepancy belong common-lines")
		.def("cml_spin_psi_now", &EMAN::Util::cml_spin_psi_now, args("data", "com", "iprj", "iw", "n_psi", "d_psi", "n_prj"), "new code common-lines\n2009-03-26 11:37:53 JB. This function spin all angle psi and evaluate the partial discrepancy belong common-lines")
		.def("cml_disc", &EMAN::Util::cml_disc, args("data", "com", "seq", "weights", "n_lines"), "new code common-lines\n2009-03-30 15:44:05 JB. Compute the discrepancy belong all common-lines")
		.def("cml_prepare_line", &EMAN::Util::cml_prepare_line, args("sino", "line", "ilf", "ihf", "pos_line", "nblines"), "new code common-lines\nThis function prepare the line from sinogram by cutting off some frequencies,\nand creating the mirror part (complexe conjugate of the first part). Then\nboth lines (mirror and without) are drop to the sinogram.\nline is in Fourrier space, ilf low frequency, ihf high frequency, nblines\nnumber of lines of the half sinogram (the non miror part), sino the sinogram,\npos_line the position of the line in the sino.")
		.def("set_line", &EMAN::Util::set_line, args("img", "posline", "line", "offset", "length"), "new code common-lines\nThis function drop a line (line) to an 2D image (img).\nThe position of the line to the image is defined by (postline).\nThe part of the line paste is defined by (offset), the begin position\nand (length) the size.")
		.def("decimate", &EMAN::Util::decimate, EMAN_Util_decimate_overloads_2_4(args("img", "x_step", "y_step", "z_step"), "Decimates the image with respect to the image center.\n(i.e) the center of the original image is kept the same\nand then the initial start pixel is calculated with respect to the\ncenter of the image.\nworks for all 3 dimensions\n \nimg - image\nx_step - x-pixel\ny_step - y-pixel(default=1)\nz_step - z-pixel(default=1)")[ return_value_policy< manage_new_object >() ])
		.def("window", &EMAN::Util::window, EMAN_Util_window_overloads_2_7(args("img", "new_nx", "new_ny", "new_nz", "x_offset", "y_offset", "z_offset"), "img - \nnew_nx - \nnew_ny - (defalt=1)\nnew_nz - (defalt=1)\nx_offset - (default=0)\ny_offset - (default=0)\nz_offset - (default=0)")[ return_value_policy< manage_new_object >() ])
		.def("pad", &EMAN::Util::pad, EMAN_Util_pad_overloads_2_8(args("img", "new_nx", "new_ny", "new_nz", "x_offset", "y_offset", "z_offset", "params"), "img - \nnew_nx - \nnew_ny - (defalt=1)\nnew_nz - (defalt=1)\nx_offset - (default=0)\ny_offset - (default=0)\nz_offset - (default=0)\nparams - (default='average')")[ return_value_policy< manage_new_object >() ])
		.def("histc", &EMAN::Util::histc, args("ref", "img", "mask"), "")
		.def("hist_comp_freq", &EMAN::Util::hist_comp_freq, args("PA", "PB", "size_img", "hist_len", "img", "ref_freq_hist", "mask", "ref_h_diff", "ref_h_min"), "")
		.def("tf", &EMAN::Util::tf, EMAN_Util_tf_overloads_2_7(args("dzz", "ak", "voltage", "cs", "wgh", "b_factor", "sign"), "The unit in the ctf function: dz: Angstrom, cs: CM  Ps: Angstrom, Voltage: Kv,dza: Angstrom, azz: degree wgh: None unit. b_factor: Angstrom^2\nThe CTF function takes form of   *sin(-quadpi*(dz*lambda*ak^2-cs*lambda^3*ak^4/2.)-wgh)*exp(-b_factor*ak^2)*sign\nsign can be set as +1 or -1 . The unit of frequency ak is 1/Angstrom\nAttention: Envelope function in power spectrum has a form of exp(-b_factor*ak^2\ndzz - \nak - \nvoltage - (default=300.0)\ncs - (default=2.0)\nwgh - (1.0)\nb_factor - (default=0.0)\nsign - (default=-1.0)"))
		.def("compress_image_mask", &EMAN::Util::compress_image_mask, return_value_policy< manage_new_object >(), args("img", "mask"), "")
		.def("reconstitute_image_mask", &EMAN::Util::reconstitute_image_mask, return_value_policy< manage_new_object >(), args("img", "mask"), "Recreates a n-d image using its compressed 1-D form and the mask")
		.def("pw_extract", &EMAN::Util::pw_extract, args("pw", "n", "iswi", "ps"), "")
		.def("eval", &EMAN::Util::eval, args("images", "img", "S", "N", "K", "size"), "")
		.def("vrdg", &EMAN::Util::vrdg, args("ph", "th"), "")
		.def("cmp1", &EMAN::Util::cmp1, args("tmp1", "tmp2"), "")
		.def("cmp2", &EMAN::Util::cmp2, args("tmp1", "tmp2"), "")
		.def("areav_", &EMAN::Util::areav_, args("k", "n", "x", "y", "z__", "list", "lptr", "lend", "ier"), "STRIDPACK USED COMMANDS FOR VORONOI")
		.def("madn_scalar", &EMAN::Util::madn_scalar, return_value_policy< manage_new_object >(), args("img", "img1", "scalar"), "out = img + scalar * img1 ")
		.def("mult_scalar", &EMAN::Util::mult_scalar, return_value_policy< manage_new_object >(), args("img", "scalar"), "out = scalar * img")
		.def("addn_img", &EMAN::Util::addn_img, return_value_policy< manage_new_object >(), args("img", "img1"), "out = img + img1")
		.def("subn_img", &EMAN::Util::subn_img, return_value_policy< manage_new_object >(), args("img", "img1"), "out = img - img1")
		.def("muln_img", &EMAN::Util::muln_img, return_value_policy< manage_new_object >(), args("img", "img1"), "out = img * img1")
		.def("divn_img", &EMAN::Util::divn_img, return_value_policy< manage_new_object >(), args("img", "img1"), "out = img / img1")
		.def("squaren_img", &EMAN::Util::squaren_img, return_value_policy< manage_new_object >(), args("img"), "out = |img|^2")
		.def("divn_filter", &EMAN::Util::divn_filter, return_value_policy< manage_new_object >(), args("img", "img1"), "img /= Re(img1) with zero check")
		.def("mad_scalar", &EMAN::Util::mad_scalar, args("img", "img1", "scalar"), "img += scalar * img1")
		.def("mul_scalar", &EMAN::Util::mul_scalar, args("img", "scalar"), "img *= scalar")
		.def("add_img", &EMAN::Util::add_img, args("img", "img1"), "img += img1")
		.def("add_img2", &EMAN::Util::add_img2, args("img", "img1"), "img += img1**2")
		.def("sub_img", &EMAN::Util::sub_img, args("img", "img1"), "img -= img1")
		.def("mul_img", &EMAN::Util::mul_img, args("img", "img1"), "img *= img1")
		.def("div_img", &EMAN::Util::div_img, args("img", "img1"), "img /= img1")
		.def("square_img", &EMAN::Util::square_img, args("img"), "img = |img|^2")
		.def("div_filter", &EMAN::Util::div_filter, args("img", "img1"), "img /= Re(img1) with zero check")
		.def("set_freq", &EMAN::Util::set_freq, args("freqvol", "freqvol"), "utility for sxlocres")
		.def("ctf_img", &EMAN::Util::ctf_img, return_value_policy< manage_new_object >(), args("nx", "ny", "nz", "dz", "ps", "voltage", "cs", "wgh", "b_factor", "dza", "azz", "sign"), "")
		.def("ctf2_rimg", &EMAN::Util::ctf2_rimg, return_value_policy< manage_new_object >(), args("nx", "ny", "nz", "dz", "ps", "voltage", "cs", "wgh", "b_factor", "dza", "azz", "sign"), "")
		.def("ctf_rimg", &EMAN::Util::ctf_rimg, return_value_policy< manage_new_object >(), args("nx", "ny", "nz", "dz", "ps", "voltage", "cs", "wgh", "b_factor", "dza", "azz", "sign"), "")
		.def("pack_complex_to_real", &EMAN::Util::pack_complex_to_real, return_value_policy< manage_new_object >(), args("img"), "pack absolute values of complex image into  real image with addition of Friedel part ")
		.def("histogram", &EMAN::Util::histogram, EMAN_Util_histogram_overloads_2_5(args("image", "mask", "nbins", "hmin", "hmax"), "image - \nmask - \nnbins - (default = 128)\nhmin - (default = 0.0)\nhmax - (default = 0.0)"))
		.def("multiref_polar_ali_2d", &EMAN::Util::multiref_polar_ali_2d, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny"), "formerly known as apmq\nDetermine shift and rotation between image and many referenceimages (crefim, weights have to be applied) quadratic\ninterpolation")
		.def("multiref_polar_ali_2d_delta", &EMAN::Util::multiref_polar_ali_2d_delta, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny"), "formerly known as apmq\nDetermine shift and rotation between image and many referenceimages (crefim, weights have to be applied) quadratic\ninterpolation")
		.def("multiref_polar_ali_2d_nom", &EMAN::Util::multiref_polar_ali_2d_nom, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny"), "formerly known as apnq DO NOT CONSIDER MIRROR\nDetermine shift and rotation between image and many reference\nimages (crefim, weights have to be applied) quadratic\ninterpolation")
		.def("multiref_polar_ali_2d_local", &EMAN::Util::multiref_polar_ali_2d_local, args("image", "crefim", "xrng", "yrng", "step", "ant", "mode", "numr", "cnx", "cny"), "formerly known as apmq\nDetermine shift and rotation between image and many reference\nimages (crefim, weights have to be applied) quadratic\ninterpolation")
		.def("shc", &EMAN::Util::shc, args("image", "crefim", "xrng", "yrng", "step", "ant", "mode", "numr", "cnx", "cny"), "")
		.def("shc_multipeaks", &EMAN::Util::shc_multipeaks, args("image", "crefim", "xrng", "yrng", "step", "ant", "mode", "numr", "cnx", "cny", "max_peaks_count"), "")
		.def("multiref_polar_ali_2d_local_psi", &EMAN::Util::multiref_polar_ali_2d_local_psi, args("image", "crefim", "xrng", "yrng", "step", "ant", "psi_max", "mode", "numr", "cnx", "cny"), "formerly known as apmq\nDetermine shift and rotation between image and many reference\nimages (crefim, weights have to be applied) quadratic\ninterpolation")
		.def("multiref_polar_ali_helical", &EMAN::Util::multiref_polar_ali_helical,EMAN_Util_multiref_polar_ali_helical_overloads_10_11(args("image", "crefim", "xrng", "yrng", "step", "psi_max", "mode", "numr", "cnx", "cny","ynumber"), "formerly known as apmq\nDetermine shift and rotation between image and many reference\nimages (crefim, weights have to be applied) quadratic\ninterpolation\nSearch for peaks only within +/-psi_max from 0 and 180 (helical)"))
		.def("multiref_polar_ali_helical_local", &EMAN::Util::multiref_polar_ali_helical_local,EMAN_Util_multiref_polar_ali_helical_local_overloads_11_13(args("image", "crefim", "xrng", "yrng", "step", "ant", "psi_max", "mode", "numr", "cnx", "cny","ynumber","yrnglocal"), "formerly known as apmq\nDetermine shift and rotation between image and many reference\nimages (crefim, weights have to be applied) quadratic\ninterpolation\nSearch for peaks only within +/-psi_max from 0 and 180 (helical)"))
		.def("multiref_polar_ali_helical_90", &EMAN::Util::multiref_polar_ali_helical_90,EMAN_Util_multiref_polar_ali_helical_90_overloads_10_11(args("image", "crefim", "xrng", "yrng", "step", "psi_max", "mode", "numr", "cnx", "cny","ynumber"), "formerly known as apmq\nDetermine shift and rotation between image and many reference with theta equals to 90 degree\nimages (crefim, weights have to be applied) quadratic\ninterpolation\nSearch for peaks only within +/-psi_max from 0 and 180 (helical)"))
		.def("multiref_polar_ali_helical_90_local", &EMAN::Util::multiref_polar_ali_helical_90_local,EMAN_Util_multiref_polar_ali_helical_90_local_overloads_11_13(args("image", "crefim", "xrng", "yrng", "step", "ant", "psi_max", "mode", "numr", "cnx", "cny","ynumber","yrnglocal"), "formerly known as apmq\nDetermine shift and rotation between image and many reference with theta equals to 90 degree\nimages (crefim, weights have to be applied) quadratic\ninterpolation\nSearch for peaks only within +/-psi_max from 0 and 180 (helical)"))
		.def("multiref_polar_ali_helicon_local", &EMAN::Util::multiref_polar_ali_helicon_local,EMAN_Util_multiref_polar_ali_helicon_local_overloads_11_13(args("image", "crefim", "xrng", "yrng", "step", "ant", "psi_max", "mode", "numr", "cnx", "cny","ynumber","yrnglocal"), "formerly known as apmq\nDetermine shift and rotation between image and many reference\nimages (crefim, weights have to be applied) quadratic\ninterpolation\nSearch for peaks only within +/-psi_max from 0 and 180 (helicon)"))
		.def("multiref_polar_ali_helicon_90_local", &EMAN::Util::multiref_polar_ali_helicon_90_local,EMAN_Util_multiref_polar_ali_helicon_90_local_overloads_11_13(args("image", "crefim", "xrng", "yrng", "step", "ant", "psi_max", "mode", "numr", "cnx", "cny","ynumber","yrnglocal"), "formerly known as apmq\nDetermine shift and rotation between image and many reference with theta equals to 90 degree\nimages (crefim, weights have to be applied) quadratic\ninterpolation\nSearch for peaks only within +/-psi_max from 0 and 180 (helicon)"))
		.def("multiref_peaks_ali2d", &EMAN::Util::multiref_peaks_ali2d, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny", "peaks", "peakm"), "Determine shift and rotation between image and one reference\nimage (crefim, weights have to be applied) using quadratic\ninterpolation, return a list of peaks  PAP  07/21/08\n\nccf1d keeps 1d ccfs stored as (maxrin, -kx-1:kx+1, -ky-1:ky+1)\nmargin is needed for peak search and both arrays are initialized with -1.0e20")
		.def("multiref_peaks_compress_ali2d", &EMAN::Util::multiref_peaks_compress_ali2d, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny", "peaks", "peakm", "peaks_compress", "peakm_compress"), "Determine shift and rotation between image and one reference\nimage (crefim, weights have to be applied) using quadratic\ninterpolation, return a list of peaks  PAP  07/21/08\n\nccf1d keeps 1d ccfs stored as (maxrin, -kx-1:kx+1, -ky-1:ky+1)\nmargin is needed for peak search and both arrays are initialized with -1.0e20")
		.def("ali2d_ccf_list", &EMAN::Util::ali2d_ccf_list, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny", "T"), "Determine shift and rotation between image and one reference\nimage (crefim, weights have to be applied) using quadratic\ninterpolation")
		.def("ali2d_ccf_list_snake", &EMAN::Util::ali2d_ccf_list_snake, args("image", "crefim", "xrng", "yrng", "step", "mode", "numr", "cnx", "cny", "T"), "Determine shift and rotation between image and one reference\nimage (crefim, weights have to be applied) using quadratic\ninterpolation")
		.def("compress_image_mask", &EMAN::Util::compress_image_mask, return_value_policy< manage_new_object >())
		.def("reconstitute_image_mask", &EMAN::Util::reconstitute_image_mask, return_value_policy< manage_new_object >())
		.def("pw_extract", &EMAN::Util::pw_extract)
		.def("eval", &EMAN::Util::eval)
		.def("local_inner_product", &EMAN::Util::local_inner_product)
		.def("vrdg", &EMAN::Util::vrdg)
		.def("cmp1", &EMAN::Util::cmp1)
		.def("cmp2", &EMAN::Util::cmp2)
		.def("areav_", &EMAN::Util::areav_)
		.def("madn_scalar", &EMAN::Util::madn_scalar, return_value_policy< manage_new_object >())
		.def("mult_scalar", &EMAN::Util::mult_scalar, return_value_policy< manage_new_object >())
		.def("addn_img", &EMAN::Util::addn_img, return_value_policy< manage_new_object >())
		.def("subn_img", &EMAN::Util::subn_img, return_value_policy< manage_new_object >())
		.def("muln_img", &EMAN::Util::muln_img, return_value_policy< manage_new_object >())
		.def("divn_img", &EMAN::Util::divn_img, return_value_policy< manage_new_object >())
		.def("squaren_img", &EMAN::Util::squaren_img, return_value_policy< manage_new_object >())
		.def("divn_filter", &EMAN::Util::divn_filter, return_value_policy< manage_new_object >())
		.def("mad_scalar", &EMAN::Util::mad_scalar)
		.def("mul_scalar", &EMAN::Util::mul_scalar)
		.def("add_img", &EMAN::Util::add_img)
		.def("add_img_abs", &EMAN::Util::add_img_abs)
		.def("add_img2", &EMAN::Util::add_img2)
		.def("sub_img", &EMAN::Util::sub_img)
		.def("mul_img", &EMAN::Util::mul_img)
		.def("div_img", &EMAN::Util::div_img)
		.def("square_img", &EMAN::Util::square_img)
		.def("div_filter", &EMAN::Util::div_filter)
		.def("set_freq", &EMAN::Util::set_freq)
		.def("ctf_img", &EMAN::Util::ctf_img, return_value_policy< manage_new_object >())
		.def("ctf2_rimg", &EMAN::Util::ctf2_rimg, return_value_policy< manage_new_object >())
		.def("ctf_rimg", &EMAN::Util::ctf_rimg, return_value_policy< manage_new_object >())
		.def("pack_complex_to_real", &EMAN::Util::pack_complex_to_real, return_value_policy< manage_new_object >())
		.def("histogram", &EMAN::Util::histogram)
		.def("multiref_polar_ali_2d", &EMAN::Util::multiref_polar_ali_2d)
		.def("multiref_polar_ali_2d_peaklist", &EMAN::Util::multiref_polar_ali_2d_peaklist)
		.def("multiref_polar_ali_2d_peaklist_local", &EMAN::Util::multiref_polar_ali_2d_peaklist_local)
		.def("assign_projangles", &EMAN::Util::assign_projangles)
		.def("nearestk_to_refdir", &EMAN::Util::nearestk_to_refdir)
		.def("nearest_ang", &EMAN::Util::nearest_ang)
		.def("assign_groups", &EMAN::Util::assign_groups)
		.def("group_proj_by_phitheta", &EMAN::Util::group_proj_by_phitheta)
		.def("multiref_polar_ali_2d_nom", &EMAN::Util::multiref_polar_ali_2d_nom)
		.def("multiref_polar_ali_2d_local", &EMAN::Util::multiref_polar_ali_2d_local)
		.def("multiref_polar_ali_2d_local_psi", &EMAN::Util::multiref_polar_ali_2d_local_psi)
		.def("multiref_polar_ali_helical", &EMAN::Util::multiref_polar_ali_helical)
		.def("multiref_polar_ali_helical_local", &EMAN::Util::multiref_polar_ali_helical_local)
		.def("multiref_polar_ali_helical_90", &EMAN::Util::multiref_polar_ali_helical_90)
		.def("multiref_polar_ali_helical_90_local", &EMAN::Util::multiref_polar_ali_helical_90_local)
		.def("multiref_polar_ali_helicon_local", &EMAN::Util::multiref_polar_ali_helicon_local)
		.def("multiref_polar_ali_helicon_90_local", &EMAN::Util::multiref_polar_ali_helicon_90_local)
		.def("multiref_peaks_ali2d", &EMAN::Util::multiref_peaks_ali2d)
		.def("multiref_peaks_compress_ali2d", &EMAN::Util::multiref_peaks_compress_ali2d)
		.def("ali2d_ccf_list", &EMAN::Util::ali2d_ccf_list)
		.def("ali2d_ccf_list_snake", &EMAN::Util::ali2d_ccf_list_snake)
		//.def("multiref_peaks_ali", &EMAN::Util::multiref_peaks_ali)
		.def("local_inner_product", &EMAN::Util::local_inner_product, args("image1", "image2", "lx", "ly", "lz", "w"), "")
		.def("move_points", &EMAN::Util::move_points, return_value_policy< manage_new_object >(), args("image", "qprob", "ri", "ro"), "")
		.def("is_file_exist", &EMAN::Util::is_file_exist,args("filename"), "check whether a file exists or not\n \nreturn True if the file exists; False if not.")
		.def("svdcmp", &EMAN::Util::svdcmp, args("data", "nvec"), "Perform singular value decomposition on a set of images\n \ndata - A List of data objects to be decomposed\nnvec - Number of basis vectors to return, 0 returns full decomposition\n \nreturn A list of images representing basis vectors in the SVD generated subspace")
		.def("sstrncmp", &EMAN::Util::sstrncmp, args("s1", "s2"), "Safe string compare. It compares 's2' with the first N\ncharacters of 's1', where N is the length of 's2'.\n \ns1 - String 1. Its first strlen(s2) characters will be used to do the comparison.\ns2 - String 2. Its whole string will be used to do the comparison.\n \nreturn True if the comparison is equal. False if not equal.")
		.def("int2str", &EMAN::Util::int2str, args("n"), "Get a string format of an integer, e.g. 123 will be '123'.\n \nn - The input integer.\n \nreturn The string format of the given integer.")
		.def("change_filename_ext", &EMAN::Util::change_filename_ext, args("old_filename", "new_ext"), "Change a file's extension and return the new filename.\nIf the given new extension is empty, the old filename is\nnot changed. If the old filename has no extension, add the\nnew extension to it.")
		.def("remove_filename_ext", &EMAN::Util::remove_filename_ext,args("filename"), "Remove a filename's extension and return the new filename.\n \nfilename The old filename whose extension is going to be removed.\n \nreturn The new filename without extension.")
		.def("save_data",(void (*)(float,float,const vector<float> &,const string &))&EMAN::Util::save_data, args("x0", "dx", "y_array", "filename"), "Save x, y data into a file. Each line of the file have the\nformat \"x1TABy1\", where x1 = x0 + dx*i; y1 = y_array[i].\n \nx0 - The starting point of x.\ndx - delta x. The increase step of x data.\ny_array - The y data array.\nfilename - The output filename.")
//		.def("save_data",(void (*)(const vector<float> &,const vector<float> &,const string &))&EMAN::Util::save_data, args("x0", "dx", "y_array", "filename"), "Save x, y data into a file. Each line of the file have the\nformat \"x1TABy1\"\nfilename - The output filename.")
		.def("get_filename_ext", &EMAN::Util::get_filename_ext, args("filename"), "Get a filename's extension.\n \nfilename - A given filename.\n \nreturn The filename's extension, or empty string if the file has no extension.")
		.def("sbasename", &EMAN::Util::sbasename, args("filename"), "Get a filename's basename. For example, the basename of\n'hello.c' is still 'hello.c'; The basename of\n'/tmp/abc/hello.c' is 'hello.c'.\n \nfilename - The given filename, full path or relative path.\n \nreturn The basename of the filename.")
		.def("set_randnum_seed", &EMAN::Util::set_randnum_seed, args("seed"), "Set the seed for Randnum class\n \nseed - the seed for current random number generator")
		.def("get_randnum_seed", &EMAN::Util::get_randnum_seed, "Get the seed for Randnum class\n \nreturn the seed for current random number generator")
		.def("get_irand", (int (*)(int, int))&EMAN::Util::get_irand, args("low", "high"), "Get an integer random number between low and high, [low, high]\n \nlow - The lower bound of the random number.\nhigh - The upper bound of the random number.\n \nreturn The random number between low and high.")
		.def("get_frand", (float (*)(int, int))&EMAN::Util::get_frand, args("low", "high"), "Get a float random number between low and high, [low, high)\n \nlow The lower bound of the random number.\nhigh The upper bound of the random number.\n \nreturn The random number between low and high.")
		.def("get_frand", (float (*)(float, float))&EMAN::Util::get_frand, args("low", "high"), "Get a float random number between low and high, [low, high)\n \nlow The lower bound of the random number.\nhigh The upper bound of the random number.\n \nreturn The random number between low and high.")
		.def("get_frand", (float (*)(double, double))&EMAN::Util::get_frand, args("low", "high"), "Get a float random number between low and high, [low, high)\n \nlow The lower bound of the random number.\nhigh The upper bound of the random number.\n \nreturn The random number between low and high.")
		.def("get_gauss_rand", &EMAN::Util::get_gauss_rand, args("mean", "sigma"), "Get a Gaussian random number.\n \nmean - The gaussian mean\nsigma - The gaussian sigma\n \nreturn the gaussian random number.")
		.def("round", (int (*)(float))&EMAN::Util::round, args("x"), "Get ceiling round of a float number x.\n \nx - Given float number.\n \nreturn Ceiling round of x.")
		.def("round", (int (*)(double))&EMAN::Util::round, args("x"), "Get ceiling round of a double number x.\n \nx - Given float number.\n \nreturn Ceiling round of x.")
		.def("bilinear_interpolate", &EMAN::Util::bilinear_interpolate, args("p1", "p2", "p3", "p4", "t", "u"), "Calculate bilinear interpolation.\n \np1 - The first number. corresponding to (x0,y0).\np2 - The second number. corresponding to (x1,y0).\np3 - The third number. corresponding to (x1,y1).\np4 - The fourth number. corresponding to (x0,y1).\nt - \nu - \n \nreturn The bilinear interpolation value.")
		.def("trilinear_interpolate", &EMAN::Util::trilinear_interpolate, args("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "t", "u", "v"), "Calculate trilinear interpolation.\n \np1 - The first number. corresponding to (x0,y0,z0).\np2 - The second number. corresponding to (x1,y0,z0).\np3 - The third number. corresponding to (x0,y1, z0).\np4 - The fourth number. corresponding to (x1,y1,z0).\np5 - The fifth number. corresponding to (x0,y0,z1).\np6 - The sixth number. corresponding to (x1,y0,z1).\np7 - The seventh number. corresponding to (x0,y1,z1).\np8 The eighth number. corresponding to (x1,y1,z1).\nt - \nu - \nv - \n \nreturn The trilinear interpolation value.")
		.def("calc_best_fft_size", &EMAN::Util::calc_best_fft_size, args("low"), "Search the best FFT size with good primes. It supports\nFFT size up to 4096 now.\n \nlow - low size the search starts with.\n \nreturn The best FFT size.")
		.def("nonconvex", &EMAN::Util::nonconvex, args("curve","first"),"Takes a 1-D curve (list of floats) and makes it nonconvex, by iteratively\nconstraining points to be less than the mean of the two surrounding points\n")
		.def("windowdot", &EMAN::Util::windowdot, args("curveA","curveB","windowsize","normalize"),"Computes a windowed dot product between curve A and curve B. Curve B is normalized, curve A can be normalized or not. It thus gives either an absolute or relative indicator of similarity between the two curves.\n")
		.def("calc_bilinear_least_square", &EMAN::Util::calc_bilinear_least_square, args("points"), "calculate bilinear least-square fit, z = a + b x + c y\nTakes a set of x,y,z vectors and produces an a,b,c vector\ndoes not accept error bars on z or return goodness of fit\n \npoints  a vector<float> of x,y,z values in (x1,y1,z1,x2,y2,z2...) sequence to fit a plane to\n \nreturn result as a Vec3f(a,b,c)")
		.def("square", (int (*)(int))&EMAN::Util::square, args("n"), "Calculate an int number's square.\n \nn - Given int number.")
		.def("square", (float (*)(float))&EMAN::Util::square, args("x"), "Calculate a float number's square.\nx - n Given float number.")
		.def("square", (float (*)(double))&EMAN::Util::square, args("x"), "Calculate a double number's square.\nx - n Given double number.")
		.def("square_sum", &EMAN::Util::square_sum, args("x", "y"), "Calcuate (x*x + y*y).\n \nx - The first number.\ny - The second number.\n \nreturn (x*x + y*y).")
		.def("hypot3", (float (*)(int, int, int))&EMAN::Util::hypot3, args("x", "y", "z"), "Euclidean distance function in 3D: f(x,y,z) = sqrt(x*x + y*y + z*z)\n \nx - The first number\ny - The second number\nz - The third number\n \nreturn sqrt(x*x + y*y + z*z)")
		.def("hypot3", (float (*)(float, float, float))&EMAN::Util::hypot3, args("x", "y", "z"), "Euclidean distance function in 3D: f(x,y,z) = sqrt(x*x + y*y + z*z)\n \nx - The first number\ny - The second number\nz - The third number\n \nreturn sqrt(x*x + y*y + z*z)")
		.def("hypot3", (double (*)(double, double, double))&EMAN::Util::hypot3, args("x", "y", "z"), "Euclidean distance function in 3D: f(x,y,z) = sqrt(x*x + y*y + z*z)\n \nx - The first number\ny - The second number\nz - The third number\n \nreturn sqrt(x*x + y*y + z*z)")
		.def("fast_floor", &EMAN::Util::fast_floor, args("x"), "A fast way to calculate a floor, which is largest integral\nvalue not greater than argument.\n \nx - A float point number.\n \nreturn floor of x.")
		.def("agauss", &EMAN::Util::agauss, args("a", "dx", "dy", "dz", "d"), "Calculate Gaussian value.  a * exp(-(dx * dx + dy * dy + dz * dz) / d)\n \na - amplitude\ndx - x center\ndy - y center\ndz - z center\nd - width of gaussian\n \nreturn The Gaussian value.")
		.def("get_min", (int (*)(int, int))&EMAN::Util::get_min, args("f1", "f2"), "Get the minimum of 2 numbers.\n \nf1 - The first number.\nf2 - The second number.\n \nreturn The minimum of 2 numbers.")
		.def("get_min", (int (*)(int, int, int))&EMAN::Util::get_min, args("f1", "f2", "f3"), "Get the minimum of 3 numbers.\n \nf1 - The first number.\nf2 - The second number.\nf3 - The third number.\n \nreturn The minimum of 3 numbers.")
		.def("get_min", (float (*)(float, float))&EMAN::Util::get_min, args("f1", "f2"), "Get the minimum of 2 numbers.\n \nf1 - The first number.\nf2 - The second number.\n \nreturn The minimum of 2 numbers.")
		.def("get_min", (float (*)(float, float, float))&EMAN::Util::get_min, args("f1", "f2", "f3"), "Get the minimum of 3 numbers.\n \nf1 - The first number.\nf2 - The second number.\nf3 - The third number.\n \nreturn The minimum of 3 numbers.")
		.def("get_min", (float (*)(float, float, float, float))&EMAN::Util::get_min, args("f1", "f2", "f3", "f4"), "Get the minimum of 4 numbers.\n \nf1 - The first number.\nf2 - The second number.\nf3 - The third number.\nf4 - The fourth number.\n \nreturn The minimum of 3 numbers.")
		.def("get_max", (float (*)(float, float))&EMAN::Util::get_max, args("f1", "f2"), "Get the maximum of 2 numbers.\n \nf1 - The first number.\nf2 - The second number.\n \nreturn The maximum of 2 numbers.")
		.def("get_max", (float (*)(float, float, float))&EMAN::Util::get_max, args("f1", "f2", "f3"), "Get the maximum of 3 numbers.\n \nf1 - The first number.\nf2 - The second number.\nf3 - The third number.\n \nreturn The maximum of 3 numbers.")
		.def("get_max", (float (*)(float, float, float, float))&EMAN::Util::get_max, args("f1", "f2", "f3", "f4"), "Get the maximum of 3 numbers.\n \nf1 - The first number.\nf2 - The second number.\nf3 - The third number.\nf4 - The fourth number.\n \nreturn The maximum of 4 numbers.")
		.def("get_stats", &EMAN::Util::get_stats, args("data"), "Get the mean, standard deviation, skewness and kurtosis of the input data\n \ndata - the vector of input data\n \nexception EmptyContainerException when the argument vector is empty")
		.def("get_stats_cstyle", &EMAN::Util::get_stats_cstyle, args("data"), "Get the mean, standard deviation, skewness and kurtosis of the input data\n \ndata - the vector of input data\n \nexception EmptyContainerException when the argument vector is empty\n\nPerforms the same calculations as in get_stats, but uses a single pass, optimized c approach\nShould perform better than get_stats")
		.def("angle_sub_2pi", &EMAN::Util::angle_sub_2pi, args("x", "y"), "Calculate the difference of 2 angles and makes the\nequivalent result to be less than Pi.\n \nx - The first angle.\ny - The second angle.\n \nreturn The difference of 2 angles.")
		.def("angle_sub_pi", &EMAN::Util::angle_sub_pi, args("x", "y"), "Calculate the difference of 2 angles and makes the\nequivalent result to be less than Pi/2.\n \nx - The first angle.\ny - The second angle.\n \nreturn The difference of 2 angles.")
#ifndef _WIN32
		.def("recv_broadcast", &EMAN::Util::recv_broadcast, args("port"), "")
#endif	//_WIN32
		.def("get_time_label", &EMAN::Util::get_time_label, "Get the current time in a string with format 'mm/dd/yyyy hh:mm'.\n \nreturn The current time string.")
		.def("eman_copysign", &EMAN::Util::eman_copysign, args("a", "b"), "copy sign of a number. return a value whose absolute value\nmatches that of 'a', but whose sign matches that of 'b'.  If 'a'\nis a NaN, then a NaN with the sign of 'b' is returned.\nIt is exactly copysign() on non-Windows system.\n \na - The first number.\nb - The second number.\n \nreturn Copy sign of a number.")
		.def("eman_erfc", &EMAN::Util::eman_erfc, args("x"), "complementary error function. It is exactly erfc() on\nnon-Windows system. On Windows, it tries to simulate erfc().\n \nThe erf() function returns the error function of x; defined as\nerf(x) = 2/sqrt(pi)* integral from 0 to x of exp(-t*t) dt\n \nThe erfc() function returns the complementary error function of x, that\nis 1.0 - erf(x).\n \nx - A float number.\n \nreturn The complementary error function of x.")
		.def("twoD_fine_ali", &EMAN::Util::twoD_fine_ali, args("image", "refim", "mask", "ang", "sxs", "sys"), "")
		.def("twoD_fine_ali_G", &EMAN::Util::twoD_fine_ali_G, args("image", "refim", "mask", "kb", "ang", "sxs", "sys"), "")
		.def("twoD_fine_ali_SD", &EMAN::Util::twoD_fine_ali_SD, args("image", "refim", "mask", "ang", "sxs", "sys"), "")
		.def("twoD_fine_ali_SD_G", &EMAN::Util::twoD_fine_ali_SD_G, args("image", "refim", "mask", "kb", "ang", "sxs", "sys"), "")
		.def("twoD_to_3D_ali", &EMAN::Util::twoD_to_3D_ali, args("volft", "kb", "refim", "mask", "phi", "theta", "psi", "sxs", "sxy"), "")
		.def("multi_align_error", (vector<float> (*)(vector<float>, vector<float>, int))&EMAN::Util::multi_align_error, args("args", "all_ali_params", "d"), "")
		.def("multi_align_error_func", &EMAN::Util::multi_align_error_func, args("args", "all_ali_params", "nima", "num_ali", "d"), "")
		.def("multi_align_error_func2", &EMAN::Util::multi_align_error_func2, args("args", "all_ali_params", "nima", "num_ali", "d"), "")
		.def("multi_align_error_dfunc", &EMAN::Util::multi_align_error_dfunc, args("args", "all_ali_params", "nima", "num_ali", "g", "d"), "")
		.def("get_biggest_cluster", &EMAN::Util::get_biggest_cluster, return_value_policy< manage_new_object >(), args("mg"), "")
		.def("get_slice", &EMAN::Util::get_slice, return_value_policy< manage_new_object >(), args("vol", "dim", "index"), "This function returns a 2-D slice from a 3-D EMData object\ndim denotes the slice is perpendicular to which dimension\n1 for x-dimension, 2 for y-dimension and 3 for z-dimension")
		.def("merge_peaks", &EMAN::Util::merge_peaks, args("peak1", "peak2", "p_size"), "")
		.def("point_is_in_triangle_2d", &EMAN::Util::point_is_in_triangle_2d, args("p1", "p2", "p3", "actual_point"), "Determines if a point is in a 2D triangle using the Barycentric method, which is\na fast way of performing the query\nTriangle points can be specified in any order\n \np1 - point one\np2 - point two\np3 - point three\nactual_point - the point which might be in the triangle described by p1,p2 and p3\n \nreturn true if the point is in the triangle, false otherwise")
		.def("point_is_in_convex_polygon_2d", &EMAN::Util::point_is_in_convex_polygon_2d, args("p1", "p2", "p3", "p4", "actual_point"), "Determines if a point is in a 2D convex polygon described by 4 points using\nthe Barycentric method, which is a fast way of performing the query.\nThe points must be ordered in the way you would encounter them if you traversed\nthe boundary of the polygon. Direction is irrelevant.\nCould be generalized for polygons with more points\n \np1 - point one\np2 - point two\np3 - point three\np4 - point three\nactual_point - the point which might be in the polygon described by p1,p2,p3 and p4\n \nreturn true if the point is in the polygon, false otherwise")
		.def("sstevd", &pysstevd, args("jobz", "n", "diag", "subdiag", "qmat", "kstep", "fwork", "lfwrk", "iwork", "liwrk"), "")
		.def("snrm2",  &pysnrm2, args("n", "a", "incx"), "")
		.def("sgemv",  &pysgemv, args("trans", "m", "n", "alpha", "a", "lda", "x", "incx", "beta", "y", "incy"), "")
		.def("saxpy",  &pysaxpy, args("n", "alpha", "x", "incx", "y", "incy"), "")
		.def("sdot",   &pysdot, args("n", "x", "incx", "y", "incy"), "")
		.def("readarray", &readarray, args("f", "x", "size"), "")
		.def("k_means_cont_table", &pyk_means_cont_table, args("group1", "group2", "stb", "s1", "s2", "flag"), "k_means_cont_table_ is locate to util_sparx.cpp\nhelper to create the contengency table for partition matching (k-means)\nflag define is the list of stable obj must be store to stb, but the size st\nmust be know before. The trick is first start wihtout the flag to get number\nof elements stable, then again with the flag to get the list. This avoid to\nhave two differents functions for the same thing.")
		.def("bb_enumerateMPI", &pybb_enumerateMPI, args("parts", "classDims", "nParts", "nClasses", "T", "nguesses", "LARGEST_CLASS","J","max_branching","stmult","branchfunc", "LIM"), "bb_enumerateMPI is locate in util_sparx.cpp\nK is the number of classes in each partition (should be the same for all partitions)\nthe first element of each class is its original index in the partition, and second is dummy var\nMPI: if nTop <= 0, then initial prune is called, and the pruned partitions are returned in a 1D array.\nThe first element is reserved for max_levels (the size of the smallest\npartition after pruning).\nif nTop > 0, then partitions are assumed to have been pruned, where only dummy variables of un-pruned partitions are set to 1, and findTopLargest is called\nto find the top weighted matches. The matches, where each match is preceded by its cost, is returned in a one dimensional vector.\nessentially the same as bb_enumerate but with the option to do mpi version.")
		.def("Normalize_ring", &EMAN::Util::Normalize_ring, args("ring", "numr"), "")
		.def("image_mutation", &EMAN::Util::image_mutation, args("img", "mutation_rate"), "")
		.def("list_mutation", &EMAN::Util::list_mutation, args("list", "rate", "min_val", "max_val", "K", "is_mirror"), "")
		.def("get_transform_params", &EMAN::Util::get_transform_params, args("image", "xform", "convention"), "")
		.def("constrained_helix_exhaustive", &EMAN::Util::constrained_helix_exhaustive, args("data", "fdata", "refproj", "rotproj", "dp_dphi_rise_delta", "nphi_phiwobble_range_ywobble_Dsym_nwx_nwy_nwxc_nwyc", "FindPsi", "psi_max", "crefim", "numr", "maxrin", "mode", "cnx", "cny"), "")
		.def("diff_between_matrix_of_3D_parameters_angles", &EMAN::Util::diff_between_matrix_of_3D_parameters_angles, args("all_params", "rotations"), "")
		.def("max_clique", &EMAN::Util::max_clique, args("edges"), "")
		.staticmethod("point_is_in_triangle_2d")
		.staticmethod("point_is_in_convex_polygon_2d")
		.staticmethod("infomask")
		.staticmethod("helixshiftali")
		.staticmethod("snakeshiftali")
		.staticmethod("curhelixshiftali")
		.staticmethod("bsplineBase")
		.staticmethod("bsplineBasedu")
		.staticmethod("convertTocubicbsplineCoeffs")
		.staticmethod("CANG")
		.staticmethod("ener")
		.staticmethod("ener_tot")
		.staticmethod("min_dist_real")
		.staticmethod("min_dist_four")
		.staticmethod("cml_weights")
		.staticmethod("cml_init_rot")
		.staticmethod("cml_update_rot")
		.staticmethod("cml_line_insino")
		.staticmethod("cml_line_insino_all")
		.staticmethod("cml_line_in3d")
		.staticmethod("cml_spin_psi")
		.staticmethod("cml_disc")
		.staticmethod("cml_spin_psi_now")
		.staticmethod("set_line")
		.staticmethod("cml_prepare_line")
		.staticmethod("sstrncmp")
		.staticmethod("int2str")
		.staticmethod("square_sum")
		.staticmethod("Polar2D")
		.staticmethod("get_time_label")
		.staticmethod("get_max")
		.staticmethod("mul_img")
		.staticmethod("div_img")
		.staticmethod("square_img")
		.staticmethod("div_filter")
		.staticmethod("set_freq")
		.staticmethod("histogram")
#ifndef _WIN32
		.staticmethod("recv_broadcast")
#endif	//_WIN32
        .staticmethod("Crosrng_e")
        .staticmethod("Crosrng_rand_e")
		.staticmethod("Crosrng_ew")
		.staticmethod("Crosrng_ms")
		.staticmethod("Crosrng_ms_delta")
		.staticmethod("Crosrng_sm_psi")
		.staticmethod("Crosrng_ns")
		.staticmethod("Crosrng_msg")
		.staticmethod("Crosrng_msg_s")
		.staticmethod("Crosrng_msg_m")
		.staticmethod("Crosrng_msg_vec_p")
		.staticmethod("alrl_ms")
		.staticmethod("window")
		.staticmethod("angle_sub_2pi")
		.staticmethod("WTF")
		.staticmethod("coveig")
		.staticmethod("coveig_for_py")
		.staticmethod("tf")
		.staticmethod("cmp1")
		.staticmethod("cmp2")
		.staticmethod("hist_comp_freq")
		.staticmethod("even_angles")
		.staticmethod("get_min")
		.staticmethod("eman_erfc")
		.staticmethod("bilinear_interpolate")
		.staticmethod("round")
		.staticmethod("square")
		.staticmethod("areav_")
		.staticmethod("TwoDTestFunc")
		.staticmethod("set_randnum_seed")
		.staticmethod("get_randnum_seed")
		.staticmethod("get_frand")
		.staticmethod("get_irand")
		.staticmethod("BPCQ")
		.staticmethod("trilinear_interpolate")
		.staticmethod("bilinear")
		.staticmethod("fftc_d")
		.staticmethod("addn_img")
		.staticmethod("pw_extract")
		.staticmethod("muln_img")
		.staticmethod("divn_img")
		.staticmethod("squaren_img")
		.staticmethod("divn_filter")
		.staticmethod("ctf_img")
		.staticmethod("ctf2_rimg")
		.staticmethod("ctf_rimg")
		.staticmethod("pack_complex_to_real")
		.staticmethod("fftc_q")
		.staticmethod("remove_filename_ext")
		.staticmethod("save_data")
		.staticmethod("quadri")
		.staticmethod("get_pixel_conv_new")
		.staticmethod("cluster_pairwise")
		.staticmethod("cluster_equalsize")
		.staticmethod("vareas")
		.staticmethod("decimate")
		.staticmethod("eman_copysign")
		.staticmethod("sbasename")
		.staticmethod("madn_scalar")
		.staticmethod("local_inner_product")
		.staticmethod("move_points")
		.staticmethod("splint")
		.staticmethod("compress_image_mask")
		.staticmethod("calc_best_fft_size")
		.staticmethod("nonconvex")
		.staticmethod("windowdot")
		.staticmethod("calc_bilinear_least_square")
		.staticmethod("angle_sub_pi")
		.staticmethod("sub_img")
		.staticmethod("reconstitute_image_mask")
		.staticmethod("vrdg")
		.staticmethod("update_fav")
		.staticmethod("change_filename_ext")
		.staticmethod("cyclicshift")
		.staticmethod("sub_fav")
		.staticmethod("get_gauss_rand")
		.staticmethod("mult_scalar")
		.staticmethod("svdcmp")
		.staticmethod("agauss")
		.staticmethod("WTM")
		.staticmethod("subn_img")
		.staticmethod("multiref_polar_ali_2d")
		.staticmethod("multiref_polar_ali_2d_peaklist")
		.staticmethod("multiref_polar_ali_2d_peaklist_local")
		.staticmethod("nearest_ang")
		.staticmethod("assign_groups")
		.staticmethod("assign_projangles")
		.staticmethod("nearestk_to_refdir")
		.staticmethod("group_proj_by_phitheta")
		.staticmethod("multiref_polar_ali_2d_delta")
		.staticmethod("multiref_polar_ali_2d_nom")
		.staticmethod("multiref_polar_ali_2d_local")
		.staticmethod("shc")
		.staticmethod("shc_multipeaks")
		.staticmethod("multiref_polar_ali_2d_local_psi")
		.staticmethod("multiref_polar_ali_helical")
		.staticmethod("multiref_polar_ali_helical_local")
		.staticmethod("multiref_polar_ali_helical_90")
		.staticmethod("multiref_polar_ali_helical_90_local")
		.staticmethod("multiref_polar_ali_helicon_local")
		.staticmethod("multiref_polar_ali_helicon_90_local")
		.staticmethod("multiref_peaks_ali2d")
		//.staticmethod("multiref_peaks_ali")
		.staticmethod("multiref_peaks_compress_ali2d")
		.staticmethod("ali2d_ccf_list")
		.staticmethod("ali2d_ccf_list_snake")
		.staticmethod("hypot3")
		.staticmethod("histc")
		.staticmethod("get_stats")
		.staticmethod("get_stats_cstyle")
		.staticmethod("fftr_d")
		.staticmethod("mad_scalar")
		.staticmethod("fftr_q")
		.staticmethod("is_file_exist")
		.staticmethod("pad")
		.staticmethod("add_img")
		.staticmethod("add_img_abs")
		.staticmethod("im_diff")
		.staticmethod("Polar2Dmi")
		.staticmethod("fast_floor")
		.staticmethod("Frngs")
		.staticmethod("Frngs_inv")
		.staticmethod("Applyws")
		.staticmethod("eval")
		.staticmethod("add_img2")
		.staticmethod("mul_scalar")
		.staticmethod("get_filename_ext")
		.staticmethod("Polar2Dm")
		.staticmethod("twoD_fine_ali")
		.staticmethod("twoD_fine_ali_G")
		.staticmethod("twoD_fine_ali_SD")
		.staticmethod("twoD_fine_ali_SD_G")
		.staticmethod("twoD_to_3D_ali")
		.staticmethod("multi_align_error")
		.staticmethod("multi_align_error_func")
		.staticmethod("multi_align_error_func2")
		.staticmethod("multi_align_error_dfunc")
		.staticmethod("get_biggest_cluster")
		.staticmethod("merge_peaks")
		.staticmethod("get_slice")
		.staticmethod("sstevd")
		.staticmethod("sgemv")
		.staticmethod("snrm2")
		.staticmethod("saxpy")
		.staticmethod("sdot")
		.staticmethod("readarray")
		.staticmethod("k_means_cont_table")
		.staticmethod("bb_enumerateMPI")
		.staticmethod("Normalize_ring")
		.staticmethod("image_mutation")
		.staticmethod("list_mutation")
		.staticmethod("get_transform_params")
		.staticmethod("constrained_helix_exhaustive")
		.staticmethod("diff_between_matrix_of_3D_parameters_angles")
		.staticmethod("max_clique")
	);

    scope* EMAN_Util_sincBlackman_scope = new scope(
    class_< EMAN::Util::sincBlackman, EMAN_Util_sincBlackman_Wrapper >("sincBlackman", init< const EMAN::Util::sincBlackman& >())
        .def(init< int, float, optional< int > >())
        .def("sBwin_tab", &EMAN::Util::sincBlackman::sBwin_tab, args("x"), "")
        .def("get_sB_size", &EMAN::Util::sincBlackman::get_sB_size, "Return the size of the kernel")
    );
    delete EMAN_Util_sincBlackman_scope;

    scope* EMAN_Util_KaiserBessel_scope = new scope(
    class_< EMAN::Util::KaiserBessel, EMAN_Util_KaiserBessel_Wrapper >("KaiserBessel",
    		"1-D Kaiser-Bessel window function class.\n"
    		"(It's a class so that the windowing parameters may be\n"
    		"instantiated and held in the instance object.)\n\n"
    		"The I0 version can be tabulated and interpolated upon\n"
    		"demand, but the max error needs to be checked.  The\n"
    		"\"vtable\" parameter corresponds to the maximum value of x\n"
    		"for which the I0 window is non-zero.  Setting \"vtable\"\n"
    		"different from \"v\" corresponds to a change in units of x.\n"
    		"In practice, it is often handy to replace x in some sort\n"
    		"of absolute units with x described in terms of grid\n"
    		"intervals.\n\n"
    		"The get_kbsinh_win and get_kbi0_win functions return\n"
    		"single-argument function objects, which is what a\n"
    		"generic routine is likely to want.\n\n"
    		"see P. A. Penczek, R. Renka, and H. Schomberg, J. Opt. Soc. Am. _21_, 449 (2004)",
    		init< const EMAN::Util::KaiserBessel& >())
        .def(init< float, int, float, float, int, optional< float, int > >())
        .def("sinhwin", &EMAN::Util::KaiserBessel::sinhwin, &EMAN_Util_KaiserBessel_Wrapper::default_sinhwin, args("x"), "Kaiser-Bessel Sinh window function")
        .def("i0win", &EMAN::Util::KaiserBessel::i0win, &EMAN_Util_KaiserBessel_Wrapper::default_i0win, args("x"), "Kaiser-Bessel I0 window function")
        .def("I0table_maxerror", &EMAN::Util::KaiserBessel::I0table_maxerror, "Compute the maximum error in the table")
        .def("dump_table", &EMAN::Util::KaiserBessel::dump_table)
        .def("i0win_tab", &EMAN::Util::KaiserBessel::i0win_tab, args("x"), "Kaiser-Bessel I0 window function (uses table lookup)")
        .def("get_window_size", &EMAN::Util::KaiserBessel::get_window_size, "Return the size of the I0 window")
        .def("get_kbsinh_win", &EMAN::Util::KaiserBessel::get_kbsinh_win, "Sinh window function object factory")
        .def("get_kbi0_win", &EMAN::Util::KaiserBessel::get_kbi0_win, "I0 window function object factory")
    );

    class_< EMAN::Util::KaiserBessel::kbsinh_win >("kbsinh_win", "Sinh window function object", init< const EMAN::Util::KaiserBessel::kbsinh_win& >())
        .def(init< EMAN::Util::KaiserBessel& >())
        .def("get_window_size", &EMAN::Util::KaiserBessel::kbsinh_win::get_window_size)
        .def("__call__", &EMAN::Util::KaiserBessel::kbsinh_win::operator ())
    ;


    class_< EMAN::Util::KaiserBessel::kbi0_win >("kbi0_win", "I0 window function object", init< const EMAN::Util::KaiserBessel::kbi0_win& >())
        .def(init< EMAN::Util::KaiserBessel& >())
        .def("get_window_size", &EMAN::Util::KaiserBessel::kbi0_win::get_window_size)
        .def("__call__", &EMAN::Util::KaiserBessel::kbi0_win::operator ())
    ;

    delete EMAN_Util_KaiserBessel_scope;


    class_< EMAN::Util::FakeKaiserBessel, bases< EMAN::Util::KaiserBessel > , EMAN_Util_FakeKaiserBessel_Wrapper >("FakeKaiserBessel", init< const EMAN::Util::FakeKaiserBessel& >())
        .def(init< float, int, float, float, int, optional< float, int > >())
        .def("sinhwin", (float (EMAN::Util::FakeKaiserBessel::*)(float) const)&EMAN::Util::FakeKaiserBessel::sinhwin, (float (EMAN_Util_FakeKaiserBessel_Wrapper::*)(float) const)&EMAN_Util_FakeKaiserBessel_Wrapper::default_sinhwin)
        .def("i0win", (float (EMAN::Util::FakeKaiserBessel::*)(float) const)&EMAN::Util::FakeKaiserBessel::i0win, (float (EMAN_Util_FakeKaiserBessel_Wrapper::*)(float) const)&EMAN_Util_FakeKaiserBessel_Wrapper::default_i0win)
        .def("build_I0table", (void (EMAN::Util::FakeKaiserBessel::*)() )&EMAN::Util::FakeKaiserBessel::build_I0table, (void (EMAN_Util_FakeKaiserBessel_Wrapper::*)())&EMAN_Util_FakeKaiserBessel_Wrapper::default_build_I0table)
    ;


    class_< EMAN::Util::Gaussian >("Gaussian",
    		"Gaussian function class.\n\n"
    		"Usage:\n"
    		"   Gaussian gauss(sigma);\n"
    		"   float g = gauss(x);",
    		init< const EMAN::Util::Gaussian& >())
        .def(init< optional< float > >())
        .def("__call__", &EMAN::Util::Gaussian::operator ())
    ;


    class_< EMAN::Util::tmpstruct >("tmpstruct", init<  >())
        .def(init< const EMAN::Util::tmpstruct& >())
        .def_readwrite("theta1", &EMAN::Util::tmpstruct::theta1)
        .def_readwrite("phi1", &EMAN::Util::tmpstruct::phi1)
        .def_readwrite("key1", &EMAN::Util::tmpstruct::key1)
    ;

    delete EMAN_Util_scope;


	
    scope* EMAN_EMUtil_scope = new scope(
    class_< EMAN::EMUtil >("EMUtil", init<  >())
        .def(init< const EMAN::EMUtil& >())
        .def("vertical_acf", &EMAN::EMUtil::vertical_acf, return_value_policy< manage_new_object >(), args("image", "maxdy"), "")
        .def("make_image_median", &EMAN::EMUtil::make_image_median, return_value_policy< manage_new_object >(), args("image_list"), "")
        .def("get_image_ext_type", &EMAN::EMUtil::get_image_ext_type, args("file_ext"), "Get an image's format type from its filename extension.\n \nfile_ext - File extension.\n \nreturn image format type.")
        .def("get_image_type", &EMAN::EMUtil::get_image_type, args("filename"), "Get an image's format type by processing the first 1K of the image.\n \nfilename - Image file name.\n \nreturn image format type.")
        .def("get_image_count", &EMAN::EMUtil::get_image_count, args("filename"), "Get the number of images in an image file.\n \nfilename Image file name.\n \nreturn Number of images in the given file.")
        .def("get_imageio", &EMAN::EMUtil::get_imageio, EMAN_EMUtil_get_imageio_overloads_2_3(args("filename", "rw_mode", "image_type"), "Get an ImageIO object. It may be a newly created\nobject. Or an object stored in the cache.\n \nfilename - Image file name.\nrw_mode - ImageIO read/write mode.\nimage_type - Image format type.(default=IMAGE_UNKNOWN)\n \nreturn An ImageIO object.")[ return_internal_reference< 1 >() ])
        .def("get_imagetype_name", &EMAN::EMUtil::get_imagetype_name, args("type"), "Give each image type a meaningful name.\n \ntype - Image format type.\n \nreturn A name for that type.")
        .def("get_datatype_string", &EMAN::EMUtil::get_datatype_string, args("type"), "Give each data type a meaningful name\n \ntype - the EMDataType\n \nreturn a name for that data type")
        .def("process_ascii_region_io", &EMAN::EMUtil::process_ascii_region_io, args("data", "file", "rw_mode", "image_index", "mode_size", "nx", "ny", "nz", "area", "has_index_line", "nitems_per_line", "outformat"), "Works for regions that are outside the image data dimension area.\nThe only function that calls this is in xplorio.cpp - that function\nthrows if the region is invalid.")
        .def("read_raw_emdata", read_raw_emdata, args("emdata", "filename", "offset", "rw_mode", "image_index", "mode",  "area"), "This function will read raw binary data from disk into an existing EMData objects (with the correct dimensions). mode: 0- float, 1- unsigned 32 bit int, 2- signed 32 bit int, 3- unsigned 16 bit int, 4- signed 16 bit int")
        .def("dump_dict", &EMAN::EMUtil::dump_dict, args("dict"), "Dump a Dict object.\n \ndict - A Dict object")
        .def("is_same_size", &EMAN::EMUtil::is_same_size, args("image1", "image2"), "Check whether two EMData images are of the same size.\n \nimage1 - The first EMData image.\nimage2 - The second EMData image.return Whether two EMData images are of the same size.")
        .def("is_same_ctf", &EMAN::EMUtil::is_same_ctf,args("image1", "image2"), "Check whether two EMData images have the same CTF parameters.\n \nimage1 - The first EMData image.\nimage2 The second EMData image.\n \nreturn whether two EMData images have the same CTF.")
        .def("is_complex_type", &EMAN::EMUtil::is_complex_type, args("datatype"), "")
        .def("is_valid_filename", &EMAN::EMUtil::is_valid_filename, args("filename"), "Ask whether or not the given filename is a valid EM image filename\nThis is the same thing as checking whether or not the return value of EMUtil.get_image_ext_type\nis IMAGE_UNKNOWN\n \nfilename - Image file name.\n \nreturn whether or not it is a valid filename")
        .def("jump_lines", &EMAN::EMUtil::jump_lines, args("file", "nlines"), "")
        .def("get_euler_names", &EMAN::EMUtil::get_euler_names, args("euler_type"), "")
        .def("get_all_attributes", &EMAN::EMUtil::get_all_attributes, args("file_name", "attr_name"), "Get an attribute from a stack of image, returned as a vector\n \nfile_name - the image file name\nattr_name - The header attribute name.\n \nreturn the vector of attribute value\n \nexception - NotExistingObjectException when access an non-existing attribute\nexception - InvalidCallException when call this function for a non-stack image")
		.def("cuda_available", &EMAN::EMUtil::cuda_available)
#ifdef EM_HDF5
		.def("read_hdf_attribute", &EMAN::EMUtil::read_hdf_attribute, EMAN_EMUtil_read_hdf_attribute_2_3(args("filename", "key", "image_index"), "Retrive a single attribute value from a HDF5 image file.\n \nfilename - HDF5 image's file name\nkey - the attribute's key name\nimage_index - the image index, default=0\n \nreturn the attribute value for the given key"))
		.def("write_hdf_attribute", &EMAN::EMUtil::write_hdf_attribute, EMAN_EMUtil_write_hdf_attribute_3_4(args("filename", "key", "value", "image_index"), "Write a single attribute value from a HDF5 image file.\n \nfilename - HDF5 image's file name\nkey - the attribute's key name\nvalue - the attribute's value\nimage_index - the image index, default=0\n \nreturn 0 for success"))
		.def("delete_hdf_attribute", &EMAN::EMUtil::delete_hdf_attribute, EMAN_EMUtil_delete_hdf_attribute_2_3(args("filename", "key", "image_index"), "Delete a single attribute from a HDF5 image file.\n \nfilename - HDF5 image's file name\nkey - the attribute's key name\nimage_index - the image index, default=0\n \nreturn 0 for success, -1 for failure."))
#endif	//EM_HDF5
        .staticmethod("cuda_available")
        .staticmethod("read_raw_emdata")
        .staticmethod("vertical_acf")
        .staticmethod("get_datatype_string")
        .staticmethod("dump_dict")
        .staticmethod("get_all_attributes")
        .staticmethod("get_imageio")
        .staticmethod("get_image_count")
        .staticmethod("get_imagetype_name")
        .staticmethod("get_image_type")
        .staticmethod("is_same_size")
        .staticmethod("is_valid_filename")
        .staticmethod("make_image_median")
        .staticmethod("jump_lines")
        .staticmethod("get_euler_names")
        .staticmethod("is_same_ctf")
        .staticmethod("get_image_ext_type")
        .staticmethod("process_ascii_region_io")
        .staticmethod("is_complex_type")
#ifdef EM_HDF5
        .staticmethod("read_hdf_attribute")
        .staticmethod("write_hdf_attribute")
        .staticmethod("delete_hdf_attribute")
#endif	//EM_HDF5
    );

    enum_< EMAN::EMUtil::EMDataType >("EMDataType")
        .value("EM_SHORT_COMPLEX", EMAN::EMUtil::EM_SHORT_COMPLEX)
        .value("EM_SHORT", EMAN::EMUtil::EM_SHORT)
        .value("EM_UCHAR", EMAN::EMUtil::EM_UCHAR)
        .value("EM_FLOAT_COMPLEX", EMAN::EMUtil::EM_FLOAT_COMPLEX)
        .value("EM_CHAR", EMAN::EMUtil::EM_CHAR)
        .value("EM_INT", EMAN::EMUtil::EM_INT)
        .value("EM_USHORT", EMAN::EMUtil::EM_USHORT)
        .value("EM_USHORT_COMPLEX", EMAN::EMUtil::EM_USHORT_COMPLEX)
        .value("EM_UNKNOWN", EMAN::EMUtil::EM_UNKNOWN)
        .value("EM_UINT", EMAN::EMUtil::EM_UINT)
        .value("EM_DOUBLE", EMAN::EMUtil::EM_DOUBLE)
        .value("EM_FLOAT", EMAN::EMUtil::EM_FLOAT)
    ;


    enum_< EMAN::EMUtil::ImageType >("ImageType")
        .value("IMAGE_FITS", EMAN::EMUtil::IMAGE_FITS)
        .value("IMAGE_ICOS", EMAN::EMUtil::IMAGE_ICOS)
        .value("IMAGE_UNKNOWN", EMAN::EMUtil::IMAGE_UNKNOWN)
        .value("IMAGE_GATAN2", EMAN::EMUtil::IMAGE_GATAN2)
        .value("IMAGE_EMIM", EMAN::EMUtil::IMAGE_EMIM)
        .value("IMAGE_VTK", EMAN::EMUtil::IMAGE_VTK)
        .value("IMAGE_MRC", EMAN::EMUtil::IMAGE_MRC)
        .value("IMAGE_LST", EMAN::EMUtil::IMAGE_LST)
        .value("IMAGE_DM3", EMAN::EMUtil::IMAGE_DM3)
        .value("IMAGE_DM4", EMAN::EMUtil::IMAGE_DM4)
        .value("IMAGE_PNG", EMAN::EMUtil::IMAGE_PNG)
        .value("IMAGE_XPLOR", EMAN::EMUtil::IMAGE_XPLOR)
        .value("IMAGE_AMIRA", EMAN::EMUtil::IMAGE_AMIRA)
        .value("IMAGE_SAL", EMAN::EMUtil::IMAGE_SAL)
        .value("IMAGE_SPIDER", EMAN::EMUtil::IMAGE_SPIDER)
        .value("IMAGE_SINGLE_SPIDER", EMAN::EMUtil::IMAGE_SINGLE_SPIDER)
        .value("IMAGE_PGM", EMAN::EMUtil::IMAGE_PGM)
        .value("IMAGE_EM", EMAN::EMUtil::IMAGE_EM)
        .value("IMAGE_TIFF", EMAN::EMUtil::IMAGE_TIFF)
        .value("IMAGE_IMAGIC", EMAN::EMUtil::IMAGE_IMAGIC)
        .value("IMAGE_HDF", EMAN::EMUtil::IMAGE_HDF)
        .value("IMAGE_JPEG", EMAN::EMUtil::IMAGE_JPEG)
        .value("IMAGE_V4L", EMAN::EMUtil::IMAGE_V4L)
        .value("IMAGE_PIF", EMAN::EMUtil::IMAGE_PIF)
        .value("IMAGE_DF3", EMAN::EMUtil::IMAGE_DF3)
    ;

    delete EMAN_EMUtil_scope;

    class_< EMAN::ImageSort >("ImageSort", init< const EMAN::ImageSort& >())
        .def(init< int >())
        .def("sort", &EMAN::ImageSort::sort)
        .def("set", &EMAN::ImageSort::set)
        .def("get_index", &EMAN::ImageSort::get_index)
        .def("get_score", &EMAN::ImageSort::get_score)
        .def("size", &EMAN::ImageSort::size)
    ;

    class_< EMAN::TestUtil >("TestUtil", "TestUtil defines function assisting testing of EMAN2.", init<  >())
        .def(init< const EMAN::TestUtil& >())
        .add_static_property("EMDATA_HEADER_EXT", make_getter(EMAN::TestUtil::EMDATA_HEADER_EXT))
        .add_static_property("EMDATA_DATA_EXT", make_getter(EMAN::TestUtil::EMDATA_DATA_EXT))
        .def("get_debug_int", &EMAN::TestUtil::get_debug_int)
        .def("get_debug_float", &EMAN::TestUtil::get_debug_float)
        .def("get_debug_string", &EMAN::TestUtil::get_debug_string)
        .def("get_debug_transform", &EMAN::TestUtil::get_debug_transform)
        .def("get_debug_image", &EMAN::TestUtil::get_debug_image)
        .def("get_golden_image", &EMAN::TestUtil::get_golden_image)
        .def("to_emobject", &EMAN::TestUtil::to_emobject)
        .def("emobject_to_py", (EMAN::EMObject (*)(bool))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(unsigned int))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(int))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(float))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(double))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(const std::string&))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(EMAN::EMData*))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(EMAN::XYData*))&EMAN::TestUtil::emobject_to_py)
		.def("emobject_to_py", (EMAN::EMObject (*)(EMAN::Transform*))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_to_py", (EMAN::EMObject (*)(EMAN::Ctf*))&EMAN::TestUtil::emobject_to_py)
        .def("emobject_farray_to_py", (EMAN::EMObject (*)())&EMAN::TestUtil::emobject_farray_to_py)
        .def("emobject_strarray_to_py", (EMAN::EMObject (*)())&EMAN::TestUtil::emobject_strarray_to_py)
        .def("emobject_transformarray_to_py", (EMAN::EMObject (*)())&EMAN::TestUtil::emobject_transformarray_to_py)
        .def("test_IntPoint", &EMAN::TestUtil::test_IntPoint)
        .def("test_FloatPoint", &EMAN::TestUtil::test_FloatPoint)
        .def("test_IntSize", &EMAN::TestUtil::test_IntSize)
        .def("test_FloatSize", &EMAN::TestUtil::test_FloatSize)
        .def("test_Vec3i", &EMAN::TestUtil::test_Vec3i)
        .def("test_Vec3f", &EMAN::TestUtil::test_Vec3f)
        .def("test_vector_int", &EMAN::TestUtil::test_vector_int)
        .def("test_vector_float", &EMAN::TestUtil::test_vector_float)
        .def("test_vector_long", &EMAN::TestUtil::test_vector_long)
        .def("test_vector_string", &EMAN::TestUtil::test_vector_string)
        .def("test_vector_emdata", &EMAN::TestUtil::test_vector_emdata)
        .def("test_vector_pixel", &EMAN::TestUtil::test_vector_pixel)
        .def("test_map_int", &EMAN::TestUtil::test_map_int)
        .def("test_map_long", &EMAN::TestUtil::test_map_long)
        .def("test_map_float", &EMAN::TestUtil::test_map_float)
        .def("test_map_string", &EMAN::TestUtil::test_map_string)
        .def("test_map_emobject", &EMAN::TestUtil::test_map_emobject)
        .def("test_map_vecstring", &EMAN::TestUtil::test_map_vecstring)
        .def("test_dict", &EMAN::TestUtil::test_dict)
        .def("dump_image_from_file", &EMAN::TestUtil::dump_image_from_file)
        .def("dump_emdata", &EMAN::TestUtil::dump_emdata)
        .def("check_image", &EMAN::TestUtil::check_image, EMAN_TestUtil_check_image_overloads_1_2())
        .def("set_progname", &EMAN::TestUtil::set_progname)
        .def("make_image_file", &EMAN::TestUtil::make_image_file, EMAN_TestUtil_make_image_file_overloads_2_6())
        .def("verify_image_file", &EMAN::TestUtil::verify_image_file, EMAN_TestUtil_verify_image_file_overloads_2_6())
        .def("make_image_file2", &EMAN::TestUtil::make_image_file2, EMAN_TestUtil_make_image_file2_overloads_2_6())
        .def("verify_image_file2", &EMAN::TestUtil::verify_image_file2, EMAN_TestUtil_verify_image_file2_overloads_2_6())
        .staticmethod("test_Vec3f")
        .staticmethod("verify_image_file2")
        .staticmethod("test_vector_float")
        .staticmethod("test_Vec3i")
        .staticmethod("test_map_int")
        .staticmethod("make_image_file")
        .staticmethod("emobject_strarray_to_py")
        .staticmethod("emobject_transformarray_to_py")
        .staticmethod("test_map_vecstring")
        .staticmethod("get_debug_int")
        .staticmethod("test_vector_long")
        .staticmethod("emobject_farray_to_py")
        .staticmethod("test_IntPoint")
        .staticmethod("dump_image_from_file")
        .staticmethod("dump_emdata")
        .staticmethod("set_progname")
        .staticmethod("test_map_float")
        .staticmethod("test_IntSize")
        .staticmethod("test_FloatSize")
        .staticmethod("verify_image_file")
        .staticmethod("test_map_long")
        .staticmethod("to_emobject")
        .staticmethod("make_image_file2")
        .staticmethod("emobject_to_py")
        .staticmethod("get_golden_image")
        .staticmethod("test_map_string")
        .staticmethod("test_vector_int")
        .staticmethod("test_vector_emdata")
        .staticmethod("test_vector_string")
        .staticmethod("get_debug_string")
        .staticmethod("get_debug_transform")
        .staticmethod("test_FloatPoint")
        .staticmethod("test_dict")
        .staticmethod("test_map_emobject")
        .staticmethod("test_vector_pixel")
        .staticmethod("get_debug_float")
        .staticmethod("get_debug_image")
        .staticmethod("check_image")
    ;
}

