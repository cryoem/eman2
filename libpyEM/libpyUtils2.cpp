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

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_im_diff_overloads_2_3, EMAN::Util::im_diff, 2, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_TwoDTestFunc_overloads_5_7, EMAN::Util::TwoDTestFunc, 5, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_even_angles_overloads_1_5, EMAN::Util::even_angles, 1, 5)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_decimate_overloads_2_4, EMAN::Util::decimate, 2, 4)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_window_overloads_2_7, EMAN::Util::window, 2, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_pad_overloads_2_8, EMAN::Util::pad, 2, 8)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_tf_overloads_2_7, EMAN::Util::tf, 2, 7)

BOOST_PYTHON_FUNCTION_OVERLOADS(EMAN_Util_ctf_img_overloads_5_12, EMAN::Util::ctf_img, 5, 12)


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
        if( datatype != 'i' )
        {
        	PyErr_SetString(PyExc_ValueError, "expected an integer PyArrayObject for get_iptr");
		return NULL;
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
    if( !PyFile_Check(f.ptr()) )
    {
        std::cout << "Error: expecting a file object" << std::endl;
        return;
    }

    FILE*  fh = PyFile_AsFile( f.ptr() );
    float* fx = get_fptr( x );

    fread( fx, sizeof(float), size, fh );
}
// Module ======================================================================
BOOST_PYTHON_MODULE(libpyUtils2)
{
	scope* EMAN_Util_scope = new scope(
		 class_< EMAN::Util>("Util", init<  >())
		.def(init< const EMAN::Util& >())
		.def("coveig", &EMAN::Util::coveig)
		.def("coveig_for_py", &EMAN::Util::coveig_for_py)
		.def("ExpMinus4YSqr", &EMAN::Util::ExpMinus4YSqr)
		.def("WTM", &EMAN::Util::WTM)
		.def("WTF", &EMAN::Util::WTF)
		.def("CANG", &EMAN::Util::CANG)
		.def("BPCQ", &EMAN::Util::BPCQ)
		.def("infomask", &EMAN::Util::infomask)
		.def("cluster_pairwise", &EMAN::Util::cluster_pairwise)
		.def("cluster_equalsize", &EMAN::Util::cluster_equalsize)
		.def("vareas", &EMAN::Util::vareas)
		.def("cyclicshift", &EMAN::Util::cyclicshift)
		.def("im_diff", &EMAN::Util::im_diff, EMAN_Util_im_diff_overloads_2_3())
		.def("TwoDTestFunc", &EMAN::Util::TwoDTestFunc, return_value_policy< manage_new_object >(), EMAN_Util_TwoDTestFunc_overloads_5_7())
		.def("splint", &EMAN::Util::splint)
		.def("even_angles", &EMAN::Util::even_angles, EMAN_Util_even_angles_overloads_1_5())
		.def("quadri", &EMAN::Util::quadri)
		.def("get_pixel_conv_new", &EMAN::Util::get_pixel_conv_new)
		.def("bilinear", &EMAN::Util::bilinear)
		.def("Polar2D", &EMAN::Util::Polar2D, return_value_policy< manage_new_object >())
		.def("Polar2Dm", &EMAN::Util::Polar2Dm, return_value_policy< manage_new_object >())
		.def("alrl_ms", &EMAN::Util::alrl_ms)
		.def("Polar2Dmi", &EMAN::Util::Polar2Dmi, return_value_policy< manage_new_object >())
		.def("fftr_q", &EMAN::Util::fftr_q)
		.def("fftr_d", &EMAN::Util::fftr_d)
		.def("fftc_q", &EMAN::Util::fftc_q)
		.def("fftc_d", &EMAN::Util::fftc_d)
		.def("Frngs", &EMAN::Util::Frngs)
		.def("Frngs_inv", &EMAN::Util::Frngs_inv)
		.def("Crosrng_e", &EMAN::Util::Crosrng_e)
		.def("Crosrng_ew", &EMAN::Util::Crosrng_ew)
		.def("Crosrng_ms", &EMAN::Util::Crosrng_ms)
		.def("Crosrng_ns", &EMAN::Util::Crosrng_ns)
		.def("Crosrng_msg", &EMAN::Util::Crosrng_msg, return_value_policy< manage_new_object >())
		.def("Crosrng_msg_s", &EMAN::Util::Crosrng_msg_s, return_value_policy< manage_new_object >())
		.def("Crosrng_msg_m", &EMAN::Util::Crosrng_msg_m, return_value_policy< manage_new_object >())
		.def("Crosrng_msg_vec_p",&EMAN::Util::Crosrng_msg_vec_p )
		.def("prb1d", &EMAN::Util::prb1d)
		.def("update_fav", &EMAN::Util::update_fav)
		.def("sub_fav", &EMAN::Util::sub_fav)
		.def("ener", &EMAN::Util::ener)
		.def("ener_tot", &EMAN::Util::ener_tot)
		.def("min_dist_real", &EMAN::Util::min_dist_real)
		.def("min_dist_four", &EMAN::Util::min_dist_four)
		.def("cml_weights", &EMAN::Util::cml_weights)
		.def("cml_init_rot", &EMAN::Util::cml_init_rot)
		.def("cml_update_rot", &EMAN::Util::cml_update_rot)
		.def("cml_line_insino", &EMAN::Util::cml_line_insino)
	        .def("cml_line_insino_all", &EMAN::Util::cml_line_insino_all)
		.def("cml_line_in3d", &EMAN::Util::cml_line_in3d)
		.def("cml_spin_psi", &EMAN::Util::cml_spin_psi)
		.def("cml_disc", &EMAN::Util::cml_disc)
		.def("cml_prepare_line", &EMAN::Util::cml_prepare_line)
		.def("set_line", &EMAN::Util::set_line)
		.def("decimate", &EMAN::Util::decimate, return_value_policy< manage_new_object >(), EMAN_Util_decimate_overloads_2_4())
		.def("window", &EMAN::Util::window, return_value_policy< manage_new_object >(), EMAN_Util_window_overloads_2_7())
		.def("pad", &EMAN::Util::pad, return_value_policy< manage_new_object >(), EMAN_Util_pad_overloads_2_8())
		.def("histc", &EMAN::Util::histc)
		.def("hist_comp_freq", &EMAN::Util::hist_comp_freq)
		.def("tf", &EMAN::Util::tf, EMAN_Util_tf_overloads_2_7())
		.def("compress_image_mask", &EMAN::Util::compress_image_mask, return_value_policy< manage_new_object >())
		.def("reconstitute_image_mask", &EMAN::Util::reconstitute_image_mask, return_value_policy< manage_new_object >())
		.def("pw_extract", &EMAN::Util::pw_extract)
		.def("eval", &EMAN::Util::eval)
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
		.def("divn_filter", &EMAN::Util::divn_filter, return_value_policy< manage_new_object >())
		.def("mad_scalar", &EMAN::Util::mad_scalar)
		.def("mul_scalar", &EMAN::Util::mul_scalar)
		.def("add_img", &EMAN::Util::add_img)
		.def("add_img2", &EMAN::Util::add_img2)
		.def("sub_img", &EMAN::Util::sub_img)
		.def("mul_img", &EMAN::Util::mul_img)
		.def("div_img", &EMAN::Util::div_img)
		.def("div_filter", &EMAN::Util::div_filter)
		.def("ctf_img", &EMAN::Util::ctf_img, return_value_policy< manage_new_object >())
		.def("pack_complex_to_real", &EMAN::Util::pack_complex_to_real, return_value_policy< manage_new_object >())
		.def("histogram", &EMAN::Util::histogram)
		.def("multiref_polar_ali_2d", &EMAN::Util::multiref_polar_ali_2d)
		.def("multiref_polar_ali_2d_nom", &EMAN::Util::multiref_polar_ali_2d_nom)
		.def("multiref_polar_ali_2d_local", &EMAN::Util::multiref_polar_ali_2d_local)
		.def("multiref_peaks_ali2d", &EMAN::Util::multiref_peaks_ali2d)
		.def("multiref_peaks_compress_ali2d", &EMAN::Util::multiref_peaks_compress_ali2d)
		.def("ali2d_ccf_list", &EMAN::Util::ali2d_ccf_list)
		//.def("multiref_peaks_ali", &EMAN::Util::multiref_peaks_ali)
		.def("move_points", &EMAN::Util::move_points, return_value_policy< manage_new_object >())
		.def("is_file_exist", &EMAN::Util::is_file_exist)
		.def("svdcmp", &EMAN::Util::svdcmp)
		.def("sstrncmp", &EMAN::Util::sstrncmp)
		.def("int2str", &EMAN::Util::int2str)
		.def("change_filename_ext", &EMAN::Util::change_filename_ext)
		.def("remove_filename_ext", &EMAN::Util::remove_filename_ext)
		.def("save_data",(void (*)(float,float,const vector<float> &,const string &))&EMAN::Util::save_data)
		.def("get_filename_ext", &EMAN::Util::get_filename_ext)
		.def("sbasename", &EMAN::Util::sbasename)
		.def("set_randnum_seed", &EMAN::Util::set_randnum_seed)
		.def("get_randnum_seed", &EMAN::Util::get_randnum_seed)
		.def("get_irand", (int (*)(int, int))&EMAN::Util::get_irand)
		.def("get_frand", (float (*)(int, int))&EMAN::Util::get_frand)
		.def("get_frand", (float (*)(float, float))&EMAN::Util::get_frand)
		.def("get_frand", (float (*)(double, double))&EMAN::Util::get_frand)
		.def("get_gauss_rand", &EMAN::Util::get_gauss_rand)
		.def("round", (int (*)(float))&EMAN::Util::round)
		.def("round", (int (*)(double))&EMAN::Util::round)
		.def("bilinear_interpolate", &EMAN::Util::bilinear_interpolate)
		.def("trilinear_interpolate", &EMAN::Util::trilinear_interpolate)
		.def("calc_best_fft_size", &EMAN::Util::calc_best_fft_size)
		.def("calc_bilinear_least_square", &EMAN::Util::calc_bilinear_least_square)
		.def("square", (int (*)(int))&EMAN::Util::square)
		.def("square", (float (*)(float))&EMAN::Util::square)
		.def("square", (float (*)(double))&EMAN::Util::square)
		.def("square_sum", &EMAN::Util::square_sum)
		.def("hypot3", (float (*)(int, int, int))&EMAN::Util::hypot3)
		.def("hypot3", (float (*)(float, float, float))&EMAN::Util::hypot3)
		.def("hypot3", (float (*)(double, double, double))&EMAN::Util::hypot3)
		.def("fast_floor", &EMAN::Util::fast_floor)
		.def("agauss", &EMAN::Util::agauss)
		.def("get_min", (int (*)(int, int))&EMAN::Util::get_min)
		.def("get_min", (int (*)(int, int, int))&EMAN::Util::get_min)
		.def("get_min", (float (*)(float, float))&EMAN::Util::get_min)
		.def("get_min", (float (*)(float, float, float))&EMAN::Util::get_min)
		.def("get_min", (float (*)(float, float, float, float))&EMAN::Util::get_min)
		.def("get_max", (float (*)(float, float))&EMAN::Util::get_max)
		.def("get_max", (float (*)(float, float, float))&EMAN::Util::get_max)
		.def("get_max", (float (*)(float, float, float, float))&EMAN::Util::get_max)
		.def("get_stats", &EMAN::Util::get_stats)
		.def("get_stats_cstyle", &EMAN::Util::get_stats_cstyle)
		.def("angle_sub_2pi", &EMAN::Util::angle_sub_2pi)
		.def("angle_sub_pi", &EMAN::Util::angle_sub_pi)
		.def("get_time_label", &EMAN::Util::get_time_label)
		.def("eman_copysign", &EMAN::Util::eman_copysign)
		.def("eman_erfc", &EMAN::Util::eman_erfc)
		.def("twoD_fine_ali", &EMAN::Util::twoD_fine_ali)
		.def("twoD_fine_ali_G", &EMAN::Util::twoD_fine_ali_G)
		.def("twoD_fine_ali_SD", &EMAN::Util::twoD_fine_ali_SD)
		.def("twoD_fine_ali_SD_G", &EMAN::Util::twoD_fine_ali_SD_G)
		.def("twoD_to_3D_ali", &EMAN::Util::twoD_to_3D_ali)
		.def("get_biggest_cluster", &EMAN::Util::get_biggest_cluster, return_value_policy< manage_new_object >())
		.def("get_slice", &EMAN::Util::get_slice, return_value_policy< manage_new_object >())
		.def("merge_peaks", &EMAN::Util::merge_peaks )
		.def("point_is_in_triangle_2d", &EMAN::Util::point_is_in_triangle_2d )
		.def("point_is_in_convex_polygon_2d", &EMAN::Util::point_is_in_convex_polygon_2d )
		.def("sstevd", &pysstevd )
		.def("snrm2",  &pysnrm2 )
		.def("sgemv",  &pysgemv )
		.def("saxpy",  &pysaxpy )
		.def("sdot",   &pysdot  )
		.def("readarray", &readarray )
		.staticmethod("point_is_in_triangle_2d")
		.staticmethod("point_is_in_convex_polygon_2d")
		.staticmethod("infomask")
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
		.staticmethod("set_line")
                .staticmethod("cml_prepare_line")
		.staticmethod("sstrncmp")
		.staticmethod("int2str")
		.staticmethod("square_sum")
		.staticmethod("prb1d")
		.staticmethod("Polar2D")
		.staticmethod("get_time_label")
		.staticmethod("get_max")
		.staticmethod("mul_img")
		.staticmethod("div_img")
		.staticmethod("div_filter")
		.staticmethod("histogram")
		.staticmethod("Crosrng_ms")
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
		.staticmethod("ExpMinus4YSqr")
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
		.staticmethod("divn_filter")
		.staticmethod("ctf_img")
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
		.staticmethod("move_points")
		.staticmethod("splint")
		.staticmethod("compress_image_mask")
		.staticmethod("calc_best_fft_size")
		.staticmethod("calc_bilinear_least_square")
		.staticmethod("angle_sub_pi")
		.staticmethod("sub_img")
		.staticmethod("reconstitute_image_mask")
		.staticmethod("vrdg")
		.staticmethod("update_fav")
		.staticmethod("change_filename_ext")
		.staticmethod("cyclicshift")
		.staticmethod("sub_fav")
		.staticmethod("Crosrng_ew")
		.staticmethod("get_gauss_rand")
		.staticmethod("mult_scalar")
		.staticmethod("svdcmp")
		.staticmethod("agauss")
		.staticmethod("WTM")
		.staticmethod("subn_img")
		.staticmethod("multiref_polar_ali_2d")
		.staticmethod("multiref_polar_ali_2d_nom")
		.staticmethod("multiref_polar_ali_2d_local")
		.staticmethod("multiref_peaks_ali2d")
		//.staticmethod("multiref_peaks_ali")
		.staticmethod("multiref_peaks_compress_ali2d")
		.staticmethod("ali2d_ccf_list")
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
		.staticmethod("im_diff")
		.staticmethod("Polar2Dmi")
		.staticmethod("Crosrng_e")
		.staticmethod("fast_floor")
		.staticmethod("Frngs")
		.staticmethod("Frngs_inv")
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
		.staticmethod("get_biggest_cluster")
		.staticmethod("merge_peaks")
		.staticmethod("get_slice")
		.staticmethod("sstevd")
		.staticmethod("sgemv")
		.staticmethod("snrm2")
		.staticmethod("saxpy")
		.staticmethod("sdot")
		.staticmethod("readarray")
	);

    scope* EMAN_Util_sincBlackman_scope = new scope(
    class_< EMAN::Util::sincBlackman, EMAN_Util_sincBlackman_Wrapper >("sincBlackman", init< const EMAN::Util::sincBlackman& >())
        .def(init< int, float, optional< int > >())
        .def("sBwin_tab", &EMAN::Util::sincBlackman::sBwin_tab)
        .def("get_sB_size", &EMAN::Util::sincBlackman::get_sB_size)
    );
    delete EMAN_Util_sincBlackman_scope;

    scope* EMAN_Util_KaiserBessel_scope = new scope(
    class_< EMAN::Util::KaiserBessel, EMAN_Util_KaiserBessel_Wrapper >("KaiserBessel", init< const EMAN::Util::KaiserBessel& >())
        .def(init< float, int, float, float, int, optional< float, int > >())
        .def("sinhwin", &EMAN::Util::KaiserBessel::sinhwin, &EMAN_Util_KaiserBessel_Wrapper::default_sinhwin)
        .def("i0win", &EMAN::Util::KaiserBessel::i0win, &EMAN_Util_KaiserBessel_Wrapper::default_i0win)
        .def("I0table_maxerror", &EMAN::Util::KaiserBessel::I0table_maxerror)
        .def("dump_table", &EMAN::Util::KaiserBessel::dump_table)
        .def("i0win_tab", &EMAN::Util::KaiserBessel::i0win_tab)
        .def("get_window_size", &EMAN::Util::KaiserBessel::get_window_size)
        .def("get_kbsinh_win", &EMAN::Util::KaiserBessel::get_kbsinh_win)
        .def("get_kbi0_win", &EMAN::Util::KaiserBessel::get_kbi0_win)
    );

    class_< EMAN::Util::KaiserBessel::kbsinh_win >("kbsinh_win", init< const EMAN::Util::KaiserBessel::kbsinh_win& >())
        .def(init< EMAN::Util::KaiserBessel& >())
        .def("get_window_size", &EMAN::Util::KaiserBessel::kbsinh_win::get_window_size)
        .def("__call__", &EMAN::Util::KaiserBessel::kbsinh_win::operator ())
    ;


    class_< EMAN::Util::KaiserBessel::kbi0_win >("kbi0_win", init< const EMAN::Util::KaiserBessel::kbi0_win& >())
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


    class_< EMAN::Util::Gaussian >("Gaussian", init< const EMAN::Util::Gaussian& >())
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
        .def("vertical_acf", &EMAN::EMUtil::vertical_acf, return_value_policy< manage_new_object >())
        .def("make_image_median", &EMAN::EMUtil::make_image_median, return_value_policy< manage_new_object >())
        .def("get_image_ext_type", &EMAN::EMUtil::get_image_ext_type)
        .def("get_image_type", &EMAN::EMUtil::get_image_type)
        .def("get_image_count", &EMAN::EMUtil::get_image_count)
        .def("get_imageio", &EMAN::EMUtil::get_imageio, return_internal_reference< 1 >(), EMAN_EMUtil_get_imageio_overloads_2_3())
        .def("get_imagetype_name", &EMAN::EMUtil::get_imagetype_name)
        .def("get_datatype_string", &EMAN::EMUtil::get_datatype_string)
        .def("process_ascii_region_io", &EMAN::EMUtil::process_ascii_region_io)
        .def("dump_dict", &EMAN::EMUtil::dump_dict)
        .def("is_same_size", &EMAN::EMUtil::is_same_size)
        .def("is_same_ctf", &EMAN::EMUtil::is_same_ctf)
        .def("is_complex_type", &EMAN::EMUtil::is_complex_type)
        .def("is_valid_filename", &EMAN::EMUtil::is_valid_filename)
        .def("jump_lines", &EMAN::EMUtil::jump_lines)
        .def("get_euler_names", &EMAN::EMUtil::get_euler_names)
        .def("get_all_attributes", &EMAN::EMUtil::get_all_attributes)
		.def("cuda_available", &EMAN::EMUtil::cuda_available)
		.staticmethod("cuda_available")
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

    class_< EMAN::TestUtil >("TestUtil", init<  >())
        .def(init< const EMAN::TestUtil& >())
        .def_readonly("EMDATA_HEADER_EXT", &EMAN::TestUtil::EMDATA_HEADER_EXT)
        .def_readonly("EMDATA_DATA_EXT", &EMAN::TestUtil::EMDATA_DATA_EXT)
        .def("get_debug_int", &EMAN::TestUtil::get_debug_int)
        .def("get_debug_float", &EMAN::TestUtil::get_debug_float)
        .def("get_debug_string", &EMAN::TestUtil::get_debug_string)
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
        .def("emobject_farray_to_py", &EMAN::TestUtil::emobject_farray_to_py)
        .def("emobject_strarray_to_py", &EMAN::TestUtil::emobject_strarray_to_py)
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
        .staticmethod("test_FloatPoint")
        .staticmethod("test_dict")
        .staticmethod("test_map_emobject")
        .staticmethod("test_vector_pixel")
        .staticmethod("get_debug_float")
        .staticmethod("get_debug_image")
        .staticmethod("check_image")
    ;
    ;
}

