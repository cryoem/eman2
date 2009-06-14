/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

#include "reconstructor.h"
#include "ctf.h"
#include "emassert.h"
#include "symmetry.h"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <boost/bind.hpp>
#include <boost/format.hpp>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>

using namespace EMAN;
using std::complex;


#include <iostream>
using std::cerr;
using std::endl;
using std::cout; // for debug

#include <algorithm>
// find, for_each

#include <iomanip>
using std::setprecision;

template < typename T > void checked_delete( T*& x )
{
    typedef char type_must_be_complete[ sizeof(T)? 1: -1 ];
    (void) sizeof(type_must_be_complete);
    delete x;
    x = NULL;
}


template <> Factory < Reconstructor >::Factory()
{
	force_add(&FourierReconstructor::NEW);
	force_add(&FourierReconstructorSimple2D::NEW);
	force_add(&BaldwinWoolfordReconstructor::NEW);
	force_add(&WienerFourierReconstructor::NEW);
	force_add(&BackProjectionReconstructor::NEW);
	force_add(&nn4Reconstructor::NEW);
	force_add(&nnSSNR_Reconstructor::NEW);
	force_add(&nn4_ctfReconstructor::NEW);
	force_add(&nnSSNR_ctfReconstructor::NEW);
}

void FourierReconstructorSimple2D::setup()
{
	nx = params.set_default("nx",0);

	if ( nx < 0 ) throw InvalidValueException(nx, "nx must be positive");

	bool is_fftodd = (nx % 2 == 1);

	ny = nx;
	nx += 2-is_fftodd;

	image = new EMData();
	image->set_size(nx, ny);
	image->set_complex(true);
	image->set_fftodd(is_fftodd);
	image->set_ri(true);

	tmp_data = new EMData();
	tmp_data->set_size(nx/2, nx);
}

int FourierReconstructorSimple2D::insert_slice(const EMData* const slice, const Transform & euler)
{

	// Are these exceptions really necessary? (d.woolford)
	if (!slice) throw NullPointerException("EMData pointer (input image) is NULL");

	if ( slice->get_ndim() != 1 ) throw ImageDimensionException("Image dimension must be 1");

	// I could also test to make sure the image is the right dimensions...
	if (slice->is_complex()) throw ImageFormatException("The image is complex, expecting real");

	EMData* working_slice = slice->process("xform.phaseorigin.tocorner");

	// Fourier transform the slice
	working_slice->do_fft_inplace();

	float* rdata = image->get_data();
	float* norm = tmp_data->get_data();
	float* dat = working_slice->get_data();

	float g[4];
	int offset[4];
	float dt[2];
	offset[0] = 0; offset[1] = 2; offset[2] = nx; offset[3] = nx+2;

	float alt = -((float)(euler.get_rotation("2d"))["alpha"])*M_PI/180.0f;
	for (int x = 0; x < working_slice->get_xsize() / 2; x++) {

		float rx = (float) x;

		float xx = rx*cos(alt);
		float yy = rx*sin(alt);
		float cc = 1.0;

		if (xx < 0) {
			xx = -xx;
			yy = -yy;
			cc = -1.0;
		}

		yy += ny / 2;


		dt[0] = dat[2*x];
		dt[1] = cc * dat[2*x+1];

		// PHIL IS INTERESTED FROM HERE DOWN
		int x0 = (int) floor(xx);
		int y0 = (int) floor(yy);

		int i = 2*x0 + y0*nx;

		float dx = xx - x0;
		float dy = yy - y0;

		g[0] = Util::agauss(1, dx, dy, 0, EMConsts::I2G);
		g[1] = Util::agauss(1, 1 - dx, dy, 0, EMConsts::I2G);
		g[2] = Util::agauss(1, dx, 1 - dy, 0, EMConsts::I2G);
		g[3] = Util::agauss(1, 1 - dx, 1 - dy, 0, EMConsts::I2G);

		// At the extreme we can only do some much...
		if ( x0 == nx-2 ) {
			int k = i + offset[0];
			rdata[k] += g[0] * dt[0];
			rdata[k + 1] += g[0] * dt[1];
			norm[k/2] += g[0];

			k = i + offset[2];
			rdata[k] += g[2] * dt[0];
			rdata[k + 1] += g[2] * dt[1];
			norm[k/2] += g[2];
			continue;

		}
		// capture and accommodate for periodic boundary conditions in the x direction
		if ( x0 > nx-2 ) {
			int dif = x0 - (nx-2);
			x0 -= dif;
		}
		// At the extreme we can only do some much...
		if ( y0 == ny -1 ) {
			int k = i + offset[0];
			rdata[k] += g[0] * dt[0];
			rdata[k + 1] += g[0] * dt[1];
			norm[k/2] += g[0];

			k = i + offset[1];
			rdata[k] += g[1] * dt[0];
			rdata[k + 1] += g[1] * dt[1];
			norm[k/2] += g[1];
			continue;
		}
		// capture and accommodate for periodic boundary conditions in the y direction
		if ( y0 > ny-1) {
			int dif = y0 - (ny-1);
			y0 -= dif;
		}

		if (x0 >= nx - 2 || y0 >= ny - 1) continue;




		for (int j = 0; j < 4; j++)
		{
			int k = i + offset[j];
			rdata[k] += g[j] * dt[0];
			rdata[k + 1] += g[j] * dt[1];
			norm[k/2] += g[j];

		}
	}

	return 0;

}

EMData *FourierReconstructorSimple2D::finish()
{
	normalize_threed();

	image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
	image->process_inplace("xform.phaseorigin.tocenter");

	return  	image;
}

void ReconstructorVolumeData::normalize_threed()
{
	float* norm = tmp_data->get_data();
	float* rdata = image->get_data();

	// FIXME should throw a sensible error
	if ( 0 == norm ) throw NullPointerException("The normalization volume was null!");
	if ( 0 == rdata ) throw NullPointerException("The complex reconstruction volume was null!");

	for (int i = 0; i < nx * ny * nz; i += 2) {
		float d = norm[i/2];
		if (d == 0) {
			rdata[i] = 0;
			rdata[i + 1] = 0;
		}
		else {
			rdata[i] /= d;
			rdata[i + 1] /= d;
		}
	}
}

void FourierReconstructor::free_memory()
{
	if ( inserter != 0 )
	{
		delete inserter;
		inserter = 0;
	}
	if ( interp_FRC_calculator != 0 )
	{
		delete interp_FRC_calculator;
		interp_FRC_calculator = 0;
	}
}

#include <sstream>

void FourierReconstructor::load_inserter()
{
	// ints get converted to strings byt the Dict object here

	stringstream ss;
	ss << (int)params["mode"];
	string mode;
	ss >> mode;
//	ss
//	string mode = (string)params["mode"];
	Dict parms;
	parms["rdata"] = image->get_data();
	parms["norm"] = tmp_data->get_data();
	parms["nx"] = nx;
	parms["ny"] = ny;
	parms["nz"] = nz;

	if ( inserter != 0 )
	{
		delete inserter;
	}

	inserter = Factory<FourierPixelInserter3D>::get(mode, parms);
	inserter->init();
}



void FourierReconstructor::load_interp_FRC_calculator()
{
	Dict init_parms;
	init_parms["rdata"] = image->get_data();
	init_parms["norm"] = tmp_data->get_data();
	init_parms["nx"] = nx;
	init_parms["ny"] = ny;
	init_parms["nz"] = nz;

	if ( x_scale_factor != 0 ) init_parms["x_scale"] = 1.0/x_scale_factor;
	if ( y_scale_factor != 0 ) init_parms["y_scale"] = 1.0/y_scale_factor;
	if ( z_scale_factor != 0 ) init_parms["z_scale"] = 1.0/z_scale_factor;

	if ( interp_FRC_calculator != 0 )
	{
		delete interp_FRC_calculator;
	}

	stringstream ss;
	ss << (int)params["mode"];
	string mode;
	ss >> mode;

	interp_FRC_calculator = Factory<InterpolatedFRC>::get(mode, init_parms);
}

void FourierReconstructor::setup()
{
	// default setting behavior - does not override if the parameter is already set
	params.set_default("mode",2);

	int x_size = params["x_in"];
	if ( x_size < 0 ) throw InvalidValueException(x_size, "x size of images must be greater than 0");
	int y_size = params["y_in"];
	if ( y_size < 0 ) throw InvalidValueException(y_size, "y size of images must be greater than 0");

	if ( x_size > y_size ) max_input_dim = x_size;
	else max_input_dim = y_size;

	// This is a helical adaptation - FIXME explain
	bool helical_special_behavior = false;
	if ( helical_special_behavior )
	{
		if ( x_size > y_size )
		{
			if ( (int) params["xsample"] == 0 ) params["xsample"] = y_size;
			if ( (int) params["zsample"] == 0 ) params["zsample"] = x_size;
		}
		if ( y_size > x_size )
		{
			if ( (int) params["ysample"] == 0 ) params["ysample"] = x_size;
			if ( (int) params["zsample"] == 0 ) params["zsample"] = y_size;
		}
	}

	if ( (int) params["xsample"] != 0  ) output_x = (int) params["xsample"];
	else output_x = x_size;

	if ( (int) params["ysample"] != 0  ) output_y = (int) params["ysample"];
	else output_y = y_size;

	if ( (int) params["zsample"] != 0  ) output_z = (int) params["zsample"];
	else
	{
		if ( x_size == y_size ) output_z = x_size;
		else if ( x_size > y_size ) output_z = x_size;
		else output_z = y_size;
	}

	nx = output_x;
	ny = output_y;
	nz = output_z;

	// Adjust nx if for Fourier transform even odd issues
	bool is_fftodd = (nx % 2 == 1);
	// The Fourier transform requires one extra pixel in the x direction,
	// which is two spaces in memory, one each for the complex and the
	// real components
	nx += 2-is_fftodd;

	if ( (int) params["zsample"] == 0  ) nz = max_input_dim;

	if ( nz < max_input_dim ) z_scale_factor = (float) nz / (float) max_input_dim;
	if ( ny < max_input_dim ) y_scale_factor = (float) ny / (float) max_input_dim;
	if ( nx < max_input_dim ) x_scale_factor = (float) nx / (float) max_input_dim;


	// Odd dimension support is here atm, but not above.
	image = new EMData();
	image->set_size(nx, ny, nz);
	image->set_complex(true);
	image->set_fftodd(is_fftodd);
	image->set_ri(true);

	tmp_data = new EMData();
	tmp_data->set_size(nx/2, ny, nz);
	tmp_data->to_zero();
	tmp_data->update();

	load_inserter();
	load_interp_FRC_calculator();

	if ( (bool) params["quiet"] == false )
	{
		cout << "3D Fourier dimensions are " << nx << " " << ny << " " << nz << endl;
		cout << "You will require approximately " << setprecision(3) << (nx*ny*nz*sizeof(float)*1.5)/1000000000.0 << "GB of memory to reconstruct this volume" << endl;
		cout << "Scale factors are " << x_scale_factor << " " << y_scale_factor << " " << z_scale_factor << endl;
		cout << "Max input dim is " << max_input_dim << endl;
	}
}

EMData* FourierReconstructor::preprocess_slice( const EMData* const slice,  const Transform& t )
{
	// Shift the image pixels so the real space origin is now located at the phase origin (at the bottom left of the image)
	EMData* return_slice = 0;
	Transform tmp(t);
	tmp.set_rotation(Dict("type","eman")); // resets the rotation to 0 implicitly

	Vec2f trans = tmp.get_trans_2d();
	float scale = tmp.get_scale();
	bool mirror = tmp.get_mirror();
	if (trans[0] != 0 || trans[1] != 0 || scale != 1.0 ) {
		return_slice = slice->process("math.transform",Dict("transform",&tmp));
	} else if ( mirror == true ) {
		return_slice = slice->process("xform.flip",Dict("axis","x"));
	}

	if (return_slice == 0) {
		return_slice = slice->process("xform.phaseorigin.tocorner");
	} else {
		return_slice->process_inplace("xform.phaseorigin.tocorner");
	}

// 	EMData* return_slice = slice->process("normalize.edgemean");
// 	return_slice->process_inplace("xform.phaseorigin.tocorner");

	// Fourier transform the slice
	return_slice->do_fft_inplace();

	return_slice->mult((float)sqrt(1.0f/(return_slice->get_ysize())*return_slice->get_xsize()));

	// Shift the Fourier transform so that it's origin is in the center (bottom) of the image.
	return_slice->process_inplace("xform.fourierorigin.tocenter");

	return return_slice;
}

void FourierReconstructor::zero_memory()
{
	if (tmp_data != 0 ) tmp_data->to_zero();
	if (image != 0 ) image->to_zero();

	EMData* start_model = params.set_default("start_model",(EMData*)0);

	if (start_model != 0) {
		if (!start_model->is_complex()) {
			start_model->do_fft_inplace();
			start_model->process_inplace("xform.phaseorigin.tocenter");
		}

		if (start_model->get_xsize() != nx || start_model->get_ysize() != ny || start_model->get_zsize() != nz ) {
			throw ImageDimensionException("The dimensions of the start_model are incorrect");
		}

		if (!start_model->is_shuffled()) {
			start_model->process_inplace("xform.fourierorigin.tocenter");
		}

		memcpy(image->get_data(),start_model->get_data(),nx*ny*nz*sizeof(float));

		float start_model_weight =  params.set_default("start_model_weight",1.0f);

		if (start_model_weight != 1.0) {
			image->mult(start_model_weight);
		}

		tmp_data->add(start_model_weight);
	}
}

int FourierReconstructor::insert_slice(const EMData* const input_slice, const Transform & arg)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) throw NullPointerException("EMData pointer (input image) is NULL");

	// I could also test to make sure the image is the right dimensions...
	if (input_slice->is_complex()) throw ImageFormatException("The image is complex, expecting real");

	// Get the proprecessed slice - there are some things that always happen to a slice,
	// such as as Fourier conversion and optional padding etc.
	//
	// We must use only the rotational component of the transform, scaling, translation and mirroring
	// are not implemented in Fourier space
	Transform rotation;
	if ( input_slice->has_attr("xform.projection") ) {
		rotation = *((Transform*) input_slice->get_attr("xform.projection")); // assignment operator
	} else {
		rotation = arg; // assignment operator
	}

	EMData* slice = preprocess_slice( input_slice, rotation);

	// Catch the first time a slice is inserted to do some things....
	if ( slice_insertion_flag == true )
	{
		// Zero the memory in the associated EMData objects
		zero_memory();
		// Reset the image_idx to zero
		image_idx = 0;
		// Set flags appropriately
		slice_agreement_flag = true;
		slice_insertion_flag = false;

		// should be a if verbose here
		//if ( prev_quality_scores.size() != 0 ) print_stats( prev_quality_scores );
	}

	// quality_scores.size() is zero on the first run, so this enforcement of slice quality does not take
	// place until after the first round of slice insertion
	if ( quality_scores.size() != 0 )
	{
		if ( quality_scores[image_idx].get_snr_normed_frc_integral() < (float) params["hard"] )
		{
			image_idx++;
			return 1;
		}
		slice->mult(1.f/quality_scores[image_idx].get_norm());
// 		cout << "Norm multiplied by " << 1.f/quality_scores[image_idx].get_norm() << endl;
		image_idx++;
	}

	// Finally to the pixel wise slice insertion
	rotation.set_scale(1.0);
	rotation.set_mirror(false);
	rotation.set_trans(0,0,0);
	do_insert_slice_work(slice, rotation);

	delete slice;

// 	image->update();
	return 0;
}

void FourierReconstructor::do_insert_slice_work(const EMData* const input_slice, const Transform & arg)
{
	// Reload the inserter if the mode has changed
	string mode = (string) params["mode"];
	if ( mode != inserter->get_name() )	load_inserter();

	int rl = Util::square( max_input_dim / 2);

	float dt[2];

	float y_scale = 0, x_scale = 0;

	int y_in = input_slice->get_ysize();
	int x_in = input_slice->get_xsize();
	// Adjust the dimensions to account for odd and even ffts
	if (input_slice->is_fftodd()) x_in -= 1;
	else x_in -= 2;

	if ( y_in != x_in )
	{
		if ( x_in > y_in ) y_scale = (float) x_in / (float) y_in;
		else x_scale = (float) y_in / (float) x_in;
	}

	float *dat = input_slice->get_data();
	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);
	float weight = params.set_default("weight",1.0f);

	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform t3d = arg*(*it);
		for (int y = 0; y < input_slice->get_ysize(); y++) {
			for (int x = 0; x < input_slice->get_xsize() / 2; x++) {

				float rx = (float) x;
				float ry = (float) y;

				if ( y_in != x_in )
				{
					if ( x_in > y_in ) ry *= y_scale;
					else rx *= x_scale;
				}

				if ((rx * rx + Util::square(ry - max_input_dim / 2)) > rl)
					continue;

				Vec3f coord(rx,(ry - max_input_dim / 2),0);
				coord = coord*t3d; // transpose multiplication
				float xx = coord[0];
				float yy = coord[1];
				float zz = coord[2];

				float cc = 1;

				if (xx < 0) {
					xx = -xx;
					yy = -yy;
					zz = -zz;
					cc = -1.0;
				}

				if ( z_scale_factor != 0 ) zz *= z_scale_factor;
				if ( y_scale_factor != 0 ) yy *= y_scale_factor;
				if ( x_scale_factor != 0 ) xx *= x_scale_factor;

				yy += ny / 2;
				zz += nz / 2;

				int idx = x * 2 + y * (input_slice->get_xsize());
				dt[0] = dat[idx];
				dt[1] = cc * dat[idx+1];

				inserter->insert_pixel(xx,yy,zz,dt,weight);
			}
		}
	}
}

int FourierReconstructor::determine_slice_agreement(const EMData* const input_slice, const Transform & t3d, const unsigned int num_particles_in_slice)
{
	// Are these exceptions really necessary? (d.woolford)
	if (!input_slice) {
		LOGERR("Insertion of NULL slice in FourierReconstructor::insert_slice");
		throw NullPointerException("EMData pointer (input image) is NULL");
	}

	// I could also test to make sure the image is the right dimensions...
	if (input_slice->is_complex()) {
		LOGERR("Do not Fourier transform the image before it is passed to determine_slice_agreement in FourierReconstructor, this is performed internally");
		throw ImageFormatException("The image is complex, expecting real");
	}

	// The first time determine_slice_agreement is called (in each iteration) some things need to happen...
	if ( slice_agreement_flag == true )
	{
		// Normalize the real data using the normalization values in the normalization volume
		// This is important because the InterpolatedFRC objects assumes this behaviour
		normalize_threed();
		// Reset the image_idx to index into prev_quality_scores correctly
		image_idx = 0;
		// Make a copy of the previous quality scores - the first time around this will be like copying nothing
		prev_quality_scores = quality_scores;
		// Now clear the old scores so they can be redetermined
		quality_scores.clear();
		// Reset the flags
		slice_insertion_flag = true;
		slice_agreement_flag = false;
	}

	// Get the proprecessed slice - there are some things that always happen to a slice,
	// such as as Fourier conversion and optional padding etc.
	Transform rotation;
	if ( input_slice->has_attr("xform.projection") ) {
		rotation = *((Transform*) input_slice->get_attr("xform.projection")); // assignment operator

	} else {
		rotation = t3d; // assignment operator
	}
	EMData* slice = preprocess_slice( input_slice, rotation );

	// quality_scores.size() is zero on the first run, so this enforcement of slice quality does not take
	// place until after the first round of slice insertion
	if (  prev_quality_scores.size() != 0 )
	{
		// The slice must be multiplied by the normalization value used in insert or else the calculations will
		// be inconsistent
		slice->mult(1.f/prev_quality_scores[image_idx].get_norm());
	}

	// Reset zeros the associated memory in the ifrc
	interp_FRC_calculator->reset();

	int rl = Util::square( max_input_dim / 2);

	float dt[2];

	float y_scale = 0, x_scale = 0;

	int y_in = slice->get_ysize();
	int x_in = slice->get_xsize();
	if (input_slice->is_fftodd()) x_in -= 1;
	else x_in -= 2;

	if ( y_in != x_in )
	{
		if ( x_in > y_in ) y_scale = (float) x_in / (float) y_in;
		else x_scale = (float) y_in / (float) x_in;
	}
	float *dat = slice->get_data();

	float weight = params.set_default("weight",1.0f);


	rotation.set_scale(1.0);
	rotation.set_mirror(false);
	rotation.set_trans(0,0,0);

	for (int y = 0; y < slice->get_ysize(); y++) {
		for (int x = 0; x < slice->get_xsize() / 2; x++) {

			float rx = (float) x;
			float ry = (float) y;

			if ( y_in != x_in )
			{
				if ( x_in > y_in ) ry *= y_scale;
				else rx *= x_scale;
			}

			if ((rx * rx + Util::square(ry - max_input_dim / 2)) > rl)
				continue;

// 			float xx = (float) (rx * t3d[0][0] + (ry - max_input_dim / 2) * t3d[1][0]);
// 			float yy = (float) (rx * t3d[0][1] + (ry - max_input_dim / 2) * t3d[1][1]);
// 			float zz = (float) (rx * t3d[0][2] + (ry - max_input_dim / 2) * t3d[1][2]);

			Vec3f coord(rx,(ry - max_input_dim / 2),0);
			coord = coord*rotation; // transpose multiplication
			float xx = coord[0];
			float yy = coord[1];
			float zz = coord[2];

			float cc = 1;

			if (xx < 0) {
				xx = -xx;
				yy = -yy;
				zz = -zz;
				cc = -1.0;
			}

			if ( z_scale_factor != 0 ) zz *= z_scale_factor;
			if ( y_scale_factor != 0 ) yy *= y_scale_factor;
			if ( x_scale_factor != 0 ) xx *= x_scale_factor;

			yy += ny / 2;
			zz += nz / 2;

			int idx = x * 2 + y * (slice->get_xsize());
			dt[0] = dat[idx];
			dt[1] = cc * dat[idx+1];


			if ( prev_quality_scores.size() != 0 )
			{
				// If the slice was not inserted into the 3D volume in the previous round of slice insertion
				// then the weight must be set to zero, or else the result produced by the InterpolatedFRC will be wrong.
				// This is because the InterpolatedFRC subtracts the incoming pixels from the 3D volume in order
				// to estimate quality using only data from the other slices in the volume. If the slice was not previously inserted
				// then it should not be subtracted prior to quality estimation....
				if ( prev_quality_scores[image_idx].get_snr_normed_frc_integral() < (float) params["hard"] )
				{
					weight = 0;
				}
			}
			// FIXME: this could be replaced in favor of a class implementation with no switch statement, similar to the
			// inserter in the Fourier reconstructor do_insert_slice_work method

			interp_FRC_calculator->continue_frc_calc(xx, yy, zz, dt, weight);
		}
	}

	InterpolatedFRC::QualityScores q_scores = interp_FRC_calculator->finish( num_particles_in_slice );
	// Print the quality scores here for debug information
	//q_scores.debug_print();

	if (prev_quality_scores.size() != 0 )
	{
		// If a previous normalization value has been calculated then average it with the one that
		// was just calculated, this will cause the normalization values to converge more rapidly.
		q_scores.set_norm((q_scores.get_norm() + prev_quality_scores[image_idx].get_norm())/2.0f);
		image_idx++;
	}

	quality_scores.push_back(q_scores);

	delete slice;
	return 0;
}

void FourierReconstructor::print_stats( const vector<InterpolatedFRC::QualityScores>& scores )
{
	if ( prev_quality_scores.size() == 0 )
	{
		//cout << "No quality scores present in FourierReconstructor::print_stats, nothing to print" << endl;
		return;
	}

	unsigned int size = scores.size();

	unsigned int contributing_images = 0;
	for( unsigned int i = 0; i < size; ++ i )
	{
		if (scores[i].get_snr_normed_frc_integral() < (float) params["hard"])
			continue;
		else contributing_images++;
	}

	double* norm_frc = new double[contributing_images];
	double* frc = new double[contributing_images];
	double* norm_snr = new double[contributing_images];

	unsigned int idx = 0;
	for( unsigned int i = 0; i < size; ++ i )
	{
		if (scores[i].get_snr_normed_frc_integral() < (float) params["hard"])
			continue;

		norm_frc[idx] = scores[i].get_snr_normed_frc_integral();
		frc[idx] = scores[i].get_frc_integral();
		norm_snr[idx] = scores[i].get_normed_snr_integral();

		++idx;
	}

	double mean = gsl_stats_mean(norm_frc, 1, contributing_images);
	double variance = gsl_stats_variance_m(norm_frc, 1, contributing_images, mean);

	cout << "Normalized FRC mean " << mean << " std dev " << sqrt(variance) << endl;

	mean = gsl_stats_mean(frc, 1, contributing_images);
	variance = gsl_stats_variance_m(frc, 1, contributing_images, mean);

	cout << "FRC mean " << mean << " std dev " << sqrt(variance) << endl;

	mean = gsl_stats_mean(norm_snr, 1, contributing_images);
	variance = gsl_stats_variance_m(norm_snr, 1, contributing_images, mean);
	cout << "SNR mean " << mean << " std dev " << sqrt(variance) << endl;

	double c0, c1, cov00, cov01, cov11, sumsq;
	gsl_fit_linear (norm_frc, 1, frc, 1, contributing_images, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	cout << "The correlation between frc and norm_frc is " << c0 << " + " << c1 << "x" << endl;

	gsl_fit_linear (norm_frc, 1, norm_snr, 1, contributing_images, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	cout << "The correlation between norm_snr and norm_frc is " << c0 << " + " << c1 << "x" << endl;

	gsl_fit_linear (norm_snr, 1, frc, 1, contributing_images, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	cout << "The correlation between frc and norm_snr is " << c0 << " + " << c1 << "x" << endl;

	delete [] norm_frc;
	delete [] frc;
	delete [] norm_snr;
}

EMData *FourierReconstructor::finish()
{
	float *norm = tmp_data->get_data();
	float *rdata = image->get_data();

	if (params["dlog"]) {
		size_t size = nx*ny*nz;
		for (size_t i = 0; i < size; i += 2) {
			float d = norm[i];
			if (d == 0) {
				rdata[i] = 0;
				rdata[i + 1] = 0;
			}
			else {
				float in = norm[i + 1] / norm[i];
				float cin = Util::square_sum(rdata[i], rdata[i + 1]);
				rdata[i] *= sqrt(in / cin);
				rdata[i + 1] *= sqrt(in / cin);
			}
		}
	}
	else {
		normalize_threed();
	}

// 	tmp_data->write_image("density.mrc");

	// we may as well delete the tmp data now... it saves memory and the calling program might
	// need memory after it gets the return volume.
	// If this delete didn't happen now, it would happen when the deconstructor was called,
	if ( tmp_data != 0 )
	{
		delete tmp_data;
		tmp_data = 0;
	}


	if ( params["3damp"]) {
		EMData* fftimage = image->get_fft_amplitude ();
		fftimage->write_image("threed_fft_amp.mrc");
		delete fftimage;
	}

	//
	image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
	image->process_inplace("xform.phaseorigin.tocenter");

	// If the image was padded it should be the original size, as the client would expect
	//  I blocked the rest, it is almost certainly incorrect  PAP 07/31/08
	// No, it's not incorrect. You are wrong. You have the meaning of nx mixed up. DSAW 09/23/cd
	bool is_fftodd = (nx % 2 == 1);
	if ( (nx-2*(!is_fftodd)) != output_x || ny != output_y || nz != output_z )
	{
		FloatPoint origin( (nx-output_x)/2, (ny-output_y)/2, (nz-output_z)/2 );
		FloatSize region_size( output_x, output_y, output_z);
		Region clip_region( origin, region_size );
		image->clip_inplace( clip_region );
	}

	// Should be an "if (verbose)" here or something
	//print_stats(quality_scores);

	image->update();

	return image;
}


void BaldwinWoolfordReconstructor::setup()
{
	//This is a bit of a hack - but for now it suffices
	params.set_default("mode",1);
	FourierReconstructor::setup();
	// Set up the Baldwin Kernel
	int P = (int)((1.0+0.25)*max_input_dim+1);
	float r = (float)(max_input_dim+1)/(float)P;
	dfreq = 0.2f;
	if (W != 0) delete [] W;
	int maskwidth = params.set_default("maskwidth",2);
	W = Util::getBaldwinGridWeights(maskwidth, (float)P, r,dfreq,0.5f,0.2f);
}

EMData* BaldwinWoolfordReconstructor::finish()
{
	tmp_data->write_image("density.mrc");
	image->process_inplace("xform.fourierorigin.tocorner");
	image->do_ift_inplace();
	image->depad();
	image->process_inplace("xform.phaseorigin.tocenter");

	if ( (bool) params.set_default("postmultiply", false) )
	{
		cout << "POST MULTIPLYING" << endl;
	// now undo the Fourier convolution with real space division
		float* d = image->get_data();
		float N = (float) image->get_xsize()/2.0f;
		N *= N;
		size_t rnx = image->get_xsize();
		size_t rny = image->get_ysize();
		size_t rnxy = rnx*rny;
		int cx = image->get_xsize()/2;
		int cy = image->get_ysize()/2;
		int cz = image->get_zsize()/2;
		size_t idx;
		for (int k = 0; k < image->get_zsize(); ++k ){
			for (int j = 0; j < image->get_ysize(); ++j ) {
				for (int i =0; i < image->get_xsize(); ++i ) {
					float xd = (float)(i-cx); xd *= xd;
					float yd = (float)(j-cy); yd *= yd;
					float zd = (float)(k-cz); zd *= zd;
					float weight = exp((xd+yd+zd)/N);
					idx = k*rnxy + j*rnx + i;
					d[idx] *=  weight;
				}
			}
		}
	}
	image->update();
	return  image;
}

#include <iomanip>

int BaldwinWoolfordReconstructor::insert_slice_weights(const Transform& t3d)
{
	bool fftodd = image->is_fftodd();
	int rnx = nx-2*!fftodd;

	float y_scale = 1.0, x_scale = 1.0;

	if ( ny != rnx  )
	{
		if ( rnx > ny ) y_scale = (float) rnx / (float) ny;
		else x_scale = (float) ny / (float) rnx;
	}

	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();

	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);
	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform n3d = t3d*(*it);

		for (int y = 0; y < tny; y++) {
			for (int x = 0; x < tnx; x++) {

				float rx = (float) x;
				float ry = (float) y;

				if ( ny != rnx )
				{
					if ( rnx > ny ) ry *= y_scale;
					else rx *= x_scale;
				}
// 				float xx = rx * n3d[0][0] + (ry - tny/2) * n3d[1][0];
// 				float yy = rx * n3d[0][1] + (ry - tny/2) * n3d[1][1];
// 				float zz = rx * n3d[0][2] + (ry - tny/2) * n3d[1][2];

				Vec3f coord(rx,(ry - tny/2),0);
				coord = coord*n3d; // transpose multiplication
				float xx = coord[0];
				float yy = coord[1];
				float zz = coord[2];

				if (xx < 0 ){
					xx = -xx;
					yy = -yy;
					zz = -zz;
				}

				yy += tny/2;
				zz += tnz/2;
				insert_density_at(xx,yy,zz);
			}
		}
	}

	return 0;
}

void BaldwinWoolfordReconstructor::insert_density_at(const float& x, const float& y, const float& z)
{
	int xl = Util::fast_floor(x);
	int yl = Util::fast_floor(y);
	int zl = Util::fast_floor(z);

	// w is the windowing width
	int w = params.set_default("maskwidth",2);
	float wsquared = (float) w*w;
	float dw = 1.0f/w;
	dw *= dw;

	// w minus one - this control the number of
	// pixels/voxels to the left of the main pixel
	// that will have density
	int wmox = w-1;
	int wmoy = w-1;
	int wmoz = w-1;

	// If any coordinate is incedental with a vertex, then
	// make sure there is symmetry in density accruing.
	// i.e. the window width must be equal in both directions
	if ( ((float) xl) == x ) wmox = w;
	if ( ((float) yl) == y ) wmoy = w;
	if ( ((float) zl) == z ) wmoz = w;

	float* d = tmp_data->get_data();
	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();
	size_t tnxy = tnx*tny;

	int mode = params.set_default("mode",1);

	for(int k = zl-wmoz; k <= zl+w; ++k ) {
		for(int j = yl-wmoy; j <= yl+w; ++j) {
			for( int i = xl-wmox; i <= xl+w; ++i) {
				float fac = 1.0;
				int ic = i, jc = j, kc = k;

				// Fourier space is periodic, which is enforced
				// by the next 6 if statements. These if statements
				// assume that the Fourier DC components is at
				// (0,ny/2,nz/2).
				if ( i <= 0 ) {

					if ( x != 0 && i == 0 ) fac = 1.0;
					else if ( x == 0 && i < 0) continue;
// 					if (i < 0 ) ic = -i;
					if (i < 0 ) {
						continue;
						ic = -i;
						jc = tny-jc;
						kc = tnz-kc;
					}
				}
				if ( ic >= tnx ) ic = 2*tnx-ic-1;
				if ( jc < 0 ) jc = tny+jc;
				if ( jc >= tny ) jc = jc-tny;
				if ( kc < 0 ) kc = tnz+kc;
				if ( kc >= tnz ) kc = kc-tnz;
// 				if ( ic >= tnx ) continue;
// 				if ( jc < 0 ) continue;
// 				if ( jc >= tny ) continue;
// 				if ( kc < 0 ) continue;
// 				if ( kc >= tnz ) continue;
				// This shouldn't happen
				// Debug remove later
				if ( ic < 0 ) { cout << "wo 1" << endl; }
				if ( ic >= tnx  ){ cout << "wo 2" << endl; }
				if ( jc < 0 ) { cout << "wo 3" << endl; }
				if ( jc >= tny ) { cout << "wo 4" << endl; }
				if ( kc < 0 ) { cout << "wo 5" << endl; }
				if ( kc >= tnz ) { cout << "wo 6" << endl; }


				float zd = (z-(float)k);
				float yd = (y-(float)j);
				float xd = (x-(float)i);
				zd *= zd; yd *= yd; xd *= xd;
				float distsquared = xd+yd+zd;
				// We enforce a spherical kernel
				if ( mode == 1 && distsquared > wsquared ) continue;

// 				float f = fac*exp(-dw*(distsquared));
				float f = fac*exp(-2.467f*(distsquared));
				// Debug - this error should never occur.
				if ( (kc*tnxy+jc*tnx+ic) >= tnxy*tnz ) throw OutofRangeException(0,tnxy*tnz,kc*tnxy+jc*tnx+ic, "in density insertion" );
				d[kc*tnxy+jc*tnx+ic] += f;
			}
		}
	}
}

int BaldwinWoolfordReconstructor::insert_slice(const EMData* const input_slice, const Transform & t)
{
	Transform rotation;
	if ( input_slice->has_attr("xform.projection") ) {
		rotation = *((Transform*) input_slice->get_attr("xform.projection")); // assignment operator
	} else {
		rotation = t; // assignment operator
	}
	Transform tmp(rotation);
	tmp.set_rotation(Dict("type","eman")); // resets the rotation to 0 implicitly

	Vec2f trans = tmp.get_trans_2d();
	float scale = tmp.get_scale();
	bool mirror = tmp.get_mirror();
	EMData* slice = 0;
	if (trans[0] != 0 || trans[1] != 0 || scale != 1.0 ) {
		slice = input_slice->process("math.transform",Dict("transform",&tmp));
	} else if ( mirror == true ) {
		slice = input_slice->process("xform.flip",Dict("axis","x"));
	}
	if ( slice == 0 ) {
		slice = input_slice->process("xform.phaseorigin.tocorner");
	} else {
		slice->process_inplace("xform.phaseorigin.tocorner");
	}

	slice->do_fft_inplace();
	slice->process_inplace("xform.fourierorigin.tocenter");
	float *dat = slice->get_data();
	float dt[2];

	bool fftodd = image->is_fftodd();
	int rnx = nx-2*!fftodd;

	float y_scale = 1.0, x_scale = 1.0;

	if ( ny != rnx  )
	{
		if ( rnx > ny ) y_scale = (float) rnx / (float) ny;
		else x_scale = (float) ny / (float) rnx;
	}

	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();

	vector<Transform> syms = Symmetry3D::get_symmetries((string)params["sym"]);
// 	float weight = params.set_default("weight",1.0f);

	rotation.set_scale(1.0); rotation.set_mirror(false); rotation.set_trans(0,0,0);
	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {
		Transform t3d = rotation*(*it);

		for (int y = 0; y < tny; y++) {
			for (int x = 0; x < tnx; x++) {
				float rx = (float) x;
				float ry = (float) y;

				if ( ny != rnx )
				{
					if ( rnx > ny ) ry *= y_scale;
					else rx *= x_scale;
				}

// 				float xx = rx * n3d[0][0] + (ry - tny/2) * n3d[1][0];
// 				float yy = rx * n3d[0][1] + (ry - tny/2) * n3d[1][1];
// 				float zz = rx * n3d[0][2] + (ry - tny/2) * n3d[1][2];

				Vec3f coord(rx,(ry - tny/2),0);
				coord = coord*t3d; // transpose multiplication
				float xx = coord[0];
				float yy = coord[1];
				float zz = coord[2];


				float cc = 1;
				if (xx < 0 ){
					xx = -xx;
					yy = -yy;
					zz = -zz;
					cc = -1;
				}

				yy += tny/2;
				zz += tnz/2;

				int idx = x * 2 + y * (slice->get_xsize());
				dt[0] = dat[idx];
				dt[1] = cc * dat[idx+1];

				insert_pixel(xx,yy,zz,dt);
			}
		}
	}

	delete slice;

	return 0;
}

void BaldwinWoolfordReconstructor::insert_pixel(const float& x, const float& y, const float& z, const float dt[2])
{
	int xl = Util::fast_floor(x);
	int yl = Util::fast_floor(y);
	int zl = Util::fast_floor(z);

	// w is the windowing width
	int w = params.set_default("maskwidth",2);
	float wsquared = (float) w*w;
	float dw = 1.0f/w;
	dw *= dw;

	int wmox = w-1;
	int wmoy = w-1;
	int wmoz = w-1;

	// If any coordinate is incedental with a vertex, then
	// make sure there is symmetry in density accruing.
	// i.e. the window width must be equal in both directions
	if ( ((float) xl) == x ) wmox = w;
	if ( ((float) yl) == y ) wmoy = w;
	if ( ((float) zl) == z ) wmoz = w;

	float* we = tmp_data->get_data();
	int tnx = tmp_data->get_xsize();
	int tny = tmp_data->get_ysize();
	int tnz = tmp_data->get_zsize();
	int tnxy = tnx*tny;

	int rnx = 2*tnx;
	int rnxy = 2*tnxy;

	int mode = params.set_default("mode",1);

	float* d = image->get_data();
	for(int k = zl-wmoz; k <= zl+w; ++k ) {
		for(int j = yl-wmoy; j <= yl+w; ++j) {
			for( int i = xl-wmox; i <= xl+w; ++i) {
				float fac = 1.0;
				int ic = i, jc = j, kc = k;

				// Fourier space is periodic, which is enforced
				// by the next 6 if statements. These if statements
				// assume that the Fourier DC component is at
				// (0,ny/2,nz/2).
				float negfac=1.0;
				if ( i <= 0 ) {
					if ( x != 0 && i == 0 ) fac = 1.0;
					else if ( x == 0 && i < 0) continue;
					if (i < 0 ) {
						continue;
						ic = -i;
						jc = tny-jc;
						kc = tnz-kc;
						negfac=-1.0;
					}
				}
				if ( ic >= tnx ) ic = 2*tnx-ic-1;
				if ( jc < 0 ) jc = tny+jc;
				if ( jc >= tny ) jc = jc-tny;
				if ( kc < 0 ) kc = tnz+kc;
				if ( kc >= tnz ) kc = kc-tnz;
// 				if ( ic >= tnx ) continue;
// 				if ( jc < 0 ) continue;
// 				if ( jc >= tny ) continue;
// 				if ( kc < 0 ) continue;
// 				if ( kc >= tnz ) continue;

				float zd = (z-(float)k);
				float yd = (y-(float)j);
				float xd = (x-(float)i);
				zd *= zd; yd *= yd; xd *= xd;
				float distsquared = xd+yd+zd;
// 				float f = fac*exp(-dw*(distsquared));
				float f = fac*exp(-2.467f*(distsquared));
				float weight = f/we[kc*tnxy+jc*tnx+ic];
				// debug - this error should never occur
				if ( (kc*rnxy+jc*rnx+2*ic+1) >= rnxy*tnz ) throw OutofRangeException(0,rnxy*tnz,kc*rnxy+jc*rnx+2*ic+1, "in pixel insertion" );
				size_t k = kc*rnxy+jc*rnx+2*ic;

				float factor, dist,residual;
				int sizeW,sizeWmid,idx;
				switch (mode) {
					case 0:
						d[k] += weight*f*dt[0];
						d[k+1] += negfac*weight*f*dt[1];
						cout << "hello" << endl;
					break;

					case 1:
						// We enforce a spherical kernel
						if ( distsquared > wsquared ) continue;

						sizeW = (int)(1+2*w/dfreq);
						sizeWmid = sizeW/2;

						dist = sqrtf(distsquared);
						idx = (int)(sizeWmid + dist/dfreq);
						if (idx >= sizeW) throw InvalidValueException(idx, "idx was greater than or equal to sizeW");
						residual = dist/dfreq - (int)(dist/dfreq);
						if ( fabs(residual) > 1) throw InvalidValueException(residual, "Residual was too big");

						factor = (W[idx]*(1.0f-residual)+W[idx+1]*residual)*weight;

						d[k] += dt[0]*factor;
						d[k+1] += dt[1]*factor;
					break;

					default:
						throw InvalidValueException(mode, "The mode was unsupported in BaldwinWoolfordReconstructor::insert_pixel");
					break;
				}
			}
		}
	}
}

// void BaldwinWoolfordReconstructor::insert_pixel(const float& x, const float& y, const float& z, const float dt[2])
// {
// 	int xl = Util::fast_floor(x);
// 	int yl = Util::fast_floor(y);
// 	int zl = Util::fast_floor(z);
//
// 	// w is the windowing width
// 	int w = params.set_default("maskwidth",2);
// 	float dw = 1.0/w;
// 	dw *= dw;
// // 	dw = 2;
// // 	cout << w << endl;
// 	// 	int w = 3;
// 	// w minus one - this control the number of
// 	// pixels/voxels to the left of the main pixel
// 	// that will have density
// 	int wmox = w-1;
// 	int wmoy = w-1;
// 	int wmoz = w-1;
//
// 	// If any coordinate is incedental with a vertex, then
// 	// make sure there is symmetry in density accruing.
// 	// i.e. the window width must be equal in both directions
// 	if ( ((float) xl) == x ) wmox = w;
// 	if ( ((float) yl) == y ) wmoy = w;
// 	if ( ((float) zl) == z ) wmoz = w;
//
// 	float* d = tmp_data->get_data();
// 	int tnx = tmp_data->get_xsize();
// 	int tny = tmp_data->get_ysize();
// 	int tnz = tmp_data->get_zsize();
// 	int tnxy = tnx*tny;
//
// 	float weight = 1.0;
// //
// 	for(int k = zl-wmoz; k <= zl+w; ++k ) {
// 		for(int j = yl-wmoy; j <= yl+w; ++j) {
// 			for( int i = xl-wmox; i <= xl+w; ++i) {
// 				float fac = 1.0;
// 				int ic = i, jc = j, kc = k;
//
// 				// Fourier space is periodic, which is enforced
// 				// by the next 6 if statements. These if statements
// 				// assume that the Fourier DC components is at
// 				// (0,ny/2,nz/2).
// 				if ( i <= 0 ) {
// 					if ( x != 0 && i == 0 ) fac = 1.0;
// 					else if ( x == 0 && i < 0) continue;
// // 					if (i < 0 ) ic = -i;
// 					if (i < 0 ) {
// 						ic = -i;
// 						jc = tny-jc;
// 						kc = tnz-kc;
// 					}
// 				}
// 				if ( ic >= tnx ) ic = 2*tnx-ic-1;
// 				if ( jc < 0 ) jc = tny+jc;
// 				if ( jc >= tny ) jc = jc-tny;
// 				if ( kc < 0 ) kc = tnz+kc;
// 				if ( kc >= tnz ) kc = kc-tnz;
// 				// This shouldn't happen
// 				// Debug remove later
// 				if ( ic < 0 ) { cout << "wo 1" << endl; }
// 				if ( ic >= tnx  ){ cout << "wo 2" << endl; }
// 				if ( jc < 0 ) { cout << "wo 3" << endl; }
// 				if ( jc >= tny ) { cout << "wo 4" << endl; }
// 				if ( kc < 0 ) { cout << "wo 5" << endl; }
// 				if ( kc >= tnz ) { cout << "wo 6" << endl; }
//
//
// 				float zd = (z-(float)k);
// 				float yd = (y-(float)j);
// 				float xd = (x-(float)i);
// 				zd *= zd; yd *= yd; xd *= xd;
// 				// Debug - this error should never occur.
// 				if ( (kc*tnxy+jc*tnx+ic) >= tnxy*tnz ) throw OutofRangeException(0,tnxy*tnz,kc*tnxy+jc*tnx+ic, "in weight determination insertion" );
// // 				float f = fac*exp(-dw*(xd+yd+zd)*0.5);
// 				float f = exp(-2.467*(xd+yd+zd));
// 				weight += f*(d[kc*tnxy+jc*tnx+ic] - f);
// 			}
// 		}
// 	}
// 	weight = 1.0/weight;
// 	int rnx = 2*tnx;
// 	int rnxy = 2*tnxy;
// 	d = image->get_data();
// 	for(int k = zl-wmoz; k <= zl+w; ++k ) {
// 		for(int j = yl-wmoy; j <= yl+w; ++j) {
// 			for( int i = xl-wmox; i <= xl+w; ++i) {
// 				float fac = 1.0;
// 				int ic = i, jc = j, kc = k;
//
// 				// Fourier space is periodic, which is enforced
// 				// by the next 6 if statements. These if statements
// 				// assume that the Fourier DC components is at
// 				// (0,ny/2,nz/2).
// 				float negfac=1.0;
// 				if ( i <= 0 ) {
// 					if ( x != 0 && i == 0 ) fac = 1.0;
// 					else if ( x == 0 && i < 0) continue;
// 					if (i < 0 ) {
// 						continue;
// 						ic = -i;
// 						jc = tny-jc;
// 						kc = tnz-kc;
// 						negfac=-1.0;
// 					}
// 				}
// 				if ( ic >= tnx ) ic = 2*tnx-ic-1;
// 				if ( jc < 0 ) jc = tny+jc;
// 				if ( jc >= tny ) jc = jc-tny;
// 				if ( kc < 0 ) kc = tnz+kc;
// 				if ( kc >= tnz ) kc = kc-tnz;
// 				// This shouldn't happen
// 				// Debug remove later
//
//
// 				float zd = (z-(float)k);
// 				float yd = (y-(float)j);
// 				float xd = (x-(float)i);
// 				zd *= zd; yd *= yd; xd *= xd;
// // 				float f = fac*exp(-dw*(xd+yd+zd));
// 				float f = exp(-4.934*(xd+yd+zd));
// 				// Debug - this error should never occur.
// 				if ( (kc*rnxy+jc*rnx+2*ic+1) >= rnxy*tnz ) throw OutofRangeException(0,rnxy*tnz,kc*rnxy+jc*rnx+2*ic+1, "in pixel insertion" );
//
// 				d[kc*rnxy+jc*rnx+2*ic] += weight*f*dt[0];
// 				d[kc*rnxy+jc*rnx+2*ic+1] += negfac*weight*f*dt[1];
// 			}
// 		}
// 	}
// }


void WienerFourierReconstructor::setup()
{
	int size = params["size"];
	image = new EMData();
	image->set_size(size + 2, size, size);
	image->set_complex(true);
	image->set_ri(true);

	nx = image->get_xsize();
	ny = image->get_ysize();
	nz = image->get_zsize();

	int n = nx * ny * nz;
	float *rdata = image->get_data();

	for (int i = 0; i < n; i += 2) {
		float f = Util::get_frand(0.0, 2.0 * M_PI);
		rdata[i] = 1.0e-10f * sin(f);
		rdata[i + 1] = 1.0e-10f * cos(f);
	}
	image->update();

	tmp_data = new EMData();
	tmp_data->set_size(size + 1, size, size);
}


EMData *WienerFourierReconstructor::finish()
{
	float *norm = tmp_data->get_data();
	float *rdata = image->get_data();

	size_t size = nx * ny * nz;
	for (size_t i = 0; i < size; i += 2) {
		float d = norm[i];
		if (d == 0) {
			rdata[i] = 0;
			rdata[i + 1] = 0;
		}
		else {
			float w = 1 + 1 / d;
			rdata[i] /= d * w;
			rdata[i + 1] /= d * w;
		}
	}

	if( tmp_data ) {
		delete tmp_data;
		tmp_data = 0;
	}
	image->update();
	return image;
}


int WienerFourierReconstructor::insert_slice(const EMData* const slice, const Transform3D & euler)
{
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int mode = params["mode"];
	float padratio = params["padratio"];
	vector < float >snr = params["snr"];

	if (!slice->is_complex()) {
		LOGERR("Only complex slice can be inserted.");
		return 1;
	}
	float *gimx = 0;
	if (mode == 5) {
		gimx = Interp::get_gimx();
	}

	int nxy = nx * ny;
	int off[8] = {0,0,0,0,0,0,0,0};
	if (mode == 2) {
		off[0] = 0;
		off[1] = 2;
		off[2] = nx;
		off[3] = nx + 2;
		off[4] = nxy;
		off[5] = nxy + 2;
		off[6] = nxy + nx;
		off[7] = nxy + nx + 2;
	}

	float *norm = tmp_data->get_data();
	float *dat = slice->get_data();
	float *rdata = image->get_data();

	int rl = Util::square(ny / 2 - 1);
	float dt[2];
	float g[8];

	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx / 2; x++) {
			if ((x * x + Util::square(y - ny / 2)) >= rl)
			{
				continue;
			}

#ifdef _WIN32
			int r = Util::round((float)_hypot(x, (float) y - ny / 2) * Ctf::CTFOS / padratio);
#else
			int r = Util::round((float)hypot(x, (float) y - ny / 2) * Ctf::CTFOS / padratio);
#endif
			if (r >= Ctf::CTFOS * ny / 2) {
				r = Ctf::CTFOS * ny / 2 - 1;
			}

			float weight = snr[r];

			float xx = (x * euler[0][0] + (y - ny / 2) * euler[0][1]);
			float yy = (x * euler[1][0] + (y - ny / 2) * euler[1][1]);
			float zz = (x * euler[2][0] + (y - ny / 2) * euler[2][1]);
			float cc = 1;

			if (xx < 0) {
				xx = -xx;
				yy = -yy;
				zz = -zz;
				cc = -1.0;
			}

			yy += ny / 2;
			zz += nz / 2;

			dt[0] = dat[x * 2 + y * nx] * (1 + 1.0f / weight);
			dt[1] = cc * dat[x * 2 + 1 + y * nx] * (1 + 1.0f / weight);

			int x0 = 0;
			int y0 = 0;
			int z0 = 0;
			int i = 0;
			int l = 0;
			float dx = 0;
			float dy = 0;
			float dz = 0;

			int mx0 = 0;
			int my0 = 0;
			int mz0 = 0;

			size_t idx;
			switch (mode) {
			case 1:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				rdata[x0 + y0 * nx + z0 * nxy] += weight * dt[0];
				rdata[x0 + y0 * nx + z0 * nxy + 1] += weight * dt[1];
				norm[x0 + y0 * nx + z0 * nxy] += weight;
				break;

			case 2:
				x0 = (int) floor(xx);
				y0 = (int) floor(yy);
				z0 = (int) floor(zz);

				dx = xx - x0;
				dy = yy - y0;
				dz = zz - z0;

				weight /= (float)pow((float)(EMConsts::I2G * M_PI), 1.5f);

				if (x0 > nx - 2 || y0 > ny - 1 || z0 > nz - 1) {
					break;
				}

				i = (int) (x0 * 2 + y0 * nx + z0 * nxy);


				g[0] = Util::agauss(1, dx, dy, dz, EMConsts::I2G);
				g[1] = Util::agauss(1, 1 - dx, dy, dz, EMConsts::I2G);
				g[2] = Util::agauss(1, dx, 1 - dy, dz, EMConsts::I2G);
				g[3] = Util::agauss(1, 1 - dx, 1 - dy, dz, EMConsts::I2G);
				g[4] = Util::agauss(1, dx, dy, 1 - dz, EMConsts::I2G);
				g[5] = Util::agauss(1, 1 - dx, dy, 1 - dz, EMConsts::I2G);
				g[6] = Util::agauss(1, dx, 1 - dy, 1 - dz, EMConsts::I2G);
				g[7] = Util::agauss(1, 1 - dx, 1 - dy, 1 - dz, EMConsts::I2G);

				for (int j = 0; j < 8; j++) {
					int k = i + off[j];
					rdata[k] += weight * dt[0] * g[j];
					rdata[k + 1] += weight * dt[1] * g[j];
					norm[k] += weight * g[j];
				}

				break;
			case 3:
				x0 = 2 * (int) floor(xx + 0.5f);
				y0 = (int) floor(yy + 0.5f);
				z0 = (int) floor(zz + 0.5f);

				weight /= (float)pow((float)(EMConsts::I3G * M_PI), 1.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = x0 - 2;
				if (x0 == 0) {
					l = x0;
				}

				for (int k = z0 - 1; k <= z0 + 1; k++) {
					for (int j = y0 - 1; j <= y0 + 1; j++) {
						for (int i = l; i <= x0 + 2; i += 2) {
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = exp(-r / EMConsts::I3G);

							idx = i + j * nx + k * nxy;
							rdata[idx] += weight * gg * dt[0];
							rdata[idx + 1] += weight * gg * dt[1];
							norm[idx] += weight * gg;
						}
					}
				}
				break;

			case 4:
				x0 = 2 * (int) floor(xx);
				y0 = (int) floor(yy);
				z0 = (int) floor(zz);

				weight /= (float)pow((float)(EMConsts::I4G * M_PI), 1.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
					break;
				}

				l = x0 - 2;
				if (x0 == 0) {
					l = x0;
				}

				for (int k = z0 - 1; k <= z0 + 2; ++k) {
					for (int j = y0 - 1; j <= y0 + 2; ++j) {
						for (int i = l; i <= x0 + 4; i += 2) {
							float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
							float gg = exp(-r / EMConsts::I4G);

							idx = i + j * nx + k * nxy;
							rdata[idx] += weight * gg * dt[0];
							rdata[idx + 1] += weight * gg * dt[1];
							norm[idx] += weight * gg;
						}
					}
				}
				break;

			case 5:
				x0 = (int) floor(xx + .5);
				y0 = (int) floor(yy + .5);
				z0 = (int) floor(zz + .5);

				weight /= (float)pow((float)(EMConsts::I5G * M_PI), 1.5f);

				mx0 = -(int) floor((xx - x0) * 39.0f + 0.5) - 78;
				my0 = -(int) floor((yy - y0) * 39.0f + 0.5) - 78;
				mz0 = -(int) floor((zz - z0) * 39.0f + 0.5) - 78;
				x0 *= 2;

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				if (x0 == 0) {
					l = 0;
					mx0 += 78;
				}
				else if (x0 == 2) {
					l = 0;
					mx0 += 39;
				}
				else
					l = x0 - 4;
				for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
					for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
						for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
							size_t ii = i + j * nx + k * nxy;
							float gg = weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);
					int mx0 = -(int) floor((xx - x0) * 39.0f + .5);
					x0 *= 2;

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
						for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
							for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
								size_t ii = i + j * nx + k * nxy;
								float gg =
									weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
				break;
				// mode 6 is now mode 5 without the fast interpolation
			case 6:
				x0 = 2 * (int) floor(xx + .5);
				y0 = (int) floor(yy + .5);
				z0 = (int) floor(zz + .5);

				weight /= (float)pow((float)(EMConsts::I5G * M_PI), 1.5f);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				if (x0 <= 2)
					l = 0;
				else
					l = x0 - 4;
				for (int k = z0 - 2; k <= z0 + 2; ++k) {
					for (int j = y0 - 2; j <= y0 + 2; ++j) {
						for (int i = l; i <= x0 + 4; i += 2) {
							size_t ii = i + j * nx + k * nxy;
							float r = Util::hypot3((float) i / 2 - xx, (float) j - yy,
												   (float) k - zz);
							float gg = weight * exp(-r / EMConsts::I5G);

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = 2 * (int) floor(xx + 0.5f);
					y0 = (int) floor(yy + 0.5f);
					z0 = (int) floor(zz + 0.5f);

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2; k <= z0 + 2; ++k) {
						for (int j = y0 - 2; j <= y0 + 2; ++j) {
							for (int i = 0; i <= x0 + 4; i += 2) {
								size_t ii = i + j * nx + k * nxy;
								float r = Util::hypot3((float) i / 2 - xx, (float) j - yy,
													   (float) k - zz);
								float gg = weight * exp(-r / EMConsts::I5G);

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];	// note the -, complex conj.
								norm[ii] += gg;
							}
						}
					}
				}
				break;

			case 7:
				x0 = 2 * (int) floor(xx + .5);
				y0 = (int) floor(yy + .5);
				z0 = (int) floor(zz + .5);

				if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
					break;

				if (x0 <= 2)
					l = 0;
				else
					l = x0 - 4;
				for (int k = z0 - 2; k <= z0 + 2; k++) {
					for (int j = y0 - 2; j <= y0 + 2; j++) {
						for (int i = l; i <= x0 + 4; i += 2) {
							size_t ii = i + j * nx + k * nxy;
							float r = (float)sqrt(Util::hypot3((float) i / 2 - xx,
															   (float) j - yy,
															   (float) k - zz));
							float gg = weight * Interp::hyperg(r);

							rdata[ii] += gg * dt[0];
							rdata[ii + 1] += gg * dt[1];
							norm[ii] += gg;
						}
					}
				}

				if (x0 <= 2) {
					xx = -xx;
					yy = -(yy - ny / 2) + ny / 2;
					zz = -(zz - nz / 2) + nz / 2;
					x0 = 2 * (int) floor(xx + .5);
					y0 = (int) floor(yy + .5);
					z0 = (int) floor(zz + .5);

					if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
						break;

					for (int k = z0 - 2; k <= z0 + 2; ++k) {
						for (int j = y0 - 2; j <= y0 + 2; ++j) {
							for (int i = 0; i <= x0 + 4; i += 2) {
								size_t ii = i + j * nx + k * nxy;
								float r = sqrt(Util::hypot3((float) i / 2 - xx, (float) j - yy,
															(float) k - zz));
								float gg = weight * Interp::hyperg(r);

								rdata[ii] += gg * dt[0];
								rdata[ii + 1] -= gg * dt[1];
								norm[ii] += gg;
							}
						}
					}
				}
				break;
			}

		}
	}

	image->update();
	tmp_data->update();
//	slice->update();

	return 0;
}

void BackProjectionReconstructor::setup()
{
	int size = params["size"];
	image = new EMData();
	nx = size;
	ny = size;
	if ( (int) params["zsample"] != 0 ) nz = params["zsample"];
	else nz = size;
	image->set_size(nx, ny, nz);
}

EMData* BackProjectionReconstructor::preprocess_slice(const EMData* const slice, const Transform& t)
{

	EMData* return_slice = slice->process("normalize.edgemean");
	return_slice->process_inplace("filter.linearfourier");

	Transform tmp(t);
	tmp.set_rotation(Dict("type","eman")); // resets the rotation to 0 implicitly
	Vec2f trans = tmp.get_trans_2d();
	float scale = tmp.get_scale();
	bool mirror = tmp.get_mirror();
	if (trans[0] != 0 || trans[1] != 0 || scale != 1.0 ) {
		return_slice->transform(tmp);
	} else if ( mirror == true ) {
		return_slice = slice->process("xform.flip",Dict("axis","x"));
	}

	return return_slice;
}

int BackProjectionReconstructor::insert_slice(const EMData* const input, const Transform &t)
{
	if (!input) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	if (input->get_xsize() != input->get_ysize() || input->get_xsize() != nx) {
		LOGERR("tried to insert image that was not correction dimensions");
		return 1;
	}

	Transform transform;
	if ( input->has_attr("xform.projection") ) {
		transform = *((Transform*) input->get_attr("xform.projection")); // assignment operator
	} else {
		transform = t; // assignment operator
	}
	EMData* slice = preprocess_slice(input, t);

	float weight = params["weight"];
	slice->mult(weight);

	EMData *tmp = new EMData();
	tmp->set_size(nx, ny, nz);

	float *slice_data = slice->get_data();
	float *tmp_data = tmp->get_data();

	size_t nxy = nx * ny;
	size_t nxy_size = nxy * sizeof(float);;
	for (int i = 0; i < nz; ++i) {
		memcpy(&tmp_data[nxy * i], slice_data, nxy_size);
	}

	transform.set_scale(1.0);
	transform.set_mirror(false);
	transform.set_trans(0,0,0);
	transform.invert();

	tmp->transform(transform);
	image->add(*tmp);

	delete tmp;
	delete slice;

	return 0;
}

EMData *BackProjectionReconstructor::finish()
{

	Symmetry3D* sym = Factory<Symmetry3D>::get((string)params["sym"]);
	vector<Transform> syms = sym->get_syms();

	for ( vector<Transform>::const_iterator it = syms.begin(); it != syms.end(); ++it ) {

		EMData tmpcopy(*image);
		tmpcopy.transform(*it);
		image->add(tmpcopy);
	}

	image->mult(1.0f/(float)sym->get_nsym());
	delete sym;
	return image;
}

EMData* EMAN::padfft_slice( const EMData* const slice, const Transform& t, int npad )
{
        int nx = slice->get_xsize();
	int ny = slice->get_ysize();
        int ndim = (ny==1) ? 1 : 2;

	if( ndim==2 && nx!=ny )
	{
		// FIXME: What kind of exception should we throw here?
		throw std::runtime_error("Tried to padfft a 2D slice which is not square.");
	}

	// process 2D slice or 1D line -- subtract the average outside of the circle, zero-pad, fft extend, and fft
	EMData* temp = slice->average_circ_sub();

	Assert( temp != NULL );
	EMData* zeropadded = temp->norm_pad( false, npad );
	Assert( zeropadded != NULL );
	checked_delete( temp );

        zeropadded->do_fft_inplace();
	EMData* padfftslice = zeropadded;

	// shift the projection
	Vec2f trans = t.get_trans_2d();
	float sx = -trans[0];
	float sy = -trans[1];
	if(sx != 0.0f || sy != 0.0)
		padfftslice->process_inplace("filter.shift", Dict("x_shift", sx, "y_shift", sy, "z_shift", 0.0f));

        int remove = slice->get_attr_default("remove", 0);
        padfftslice->set_attr( "remove", remove );



      	padfftslice->center_origin_fft();
	return padfftslice;
}

nn4Reconstructor::nn4Reconstructor()
{
    m_volume = NULL;
    m_wptr   = NULL;
    m_result = NULL;
}

nn4Reconstructor::nn4Reconstructor( const string& symmetry, int size, int npad )
{
    m_volume = NULL;
    m_wptr   = NULL;
    m_result = NULL;
	setup( symmetry, size, npad );
	load_default_settings();
	print_params();
}

nn4Reconstructor::~nn4Reconstructor()
{
    if( m_delete_volume )
        checked_delete(m_volume);

    if( m_delete_weight )
        checked_delete( m_wptr );

    checked_delete( m_result );
}

enum weighting_method { NONE, ESTIMATE, VORONOI };

float max2d( int kc, const vector<float>& pow_a )
{
	float max = 0.0;
	for( int i=-kc; i <= kc; ++i ) {
		for( int j=-kc; j <= kc; ++j ) {
			if( i==0 && j==0 ) continue;
			{
				int c = 2*kc+1 - std::abs(i) - std::abs(j);
				max = max + pow_a[c];
			}
		}
	}
	return max;
}

float max3d( int kc, const vector<float>& pow_a )
{
	float max = 0.0;
	for( int i=-kc; i <= kc; ++i ) {
		for( int j=-kc; j <= kc; ++j ) {
			for( int k=-kc; k <= kc; ++k ) {
				if( i==0 && j==0 && k==0 ) continue;
				// if( i!=0 )
				{
					int c = 3*kc+1 - std::abs(i) - std::abs(j) - std::abs(k);
					max = max + pow_a[c];
					// max = max + c * c;
					// max = max + c;
				}
			}
		}
	}
	return max;
}


void nn4Reconstructor::setup()
{
	int size = params["size"];
	int npad = params["npad"];


	string symmetry;
	if( params.has_key("symmetry") ) {
		symmetry = params["symmetry"].to_str();
	} else {
		symmetry = "c1";
	}

	if( params.has_key("ndim") ) {
	    m_ndim = params["ndim"];
	} else {
		m_ndim = 3;
	}

    if( params.has_key( "snr" ) ) {
        m_osnr = 1.0f/float( params["snr"] );
    } else {
        m_osnr = 0.0;
    }

	setup( symmetry, size, npad );
}

void nn4Reconstructor::setup( const string& symmetry, int size, int npad )
{
	m_weighting = ESTIMATE;
	m_wghta = 0.2f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_vnx = size;
	m_vny = size;
	m_vnz = (m_ndim==3) ? size : 1;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = (m_ndim==3) ? size*npad : 1;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = (m_ndim==3) ? m_vnzp/2 : 1;

	buildFFTVolume();
	buildNormVolume();

}


void nn4Reconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;

	if( params.has_key("fftvol") ) {
		m_volume = params["fftvol"];
		m_delete_volume = false;
	} else {
		m_volume = new EMData();
		m_delete_volume = true;
	}

	if( m_volume->get_xsize() != m_vnxp+offset &&
	    m_volume->get_ysize() != m_vnyp &&
	    m_volume->get_zsize() != m_vnzp ) {
		m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
		m_volume->to_zero();
	}
	// ----------------------------------------------------------------
	// Added by Zhengfan Yang on 03/15/07
	// Original author: please check whether my revision is correct and
	// other Reconstructor need similiar revision.
	if ( m_vnxp % 2 == 0 )  m_volume->set_fftodd(0);
	else                    m_volume->set_fftodd(1);
	// ----------------------------------------------------------------

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4Reconstructor::buildNormVolume() {

	if( params.has_key("weight") ) {
		m_wptr = params["weight"];
		m_delete_weight = false;
	} else {
		m_wptr = new EMData();
		m_delete_weight = true;
	}

	if( m_wptr->get_xsize() != m_vnxc+1 &&
		m_wptr->get_ysize() != m_vnyp &&
		m_wptr->get_zsize() != m_vnzp ) {
		m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
		m_wptr->to_zero();
	}

	m_wptr->set_array_offsets(0,1,1);
}

void printImage( const EMData* line )
{
	Assert( line->get_zsize()==1 );


	int nx = line->get_xsize();
	int ny = line->get_ysize();
	for( int j=0; j < ny; ++j ) {
		for( int i=0; i < nx; ++i )  printf( "%10.3f ", line->get_value_at(i,j) );
		printf( "\n" );
	}
}



int nn4Reconstructor::insert_slice(const EMData* const slice, const Transform& t) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

        int padffted= slice->get_attr_default( "padffted", 0 );
        if( m_ndim==3 ) {
		if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
			// FIXME: Why doesn't this throw an exception?
			LOGERR("Tried to insert a slice that is the wrong size.");
			return 1;
		}
        } else {
		Assert( m_ndim==2 );
		if( slice->get_ysize() !=1 ) {
			LOGERR( "for 2D reconstruction, a line is excepted" );
        		return 1;
		}
	}

	EMData* padfft = NULL;

	if( padffted != 0 ) padfft = new EMData(*slice);
	else                padfft = padfft_slice( slice, t,  m_npad );

	int mult= slice->get_attr_default( "mult", 1 );
	Assert( mult > 0 );

        if( m_ndim==3 ) {
		insert_padfft_slice( padfft, t, mult );
	} else {
		float alpha = padfft->get_attr( "alpha" );
		alpha = alpha/180.0f*M_PI;
		for(int i=0; i < m_vnxc+1; ++i ) {
			float xnew = i*cos(alpha);
			float ynew = -i*sin(alpha);
			float btqr = padfft->get_value_at( 2*i, 0, 0 );
			float btqi = padfft->get_value_at( 2*i+1, 0, 0 );
			if( xnew < 0.0 ) {
				xnew *= -1;
				ynew *= -1;
				btqi *= -1;
			}

			int ixn = int(xnew+0.5+m_vnxp) - m_vnxp;
			int iyn = int(ynew+0.5+m_vnyp) - m_vnyp;

			if(iyn < 0 ) iyn += m_vnyp;

			(*m_volume)( 2*ixn, iyn+1, 1 ) += btqr *float(mult);
			(*m_volume)( 2*ixn+1, iyn+1, 1 ) += btqi * float(mult);
			(*m_wptr)(ixn,iyn+1, 1) += float(mult);
		}

	}
	checked_delete( padfft );
	return 0;
}

int nn4Reconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, int mult )
{
	Assert( padfft != NULL );
	// insert slice for all symmetry related positions
	for (int isym=0; isym < m_nsym; isym++) {
		Transform tsym = t.get_sym(m_symmetry, isym);
		m_volume->nn( m_wptr, padfft, tsym, mult);
        }
	return 0;
}

#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]

void circumference( EMData* win )
{
	float *tw = win->get_data();
	//  mask and subtract circumference average
	int ix = win->get_xsize();
	int iy = win->get_ysize();
	int iz = win->get_zsize();
	int L2 = (ix/2)*(ix/2);
	int L2P = (ix/2-1)*(ix/2-1);

	int IP = ix/2+1;
	int JP = iy/2+1;
	int KP = iz/2+1;

	float  TNR = 0.0f;
	size_t m = 0;
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				size_t LR = (k-KP)*(k-KP)+(j-JP)*(j-JP)+(i-IP)*(i-IP);
				if (LR<=(size_t)L2) {
					if(LR >= (size_t)L2P && LR <= (size_t)L2) {
						TNR += tw(i,j,k);
						++m;
					}
				}
			}
		}
	}

	TNR /=float(m);
	for (int k = 1; k <= iz; ++k) {
		for (int j = 1; j <= iy; ++j) {
			for (int i = 1; i <= ix; ++i) {
				size_t LR = (k-KP)*(k-KP)+(j-JP)*(j-JP)+(i-IP)*(i-IP);
				if (LR<=(size_t)L2) tw(i,j,k) -= TNR; else tw(i,j,k) = 0.0f;
			}
		}
	}

}


EMData* nn4Reconstructor::finish()
{
        if( m_ndim==3 ) {
		m_volume->symplane0(m_wptr);
	} else {
		for( int i=1; i <= m_vnyp; ++i ) {

			if( (*m_wptr)(0, i, 1)==0.0 ) {
				int j = m_vnyp + 1 - i;
				(*m_wptr)(0, i, 1) = (*m_wptr)(0, j, 1);
				(*m_volume)(0, i, 1) = (*m_volume)(0, j, 1);
				(*m_volume)(1, i, 1) = (*m_volume)(1, j, 1);
			}
		}
	}


	int box = 7;
	int kc = (box-1)/2;
	vector< float > pow_a( m_ndim*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a.back()=0.0;

	float alpha = 0.0;
	if( m_ndim==3) {
		int vol = box*box*box;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	} else {
		int ara = box*box;
		float max = max2d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)ara ) / max;
	}

	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0) {//(*v) should be treated as complex!!
					float tmp;
					tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+m_osnr);

					if( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;

								int kcz = (m_ndim==3) ? kc : 0;
								for( int kk = -kcz; kk <= kcz; ++kk ) {
									int nbrcz = cz + kk;
									if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
							    			nbrcy = -nbrcy;
							    			nbrcz = -nbrcz;
									}
									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = m_ndim*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
						float wght = 1.0f / ( 1.0f - alpha * sum );
						tmp = tmp * wght;
					}
					(*m_volume)(2*ix,iy,iz) *= tmp;
					(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

	//if(m_ndim==2) printImage( m_volume );

	// back fft
	m_volume->do_ift_inplace();

	// EMData* win = m_volume->window_center(m_vnx);
	m_volume->depad();
	circumference( m_volume );
	m_volume->set_array_offsets( 0, 0, 0 );

    m_result = m_volume->copy();
	return m_result;
}
#undef  tw

// Added By Zhengfan Yang on 03/16/07
// Beginning of the addition
// --------------------------------------------------------------------------------

nnSSNR_Reconstructor::nnSSNR_Reconstructor()
{
	m_volume = NULL;
	m_wptr   = NULL;
	m_wptr2  = NULL;
	m_result = NULL;
}

nnSSNR_Reconstructor::nnSSNR_Reconstructor( const string& symmetry, int size, int npad)
{
	m_volume = NULL;
	m_wptr   = NULL;
	m_wptr2  = NULL;
	m_result = NULL;

	setup( symmetry, size, npad );
}

nnSSNR_Reconstructor::~nnSSNR_Reconstructor()
{
	if( m_delete_volume ) checked_delete(m_volume);

	if( m_delete_weight ) checked_delete( m_wptr );

	if( m_delete_weight2 ) checked_delete( m_wptr2 );

	checked_delete( m_result );
}

void nnSSNR_Reconstructor::setup()
{
	int size = params["size"];
	int npad = params["npad"];

	string symmetry;
	if( params.has_key("symmetry") ) symmetry = params["symmetry"].to_str();
	else				 symmetry = "c1";

	setup( symmetry, size, npad );
}

void nnSSNR_Reconstructor::setup( const string& symmetry, int size, int npad )
{

	m_weighting = ESTIMATE;
	m_wghta = 0.2f;
	m_wghtb = 0.004f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_vnx = size;
	m_vny = size;
	m_vnz = size;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = size*npad;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
	buildNorm2Volume();

}

void nnSSNR_Reconstructor::buildFFTVolume() {

	if( params.has_key("fftvol") ) {
		m_volume = params["fftvol"]; /* volume should be defined in python when PMI is turned on*/
		m_delete_volume = false;
	} else {
		m_volume = new EMData();
		m_delete_volume = true;
	}
	m_volume->set_size(m_vnxp+ 2 - m_vnxp%2,m_vnyp,m_vnzp);
	m_volume->to_zero();
	if ( m_vnxp % 2 == 0 ) m_volume->set_fftodd(0);
	else                   m_volume->set_fftodd(1);

	m_volume->set_nxc(m_vnxc);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nnSSNR_Reconstructor::buildNormVolume() {
	if( params.has_key("weight") ) {
	 	m_wptr          = params["weight"];
	 	m_delete_weight = false;
	} else {
		m_wptr = new EMData();
		m_delete_weight = true;
	}

	m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr->to_zero();

	m_wptr->set_array_offsets(0,1,1);
}

void nnSSNR_Reconstructor::buildNorm2Volume() {

	if( params.has_key("weight2") ) {
		m_wptr2          = params["weight2"];
		m_delete_weight2 = false;
	} else {
		m_wptr2 = new EMData();
		m_delete_weight2 = true;
	}
	m_wptr2->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr2->to_zero();
	m_wptr2->set_array_offsets(0,1,1);
}


int nnSSNR_Reconstructor::insert_slice(const EMData* const slice, const Transform& t) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int padffted=slice->get_attr_default( "padffted", 0 );

	if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
		// FIXME: Why doesn't this throw an exception?
		LOGERR("Tried to insert a slice that has wrong size.");
		return 1;
	}

	EMData* padfft = NULL;

	if( padffted != 0 ) padfft = new EMData(*slice);
	else		    padfft = padfft_slice( slice, t, m_npad );

	int mult = slice->get_attr_default("mult", 1);

	Assert( mult > 0 );
	insert_padfft_slice( padfft, t, mult );

	checked_delete( padfft );
	return 0;
}

int nnSSNR_Reconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, int mult )
{
	Assert( padfft != NULL );
	// insert slice for all symmetry related positions
	for (int isym=0; isym < m_nsym; isym++) {
		Transform tsym = t.get_sym(m_symmetry, isym);
		m_volume->nn_SSNR( m_wptr, m_wptr2, padfft, tsym, mult);
	}
	return 0;
}


#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]
EMData* nnSSNR_Reconstructor::finish()
{
/*
  I changed the code on 05/15 so it only returns variance.
  Lines commented out are marked by //#
  The version prior to the currect changes is r1.190
  PAP
*/
	int kz, ky;
	//#int iix, iiy, iiz;
	int box = 7;
	int kc = (box-1)/2;
	float alpha = 0.0;
	float argx, argy, argz;
	vector< float > pow_a( 3*kc+1, 1.0 );
	float w = params["w"];
	EMData* SSNR = params["SSNR"];
	//#EMData* vol_ssnr = new EMData();
	//#vol_ssnr->set_size(m_vnxp, m_vnyp, m_vnzp);
	//#vol_ssnr->to_zero();
	//#  new line follow
	EMData* vol_ssnr = new EMData();
	vol_ssnr->set_size(m_vnxp+ 2 - m_vnxp%2, m_vnyp ,m_vnzp);
	vol_ssnr->to_zero();
	if ( m_vnxp % 2 == 0 ) vol_ssnr->set_fftodd(0);
	else                   vol_ssnr->set_fftodd(1);
	vol_ssnr->set_nxc(m_vnxc);
	vol_ssnr->set_complex(true);
	vol_ssnr->set_ri(true);
	vol_ssnr->set_fftpad(false);
	//#

	float dx2 = 1.0f/float(m_vnxc)/float(m_vnxc);
	float dy2 = 1.0f/float(m_vnyc)/float(m_vnyc);
#ifdef _WIN32
	float dz2 = 1.0f/_cpp_max(float(m_vnzc),1.0f)/_cpp_max(float(m_vnzc),1.0f);
	int  inc = Util::round(float(_cpp_max(_cpp_max(m_vnxc,m_vnyc),m_vnzc))/w);
#else
	float dz2 = 1.0f/std::max(float(m_vnzc),1.0f)/std::max(float(m_vnzc),1.0f);
	int  inc = Util::round(float(std::max(std::max(m_vnxc,m_vnyc),m_vnzc))/w);
#endif	//_WIN32
	SSNR->set_size(inc+1,4,1);

	float *nom    = new float[inc+1];
	float *denom  = new float[inc+1];
	int  *nn     = new int[inc+1];
	int  *ka     = new int[inc+1];
	float wght = 1.0f;
	for (int i = 0; i <= inc; i++) {
		nom[i] = 0.0f;
		denom[i] = 0.0f;
		nn[i] = 0;
		ka[i] = 0;
	}

	m_volume->symplane1(m_wptr, m_wptr2);

	if ( m_weighting == ESTIMATE ) {
		int vol = box*box*box;
		for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
		pow_a[3*kc] = 0.0;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	}

	for (int iz = 1; iz <= m_vnzp; iz++) {
		if ( iz-1 > m_vnzc ) kz = iz-1-m_vnzp; else kz = iz-1;
		argz = float(kz*kz)*dz2;
		for (int iy = 1; iy <= m_vnyp; iy++) {
			if ( iy-1 > m_vnyc ) ky = iy-1-m_vnyp; else ky = iy-1;
			argy = argz + float(ky*ky)*dy2;
			for (int ix = 0; ix <= m_vnxc; ix++) {
				float Kn = (*m_wptr)(ix,iy,iz);
				argx = std::sqrt(argy + float(ix*ix)*dx2);
				int r = Util::round(float(inc)*argx);
				if ( r >= 0 && Kn > 4.5f ) {
					if ( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;

						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
						    for ( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
								for( int kk = -kc; kk <= kc; ++kk ) {
									int nbrcz = cz + jj;
									if ( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
										nbrcy = -nbrcy;
										nbrcz = -nbrcz;
									}
		                        				int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
						wght = 1.0f / ( 1.0f - alpha * sum );
					} // end of ( m_weighting == ESTIMATE )
					float nominator = std::norm(m_volume->cmplx(ix,iy,iz))/Kn;
					float denominator = ((*m_wptr2)(ix,iy,iz)-nominator)/(Kn-1.0f);
					// Skip Friedel related values
					if( (ix>0 || (kz>=0 && (ky>=0 || kz!=0)))) {
						if ( r <= inc ) {
							nom[r]   += nominator*wght;
							denom[r] += denominator/Kn*wght;
							nn[r]    += 2;
							ka[r]    += int(Kn);
						}
/*
#ifdef	_WIN32
						//#float  tmp = _cpp_max(nominator/denominator/Kn-1.0f,0.0f);
#else
						//#float  tmp = std::max(nominator/denominator/Kn-1.0f,0.0f);
#endif	//_WIN32
						//  Create SSNR as a 3D array (-n/2:n/2+n%2-1)
						iix = m_vnxc + ix; iiy = m_vnyc + ky; iiz = m_vnzc + kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
						// Friedel part
						iix = m_vnxc - ix; iiy = m_vnyc - ky; iiz = m_vnzc - kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
*/

					}
					(*vol_ssnr)(2*ix, iy-1, iz-1) = denominator*wght;
					//(*vol_ssnr)(2*ix, iy-1, iz-1) =  real(m_volume->cmplx(ix,iy,iz))*wght/Kn;
					//(*vol_ssnr)(2*ix+1, iy-1, iz-1) = imag(m_volume->cmplx(ix,iy,iz))*wght/Kn;
				} // end of Kn>4.5
			}
		}
	}

	for (int i = 0; i <= inc; i++)  {
		(*SSNR)(i,0,0) = nom[i];  ///(*SSNR)(i,0,0) = nom[i]/denom[i] - 1;///
		(*SSNR)(i,1,0) = denom[i];    // variance
		(*SSNR)(i,2,0) = static_cast<float>(nn[i]);
		(*SSNR)(i,3,0) = static_cast<float>(ka[i]);
	}
	vol_ssnr->update();
	return vol_ssnr;
}
#undef  tw

// -----------------------------------------------------------------------------------
// End of this addition

//####################################################################################
//** nn4 ctf reconstructor

nn4_ctfReconstructor::nn4_ctfReconstructor()
{
	m_volume  = NULL;
	m_wptr    = NULL;
	m_result  = NULL;
}

nn4_ctfReconstructor::nn4_ctfReconstructor( const string& symmetry, int size, int npad, float snr, int sign )
{
	setup( symmetry, size, npad, snr, sign );
}

nn4_ctfReconstructor::~nn4_ctfReconstructor()
{
	if( m_delete_volume ) checked_delete(m_volume);

	if( m_delete_weight ) checked_delete( m_wptr );

	checked_delete( m_result );
}

void nn4_ctfReconstructor::setup()
{
	if( ! params.has_key("size") ) throw std::logic_error("Error: image size is not given");

	int size = params["size"];
	int npad = params.has_key("npad") ? int(params["npad"]) : 4;
	// int sign = params.has_key("sign") ? int(params["sign"]) : 1;
	int sign = 1;
	string symmetry = params.has_key("symmetry")? params["symmetry"].to_str() : "c1";

	float snr = params["snr"];

    m_varsnr = params.has_key("varsnr") ? int(params["varsnr"]) : 0;
	setup( symmetry, size, npad, snr, sign );

}

void nn4_ctfReconstructor::setup( const string& symmetry, int size, int npad, float snr, int sign )
{
	m_weighting = ESTIMATE;
	if( params.has_key("weighting") ) {
		int tmp = int( params["weighting"] );
		if( tmp==0 ) m_weighting = NONE;
	}



	m_wghta = 0.2f;
	m_wghtb = 0.004f;

	m_symmetry = symmetry;
	m_npad = npad;
	m_sign = sign;
	m_nsym = Transform::get_nsym(m_symmetry);

	m_snr = snr;

	m_vnx = size;
	m_vny = size;
	m_vnz = size;

	m_vnxp = size*npad;
	m_vnyp = size*npad;
	m_vnzp = size*npad;

	m_vnxc = m_vnxp/2;
	m_vnyc = m_vnyp/2;
	m_vnzc = m_vnzp/2;

	buildFFTVolume();
	buildNormVolume();
}

void nn4_ctfReconstructor::buildFFTVolume() {
	int offset = 2 - m_vnxp%2;
	if( params.has_key("fftvol") ) {
	    m_volume = params["fftvol"];
	    m_delete_volume = false;
	} else {
	    m_volume = new EMData();
	    m_delete_volume = true;
	}

	if( m_volume->get_xsize() != m_vnxp+offset &&
	    m_volume->get_ysize() != m_vnyp &&
	    m_volume->get_zsize() != m_vnzp )
	{
	    m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
	    m_volume->to_zero();
	}

	m_volume->set_nxc(m_vnxp/2);
	m_volume->set_complex(true);
	m_volume->set_ri(true);
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}

void nn4_ctfReconstructor::buildNormVolume()
{
	if( params.has_key("weight") )
	{
	    m_wptr = params["weight"];
	    m_delete_weight = false;
	} else {
	    m_wptr = new EMData();
	    m_delete_weight = true;
	}

	if( m_wptr->get_xsize() != m_vnxc+1 &&
	    m_wptr->get_ysize() != m_vnyp &&
	    m_wptr->get_zsize() != m_vnzp )
	{
	    m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	    m_wptr->to_zero();
	}

	m_wptr->set_array_offsets(0,1,1);

}

int nn4_ctfReconstructor::insert_slice(const EMData* const slice, const Transform& t)
{
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}

	int buffed = slice->get_attr_default( "buffed", 0 );
        if( buffed > 0 )
        {
            int mult = slice->get_attr_default( "mult", 1 );
            insert_buffed_slice( slice, mult );
            return 0;
        }

	int padffted= slice->get_attr_default("padffted", 0);
	if( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  )
        {
		// FIXME: Why doesn't this throw an exception?
		LOGERR("Tried to insert a slice that is the wrong size.");
		return 1;
	}

	EMData* padfft = NULL;

	if( padffted != 0 ) padfft = new EMData(*slice);
	else                padfft = padfft_slice( slice, t, m_npad );

	int mult= slice->get_attr_default("mult", 1);

	Assert( mult > 0 );
	insert_padfft_slice( padfft, t, mult );

	checked_delete( padfft );

	return 0;
}

int nn4_ctfReconstructor::insert_buffed_slice( const EMData* buffed, int mult )
{
	const float* bufdata = buffed->get_data();
	float* cdata = m_volume->get_data();
	float* wdata = m_wptr->get_data();

	int npoint = buffed->get_xsize()/4;
	for( int i=0; i < npoint; ++i ) {

               int pos2 = int( bufdata[4*i] );
               int pos1 = pos2 * 2;
               cdata[pos1  ] += bufdata[4*i+1]*mult;
               cdata[pos1+1] += bufdata[4*i+2]*mult;
               wdata[pos2  ] += bufdata[4*i+3]*mult;
/*
        std::cout << "pos1, pos2, ctfv1, ctfv2, ctf2: ";
        std::cout << pos1 << " " << bufdata[5*i+1] << " " << bufdata[5*i+2] << " ";
        std::cout << pos2 << " " << bufdata[5*i+4] << std::endl;
 */
	}
	return 0;
}

int nn4_ctfReconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, int mult )
{
	Assert( padfft != NULL );
	float tmp = padfft->get_attr("ctf_applied");
	int   ctf_applied = (int) tmp;

	for( int isym=0; isym < m_nsym; isym++) {
		Transform tsym = t.get_sym( m_symmetry, isym );

		if(ctf_applied) m_volume->nn_ctf_applied(m_wptr, padfft, tsym, mult);
		else            m_volume->nn_ctf(m_wptr, padfft, tsym, mult);
	}

	return 0;

}

#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]
EMData* nn4_ctfReconstructor::finish()
{
	m_volume->set_array_offsets(0, 1, 1);
	m_wptr->set_array_offsets(0, 1, 1);
	m_volume->symplane0_ctf(m_wptr);

	int box = 7;
	int vol = box*box*box;
	int kc = (box-1)/2;
	vector< float > pow_a( 3*kc+1, 1.0 );
	for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
	pow_a[3*kc]=0.0;


	float max = max3d( kc, pow_a );
	float alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	float osnr = 1.0f/m_snr;

	// normalize
	int ix,iy,iz;
	for (iz = 1; iz <= m_vnzp; iz++) {
		for (iy = 1; iy <= m_vnyp; iy++) {
			for (ix = 0; ix <= m_vnxc; ix++) {
				if ( (*m_wptr)(ix,iy,iz) > 0.0f) {//(*v) should be treated as complex!!
                    int iyp = (iy<=m_vnyc) ? iy - 1 : iy-m_vnyp-1;
                    int izp = (iz<=m_vnzc) ? iz - 1 : iz-m_vnzp-1;
                    float tmp=0.0;
                    if( m_varsnr )
                    {
					    float freq = sqrt( (float)(ix*ix+iyp*iyp+izp*izp) );
                        tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+freq*osnr)*m_sign;
                    }
                    else
                    {
                        tmp = (-2*((ix+iy+iz)%2)+1)/((*m_wptr)(ix,iy,iz)+osnr)*m_sign;
                    }

					if( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
								for( int kk = -kc; kk <= kc; ++kk ) {
									int nbrcz = cz + jj;
									if( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
										nbrcy = -nbrcy;
										nbrcz = -nbrcz;
									}

									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0.0 ) {
										int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
							          		  // if(ix%20==0 && iy%20==0 && iz%20==0)
							           		 //   std::cout << boost::format( "%4d %4d %4d %4d %10.3f" ) % nbrix % nbriy % nbriz % c % sum << std::endl;
									}
								}
							}
						}
						float wght = 1.0f / ( 1.0f - alpha * sum );
/*
                        if(ix%10==0 && iy%10==0)
                        {
                            std::cout << boost::format( "%4d %4d %4d " ) % ix % iy %iz;
                            std::cout << boost::format( "%10.3f %10.3f %10.3f " )  % tmp % wght % sum;
                            std::  << boost::format( "%10.3f %10.3e " ) % pow_b[r] % alpha;
                            std::cout << std::endl;
                        }
 */
						tmp = tmp * wght;
					}
					(*m_volume)(2*ix,iy,iz) *= tmp;
					(*m_volume)(2*ix+1,iy,iz) *= tmp;
				}
			}
		}
	}

    // back fft
    m_volume->do_ift_inplace();
    m_volume->depad();
    circumference( m_volume );
    m_volume->set_array_offsets( 0, 0, 0 );
    m_result = m_volume->copy();
    return m_result;
}
#undef  tw




// Added By Zhengfan Yang on 04/11/07
// Beginning of the addition
// --------------------------------------------------------------------------------

nnSSNR_ctfReconstructor::nnSSNR_ctfReconstructor()
{
    m_volume  = NULL;
    m_wptr    = NULL;
    m_wptr2   = NULL;
    m_wptr3   = NULL;
    m_result  = NULL;
}

nnSSNR_ctfReconstructor::nnSSNR_ctfReconstructor( const string& symmetry, int size, int npad, float snr, int sign)
{
    m_volume  = NULL;
    m_wptr    = NULL;
    m_wptr2   = NULL;
    m_wptr3   = NULL;
    m_result  = NULL;

    setup( symmetry, size, npad, snr, sign );
}

nnSSNR_ctfReconstructor::~nnSSNR_ctfReconstructor()
{

   if( m_delete_volume )
        checked_delete(m_volume);
   if( m_delete_weight )
        checked_delete( m_wptr );
   if ( m_delete_weight2 )
        checked_delete( m_wptr2 );
   if ( m_delete_weight3 )
        checked_delete( m_wptr3 );
   checked_delete( m_result );
}

void nnSSNR_ctfReconstructor::setup()
{
    int  size = params["size"];
    int  npad = params["npad"];
    int  sign = params["sign"];
    float snr = params["snr"];
    string symmetry;
    if( params.has_key("symmetry") ) {
		symmetry = params["symmetry"].to_str();
    } else {
		symmetry = "c1";
    }
    setup( symmetry, size, npad, snr, sign );
}
void nnSSNR_ctfReconstructor::setup( const string& symmetry, int size, int npad, float snr, int sign )
{

    m_weighting = ESTIMATE;
    m_wghta     = 0.2f;
    m_wghtb     = 0.004f;
    wiener      = 1;

    m_symmetry  = symmetry;
    m_npad      = npad;
    m_nsym      = Transform::get_nsym(m_symmetry);

    m_sign      = sign;
    m_snr       = snr;

    m_vnx       = size;
    m_vny       = size;
    m_vnz       = size;

    m_vnxp      = size*npad;
    m_vnyp      = size*npad;
    m_vnzp      = size*npad;

    m_vnxc      = m_vnxp/2;
    m_vnyc      = m_vnyp/2;
    m_vnzc      = m_vnzp/2;

    buildFFTVolume();
    buildNormVolume();
    buildNorm2Volume();
    buildNorm3Volume();
}

void nnSSNR_ctfReconstructor::buildFFTVolume() {

	int offset = 2 - m_vnxp%2;
	if( params.has_key("fftvol") ) {
		m_volume = params["fftvol"]; /* volume should be defined in python when PMI is turned on*/
		m_delete_volume = false;
    } else {
		m_volume = new EMData();
		m_delete_volume = true;
	}

	m_volume->set_size(m_vnxp+offset,m_vnyp,m_vnzp);
	m_volume->to_zero();

	if ( m_vnxp % 2 == 0 ) {
		m_volume->set_fftodd(0);
	} else {
		m_volume->set_fftodd(1);
	}

	m_volume->set_nxc(m_vnxc);
	m_volume->set_complex(true);
	m_volume->set_ri(true); //(real, imaginary) instead of polar coordinate
	m_volume->set_fftpad(true);
	m_volume->set_attr("npad", m_npad);
	m_volume->set_array_offsets(0,1,1);
}


void nnSSNR_ctfReconstructor::buildNormVolume()
{
	if( params.has_key("weight") ) {
		 m_wptr          = params["weight"];
		 m_delete_weight = false;
	} else {
		m_wptr = new EMData();
		m_delete_weight = true;
	}
	m_wptr->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr->to_zero();
	m_wptr->set_array_offsets(0,1,1);
}

void nnSSNR_ctfReconstructor::buildNorm2Volume() {

	if( params.has_key("weight2") ) {
		m_wptr2          = params["weight2"];
		m_delete_weight2 = false;
	} else {
		m_wptr2 = new EMData();
		m_delete_weight2 = true;
	}
	m_wptr2->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr2->to_zero();
	m_wptr2->set_array_offsets(0,1,1);
}

void nnSSNR_ctfReconstructor::buildNorm3Volume() {

	if( params.has_key("weight3") ) {
		m_wptr3          = params["weight3"];
		m_delete_weight3 = false;
	} else {
		m_wptr3 = new EMData();
		m_delete_weight3 = true;
	}
	m_wptr3->set_size(m_vnxc+1,m_vnyp,m_vnzp);
	m_wptr3->to_zero();
	m_wptr3->set_array_offsets(0,1,1);
}

int nnSSNR_ctfReconstructor::insert_slice(const EMData *const  slice, const Transform& t) {
	// sanity checks
	if (!slice) {
		LOGERR("try to insert NULL slice");
		return 1;
	}
	int padffted= slice->get_attr_default("padffted", 0);
	if ( padffted==0 && (slice->get_xsize()!=slice->get_ysize() || slice->get_xsize()!=m_vnx)  ) {
		// FIXME: Why doesn't this throw an exception?
		LOGERR("Tried to insert a slice that is the wrong size.");
		return 1;
	}
    EMData* padfft = NULL;

    if( padffted != 0 ) {
        padfft = new EMData(*slice);
    } else {
        padfft = padfft_slice( slice, t, m_npad );
    }

    int mult= slice->get_attr_default("mult", 1);

	Assert( mult > 0 );
	insert_padfft_slice( padfft, t, mult );

	checked_delete( padfft );
	return 0;
}
int nnSSNR_ctfReconstructor::insert_padfft_slice( EMData* padfft, const Transform& t, int mult )
{
	Assert( padfft != NULL );

	// insert slice for all symmetry related positions
	for (int isym=0; isym < m_nsym; isym++) {
		Transform tsym = t.get_sym(m_symmetry, isym);
		m_volume->nn_SSNR_ctf(m_wptr, m_wptr2, m_wptr3, padfft, tsym, mult);
	}

	return 0;
}

#define  tw(i,j,k)      tw[ i-1 + (j-1+(k-1)*iy)*ix ]
EMData* nnSSNR_ctfReconstructor::finish()
{
/*
  I changed the code on 05/15 so it only returns variance.
  Lines commented out are marked by //#
  The version prior to the currect changes is r1.190
  PAP
*/
	/***
	    m_volume ctf*(P^2D->3D(F^3D))
	    m_wptr   ctf^2
	    m_wptr2  |P^2D->3D(F^3D)|^2
	    m_wptr3  Kn
	    nominator = sum_rot [ wght*signal ]
	    denominator  = sum_rot[ wght*variance ]
						      ***/
	int kz, ky;
	int box = 7;
	int kc  = (box-1)/2;
	float alpha = 0.0;
	float argx, argy, argz;
	vector< float > pow_a( 3*kc+1, 1.0 );
	float w = params["w"];
	float dx2 = 1.0f/float(m_vnxc)/float(m_vnxc);
	float dy2 = 1.0f/float(m_vnyc)/float(m_vnyc);
#ifdef	_WIN32
	float dz2 = 1.0f/_cpp_max(float(m_vnzc),1.0f)/_cpp_max(float(m_vnzc),1.0f);
	int inc = Util::round(float(_cpp_max(_cpp_max(m_vnxc,m_vnyc),m_vnzc))/w);
#else
	float dz2 = 1.0f/std::max(float(m_vnzc),1.0f)/std::max(float(m_vnzc),1.0f);
	int inc = Util::round(float(std::max(std::max(m_vnxc,m_vnyc),m_vnzc))/w);
#endif	//_WIN32

	EMData* SSNR = params["SSNR"];
	SSNR->set_size(inc+1,4,1);
	//#EMData* vol_ssnr = new EMData();
	//#vol_ssnr->set_size(m_vnxp, m_vnyp, m_vnzp);
	//#vol_ssnr->to_zero();
	//#  new linea follow
	EMData* vol_ssnr = new EMData();
	vol_ssnr->set_size(m_vnxp+ 2 - m_vnxp%2, m_vnyp ,m_vnzp);
	vol_ssnr->to_zero();
	if ( m_vnxp % 2 == 0 ) vol_ssnr->set_fftodd(0);
	else                   vol_ssnr->set_fftodd(1);
	vol_ssnr->set_nxc(m_vnxc);
	vol_ssnr->set_complex(true);
	vol_ssnr->set_ri(true);
	vol_ssnr->set_fftpad(false);
	//#
	float *nom    = new float[inc+1];
	float *denom  = new float[inc+1];
	int  *ka     = new int[inc+1];
	int  *nn     = new int[inc+1];
	float wght = 1.f;
	for (int i = 0; i <= inc; i++) {
		nom[i]   = 0.0f;
		denom[i] = 0.0f;
		nn[i]	 = 0;
		ka[i]	 = 0;
	}
	m_volume->symplane2(m_wptr, m_wptr2, m_wptr3);
	if ( m_weighting == ESTIMATE ) {
		int vol = box*box*box;
		for( unsigned int i=1; i < pow_a.size(); ++i ) pow_a[i] = pow_a[i-1] * exp(m_wghta);
		pow_a[3*kc] = 0.0;
		float max = max3d( kc, pow_a );
		alpha = ( 1.0f - 1.0f/(float)vol ) / max;
	}
	for (int iz = 1; iz <= m_vnzp; iz++) {
		if ( iz-1 > m_vnzc ) kz = iz-1-m_vnzp; else kz = iz-1;
		argz = float(kz*kz)*dz2;
		for (int iy = 1; iy <= m_vnyp; iy++) {
			if ( iy-1 > m_vnyc ) ky = iy-1-m_vnyp; else ky = iy-1;
			argy = argz + float(ky*ky)*dy2;
			for (int ix = 0; ix <= m_vnxc; ix++) {
				float Kn = (*m_wptr3)(ix,iy,iz);
				argx = std::sqrt(argy + float(ix*ix)*dx2);
				int r = Util::round(float(inc)*argx);
				if ( r >= 0 && Kn > 4.5f ) {
					if ( m_weighting == ESTIMATE ) {
						int cx = ix;
						int cy = (iy<=m_vnyc) ? iy - 1 : iy - 1 - m_vnyp;
						int cz = (iz<=m_vnzc) ? iz - 1 : iz - 1 - m_vnzp;
						float sum = 0.0;
						for( int ii = -kc; ii <= kc; ++ii ) {
							int nbrcx = cx + ii;
							if( nbrcx >= m_vnxc ) continue;
							for ( int jj= -kc; jj <= kc; ++jj ) {
								int nbrcy = cy + jj;
								if( nbrcy <= -m_vnyc || nbrcy >= m_vnyc ) continue;
								for( int kk = -kc; kk <= kc; ++kk ) {
									int nbrcz = cz + jj;
									if ( nbrcz <= -m_vnyc || nbrcz >= m_vnyc ) continue;
									if( nbrcx < 0 ) {
										nbrcx = -nbrcx;
										nbrcy = -nbrcy;
										nbrcz = -nbrcz;
									}
									int nbrix = nbrcx;
									int nbriy = nbrcy >= 0 ? nbrcy + 1 : nbrcy + 1 + m_vnyp;
									int nbriz = nbrcz >= 0 ? nbrcz + 1 : nbrcz + 1 + m_vnzp;
									if( (*m_wptr)( nbrix, nbriy, nbriz ) == 0 ) {
										int c = 3*kc+1 - std::abs(ii) - std::abs(jj) - std::abs(kk);
										sum = sum + pow_a[c];
									}
								}
							}
						}
// 						int r = std::abs(cx) + std::abs(cy) + std::abs(cz);
						wght = 1.0f / ( 1.0f - alpha * sum );
					} // end of ( m_weighting == ESTIMATE )
					float nominator   = std::norm(m_volume->cmplx(ix,iy,iz))/(*m_wptr)(ix,iy,iz);
					float denominator = ((*m_wptr2)(ix,iy,iz)-std::norm(m_volume->cmplx(ix,iy,iz))/(*m_wptr)(ix,iy,iz))/(Kn-1.0f);
					// Skip Friedel related values
					if( (ix>0 || (kz>=0 && (ky>=0 || kz!=0)))) {
						if ( r <= inc ) {
							nom[r]   += nominator*wght;
							denom[r] += denominator/(*m_wptr)(ix,iy,iz)*wght;
							nn[r]	 += 2;
							ka[r]	 += int(Kn);
						}
/*
#ifdef	_WIN32
						float  tmp = _cpp_max(nominator/denominator/(*m_wptr)(ix,iy,iz)-1.0f,0.0f);
#else
						float  tmp = std::max(nominator/denominator/(*m_wptr)(ix,iy,iz)-1.0f,0.0f);
#endif	//_WIN32
						//  Create SSNR as a 3D array (-n/2:n/2+n%2-1)
						int iix = m_vnxc + ix; int iiy = m_vnyc + ky; int iiz = m_vnzc + kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
						// Friedel part
						iix = m_vnxc - ix; iiy = m_vnyc - ky; iiz = m_vnzc - kz;
						if( iix >= 0 && iix < m_vnxp && iiy >= 0 && iiy < m_vnyp && iiz >= 0 && iiz < m_vnzp )
							(*vol_ssnr)(iix, iiy, iiz) = tmp;
*/
					}
					(*vol_ssnr)(2*ix, iy-1, iz-1) = denominator*wght;
				} // end of Kn>4.5 or whatever
			}
		}
	}
	for (int i = 0; i <= inc; i++) {
		(*SSNR)(i,0,0) = nom[i];
		(*SSNR)(i,1,0) = denom[i];
		(*SSNR)(i,2,0) = static_cast<float>(nn[i]);
		(*SSNR)(i,3,0) = static_cast<float>(ka[i]);
	}
	vol_ssnr->update();
	return vol_ssnr;
//	}
}
#undef  tw
// -----------------------------------------------------------------------------------
// End of this addition


void EMAN::dump_reconstructors()
{
	dump_factory < Reconstructor > ();
}

map<string, vector<string> > EMAN::dump_reconstructors_list()
{
	return dump_factory_list < Reconstructor > ();
}





using std::ofstream;
using std::ifstream;


newfile_store::newfile_store( const string& filename, int npad, bool ctf )
    : m_bin_file( filename + ".bin" ),
      m_txt_file( filename + ".txt" )
{
    m_npad = npad;
    m_ctf = ctf;
}

newfile_store::~newfile_store( )
{
}

void newfile_store::add_image( EMData* emdata, const Transform& tf )
{
    if( m_bin_of == NULL )
    {
        m_bin_of = shared_ptr<ofstream>( new ofstream(m_bin_file.c_str(), std::ios::out|std::ios::binary) );
        m_txt_of = shared_ptr<ofstream>( new ofstream(m_txt_file.c_str()) );
    }


    EMData* padfft = padfft_slice( emdata, tf, m_npad );

    int nx = padfft->get_xsize();
    int ny = padfft->get_ysize();
    int n2 = ny / 2;
    int n = ny;

    float voltage=0.0f, pixel=0.0f, Cs=0.0f, ampcont=0.0f, bfactor=0.0f, defocus=0.0f;

    if( m_ctf )
    {
        Ctf* ctf = emdata->get_attr( "ctf" );
        Dict params = ctf->to_dict();
        voltage = params["voltage"];
        pixel   = params["apix"];
        Cs      = params["cs"];
        ampcont = params["ampcont"];
        bfactor = params["bfactor"];
        defocus = params["defocus"];
    }

    vector<point_t> points;
    for( int j=-ny/2+1; j <= ny/2; j++ )
    {
        int jp = (j>=0) ? j+1 : ny+j+1;
        for( int i=0; i <= n2; ++i )
        {
            int r2 = i*i + j*j;
            if( (r2<ny*ny/4) && !( (i==0) && (j<0) ) )
            {
                float ctf;
                if( m_ctf )
                {
		    float ak = std::sqrt( r2/float(ny*ny) )/pixel;
                    ctf = Util::tf( defocus, ak, voltage, Cs, ampcont, bfactor, 1);
                }
                else
                {
                    ctf = 1.0;
                }

                float xnew = i*tf[0][0] + j*tf[1][0];
                float ynew = i*tf[0][1] + j*tf[1][1];
                float znew = i*tf[0][2] + j*tf[1][2];
		std::complex<float> btq;
		if (xnew < 0.)
                {
                    xnew = -xnew;
                    ynew = -ynew;
                    znew = -znew;
                    btq = conj(padfft->cmplx(i,jp-1));
                }
                else
                {
                    btq = padfft->cmplx(i,jp-1);
                }

                int ixn = int(xnew + 0.5 + n) - n;
                int iyn = int(ynew + 0.5 + n) - n;
                int izn = int(znew + 0.5 + n) - n;
                if ((ixn <= n2) && (iyn >= -n2) && (iyn <= n2) && (izn >= -n2) && (izn <= n2))
                {
                    int ixf, iyf, izf;
                    if (ixn >= 0)
                    {
                        int iza, iya;
                        if (izn >= 0)
                            iza = izn + 1;
                        else
                            iza = n + izn + 1;

                        if (iyn >= 0)
                            iya = iyn + 1;
                        else
                            iya = n + iyn + 1;

                        ixf = ixn;
                        iyf = iya;
                        izf = iza;
                    }
                    else
                    {
                        int izt, iyt;
                        if (izn > 0)
                            izt = n - izn + 1;
                        else
                            izt = -izn + 1;

                        if (iyn > 0)
                            iyt = n - iyn + 1;
                        else
                            iyt = -iyn + 1;

                        ixf = -ixn;
                        iyf = iyt;
                        izf = izt;
                    }


                    int pos2 = ixf + (iyf-1)*nx/2 + (izf-1)*ny*nx/2;
                    float ctfv1 = btq.real() * ctf;
                    float ctfv2 = btq.imag() * ctf;
                    float ctf2 = ctf*ctf;

                    point_t p;
                    p.pos2 = pos2;
                    p.real = ctfv1;
                    p.imag = ctfv2;
                    p.ctf2 = ctf2;

                    points.push_back( p );
                }
	    }
        }
    }


    int npoint = points.size();
    std::istream::off_type offset = (m_offsets.size()==0) ? 0 : m_offsets.back();
    offset += npoint*sizeof(point_t);
    m_offsets.push_back( offset );

    *m_txt_of << m_offsets.back() << std::endl;
    m_bin_of->write( (char*)(&points[0]), sizeof(point_t)*npoint );
    checked_delete( padfft );
}

void newfile_store::get_image( int id, EMData* buf )
{
    if( m_offsets.size()==0 )
    {
        ifstream is( m_txt_file.c_str() );
        std::istream::off_type off;
        while( is >> off )
        {
            m_offsets.push_back( off );
        }

        m_bin_if = shared_ptr<std::ifstream>( new ifstream(m_bin_file.c_str(), std::ios::in|std::ios::binary) );
    }

    Assert( m_bin_if != NULL );

    std::istream::off_type offset = (id==0) ? 0 : m_offsets[id-1];
    Assert( offset >= 0 );
    m_bin_if->seekg( offset, std::ios::beg );


    if( m_bin_if->bad() || m_bin_if->fail() || m_bin_if->eof() )
    {
        std::cout << "bad or fail or eof while fetching id, offset: " << id << " " << offset << std::endl;
        throw std::logic_error( "bad happen" );
    }

    int bufsize = (m_offsets[id] - offset)/sizeof(float);
    if( buf->get_xsize() != bufsize )
    {
        buf->set_size( bufsize, 1, 1 );
    }

    char* data = (char*)(buf->get_data());
    m_bin_if->read( data, sizeof(float)*bufsize );
    buf->update();
}

void newfile_store::read( int nprj )
{
    if( m_offsets.size()==0 )
    {
        ifstream is( m_txt_file.c_str() );
        std::istream::off_type off;
        while( is >> off )
        {
            m_offsets.push_back( off );
        }
    }

    if( m_bin_if==NULL )
    {
        m_bin_if = shared_ptr< ifstream>( new ifstream(m_bin_file.c_str(), std::ios::in|std::ios::binary) );
    }


    int npoint = m_offsets[0]/sizeof(point_t);
    std::ios::off_type prjsize = m_offsets[0];

    try
    {
        m_points.resize(nprj * npoint);
    }
    catch( std::exception& e )
    {
        std::cout << "Error: " << e.what() << std::endl;
    }

    int ip = 0;
    for( int i=0; i < nprj; ++i )
    {
        m_bin_if->read( (char*)(&m_points[ip]), prjsize );
	    if( m_bin_if->bad() || m_bin_if->fail() || m_bin_if->eof() )
        {
            std::cout << "Error: file hander bad or fail or eof" << std::endl;
            return;
        }
        ip += npoint;
    }
}

void newfile_store::add_tovol( EMData* fftvol, EMData* wgtvol, const vector<int>& mults, int pbegin, int pend )
{
    float* vdata = fftvol->get_data();
    float* wdata = wgtvol->get_data();

    int npoint = m_offsets[0]/sizeof(point_t);
//    Assert( int(mults.size())==nprj );
    Assert( int(m_points.size())== (pend - pbegin)*npoint );

    for( int iprj=pbegin; iprj < pend; ++iprj )
    {
        int m = mults[iprj];
        if( m==0 ) continue;

        int ipt = (iprj-pbegin)*npoint;
        for( int i=0; i < npoint; ++i )
        {
            int pos2 = m_points[ipt].pos2;
            int pos1 = pos2*2;

            wdata[pos2] += m_points[ipt].ctf2*m;
            vdata[pos1] += m_points[ipt].real*m;
            vdata[pos1+1]+= m_points[ipt].imag*m;
            ++ipt;
        }
    }
}

void newfile_store::restart()
{
    m_bin_if = shared_ptr< ifstream>( new ifstream(m_bin_file.c_str(), std::ios::in|std::ios::binary) );
}

file_store::file_store(const string& filename, int npad, int write, bool ctf)
    : m_bin_file(filename + ".bin"),
      m_txt_file(filename + ".txt")
{
    m_ctf = ctf;
    m_prev = -1;
    m_npad = npad;
    m_write = write;
}

file_store::~file_store()
{
}

void file_store::add_image( EMData* emdata, const Transform& tf )
{

    EMData* padfft = padfft_slice( emdata, tf, m_npad );

    float* data = padfft->get_data();

    if( m_write && m_bin_ohandle == NULL )
    {
        m_bin_ohandle = shared_ptr< ofstream >( new ofstream(m_bin_file.c_str(), std::ios::out | std::ios::binary) );
        m_txt_ohandle = shared_ptr< ofstream >( new ofstream(m_txt_file.c_str() ) );
        if( m_ctf )
		 *m_txt_ohandle << "Cs pixel voltage ctf_applied amp_contrast defocus ";

	*m_txt_ohandle << "phi theta psi" << std::endl;
    }

    m_x_out = padfft->get_xsize();
    m_y_out = padfft->get_ysize();
    m_z_out = padfft->get_zsize();
    m_totsize = m_x_out*m_y_out*m_z_out;

    if( m_ctf )
    {
        Ctf* ctf = padfft->get_attr( "ctf" );
        Dict ctf_params = ctf->to_dict();

        m_ctf_applied = padfft->get_attr( "ctf_applied" );

        m_Cs = ctf_params["cs"];
        m_pixel = ctf_params["apix"];
        m_voltage = ctf_params["voltage"];
        m_amp_contrast = ctf_params["ampcont"];
        m_defocuses.push_back( ctf_params["defocus"] );
    }

    Dict params = tf.get_rotation( "spider" );
    float phi = params.get( "phi" );
    float tht = params.get( "theta" );
    float psi = params.get( "psi" );


    m_phis.push_back( phi );
    m_thetas.push_back( tht );
    m_psis.push_back( psi );

    if( m_write )
    {
        m_bin_ohandle->write( (char*)data, sizeof(float)*m_totsize );

        if( m_ctf )
        {
            *m_txt_ohandle << m_Cs << " ";
            *m_txt_ohandle << m_pixel << " ";
            *m_txt_ohandle << m_voltage << " ";
            *m_txt_ohandle << m_ctf_applied << " ";
            *m_txt_ohandle << m_amp_contrast << " ";
            *m_txt_ohandle << m_defocuses.back() << " ";
        }
        *m_txt_ohandle << m_phis.back() << " ";
        *m_txt_ohandle << m_thetas.back() << " ";
        *m_txt_ohandle << m_psis.back() << " ";
        *m_txt_ohandle << m_x_out << " ";
        *m_txt_ohandle << m_y_out << " ";
        *m_txt_ohandle << m_z_out << " ";
        *m_txt_ohandle << m_totsize << std::endl;
    }

    checked_delete(padfft);

}

void file_store::get_image( int id, EMData* padfft )
{

    if( m_phis.size() == 0 ) {
        ifstream m_txt_ifs( m_txt_file.c_str() );

	if( !m_txt_ifs )
	{
            std::cerr << "Error: file " << m_txt_file << " does not exist" << std::endl;
        }

	string line;
	std::getline( m_txt_ifs, line );

        float first, defocus, phi, theta, psi;



        while( m_txt_ifs >> first ) {

            if( m_ctf )
            {
                m_Cs = first;
                m_txt_ifs >> m_pixel >> m_voltage;
                m_txt_ifs >> m_ctf_applied >> m_amp_contrast;
                m_txt_ifs >> defocus >> phi >> theta >> psi;
                m_defocuses.push_back( defocus );
            }
            else
            {
                phi = first;
                m_txt_ifs >> theta >> psi;
            }

            m_txt_ifs >> m_x_out >> m_y_out >> m_z_out >> m_totsize;
            m_phis.push_back( phi );
            m_thetas.push_back( theta );
            m_psis.push_back( psi );
        }
    }

    Assert( m_ihandle != NULL );

    std::istream::off_type offset = id*sizeof(float)*m_totsize;
    Assert( offset >= 0 );

    if( offset > 0 )
    {
        m_ihandle->seekg(offset, std::ios::beg);
    }

    if( m_ihandle->bad() )
    {
        std::cout << "bad while fetching id, offset: " << id << " " << offset << std::endl;
        throw std::logic_error( "bad happen" );
    }

    if( m_ihandle->fail() )
    {
        std::cout << "fail while fetching id, offset, curoff: " << id << " " << offset << std::endl;
        throw std::logic_error( "fail happen" );
    }

    if( m_ihandle->eof() )
    {
        std::cout << "eof while fetching id, offset: " << id << " " << offset << std::endl;
        throw std::logic_error( "eof happen" );
    }

    if( padfft->get_xsize() != m_x_out ||
        padfft->get_ysize() != m_y_out ||
        padfft->get_zsize() != m_z_out )
    {
        padfft->set_size(m_x_out, m_y_out, m_z_out);
    }

    char* data = (char*)(padfft->get_data());
    m_ihandle->read( data, sizeof(float)*m_totsize );
    padfft->update();

    if( m_ctf )
    {
        padfft->set_attr( "Cs", m_Cs );
        padfft->set_attr( "Pixel_size", m_pixel );
        padfft->set_attr( "voltage", m_voltage );
        padfft->set_attr( "ctf_applied", m_ctf_applied );
        padfft->set_attr( "amp_contrast", m_amp_contrast );
        padfft->set_attr( "defocus", m_defocuses[id] );
    }

    padfft->set_attr( "padffted", 1 );
    padfft->set_attr( "phi", m_phis[id] );
    padfft->set_attr( "theta", m_thetas[id] );
    padfft->set_attr( "psi", m_psis[id] );

}

void file_store::restart( )
{
    if( m_ihandle == NULL )
    {
        m_ihandle = shared_ptr< ifstream >( new ifstream(m_bin_file.c_str(), std::ios::in | std::ios::binary) );
    }

    if( m_ihandle->bad() || m_ihandle->fail() || m_ihandle->eof() )
    {
        m_ihandle->open( m_bin_file.c_str(), std::ios::binary );
    }

    m_ihandle->seekg( 0, std::ios::beg );
}

/* vim: set ts=4 noet: */
