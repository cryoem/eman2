/**
 * $Id$
 */

/*
 * Author: David Woolford, 07/25/2007 (woolford@bcm.edu)
 * Copyright (c) 2000-2007 Baylor College of Medicine
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

#include <cstring>
#include "reconstructor_tools.h"

using namespace EMAN;

// Static init
float FourierPixelInserter3D::tolerance = 0.0001f;

template <> Factory < FourierPixelInserter3D >::Factory()
{
	force_add(&FourierInserter3DMode1::NEW);
	force_add(&FourierInserter3DMode2::NEW);
	force_add(&FourierInserter3DMode3::NEW);
	force_add(&FourierInserter3DMode4::NEW);
	force_add(&FourierInserter3DMode5::NEW);
	force_add(&FourierInserter3DMode6::NEW);
	force_add(&FourierInserter3DMode7::NEW);
	force_add(&FourierInserter3DMode8::NEW);
}


template <> Factory < InterpolatedFRC >::Factory()
{
	force_add(&InterpolatedFRCMode7::NEW);
	force_add(&InterpolatedFRCMode6::NEW);
	force_add(&InterpolatedFRCMode5::NEW);
	force_add(&InterpolatedFRCMode4::NEW);
	force_add(&InterpolatedFRCMode3::NEW);
	force_add(&InterpolatedFRCMode2::NEW);
	force_add(&InterpolatedFRCMode1::NEW);
}

InterpolatedFRC::InterpolatedFRC() :
		threed_rdata(0), norm_data(0), nx(0), ny(0), nz(0), nxy(0), x_scale(1.0), y_scale(1.0), z_scale(1.0), bin(1.0), r(0), rn(0), frc(0), frc_norm_rdata(0),frc_norm_dt(0)
{
}

void  InterpolatedFRC::set_params(const Dict & new_params) {

	TypeDict permissable_params = get_param_types();
	for ( Dict::const_iterator it = new_params.begin(); it != new_params.end(); ++it )
	{

		if ( !permissable_params.find_type(it->first) )
		{
			throw InvalidParameterException(it->first);
		}
		params[it->first] = it->second;
	}
	init();
}


void InterpolatedFRC::free_memory()
{
	if ( frc != 0 )
	{
		delete [] frc;
		frc = 0;
	}
	if ( frc_norm_rdata != 0 )
	{
		delete [] frc_norm_rdata;
		frc_norm_rdata = 0;
	}
	if ( frc_norm_dt != 0 )
	{
		delete [] frc_norm_dt;
		frc_norm_dt = 0;
	}
}

void InterpolatedFRC::init()
{
	if ( params.find("rdata") !=  params.end()  )
	{
		threed_rdata = params["rdata"];
		if ( threed_rdata == 0 )
			throw NotExistingObjectException("rdata", "error the rdata pointer was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("rdata", "the rdata pointer was not defined in FourierPixelInserter3D::init");

	if ( params.find("norm") !=  params.end() )
	{
		norm_data = params["norm"];
		if ( norm_data == 0 )
			throw NotExistingObjectException("norm", "error the norm pointer was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("norm", "the norm pointer was not defined in FourierPixelInserter3D::init");

	if ( params.find("nx") != params.end() )
	{
		nx = params["nx"];
		if ( nx == 0 )
			throw NotExistingObjectException("nx", "error nx was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("nx", "nx was not defined in FourierPixelInserter3D::init");

	if ( params.find("ny") != params.end() )
	{
		ny = params["ny"];
		if ( ny == 0 )
			throw NotExistingObjectException("ny", "error ny was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("ny", "ny was not defined in FourierPixelInserter3D::init");

	if ( params.find("nz") != params.end() )
	{
		nz = params["nz"];
		if ( nz == 0 )
			throw NotExistingObjectException("nz", "error nz was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("nz", "nz was not defined in FourierPixelInserter3D::init");

	if ( params.has_key("sampling") ) bin = params["sampling"];

	if ( params.has_key("x_scale") ) x_scale = params["x_scale"];
	if ( params.has_key("y_scale") ) y_scale = params["y_scale"];
	if ( params.has_key("z_scale") ) z_scale = params["z_scale"];

	nxy = nx*ny;


	if ( bin <= 0 ) throw InvalidValueException(bin, "Error: sampling must be greater than 0");

	int max_x = nx, max_y = ny, max_z = nz;

	if ( x_scale != 1.0 ) max_x = (int) x_scale*nx;
	if ( y_scale != 1.0 ) max_y = (int) y_scale*ny;
	if ( z_scale != 1.0 ) max_z = (int) z_scale*nz;

	int max = max_x;
	if ( max_y > max ) max = max_y;
	if ( max_z > max ) max = max_z;

	pixel_radius_max = max/2;
	pixel_radius_max_square = Util::square( (int) pixel_radius_max );

	size = static_cast<int>(pixel_radius_max*bin);
	frc = new float[size];
	frc_norm_rdata = new float[size];
	frc_norm_dt = new float[size];

	off[0] = 0;
	off[1] = 2;
	off[2] = nx;
	off[3] = nx + 2;
	off[4] = nxy;
	off[5] = nxy + 2;
	off[6] = nxy + nx;
	off[7] = nxy + nx + 2;
}

void InterpolatedFRC::reset()
{
	memset(frc, 0, size*sizeof(float));
	memset(frc_norm_rdata, 0, size*sizeof(float));
	memset(frc_norm_dt, 0, size*sizeof(float));

	r = 0.0;
	rn = 0.0;
}


bool InterpolatedFRC::continue_frc_calc_functoid(const float& xx, const float& yy, const float& zz, const float dt[], const InterpolationFunctoid& functoid,  const float& weight )
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	// Have to get radial coordinates - x is fine as it is but the other two need translation
	int yt = y0 - ny/2;
	int zt = z0 - nz/2;

	int radius = (int) (x_scale*x_scale*floor(xx)*floor(xx) + y_scale*y_scale*yt*yt + z_scale*z_scale*zt*zt);
	radius = static_cast<int>(sqrt((float)radius)*bin);

	if ( radius > (size-1) )
	{
		return false;
	}

	int l = x0 - 2;
	if (x0 == 0) {
		l = x0;
	}

	// The reverse interpolated point minus this pixel (mtp)
	float interp_real_mtp = 0.0, interp_comp_mtp = 0.0;

	// The reverse interpolated point
	float interp_real = 0.0, interp_comp = 0.0;

	float weight_sum = 0.0;

	size_t idx;
	float r, gg;
	for (int k = z0 - 1; k <= z0 + 2; k++) {
		for (int j = y0 - 1; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				gg = functoid.operate(r);

				idx = i + j * nx + k * nxy;

				if ( norm_data[idx/2] == 0 )
					continue;
				float norm = 1/norm_data[idx/2];

				interp_real_mtp += (threed_rdata[idx]- weight * gg * dt[0] * norm) * gg;
				interp_comp_mtp += (threed_rdata[idx + 1]- weight * gg * dt[1] * norm) * gg;

				interp_real += threed_rdata[idx] * gg;
				interp_comp += threed_rdata[idx + 1] * gg;

				weight_sum += gg;
			}
		}
	}

	interp_real_mtp /= weight_sum;
	interp_comp_mtp /= weight_sum;

	frc[radius] += interp_real_mtp*dt[0] + interp_comp_mtp*dt[1];

	frc_norm_rdata[radius] += interp_real_mtp*interp_real_mtp + interp_comp_mtp*interp_comp_mtp;

	frc_norm_dt[radius] += dt[0] * dt[0] + dt[1] * dt[1];

#ifdef	_WIN32
	r += (float)_hypot(dt[0], dt[1]);
	rn += (float)_hypot(interp_real/weight_sum, interp_comp/weight_sum);
#else
	r += (float)hypot(dt[0], dt[1]);
	rn += (float)hypot(interp_real/weight_sum, interp_comp/weight_sum);
#endif	//_WIN32

	return true;
}

bool InterpolatedFRCMode7::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	return continue_frc_calc_functoid(xx,yy,zz,dt,InterpolationFunctoidMode7(),weight);
}

bool InterpolatedFRCMode6::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	return continue_frc_calc_functoid(xx,yy,zz,dt,InterpolationFunctoidMode6(),weight);
}

bool InterpolatedFRCMode5::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx/2 - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		return true;
	}

	// Have to get radial coordinates - x is fine as it is but the other two need translation
	int yt = y0 - ny/2;
	int zt = z0 - nz/2;

	int radius = (int) (x_scale*x_scale*x0*x0 + y_scale*y_scale*yt*yt + z_scale*z_scale*zt*zt);
	radius = static_cast<int>(sqrt((float)radius)*bin);

	if ( radius > (size-1) )
	{
		return true;
	}

	x0 *= 2;

	int mx0 = -(int) floor((xx - x0) * 39.0f + 0.5f) - 78;
	int my0 = -(int) floor((yy - y0) * 39.0f + 0.5f) - 78;
	int mz0 = -(int) floor((zz - z0) * 39.0f + 0.5f) - 78;

	int l = 0;
	if (x0 == 0) {
		mx0 += 78;
	}
	else if (x0 == 2) {
		mx0 += 39;
	}
	else {
		l = x0 - 4;
	}
	// The reverse interpolated point minus this pixel (mtp)
	float interp_real_mtp = 0.0, interp_comp_mtp = 0.0;

	// The reverse interpolated point
	float interp_real = 0.0, interp_comp = 0.0;

	float weight_sum = 0.0;

	size_t idx;
	float gg;
	for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
		for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
			for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
				gg = InterpolationFunctoidMode5().operate(mmx,mmy,mmz);

				idx = i + j * nx + k * nxy;

				if ( norm_data[idx/2] == 0 )
					continue;
				float norm = 1/norm_data[idx/2];

				interp_real_mtp += (threed_rdata[idx] - weight * dt[0] * gg * norm) * gg;
				interp_comp_mtp += (threed_rdata[idx+1] - weight * dt[1] * gg * norm ) * gg;

				interp_real += threed_rdata[idx] *gg;
				interp_comp += threed_rdata[idx + 1] * gg;

				weight_sum += gg;
			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;

		x0 = (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		int mx0 = -(int) floor((xx_b - x0) * 39.0f + 0.5f);
		x0 *= 2;

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return false;

		size_t idx;
		float gg, norm;
		for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
			for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
				for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
					gg = InterpolationFunctoidMode5().operate(mmx,mmy,mmz);

					idx = i + j * nx + k * nxy;
					if ( norm_data[idx/2] == 0 )
						continue;
					norm = 1/norm_data[idx/2];

					interp_real_mtp += (threed_rdata[idx] - weight * dt[0] * gg * norm ) * gg;
					interp_comp_mtp += (threed_rdata[idx+1] + weight * dt[1] * gg * norm ) * gg; // note the +, complex conj.

					interp_real += threed_rdata[idx] *gg;
					interp_comp += threed_rdata[idx + 1] * gg;

					weight_sum += gg;
				}
			}
		}
	}

	interp_real_mtp /= weight_sum;
	interp_comp_mtp /= weight_sum;

	frc[radius] += interp_real_mtp*dt[0] + interp_comp_mtp*dt[1];

	frc_norm_rdata[radius] += interp_real_mtp*interp_real_mtp + interp_comp_mtp*interp_comp_mtp;

	frc_norm_dt[radius] +=  dt[0] * dt[0] + dt[1] * dt[1];

#ifdef	_WIN32
	r += (float)_hypot(dt[0], dt[1]);
	rn += (float)_hypot(interp_real/weight_sum, interp_comp/weight_sum);
#else
	r += (float)hypot(dt[0], dt[1]);
	rn += (float)hypot(interp_real/weight_sum, interp_comp/weight_sum);
#endif	//_WIN32

	return true;
}



bool InterpolatedFRCMode4::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	return continue_frc_calc_functoid(xx,yy,zz,dt,InterpolationFunctoidMode4(),weight);
}

bool InterpolatedFRCMode3::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	return continue_frc_calc_functoid(xx,yy,zz,dt,InterpolationFunctoidMode3(),weight);
}

bool InterpolatedFRCMode2::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	if (x0 >= nx/2 - 2 || y0 >= ny - 1 || z0 >= nz - 2) {
		return false;
	}


	// Have to get radial coordinates - x is fine as it is but the other two need translation
	int yt = y0 - ny/2;
	int zt = z0 - nz/2;

	int radius = (int) (x_scale*x_scale*x0*x0 + y_scale*y_scale*yt*yt + z_scale*z_scale*zt*zt);
// 	int radius =  x0*x0 + y_scale*y_scale*yt*yt + z_scale*z_scale*zt*zt;
	radius = static_cast<int>(sqrt((float)radius)*bin);

	if ( radius > (size-1) )
	{
		return false;
	}

	float dx = xx - x0;
	float dy = yy - y0;
	float dz = zz - z0;

	int i = (int) (x0 * 2 + y0 * nx + z0 * nxy);

	g[0] = Util::agauss(1, dx, dy, dz, EMConsts::I2G);
	g[1] = Util::agauss(1, 1 - dx, dy, dz, EMConsts::I2G);
	g[2] = Util::agauss(1, dx, 1 - dy, dz, EMConsts::I2G);
	g[3] = Util::agauss(1, 1 - dx, 1 - dy, dz, EMConsts::I2G);
	g[4] = Util::agauss(1, dx, dy, 1 - dz, EMConsts::I2G);
	g[5] = Util::agauss(1, 1 - dx, dy, 1 - dz, EMConsts::I2G);
	g[6] = Util::agauss(1, dx, 1 - dy, 1 - dz, EMConsts::I2G);
	g[7] = Util::agauss(1, 1 - dx, 1 - dy, 1 - dz, EMConsts::I2G);

	// The reverse interpolated point minus this pixel (mtp)
	float interp_real_mtp = 0.0, interp_comp_mtp = 0.0;

	// The reverse interpolated point
	float interp_real = 0.0, interp_comp = 0.0;
// 	cout << weight << endl;
	float weight_sum = 0.0;
	for (int j = 0; j < 8; j++) {
		int k = i + off[j];
		if ( norm_data[k/2] == 0 )
			continue;
		float norm = 1/norm_data[k/2];

		interp_real_mtp += (threed_rdata[k] - weight * dt[0] * g[j] * norm ) * g[j];
		interp_comp_mtp += (threed_rdata[k+1] - weight * dt[1] * g[j] * norm ) * g[j];

		interp_real += threed_rdata[k] * g[j];
		interp_comp += threed_rdata[k + 1] * g[j];

		weight_sum += g[j];

	}

	if ( weight_sum == 0 ) return false;

	interp_real_mtp /= weight_sum;
	interp_comp_mtp /= weight_sum;

	frc[radius] += interp_real_mtp*dt[0] + interp_comp_mtp*dt[1];

	frc_norm_rdata[radius] += interp_real_mtp*interp_real_mtp + interp_comp_mtp*interp_comp_mtp;

	frc_norm_dt[radius] +=  dt[0] * dt[0] + dt[1] * dt[1];

#ifdef	_WIN32
	r += (float)_hypot(dt[0], dt[1]);
	rn += (float)_hypot(interp_real/weight_sum, interp_comp/weight_sum);
#else
	r += (float)hypot(dt[0], dt[1]);
	rn += (float)hypot(interp_real/weight_sum, interp_comp/weight_sum);
#endif	//_WIN32

	return true;
}


bool InterpolatedFRCMode1::continue_frc_calc(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	int idx = x0 + y0 * nx + z0 * nxy;

	// Have to get radial coordinates - x is fine as it is but the other two need translation
	int yt = y0 - ny/2;
	int zt = z0 - nz/2;

	int radius = (int) (x_scale*x_scale* floor(xx)*floor(xx) + y_scale*y_scale*yt*yt + z_scale*z_scale*zt*zt);
	radius = static_cast<int>(sqrt((float)radius)*bin);

	if ( radius > (size-1) ) return false;

	// The reverse interpolated point minus this pixel (mtp)
	float interp_real_mtp = 0.0, interp_comp_mtp = 0.0;

	// The reverse interpolated point
	float interp_real = 0.0, interp_comp = 0.0;

	if ( norm_data[idx/2] == 0 )
		return false;
	float norm = 1/norm_data[idx/2];

	interp_real_mtp = threed_rdata[idx] - weight * dt[0] * norm;
	interp_comp_mtp = threed_rdata[idx+1] - weight * dt[1] * norm;

	interp_real += threed_rdata[idx];
	interp_comp += threed_rdata[idx + 1];

	frc[radius] += interp_real_mtp*dt[0] + interp_comp_mtp*dt[1];

	frc_norm_rdata[radius] += interp_real_mtp*interp_real_mtp + interp_comp_mtp*interp_comp_mtp;

	frc_norm_dt[radius] +=  dt[0] * dt[0] + dt[1] * dt[1];

#ifdef	_WIN32
	r += (float)_hypot(dt[0], dt[1]);
	rn += (float)_hypot(interp_real, interp_comp);
#else
	r += (float)hypot(dt[0], dt[1]);
	rn += (float)hypot(interp_real, interp_comp);
#endif	//_WIN32

	return true;
}

InterpolatedFRC::QualityScores& InterpolatedFRC::QualityScores::operator=( const QualityScores& that )
{
	if ( &that != this )
	{
		frc_integral = that.frc_integral;
	// 					snr_normed_frc_intergral = that.snr_normed_frc_intergral;
		normed_snr_integral  = that.normed_snr_integral;
		norm = that.norm;
	}
	return *this;
}

InterpolatedFRC::QualityScores InterpolatedFRC::finish(const unsigned int num_particles)
{
	float frc_integral = 0, snr_normed_frc_intergral = 0, normed_snr_integral = 0;

//	int contrib = 0;

// 	float contrib_thresh = 0.01;
	int cutoff = size-1;
	for( ; cutoff >= 0; --cutoff )
	{
		if ( frc[cutoff] !=0 ) break;
	}
	cutoff += 1;

	for( int i = 0; i < cutoff; ++i )
	{
		if ( frc_norm_rdata[i] == 0 || frc_norm_dt[i] == 0 )
			frc[i] = 0;
		else
			frc[i] /= sqrt(frc_norm_rdata[i]*frc_norm_dt[i]);

		frc_integral += frc[i];

		float tmp = frc[i]*frc[i];

		if ( tmp > .999f && tmp < 1.001f )
		{
			tmp = 0.999f;
		}

		// This shouldn't happen and at the moment is for debug only
		if ( tmp > 1.0f )
		{
			cout << " tmp " << tmp << " div by " << (1.0f-tmp) << " equals " << (tmp/(1.0f-tmp));
		}

		float adjusted_ssnr = tmp/((1.0f-tmp)*num_particles);
		normed_snr_integral += adjusted_ssnr;
		snr_normed_frc_intergral += sqrt(adjusted_ssnr/( 1.0f + adjusted_ssnr ));
	}
	frc_integral /= size;
	snr_normed_frc_intergral /= size;
	normed_snr_integral /= size;

	InterpolatedFRC::QualityScores quality_scores;
	quality_scores.set_frc_integral( frc_integral );
	quality_scores.set_snr_normed_frc_integral( snr_normed_frc_intergral );
	quality_scores.set_normed_snr_integral( normed_snr_integral );

	if (rn!=0)
	{
		r = r/rn;
	}
	else r=1.0;

	quality_scores.set_norm(r);

	return quality_scores;
}

void FourierPixelInserter3D::init()
{
	if ( params.find("rdata") !=  params.end()  )
	{
		rdata = params["rdata"];
		if ( rdata == 0 )
			throw NotExistingObjectException("rdata", "error the rdata pointer was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("rdata", "the rdata pointer was not defined in FourierPixelInserter3D::init");

	if ( params.find("norm") !=  params.end() )
	{
		norm = params["norm"];
		if ( norm == 0 )
			throw NotExistingObjectException("norm", "error the norm pointer was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("norm", "the norm pointer was not defined in FourierPixelInserter3D::init");

	if ( params.find("nx") != params.end() )
	{
		nx = params["nx"];
		if ( nx == 0 )
			throw NotExistingObjectException("nx", "error nx was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("nx", "nx was not defined in FourierPixelInserter3D::init");

	if ( params.find("ny") != params.end() )
	{
		ny = params["ny"];
		if ( ny == 0 )
			throw NotExistingObjectException("ny", "error ny was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("ny", "ny was not defined in FourierPixelInserter3D::init");

	if ( params.find("nz") != params.end() )
	{
		nz = params["nz"];
		if ( nz == 0 )
			throw NotExistingObjectException("nz", "error nz was 0 in FourierPixelInserter3D::init");
	}
	else throw NotExistingObjectException("nz", "nz was not defined in FourierPixelInserter3D::init");

	nxy = nx*ny;
}

bool FourierInserter3DMode1::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	int idx = x0 + y0 * nx + z0 * nxy;

	rdata[idx] += dt[0]*weight;
	rdata[idx + 1] +=  dt[1]*weight;
	norm[idx/2] += weight;


	return true;
}

bool FourierInserter3DMode1::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if ( fabs(rdata[x0 + y0 * nx + z0 * nxy]) < tolerance && fabs(rdata[x0 + y0 * nx + z0 * nxy+1] ) < tolerance )
	{
		return true;
	}
	else
	{
		cout << "the voxel value is non zero it is " << rdata[x0 + y0 * nx + z0 * nxy] << " " <<  rdata[x0 + y0 * nx + z0 * nxy+1] << endl;
		return false;
	}
}

void FourierInserter3DMode2::init()
{
	FourierPixelInserter3D::init();
	off[0] = 0;
	off[1] = 2;
	off[2] = nx;
	off[3] = nx + 2;
	off[4] = nxy;
	off[5] = nxy + 2;
	off[6] = nxy + nx;
	off[7] = nxy + nx + 2;
}

bool FourierInserter3DMode2::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[],const float& weight)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	float dx = xx - x0;
	float dy = yy - y0;
	float dz = zz - z0;

	//
	if (x0 >= nx/2 - 2 || y0 >= ny - 1 || z0 >= nz - 1) return false;

	int i = (int) (x0 * 2 + y0 * nx + z0 * nxy);

	g[0] = Util::agauss(1, dx, dy, dz, EMConsts::I2G);
	g[1] = Util::agauss(1, 1 - dx, dy, dz, EMConsts::I2G);
	g[2] = Util::agauss(1, dx, 1 - dy, dz, EMConsts::I2G);
	g[3] = Util::agauss(1, 1 - dx, 1 - dy, dz, EMConsts::I2G);
	g[4] = Util::agauss(1, dx, dy, 1 - dz, EMConsts::I2G);
	g[5] = Util::agauss(1, 1 - dx, dy, 1 - dz, EMConsts::I2G);
	g[6] = Util::agauss(1, dx, 1 - dy, 1 - dz, EMConsts::I2G);
	g[7] = Util::agauss(1, 1 - dx, 1 - dy, 1 - dz, EMConsts::I2G);

	for (int j = 0; j < 8; j++)
	{
		int k = i + off[j];
		float gg = g[j]*weight;
		rdata[k] += gg * dt[0];
		rdata[k + 1] += gg * dt[1];
		norm[k/2] += gg;

	}

	return true;
}

bool FourierInserter3DMode2::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	int i = (int) (x0 * 2 + y0 * nx + z0 * nxy);

	for (int j = 0; j < 8; j++)
	{
		int k = i + off[j];
		if ( fabs(rdata[k]) > tolerance  ) return false;
		if ( fabs(rdata[k + 1]) > tolerance ) return false;
	}
	return true;
}

bool FourierInserter3DMode3::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[],const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		return false;
	}

	int l = x0 - 2;
	if (x0 == 0) {
		l = x0;
	}

	size_t idx;
	float r, gg;
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {
				r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				gg = exp(-r / EMConsts::I3G)*weight;

				idx = i + j * nx + k * nxy;

				rdata[idx] += gg * dt[0];
				rdata[idx + 1] += gg * dt[1];
				norm[idx/2] += gg;

			}
		}
	}
	return true;
}

bool FourierInserter3DMode3::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		// Just continue if this is the case - pixels of in these ranges are no inserted by mode 4
		return true;
	}

	int l = x0 - 2;
	if (x0 == 0) {
		l = x0;
	}

	size_t idx;
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {
				idx = i + j * nx + k * nxy;
				if ( rdata[idx] > tolerance ) return false;
				if ( rdata[idx + 1] > tolerance ) return false;

			}
		}
	}

	return true;
}

bool FourierInserter3DMode4::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[],const float& weight)
{
	int x0 = 2 * (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		return false;
	}

	int l = x0 - 2;
	if (x0 == 0) {
		l = x0;
	}

	float r, gg;
	size_t idx;
	for (int k = z0 - 1; k <= z0 + 2; k++) {
		for (int j = y0 - 1; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				gg = exp(-r / EMConsts::I4G)*weight;

				idx = i + j * nx + k * nxy;

				rdata[idx] += gg * dt[0];
				rdata[idx + 1] +=  gg * dt[1];
				norm[idx/2] += gg;

			}
		}
	}

	return true;
}

bool FourierInserter3DMode4::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		// Just continue if this is the case - pixels of in these ranges are no inserted by mode 4
		return true;
	}

	int l = x0 - 2;
	if (x0 == 0) {
		l = x0;
	}

	size_t idx;
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {
				idx = i + j * nx + k * nxy;
				if ( rdata[idx] > tolerance ) return false;
				if ( rdata[idx + 1] > tolerance ) return false;

			}
		}
	}

	return true;
}

bool FourierInserter3DMode5::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[],const float& weight)
{
	int x0 = (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	int mx0 = -(int) floor((xx - x0) * 39.0f + 0.5f) - 78;
	int my0 = -(int) floor((yy - y0) * 39.0f + 0.5f) - 78;
	int mz0 = -(int) floor((zz - z0) * 39.0f + 0.5f) - 78;

	x0 *= 2;

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		return false;
	}

	int l = 0;
	if (x0 == 0) {
		mx0 += 78;
	}
	else if (x0 == 2) {
		mx0 += 39;
	}
	else {
		l = x0 - 4;
	}

	size_t ii;
	float gg;
	for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
		for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
			for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
				ii = i + j * nx + k * nxy;
				gg = gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000]*weight;

				rdata[ii] += gg * dt[0];
				rdata[ii + 1] += gg * dt[1];
				norm[ii/2] += gg;

			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;

		x0 = (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		int mx0 = -(int) floor((xx_b - x0) * 39.0f + 0.5f);
		x0 *= 2;

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return false;

		size_t ii;
		float gg;
		for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
			for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
				for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
					ii = i + j * nx + k * nxy;
					gg = gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

					rdata[ii] += gg * dt[0];
					rdata[ii + 1] -= gg * dt[1]; // note the -, complex conj.
					norm[ii/2] += gg;

				}
			}
		}
	}

	return true;
}

bool FourierInserter3DMode5::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int x0 = (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	int mx0 = -(int) floor((xx - x0) * 39.0f + 0.5f) - 78;
	int my0 = -(int) floor((yy - y0) * 39.0f + 0.5f) - 78;
	int mz0 = -(int) floor((zz - z0) * 39.0f + 0.5f) - 78;

	x0 *= 2;

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2) {
		return true;
	}

	int l = 0;
	if (x0 == 0) {
		mx0 += 78;
	}
	else if (x0 == 2) {
		mx0 += 39;
	}
	else {
		l = x0 - 4;
	}

	size_t ii;
	for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
		for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
			for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
				ii = i + j * nx + k * nxy;

				if ( fabs(rdata[ii]) > tolerance ) return false;
				if ( fabs(rdata[ii + 1]) > tolerance ) return false;

			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;

		x0 = (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		int mx0 = -(int) floor((xx_b - x0) * 39.0f + 0.5f);
		x0 *= 2;

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return true;

		size_t ii;
		for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
			for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
				for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
					ii = i + j * nx + k * nxy;

					if ( fabs(rdata[ii]) > tolerance ) return false;
					if ( fabs(rdata[ii + 1]) > tolerance ) return false;
				}
			}
		}
	}

	return true;
}

bool FourierInserter3DMode6::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[],const float& weight)
{
	int	x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
		return false;

	int l = x0 - 4;
	if (x0 <= 2) {
		l = 0;
	}

	size_t ii;
	float r, gg;
	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				ii = i + j * nx + k * nxy;
				r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				gg = exp(-r / EMConsts::I5G)*weight;

				rdata[ii] += gg * dt[0];
				rdata[ii + 1] += gg * dt[1];
				norm[ii/2] += gg;

			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;

		x0 = 2 * (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return false;

		size_t ii;
		float r, gg;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					ii = i + j * nx + k * nxy;
					r = Util::hypot3((float) i / 2 - xx_b, (float) j - yy_b,
											(float) k - zz_b);
					gg = exp(-r / EMConsts::I5G)*weight;

					rdata[ii] += gg * dt[0];
					rdata[ii + 1] -= gg * dt[1]; // note the -, complex conj
					norm[ii/2] += gg;

				}
			}
		}
	}

	return true;
}

bool FourierInserter3DMode6::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int	x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
		return true;

	int l = x0 - 4;
	if (x0 <= 2) {
		l = 0;
	}

	size_t ii;
	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				ii = i + j * nx + k * nxy;

				if ( fabs(rdata[ii]) > tolerance ) return false;
				if ( fabs(rdata[ii + 1]) > tolerance ) return false;
			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;

		x0 = 2 * (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return true;

		size_t ii;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					ii = i + j * nx + k * nxy;

					if ( fabs(rdata[ii]) > tolerance ) return false;
					if ( fabs(rdata[ii + 1]) > tolerance ) return false;
				}
			}
		}
	}

	return true;
}

bool FourierInserter3DMode7::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[],const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
		return false;

	int l = x0 - 4;
	if (x0 <= 2)
		l = 0;

	size_t ii;
	float r, gg;
	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				ii = i + j * nx + k * nxy;
				r =	sqrt(Util::hypot3((float) i / 2 - xx, (float) j - yy, (float) k - zz));
				gg = Interp::hyperg(r)*weight;

				rdata[ii] += gg * dt[0];
				rdata[ii + 1] += gg * dt[1];
				norm[ii/2] += gg;

			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;
		x0 = 2 * (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return false;

		size_t ii;
		float r, gg;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					ii = i + j * nx + k * nxy;
					r = sqrt(Util::hypot3((float) i / 2 - xx_b, (float) j - yy_b,
								   (float) k - zz_b));
					gg = Interp::hyperg(r)*weight;

					rdata[ii] += gg * dt[0];
					rdata[ii + 1] -= gg * dt[1];// note the -, complex conj
					norm[ii/2] += gg;

				}
			}
		}
	}

	return true;
}

bool FourierInserter3DMode7::effected_pixels_are_zero(const float& xx, const float& yy, const float& zz)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
		return true;

	int l = x0 - 4;
	if (x0 <= 2)
		l = 0;

	size_t ii;
	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				ii = i + j * nx + k * nxy;

				if ( fabs(rdata[ii]) > tolerance ) return false;
				if ( fabs(rdata[ii + 1]) > tolerance ) return false;
			}
		}
	}

	if (x0 <= 2) {
		float xx_b = -xx;
		float yy_b = -(yy - ny / 2) + ny / 2;
		float zz_b = -(zz - nz / 2) + nz / 2;
		x0 = 2 * (int) floor(xx_b + 0.5f);
		y0 = (int) floor(yy_b + 0.5f);
		z0 = (int) floor(zz_b + 0.5f);

		if (y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
			return true;

		size_t ii;
		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					ii = i + j * nx + k * nxy;

					if ( fabs(rdata[ii]) > tolerance ) return false;
					if ( fabs(rdata[ii + 1]) > tolerance ) return false;
				}
			}
		}
	}

	return true;
}

void FourierInserter3DMode8::init()
{
	FourierPixelInserter3D::init();
	int P = (int)((1.0+0.25)*nx+1);
	float r = (float)(nx+1)/(float)P;
	mFreqCutoff = 2;
	mDFreq = 0.2f;
	if (W != 0) delete [] W;
	W = Util::getBaldwinGridWeights(mFreqCutoff, (float)P, r,mDFreq,0.5f,0.2f);

}
bool FourierInserter3DMode8::insert_pixel(const float& qx, const float& qy, const float& qz, const float fq[],const float& weight)
{
	int x0 = (int) floor(qx);
	int y0 = (int) floor(qy);
	int z0 = (int) floor(qz);

	int sizeW = (int)(1+2*mFreqCutoff/mDFreq);
	int sizeWmid = sizeW/2;

	for (int z = z0-mFreqCutoff; z < z0+mFreqCutoff; ++z){
		for (int y = y0-mFreqCutoff; y < y0+mFreqCutoff; ++y){
			for (int x = x0-mFreqCutoff; x < x0+mFreqCutoff; ++x){
				if ( x < 0 || x >= nx ) continue;
				if ( y < 0 || y >= ny ) continue;
				if ( z < 0 || z >= nz ) continue;
				float dist = (float)((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
				dist = sqrtf(dist);
				// We enforce a spherical kernel
				if ( dist > mFreqCutoff ) continue;
				int idx = (int)(sizeWmid + dist/mDFreq);
				if (idx >= sizeW) throw;
				float residual = dist/mDFreq - (int)(dist/mDFreq);
				if ( fabs(residual) > 1) throw;

				float factor = W[idx]*(1.0f-residual)+W[idx+1]*residual*weight;

				size_t k = z*nxy + y*nx + 2*x;

// 				float c = Util::agauss(1, x-x0,y-y0,z-z0, EMConsts::I2G);
				rdata[k] += fq[0]*factor;
				rdata[k+1] += fq[1]*factor;


				norm[k/2] += weight;

			}
		}
	}

	return true;
}
