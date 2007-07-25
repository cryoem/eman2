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

#include "reconstructor_tools.h"

using namespace EMAN;

// Static init
float FourierPixelInserter3D::tolerance = 0.0001;

template <> Factory < FourierPixelInserter3D >::Factory()
{
	force_add(&FourierInserter3DMode1::NEW);
	force_add(&FourierInserter3DMode2::NEW);
	force_add(&FourierInserter3DMode3::NEW);
	force_add(&FourierInserter3DMode4::NEW);
	force_add(&FourierInserter3DMode5::NEW);
	force_add(&FourierInserter3DMode6::NEW);
	force_add(&FourierInserter3DMode7::NEW);
}

InterpolatedFRC::InterpolatedFRC(float* const rdata, const int xsize, const int ysize, const int zsize, const float& sampling ) :
		threed_rdata(rdata), nx(xsize), ny(ysize), nz(zsize), nxy(xsize*ysize), bin(sampling)
{
	if ( sampling <= 0 )
	{	
		throw InvalidValueException(sampling, "Error: sampling must be greater than 0");
	}

	pixel_radius_max = ny/2;
	pixel_radius_max_square = Util::square( (int) pixel_radius_max );

	size = static_cast<int>(pixel_radius_max*bin);
	// The parentheses effectively cause initialization by 0 (apparently) - this should be tested because without 0 initialization this objects behavior will be unexpected. 
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

InterpolatedFRC::InterpolatedFRC( const InterpolatedFRC& that ) :
		threed_rdata(that.threed_rdata), nx(that.nx), ny(that.ny), nz(that.nz), nxy(that.nxy), bin(that.bin),
					 size(that.size), pixel_radius_max(that.pixel_radius_max), pixel_radius_max_square(that.pixel_radius_max_square)
{
	frc = new float[size];
	frc_norm_rdata = new float[size];
	frc_norm_dt = new float[size];
	
	// Now copy the data in that
	memcpy(frc,that.frc, size*sizeof(float));
	memcpy(frc_norm_rdata,that.frc_norm_rdata, size*sizeof(float));
	memcpy(frc_norm_dt,that.frc_norm_dt, size*sizeof(float));

	off[0] = 0;
	off[1] = 2;
	off[2] = nx;
	off[3] = nx + 2;
	off[4] = nxy;
	off[5] = nxy + 2;
	off[6] = nxy + nx;
	off[7] = nxy + nx + 2;
}

InterpolatedFRC& InterpolatedFRC::operator=( const InterpolatedFRC& that)
{
	if (this != &that)
	{
		threed_rdata = that.threed_rdata;
		nx = that.nx; ny = that.ny; nz = that.nz; nxy = that.nxy; bin = that.bin;
		size = that.size; pixel_radius_max = that.pixel_radius_max; pixel_radius_max_square = that.pixel_radius_max_square;

		free_memory();

		frc = new float[size];
		frc_norm_rdata = new float[size];
		frc_norm_dt = new float[size];
		
		// Now copy the data in that
		memcpy(frc,that.frc, size*sizeof(float));
		memcpy(frc_norm_rdata,that.frc_norm_rdata, size*sizeof(float));
		memcpy(frc_norm_dt,that.frc_norm_dt, size*sizeof(float));
		
		off[0] = 0;
		off[1] = 2;
		off[2] = nx;
		off[3] = nx + 2;
		off[4] = nxy;
		off[5] = nxy + 2;
		off[6] = nxy + nx;
		off[7] = nxy + nx + 2;
	}

	return *this;
}


void InterpolatedFRC::reset()
{
	memset(frc, 0, size*sizeof(float));
	memset(frc_norm_rdata, 0, size*sizeof(float));
	memset(frc_norm_dt, 0, size*sizeof(float));
}

bool InterpolatedFRC::continue_frc_calc3(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);
	
	int offset = (ny%2 == 0? 1:0);
	int yt = (int) floor(yy) - ny/2 - offset;
	int zt = (int) floor(zz) - nz/2;
	
	int radius =  (int) floor(xx)*(int) floor(xx) + yt*yt + zt*zt;
	radius = static_cast<int>(sqrtf(radius)*bin);
	
	// debug
	if ( radius > (size-1) )
	{
		//cout is debug
		//cout << "radius " << radius << " was greater than or equal to size " << size  << endl;
		return false;
	}
	
	int l = x0 - 2;
	if (x0 == 0) {
		l = x0;
	}
	
	// The reverse interpolated point
	float interp_real = 0, interp_comp = 0;
	
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {
				float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				float gg = exp(-r / EMConsts::I3G);

				interp_real += (threed_rdata[i + j * nx + k * nxy]- weight * gg * dt[0])*gg;
				interp_comp += (threed_rdata[i + j * nx + k * nxy + 1]- weight * gg * dt[1])*gg;
			}
		}
	}
	
	frc[radius] += interp_real*dt[0] + interp_comp*dt[1];

	frc_norm_rdata[radius] += interp_real*interp_real + interp_comp*interp_comp;
	
	frc_norm_dt[radius] +=  dt[0] * dt[0] + dt[1] * dt[1];
	
	return true;
}

bool InterpolatedFRC::continue_frc_calc2(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);
 
	if (x0 > nx - 2 || y0 > ny - 1 || z0 > nz - 1)
	{
		return false;
	}
	
	// Have to get radial coordinates - x is fine as it is but the other two need translation ( after testing it seems like z does not need translation - more testing required)
	int offset = (ny%2 == 0? 1:0);
	int yt = y0 - ny/2 - offset;
	int zt = z0 - nz/2;

	int radius =  x0*x0 + yt*yt + zt*zt;
	radius = static_cast<int>(sqrtf(radius)*bin);

	// debug
	if ( radius > (size-1) )
	{
		//cout is debug
		//cout << "radius " << radius << " was greater than or equal to size " << size  << endl;
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
	
	// The reverse interpolated point
	float interp_real = 0, interp_comp = 0;
	
	for (int j = 0; j < 8; j++) {
		int k = i + off[j];
		interp_real += (threed_rdata[k] - weight * dt[0] * g[j]) * g[j];
		interp_comp += (threed_rdata[k+1] - weight * dt[1] * g[j]) * g[j];
	}
	
// 	if ( radius == 0 )
	
//	cout << "interp was " << interp_real << " " << interp_comp << " actual was " << dt[0] << " " << dt[1] << endl;
	
	frc[radius] += interp_real*dt[0] + interp_comp*dt[1];

	frc_norm_rdata[radius] += interp_real*interp_real + interp_comp*interp_comp;
	
	frc_norm_dt[radius] +=  dt[0] * dt[0] + dt[1] * dt[1];
	
//	cout << "Current values are " << frc[radius] << " " << frc_norm_rdata[radius] << " " << frc_norm_dt[radius] << endl;
	
	return true;
}


bool InterpolatedFRC::continue_frc_calc1(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);
	
	// Have to get radial coordinates - x is fine as it is but the other two need translation ( after testing it seems like z does not need translation - more testing required)
	int offset = (ny%2 == 0? 1:0);
	int yt = (int) floor(yy) - ny/2 - offset;
	int zt = (int) floor(zz) - nz/2;
	
	int radius =  (int) floor(xx)*(int) floor(xx) + yt*yt + zt*zt;
	radius = static_cast<int>(sqrtf(radius)*bin);
	
	// debug
	if ( radius > (size-1) )
	{
		//cout is debug
		//cout << "radius " << radius << " was greater than or equal to size " << size  << endl;
		return false;
	}
	
	// The reverse interpolated point
	float interp_real = 0, interp_comp = 0;
	
	interp_real = threed_rdata[x0 + y0 * nx + z0 * nxy] - weight * dt[0];
	interp_comp = threed_rdata[x0 + y0 * nx + z0 * nxy + 1] - weight * dt[1];
	
	frc[radius] += interp_real*dt[0] + interp_comp*dt[1];

	frc_norm_rdata[radius] += interp_real*interp_real + interp_comp*interp_comp;
	
	frc_norm_dt[radius] +=  dt[0] * dt[0] + dt[1] * dt[1];
	
	return true;
}


QualityScores InterpolatedFRC::finish(const unsigned int num_particles)
{
	float frc_integral = 0, snr_normed_frc_intergral = 0, normed_snr_integral = 0;

	int contrib = 0;

// 	float contrib_thresh = 0.01;

	for( int i = 0; i < size; ++i )
	{
		if ( frc_norm_rdata[i] == 0 || frc_norm_dt[i] == 0 )
			frc[i] = 0;
		else
			frc[i] /= sqrtf(frc_norm_rdata[i]*frc_norm_dt[i]);

// 		if ( frc[i] < contrib_thresh ) continue;
		contrib++;
		// Accumulate the frc integral - atm this is for testing purposes but could change
		frc_integral += frc[i];
		
		float tmp = frc[i]*frc[i];

		if ( tmp > .999 && tmp < 1.001 )
		{
			tmp = 0.999;
		}
		
		// This shouldn't happen and at the moment is for debug only
		if ( tmp > 1 )
		{
			cout << " tmp " << tmp << " div by " << (1.0-tmp) << " equals " << (tmp/(1.0-tmp));
		}
		
		float adjusted_ssnr = tmp/((1.0-tmp)*num_particles);
		normed_snr_integral += adjusted_ssnr;
		snr_normed_frc_intergral += sqrtf(adjusted_ssnr/( 1.0 + adjusted_ssnr ));
	}

	frc_integral /= contrib;
	snr_normed_frc_intergral /= contrib;
	normed_snr_integral /= contrib;

	QualityScores quality_scores;
	quality_scores.set_frc_integral( frc_integral );
	quality_scores.set_snr_normed_frc_integral( snr_normed_frc_intergral );
	quality_scores.set_normed_snr_integral( normed_snr_integral );
	
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
	
	(*pixel_operation)(rdata + x0 + y0 * nx + z0 * nxy, weight * dt[0]);
	(*pixel_operation)(rdata + x0 + y0 * nx + z0 * nxy + 1, weight  * dt[1]);
	(*pixel_operation)(norm + x0 + y0 * nx + z0 * nxy, weight);

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
		cout << "bad " << rdata[x0 + y0 * nx + z0 * nxy] << " " <<  rdata[x0 + y0 * nx + z0 * nxy+1] << endl;
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

bool FourierInserter3DMode2::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = (int) floor(xx);
	int y0 = (int) floor(yy);
	int z0 = (int) floor(zz);

	float dx = xx - x0;
	float dy = yy - y0;
	float dz = zz - z0;

	if (x0 > nx - 2 || y0 > ny - 1 || z0 > nz - 1) {
		return false;
	}

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
		(*pixel_operation)(rdata + k, weight * g[j] * dt[0]);
		(*pixel_operation)(rdata + k + 1, weight * g[j] * dt[1]);
		(*pixel_operation)(norm + k, weight * g[j]);
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

bool FourierInserter3DMode3::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
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

	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {
				float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				float gg = exp(-r / EMConsts::I3G);

				(*pixel_operation)(rdata + i + j * nx + k * nxy, weight * gg * dt[0]);
				(*pixel_operation)(rdata + i + j * nx + k * nxy + 1, weight * gg * dt[1]);
				(*pixel_operation)(norm + i + j * nx + k * nxy, weight * gg);
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
	
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {

				if ( rdata[i + j * nx + k * nxy] > tolerance ) return false;
				if ( rdata[i + j * nx + k * nxy + 1] > tolerance ) return false;

			}
		}
	}
	
	return true;
}

bool FourierInserter3DMode4::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
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

	for (int k = z0 - 1; k <= z0 + 2; k++) {
		for (int j = y0 - 1; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				float gg = exp(-r / EMConsts::I4G);

				(*pixel_operation)(rdata + i + j * nx + k * nxy, weight * gg * dt[0]);
				(*pixel_operation)(rdata + i + j * nx + k * nxy + 1, weight * gg * dt[1]);
				(*pixel_operation)(norm + i + j * nx + k * nxy, weight * gg);
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
	
	for (int k = z0 - 1; k <= z0 + 1; k++) {
		for (int j = y0 - 1; j <= y0 + 1; j++) {
			for (int i = l; i <= x0 + 2; i += 2) {

				if ( rdata[i + j * nx + k * nxy] > tolerance ) return false;
				if ( rdata[i + j * nx + k * nxy + 1] > tolerance ) return false;

			}
		}
	}
	
	return true;
}

bool FourierInserter3DMode5::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
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

	for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
		for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
			for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
				int ii = i + j * nx + k * nxy;
				float gg = weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];

				(*pixel_operation)(rdata + ii, gg * dt[0]);
				(*pixel_operation)(rdata + ii + 1, gg * dt[1]);
				(*pixel_operation)(norm + ii, gg);
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

		for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
			for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
				for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
					int ii = i + j * nx + k * nxy;
					float gg =
							weight * gimx[abs(mmx) + abs(mmy) * 100 + abs(mmz) * 10000];
					(*pixel_operation)(rdata + ii, gg * dt[0]);
					(*other_pixel_operation)(rdata+ii + 1, gg * dt[1]); // note the -, complex conj.
					(*pixel_operation)(norm + ii, gg);
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

	for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
		for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
			for (int i = l, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
				int ii = i + j * nx + k * nxy;

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

		for (int k = z0 - 2, mmz = mz0; k <= z0 + 2; k++, mmz += 39) {
			for (int j = y0 - 2, mmy = my0; j <= y0 + 2; j++, mmy += 39) {
				for (int i = 0, mmx = mx0; i <= x0 + 4; i += 2, mmx += 39) {
					int ii = i + j * nx + k * nxy;
					
					if ( fabs(rdata[ii]) > tolerance ) return false;
					if ( fabs(rdata[ii + 1]) > tolerance ) return false;
				}
			}
		}
	}
	
	return true;
}

bool FourierInserter3DMode6::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
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

	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				int ii = i + j * nx + k * nxy;
				float r = Util::hypot3((float) i / 2 - xx, j - yy, k - zz);
				float gg = weight * exp(-r / EMConsts::I5G);

				(*pixel_operation)(rdata + ii, gg * dt[0]);
				(*pixel_operation)(rdata + ii + 1, gg * dt[1]);
				(*pixel_operation)(norm + ii, gg);
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

		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					int ii = i + j * nx + k * nxy;
					float r = Util::hypot3((float) i / 2 - xx_b, (float) j - yy_b,
											(float) k - zz_b);
					float gg = weight * exp(-r / EMConsts::I5G);

					(*pixel_operation)(rdata + ii, gg * dt[0]);
					(*other_pixel_operation)(rdata+ii + 1, gg * dt[1]);// note the -, complex conj
					(*pixel_operation)(norm + ii, gg);
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

	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				int ii = i + j * nx + k * nxy;

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

		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					int ii = i + j * nx + k * nxy;
					
					if ( fabs(rdata[ii]) > tolerance ) return false;
					if ( fabs(rdata[ii + 1]) > tolerance ) return false;
				}
			}
		}
	}
	
	return true;
}

bool FourierInserter3DMode7::insert_pixel(const float& xx, const float& yy, const float& zz, const float dt[], const float& weight)
{
	int x0 = 2 * (int) floor(xx + 0.5f);
	int y0 = (int) floor(yy + 0.5f);
	int z0 = (int) floor(zz + 0.5f);

	if (x0 >= nx - 4 || y0 > ny - 3 || z0 > nz - 3 || y0 < 2 || z0 < 2)
		return false;

	int l = x0 - 4;
	if (x0 <= 2)
		l = 0;

	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				int ii = i + j * nx + k * nxy;
				float r =
					sqrt(Util::
						hypot3((float) i / 2 - xx, (float) j - yy, (float) k - zz));
				float gg = weight * Interp::hyperg(r);

				(*pixel_operation)(rdata + ii, gg * dt[0]);
				(*pixel_operation)(rdata + ii + 1, gg * dt[1]);
				(*pixel_operation)(norm + ii, gg);
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

		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					int ii = i + j * nx + k * nxy;
					float r = sqrt(Util::hypot3((float) i / 2 - xx_b, (float) j - yy_b,
								   (float) k - zz_b));
					float gg = weight * Interp::hyperg(r);

					(*pixel_operation)(rdata + ii, gg * dt[0]);
					(*other_pixel_operation)(rdata+ii + 1, gg * dt[1]);// note the -, complex conj
					(*pixel_operation)(norm + ii, gg);
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

	for (int k = z0 - 2; k <= z0 + 2; k++) {
		for (int j = y0 - 2; j <= y0 + 2; j++) {
			for (int i = l; i <= x0 + 4; i += 2) {
				int ii = i + j * nx + k * nxy;
				
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

		for (int k = z0 - 2; k <= z0 + 2; k++) {
			for (int j = y0 - 2; j <= y0 + 2; j++) {
				for (int i = 0; i <= x0 + 4; i += 2) {
					int ii = i + j * nx + k * nxy;
					
					if ( fabs(rdata[ii]) > tolerance ) return false;
					if ( fabs(rdata[ii + 1]) > tolerance ) return false;
				}
			}
		}
	}
	
	return true;
}
