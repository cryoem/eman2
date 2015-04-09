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

#include "emdata.h"
#include "ctf.h"
#include "portable_fileio.h"
#include "imageio.h"

#include <cstring>
#include <sstream>
using std::stringstream;

#include <iomanip>
using std::setprecision;
using namespace EMAN;

#ifdef EMAN2_USING_CUDA
#include "cuda/cuda_processor.h"
#endif
// EMAN2_USING_CUDA

EMData* EMData::get_fft_amplitude2D()
{
	ENTERFUNC;

//	int ndim = get_ndim();
	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}
	if (nz>1) {
		LOGERR("2D image expected. Input image is 3D");
		throw ImageFormatException("2D odd square complex image"
			" expected Input image is 3D.");
	}

	int nx2 = nx/2;

	EMData *dat = copy_head();

	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float temp=0;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx2; i++) {
			temp = (*this)(2*i,j)*(*this)(2*i,j);
			temp += (*this)(2*i+1,j)*(*this)(2*i+1,j);
			(*dat)(i,j) = std::sqrt(temp);
		}
	}

	dat->update();
	dat->set_complex(false);
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}


EMData* EMData::get_fft_amplitude()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}

	ri2ap();

	int nx2 = nx - 2;
	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float *d = dat->get_data();
	float *data = get_data();
	int ndim = get_ndim();

	size_t idx1, idx2, idx3;
	if (ndim == 3) {
		for (int k = 1; k < nz; ++k) {
			for (int j = 1; j < ny; ++j) {
				for (int i = 0; i < nx2/2; ++i) {
					idx1 = (size_t)k*nx2*ny+j*nx2+nx2/2+i;
					idx2 = (size_t)k*nx*ny+j*nx+2*i;
					idx3 = (size_t)(nz-k)*nx2*ny+(ny-j)*nx2+nx2/2-i;
					d[idx1] = data[idx2];
					d[idx3] = data[idx2];
				}
			}
		}
	}
	else if (ndim == 2) {
		for (int j = 1; j < ny; ++j) {
			for (int i = 0; i < nx2/2; ++i) {
				d[j*nx2+nx2/2+i] = data[j*nx+2*i];
				d[(ny-j)*nx2+nx2/2-i] = data[j*nx+2*i];
			}
		}
	}

	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}


EMData* EMData::get_fft_phase()
{
	ENTERFUNC;

	if (!is_complex()) {
		LOGERR("complex image expected. Input image is real image.");
		throw ImageFormatException("complex image expected. Input image is a real image.");
	}

	ri2ap();

	int nx2 = nx - 2;
	EMData *dat = copy_head();
	dat->set_size(nx2, ny, nz);
	dat->to_zero();

	float *d = dat->get_data();
	float * data = get_data();

	int ndim = get_ndim();
	size_t idx1, idx2, idx3;
	if (ndim == 3) {
		for (int k = 1; k < nz; ++k) {
			for (int j = 1; j < ny; ++j) {
				for (int i = 0; i < nx2/2; ++i) {
					idx1 = (size_t)k*nx2*ny+j*nx2+nx2/2+i;
					idx2 = (size_t)k*nx*ny+j*nx+2*i+1;
					idx3 = (size_t)(nz-k)*nx2*ny+(ny-j)*nx2+nx2/2-i;
					d[idx1] = data[idx2];
					d[idx3] = -data[idx2];
				}
			}
		}
	}
	else if (ndim == 2) {
		for (int j = 1; j < ny; ++j) {
			for (int i = 0; i < nx2/2; ++i) {
				d[j*nx2+nx2/2+i] = data[j*nx+2*i+1];
				d[(ny-j)*nx2+nx2/2-i] = -data[j*nx+2*i+1];
			}
		}
	}

	dat->update();
	dat->set_complex(false);
	if(dat->get_ysize()==1 && dat->get_zsize()==1) {
		dat->set_complex_x(false);
	}
	dat->set_ri(false);

	EXITFUNC;
	return dat;
}

#include <sys/stat.h>

void EMData::write_data(string fsp,size_t loc,const Region* area,const int file_nx, const int file_ny, const int file_nz) {

	if (area) {
		struct stat fileinfo;
		if ( stat(fsp.c_str(),&fileinfo) != 0 ) throw UnexpectedBehaviorException("To write an image using a region the file must already exist and be the correct dimensions");
	}


	FILE *f = 0;
	f=fopen(fsp.c_str(), "rb+");
	if (!f) f=fopen(fsp.c_str(), "wb");
	if (!f) throw FileAccessException(fsp);
	portable_fseek(f,loc,SEEK_SET);
	if (!area) {
		if (fwrite(get_data(),nx*ny,nz*4,f)!=(size_t)(nz*4)) throw FileAccessException(fsp);
	} else {
		int fnx = nx;
		if (file_nx != 0) fnx = file_nx;
		int fny = ny;
		if (file_ny != 0) fny = file_ny;
		int fnz = nz;
		if (file_nz != 0) fnz = file_nz;

		EMUtil::process_region_io(get_data(), f, ImageIO::READ_WRITE,
								  0, 4,fnx,fny,fnz,area);
	}
	fclose(f);
//	printf("WROTE %s %ld %ld\n",fsp.c_str(),loc,nx*ny*nz);
}

void EMData::read_data(string fsp,size_t loc,const Region* area, const int file_nx, const int file_ny, const int file_nz) {
	FILE *f = 0;
	f=fopen(fsp.c_str(), "rb");
	if (!f) throw FileAccessException(fsp);
	int fnx = nx;
	if (file_nx != 0) fnx = file_nx;
	int fny = ny;
	if (file_ny != 0) fny = file_ny;
	int fnz = nz;
	if (file_nz != 0) fnz = file_nz;

	portable_fseek(f,loc,SEEK_SET);
	EMUtil::process_region_io(get_data(), f, ImageIO::READ_ONLY,
							  0, 4,fnx,fny,fnz,area);
//	portable_fseek(f,loc,SEEK_SET);
//	if (fread(get_data(),nx*ny,nz*4,f)!=(size_t)(nz*4)) throw FileAccessException(fsp);
	fclose(f);
}

float EMData::calc_center_density()
{
	ENTERFUNC;

	float center = get_attr("mean");
	float sigma = get_attr("sigma");
	float ds = sigma / 2;
	size_t size = (size_t)nx * ny * nz;
	float *d = get_data();
	float sigma1 = sigma / 20;
	float sigma2 = sigma / 1000;

	while (ds > sigma1) {
		double sum = 0;
		int norm = 0;

		for (size_t i = 0; i < size; ++i) {
			if (fabs(d[i] - center) < ds) {
				sum += d[i];
				norm++;
			}
		}
		if (!norm) {
			break;
		}
		float mean = (float)(sum / norm);
		if (fabs(mean - center) < sigma2) {
			ds *= 0.5f;
		}
		center = mean;
	}
	EXITFUNC;

	return center;
}


float EMData::calc_sigma_diff()
{
	ENTERFUNC;

	float *d = get_data();
	float mean = get_attr("mean");
	float sigma = get_attr("sigma");

	double sum_up = 0;
	double sum_down = 0;
	int nup = 0;
	int ndown = 0;

	size_t size = (size_t)nx * ny * nz;

	for (size_t i = 0; i < size; ++i) {
		if (d[i] > mean) {
			sum_up += Util::square(d[i] - mean);
			nup++;
		}
		else {
			sum_down += Util::square(mean - d[i]);
			ndown++;
		}
	}

	float sigup = std::sqrt((float)sum_up / nup);
	float sigdown = std::sqrt((float)sum_down / ndown);
	float sig_diff = fabs(sigup - sigdown) / sigma;


	EXITFUNC;
	return sig_diff;

}


IntPoint EMData::calc_min_location() const
{
	ENTERFUNC;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	float min = FLT_MAX;
	int min_x = 0;
	int min_y = 0;
	int min_z = 0;
	int nxy = nx * ny;
	float * data = get_data();

	for (int j = 0; j < nz; ++j) {
		size_t cur_z = (size_t)j * nxy;

		for (int k = 0; k < ny; ++k) {
			size_t cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float t = data[l + cur_y];
				if (t < min) {
					min_x = l;
					min_y = k;
					min_z = j;
					min = t;
				}
			}
		}
	}

	return IntPoint(min_x, min_y, min_z);
}


IntPoint EMData::calc_max_location() const
{
	ENTERFUNC;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	float max = -FLT_MAX;
	int max_x = 0;
	int max_y = 0;
	int max_z = 0;
	int nxy = nx * ny;
	float * data = get_data();

	for (int j = 0; j < nz; ++j) {
		size_t cur_z = (size_t)j * nxy;

		for (int k = 0; k < ny; ++k) {
			size_t cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float t = data[l + cur_y];
				if (t > max) {
					max_x = l;
					max_y = k;
					max_z = j;
					max = t;
				}
			}
		}
	}

	EXITFUNC;
	return IntPoint(max_x, max_y, max_z);
}


IntPoint EMData::calc_max_location_wrap(const int maxdx, const int maxdy, const int maxdz)
{
	int maxshiftx = maxdx, maxshifty = maxdy, maxshiftz = maxdz;
	if (maxdx == -1) maxshiftx = get_xsize()/4;
	if (maxdy == -1) maxshifty = get_ysize()/4;
	if (maxdz == -1) maxshiftz = get_zsize()/4;

	float max_value = -FLT_MAX;

	IntPoint peak(0,0,0);
	
#ifdef EMAN2_USING_CUDA //CUDA
	if(EMData::usecuda == 1 && cudarwdata){
	
		CudaPeakInfo* soln = calc_max_location_wrap_cuda(cudarwdata, nx, ny, nz, maxdx, maxdy, maxdz);
		
		peak[0] = soln->px;
		peak[1] = soln->py;
		peak[2] = soln->pz;
		free(soln);
		
//		cout << "x " << peak[0] << " y " << peak[1] << " z " << peak[2] << endl;
		return peak;
	}
#endif
	for (int k = -maxshiftz; k <= maxshiftz; k++) {
		for (int j = -maxshifty; j <= maxshifty; j++) {
			for (int i = -maxshiftx; i <= maxshiftx; i++) {

				float value = get_value_at_wrap(i,j,k);

				if (value > max_value) {
					max_value = value;
					peak[0] = i;
					peak[1] = j;
					peak[2] = k;
				}
			}
		}
	}

	return peak;
}

vector<float> EMData::calc_max_location_wrap_intp(const int maxdx, const int maxdy, const int maxdz)
{
	int maxshiftx = maxdx, maxshifty = maxdy, maxshiftz = maxdz;
	if (maxdx == -1) maxshiftx = get_xsize()/4;
	if (maxdy == -1) maxshifty = get_ysize()/4;
	if (maxdz == -1) maxshiftz = get_zsize()/4;

	float max_value = -FLT_MAX;

	IntPoint peak(0,0,0);

// NOT yet working......
/**
#ifdef EMAN2_USING_CUDA //CUDA
	if(cudarwdata){
	
		CudaPeakInfo* soln = calc_max_location_wrap_cuda(cudarwdata, nx, ny, nz, maxdx, maxdy, maxdz);
		
		peak[0] = soln->px;
		peak[1] = soln->py;
		peak[2] = soln->pz;
		free(soln);
		
//		cout << "x " << peak[0] << " y " << peak[1] << " z " << peak[2] << endl;
		return peak;
	}
#endif
**/
	for (int k = -maxshiftz; k <= maxshiftz; k++) {
		for (int j = -maxshifty; j <= maxshifty; j++) {
			for (int i = -maxshiftx; i <= maxshiftx; i++) {

				float value = get_value_at_wrap(i,j,k);

				if (value > max_value) {
					max_value = value;
					peak[0] = i;
					peak[1] = j;
					peak[2] = k;
				}
			}
		}
	}
	
	// compute the center of mass
	float cmx = 0.0; float cmy = 0.0f; float cmz = 0.0f;
	float sval = 0.0f;
	for (float x = float(peak[0])-2.0f; x <= float(peak[0])+2.0f; x++) {
		for (float y = float(peak[1])-2.0f; y <= float(peak[1])+2.0f; y++) {
			for (float z = float(peak[2])-2.0f; z <= float(peak[2])+2.0f; z++) {
				//Compute center of mass
				float val = get_value_at_wrap(x,y,z);
				cmx += x*val;
				cmy += y*val;
				cmz += z*val;
				sval += val;
			}
		}
	}
	cmx /= sval;
	cmy /= sval;
	cmz /= sval;

	vector<float> mydata;
	mydata.push_back(cmx);
	mydata.push_back(cmy);
	mydata.push_back(cmz);
	mydata.push_back(max_value);

/**
	// I guess I could use GSL, but this is faster....
	float x1 = float(peak[0]) - 1.0f;
	float x2 = float(peak[0]);
	float x3 = float(peak[0]) + 1.0f;
	float y1 = float(peak[1]) - 1.0f;
	float y2 = float(peak[1]);
	float y3 = float(peak[1]) + 1.0f;
	float z1 = float(peak[2]) - 1.0f;
	float z2 = float(peak[2]);
	float z3 = float(peak[2]) + 1.0f;
	
	float yx1 = get_value_at_wrap(x1,y2,z2);
	float yx2 = get_value_at_wrap(x2,y2,z2);
	float yx3 = get_value_at_wrap(x3,y2,z2);
	float yy1 = get_value_at_wrap(x2,y1,z2);
	float yy2 = get_value_at_wrap(x2,y2,z2);
	float yy3 = get_value_at_wrap(x2,y3,z2);
	float yz1 = get_value_at_wrap(x2,y2,z1);
	float yz2 = get_value_at_wrap(x2,y2,z2);
	float yz3 = get_value_at_wrap(x2,y2,z3);

	// Fit peak in X to y = ax^2 + bx +c
	float bx = ((yx1 - yx2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2) - (yx2-yx3))/(-(x2 - x3) + (x1 - x2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2));
	float ax = ((yx1 - yx2) - bx*(x1 - x2))/(x1*x1 - x2*x2);
	//Find minima
	float xintp = -bx/(2*ax);

	// Fit peak in X to y = ax^2 + bx +c
	float by = ((yy1 - yy2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2) - (yy2-yy3))/(-(x2 - x3) + (x1 - x2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2));
	float ay = ((yy1 - yy2) - by*(x1 - x2))/(x1*x1 - x2*x2);
	//Find minima
	float yintp = -by/(2*ay);
	
	// Fit peak in X to y = ax^2 + bx +c
	float bz = ((yz1 - yz2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2) - (yz2-yz3))/(-(x2 - x3) + (x1 - x2)*(x2*x2 - x3*x3)/(x1*x1 - x2*x2));
	float az = ((yz1 - yz2) - bz*(x1 - x2))/(x1*x1 - x2*x2);
	//Find minima
	float zintp = -bz/(2*az);
	
	vector<float> mydata;
	mydata.push_back(xintp);
	mydata.push_back(yintp);
	mydata.push_back(zintp);
	mydata.push_back(max_value);
**/
	return mydata;
}

FloatPoint EMData::calc_center_of_mass(float threshold)
{
	float *data = get_data();

	//float sigma = get_attr("sigma");
	//float mean = get_attr("mean");
	float m = 0.0;

	FloatPoint com(0,0,0);
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			int j2 = nx * j;
			for (int k = 0; k < nz; ++k) {
				size_t l = i + j2 + (size_t)k * nxy;
				if (data[l] >= threshold) {		// threshold out noise (and negative density)
					com[0] += i * data[l];
					com[1] += j * data[l];
					com[2] += k * data[l];
					m += data[l];
				}
			}
		}
	}

	com[0] /= m;
	com[1] /= m;
	com[2] /= m;

	return com;
}


size_t EMData::calc_min_index() const
{
	IntPoint min_location = calc_min_location();
	size_t i = min_location[0] + min_location[1] * nx + (size_t)min_location[2] * nx * ny;
	return i;
}


size_t EMData::calc_max_index() const
{
	IntPoint max_location = calc_max_location();
	size_t i = max_location[0] + max_location[1] * nx + (size_t)max_location[2] * nx * ny;
	return i;
}


vector<Pixel> EMData::calc_highest_locations(float threshold) const 
{
	ENTERFUNC;

	vector<Pixel> result;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	int nxy = nx * ny;
	float * data = get_data();

	for (int j = 0; j < nz; ++j) {
		size_t cur_z = (size_t)j * nxy;

		for (int k = 0; k < ny; ++k) {
			size_t cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float v =data[l + cur_y];
				if (v > threshold) {
					result.push_back(Pixel(l, k, j, v));
				}
			}
		}
	}

	std::sort(result.begin(), result.end());

	EXITFUNC;
	return result;
}

vector<Pixel> EMData::calc_n_highest_locations(int n)
{
	ENTERFUNC;

	vector<Pixel> result;

	int di = 1;
	if (is_complex() && !is_ri()) {
		di = 2;
	}

	// initialize with n elements
	float * data = get_data();
	for ( int i=0; i<n; i++) result.push_back(Pixel(0,0,0,-data[0]));

	int nxy = nx * ny;

	for (int j = 0; j < nz; ++j) {
		size_t cur_z = (size_t)j * nxy;

		for (int k = 0; k < ny; ++k) {
			size_t cur_y = k * nx + cur_z;

			for (int l = 0; l < nx; l += di) {
				float v =data[l + cur_y];
				if (v<result[n-1].value) continue;
				for (vector<Pixel>::iterator i=result.begin(); i<result.end(); i++) {
					if (v>(*i).value) { result.insert(i,Pixel(l, k, j, v)); result.pop_back(); break; }
				}
			}
		}
	}

	EXITFUNC;
	return result;
}

vector<Pixel> EMData::find_pixels_with_value(float val) 
{
	ENTERFUNC;
	
	if ( is_complex() ) throw ImageFormatException("Error - find_pixels_with_value real only");

	vector<Pixel> result;

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				if (get_value_at(i,j,k)==val) result.push_back(Pixel(i,j,k,val));
			}
		}
	}

	EXITFUNC;
	return result;
}

float EMData::get_edge_mean() const
{
	ENTERFUNC;
#ifdef EMAN2_USING_CUDA
	if(EMData::usecuda == 1 && cudarwdata){
		
			return get_edgemean_cuda(cudarwdata, nx, ny, nz);
			
	}
#endif
	int di = 0;
	double edge_sum = 0;
	float edge_mean = 0;
	size_t nxy = nx * ny;
	float * data = get_data();
	if (nz == 1) {
		for (int i = 0, j = (ny - 1) * nx; i < nx; ++i, ++j) {
			edge_sum += data[i] + data[j];
		}
		for (size_t i = 0, j = nx - 1; i < nxy; i += nx, j += nx) {
			edge_sum += data[i] + data[j];
		}
		edge_mean = (float)edge_sum / (nx * 2 + ny * 2);
	}
	else {
		if (nx == ny && nx == nz * 2 - 1) {
			for (size_t j = (nxy * (nz - 1)); j < nxy * nz; ++j, ++di) {
				edge_sum += data[j];
			}
		}
		else {
			for (size_t i = 0, j = (nxy * (nz - 1)); i < nxy; ++i, ++j, ++di) {
				edge_sum += data[i] + data[j];
			}
		}

		int nxy2 = nx * (ny - 1);
		for (int k = 1; k < nz - 1; ++k) {
			size_t k2 = k * nxy;
			size_t k3 = k2 + nxy2;
			for (int i = 0; i < nx; ++i, ++di) {
				edge_sum += data[i + k2] + data[i + k3];
			}
		}
		for (int k = 1; k < nz - 1; ++k) {
			size_t k2 = k * nxy;
			size_t k3 = nx - 1 + k2;
			for (int i = 1; i < ny - 1; ++i, ++di) {
				edge_sum += data[i * nx + k2] + data[i * nx + k3];
			}
		}

		edge_mean = (float)edge_sum / (di * 2);
	}
	EXITFUNC;

	return edge_mean;
}


float EMData::get_circle_mean()
{
	ENTERFUNC;

	static bool busy = false;
	static EMData *mask = 0;

	while (busy);
	busy = true;

	if (!mask || !EMUtil::is_same_size(this, mask)) {
		if (!mask) {
			mask = new EMData();
		}
		mask->set_size(nx, ny, nz);
		mask->to_one();

		float radius = (float)(ny / 2 - 2);
		mask->process_inplace("mask.sharp", Dict("inner_radius", radius - 1,
									   "outer_radius", radius + 1));

	}
	double n = 0,s=0;
	float *d = mask->get_data();
	float * data = get_data();
	size_t size = (size_t)nx*ny*nz;
	for (size_t i = 0; i < size; ++i) {
		if (d[i]) { n+=1.0; s+=data[i]; }
	}


	float result = (float)(s/n);
	busy = false;

	EXITFUNC;
	return result;
}


void EMData::set_ctf(Ctf * new_ctf)
{
	ENTERFUNC;

	vector<float> vctf = new_ctf->to_vector();
	attr_dict["ctf"] = vctf;

	EXITFUNC;
}

Ctf * EMData::get_ctf() const
{
	if(attr_dict.has_key("ctf")) {
		EMAN1Ctf * ctf = new EMAN1Ctf();
		ctf->from_vector(attr_dict["ctf"]);

		return dynamic_cast<Ctf *>(ctf);
	}
	else {
		return 0;
	}
}

#include <iostream>
using std::cout;
using std::endl;

void EMData::set_size(int x, int y, int z, bool noalloc)
{
	ENTERFUNC;

	if (x <= 0) {
		throw InvalidValueException(x, "x size <= 0");
	}
	else if (y <= 0) {
		throw InvalidValueException(y, "y size <= 0");
	}
	else if (z <= 0) {
		throw InvalidValueException(z, "z size <= 0");
	}

#ifdef MEMDEBUG2
	printf("EMDATA sz  %4d    %p (%d,%d,%d)\n",EMData::totalalloc,this,x,y,z);
#endif


	int old_nx = nx;

	size_t size = (size_t)x*y*z*sizeof(float);
	
	if (noalloc) {
		nx = x;
		ny = y;
		nz = z;
		nxy = nx*ny;
		nxyz = (size_t)nx*ny*nz;
		return;
	}
	
	if (rdata != 0) {
		rdata = (float*)EMUtil::em_realloc(rdata,size);
	} else {
		// Just pass on this for a while....see what happens
		rdata = (float*)EMUtil::em_malloc(size);
	}
// 	rdata = static_cast < float *>(realloc(rdata, size));
	if ( rdata == 0 )
	{
		stringstream ss;
		string gigs;
		ss << (float) size/1000000000.0;
		ss >> gigs;
		string message = "Cannot allocate " + gigs + " GB - not enough memory.";
		throw BadAllocException(message);
	}

	nx = x;
	ny = y;
	nz = z;
	nxy = nx*ny;
	nxyz = (size_t)nx*ny*nz;

// once the object is resized, the CUDA need to be updated
#ifdef EMAN2_USING_CUDA
	if(cudarwdata) {
		//cout << "rw free on set size" << endl;
		rw_free();
		rw_alloc();
	}
	if(cudarodata) {
		ro_free();
		ro_alloc();
	}
#endif
// EMAN2_USING_CUDA

	if (old_nx == 0) {
		EMUtil::em_memset(get_data(),0,size);
	}

	if (supp) {
		EMUtil::em_free(supp);
		supp = 0;
	}

	update();
	EXITFUNC;
}

#ifdef EMAN2_USING_CUDA

void EMData::set_size_cuda(int x, int y, int z)
{
	ENTERFUNC;

	if (x <= 0) {
		throw InvalidValueException(x, "x size <= 0");
	}
	else if (y <= 0) {
		throw InvalidValueException(y, "y size <= 0");
	}
	else if (z <= 0) {
		throw InvalidValueException(z, "z size <= 0");
	}

	nx = x;
	ny = y;
	nz = z;

	nxy = nx*ny;
	nxyz = (size_t)nx*ny*nz;

//	gpu_update();

	EXITFUNC;
}

#endif
// EMAN2_USING_CUDA

MArray2D EMData::get_2dview() const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx, ny}};
	MArray2D marray(get_data(), dims, boost::fortran_storage_order());
	return marray;
}


MArray3D EMData::get_3dview() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx, ny, nz}};
	MArray3D marray(get_data(), dims, boost::fortran_storage_order());
	return marray;
}


MCArray2D EMData::get_2dcview() const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx/2, ny}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(get_data());
	MCArray2D marray(cdata, dims, boost::fortran_storage_order());
	return marray;
}


MCArray3D EMData::get_3dcview() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(get_data());
	MCArray3D marray(cdata, dims, boost::fortran_storage_order());
	return marray;
}


MCArray3D* EMData::get_3dcviewptr() const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(get_data());
	MCArray3D* marray = new MCArray3D(cdata, dims,
									  boost::fortran_storage_order());
	return marray;
}


MArray2D EMData::get_2dview(int x0, int y0) const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx, ny}};
	MArray2D marray(get_data(), dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0}};
	marray.reindex(bases);
	return marray;
}


MArray3D EMData::get_3dview(int x0, int y0, int z0) const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx, ny, nz}};
	MArray3D marray(get_data(), dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0, z0}};
	marray.reindex(bases);
	return marray;
}


MCArray2D EMData::get_2dcview(int x0, int y0) const
{
	const int ndims = 2;
	if (get_ndim() != ndims) {
		throw ImageDimensionException("2D only");
	}
	boost::array<std::size_t,ndims> dims = {{nx/2, ny}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(get_data());
	MCArray2D marray(cdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0}};
	marray.reindex(bases);
	return marray;
}


MCArray3D EMData::get_3dcview(int x0, int y0, int z0) const
{
	const int ndims = 3;
	boost::array<std::size_t,ndims> dims = {{nx/2, ny, nz}};
	std::complex<float>* cdata = reinterpret_cast<std::complex<float>*>(get_data());
	MCArray3D marray(cdata, dims, boost::fortran_storage_order());
	boost::array<std::size_t,ndims> bases={{x0, y0, z0}};
	marray.reindex(bases);
	return marray;
}

int greaterthan( const void* p1, const void* p2 )
{
	float*  v1 = (float*) p1;
	float*  v2 = (float*) p2;

	if ( *v1 < *v2 )  return 0;
	else return 1;
}


EMObject EMData::get_attr(const string & key) const
{
	ENTERFUNC;
	
	if ((flags & EMDATA_NEEDUPD) && (key != "is_fftpad") && (key != "xform.align2d")){update_stat();} //this gives a spped up of 7.3% according to e2speedtest
		
	size_t size = (size_t)nx * ny * nz;
	if (key == "kurtosis") {
		float mean = attr_dict["mean"];
		float sigma = attr_dict["sigma"];

		float *data = get_data();
		double kurtosis_sum = 0;

		for (size_t k = 0; k < size; ++k) {
			float t = (data[k] - mean) / sigma;
			float tt = t * t;
			kurtosis_sum += tt * tt;
		}

		float kurtosis = (float)(kurtosis_sum / size - 3.0);
		return kurtosis;
	}
	else if (key == "skewness") {
		float mean = attr_dict["mean"];
		float sigma = attr_dict["sigma"];

		float *data = get_data();
		double skewness_sum = 0;
		for (size_t k = 0; k < size; ++k) {
			float t = (data[k] - mean) / sigma;
			skewness_sum +=  t * t * t;
		}
		float skewness = (float)(skewness_sum / size);
		return skewness;
	}
	else if (key == "median")
	{
		if ( is_complex() ) throw ImageFormatException("Error - can not calculate the median of a complex image");
		size_t n = size;
		float* tmp = new float[n];
		float* d = get_data();
		if (tmp == 0 ) throw BadAllocException("Error - could not create deep copy of image data");
//		for(size_t i=0; i < n; ++i) tmp[i] = d[i]; // should just be a memcpy
		std::copy(d, d+n, tmp);
		qsort(tmp, n, sizeof(float), &greaterthan);
		float median;
		if (n%2==1) median = tmp[n/2];
		else median = (tmp[n/2-1]+tmp[n/2])/2.0f;
		delete [] tmp;
		return median;
	}
	else if (key == "nonzero_median")
	{
		if ( is_complex() ) throw ImageFormatException("Error - can not calculate the median of a complex image");
		vector<float> tmp;
		size_t n = size;
		float* d = get_data();
		for( size_t i = 0; i < n; ++i ) {
			if ( d[i] != 0 ) tmp.push_back(d[i]);
		}
		sort(tmp.begin(), tmp.end());
		unsigned int vsize = tmp.size();
		float median;
		if (vsize%2==1) median = tmp[vsize/2];
		else median = (tmp[vsize/2-1]+tmp[vsize/2])/2.0f;
		return median;
	}
	else if (key == "changecount") return EMObject(changecount);
	else if (key == "nx") {
		return nx;
	}
	else if (key == "ny") {
		return ny;
	}
	else if (key == "nz") {
		return nz;
	}

	if(attr_dict.has_key(key)) {
		return attr_dict[key];
	}
	else {
		throw NotExistingObjectException(key, "The requested key does not exist");
	}

	EXITFUNC;
}

EMObject EMData::get_attr_default(const string & key, const EMObject & em_obj) const
{
	ENTERFUNC;

	if(attr_dict.has_key(key)) {
		return get_attr(key);
	}
	else {
		return em_obj;
	}

	EXITFUNC;
}

Dict EMData::get_attr_dict() const
{
	if(rdata) {
		update_stat();
	}

	Dict tmp=Dict(attr_dict);
	tmp["nx"]=nx;
	tmp["ny"]=ny;
	tmp["nz"]=nz;
	tmp["changecount"]=changecount;

	return tmp;
}

void EMData::set_attr_dict(const Dict & new_dict)
{
	/*set nx, ny nz may resize the image*/
	// This wasn't supposed to 'clip' the image, but just redefine the size --steve
	if( new_dict.has_key("nx") || new_dict.has_key("ny") || new_dict.has_key("nz") ) {
		LOGWARN("Warning: Ignored setting dimension size by modifying attribute!!!");
		const_cast<Dict&>(new_dict).erase("nx");
		const_cast<Dict&>(new_dict).erase("ny");
		const_cast<Dict&>(new_dict).erase("nz");
	}

	vector<string> new_keys = new_dict.keys();
	vector<string>::const_iterator it;
	for(it = new_keys.begin(); it!=new_keys.end(); ++it) {
		this->set_attr(*it, new_dict[*it]);
	}
}

void EMData::set_attr_dict_explicit(const Dict & new_dict)
{
	attr_dict = new_dict;
}

void EMData::del_attr(const string & attr_name)
{
	attr_dict.erase(attr_name);
}

void EMData::del_attr_dict(const vector<string> & del_keys)
{
	vector<string>::const_iterator it;
	for(it=del_keys.begin(); it!=del_keys.end(); ++it) {
		this->del_attr(*it);
	}
}

void EMData::set_attr(const string & key, EMObject val)
{
	/* Ignore dimension attribute. */
	if(key == "nx" || key == "ny" || key == "nz")
	{
		printf("Ignore setting dimension attribute %s. Use set_size if you need resize this EMData object.", key.c_str());
		return;
	}

	if(rdata) {	//skip following for header only image
		/* Ignore 'read only' attribute. */
		if(key == "sigma" ||
			key == "sigma_nonzero" ||
			key == "square_sum" ||
			key == "maximum" ||
			key == "minimum" ||
			key == "mean" ||
			key == "mean_nonzero" )
		{
			LOGWARN("Ignore setting read only attribute %s", key.c_str());
			return;
		}
	}

	attr_dict[key] = val;
}

void EMData::set_attr_python(const string & key, EMObject val)
{
	/* Ignore dimension attribute. */
	if(key == "nx" || key == "ny" || key == "nz")
	{
		printf("Ignore setting dimension attribute %s. Use set_size if you need resize this EMData object.", key.c_str());
		return;
	}

	/* Ignore 'read only' attribute. */
	if(key == "sigma" ||
		  key == "sigma_nonzero" ||
		  key == "square_sum" ||
		  key == "maximum" ||
		  key == "minimum" ||
		  key == "mean" ||
		  key == "mean_nonzero" )
	{
		LOGWARN("Ignore setting read only attribute %s", key.c_str());
		return;
	}

	EMObject::ObjectType argtype = val.get_type();
	if (argtype == EMObject::EMDATA) {
		EMData* e = (EMData*) val;
		e = e->copy();
		EMObject v(e);
		attr_dict[key] = v;
	}
	else if (argtype == EMObject::TRANSFORM) {
		Transform* t = new Transform(*((Transform*) val));
		EMObject v(t);
		attr_dict[key] = v;
		delete t; t=0;
	} else {
		attr_dict[key] = val;
	}

}

void EMData::scale_pixel(float scale) const
{
	attr_dict["apix_x"] = ((float) attr_dict["apix_x"]) * scale;
	attr_dict["apix_y"] = ((float) attr_dict["apix_y"]) * scale;
	attr_dict["apix_z"] = ((float) attr_dict["apix_z"]) * scale;
	if (attr_dict.has_key("ctf")) {
		Ctf *ctf=(Ctf *)attr_dict["ctf"];
		ctf->apix*=scale;
		attr_dict["ctf"]=ctf;
		if(ctf) {delete ctf; ctf=0;}
	}
}

//vector<float> EMData::get_data_pickle() const
std::string EMData::get_data_pickle() const
{
//	vector<float> vf;
//	vf.resize(nx*ny*nz);
//	std::copy(rdata, rdata+nx*ny*nz, vf.begin());

	std::string vf((const char *)get_data(),nx*ny*nz*sizeof(float));

	return vf;
}

//void EMData::set_data_pickle(const vector<float>& vf)
void EMData::set_data_pickle(std::string vf)
{
//	if (rdata) printf("rdata exists\n");
//	rdata = (float *)malloc(nx*ny*nz*sizeof(float));
//	std::copy(vf.begin(), vf.end(), rdata);
	EMUtil::em_memcpy(get_data(),vf.data(),(size_t)nx*ny*nz*sizeof(float));

}

int EMData::get_supp_pickle() const
{
	return 0;
}

void EMData::set_supp_pickle(int)
{
	this->supp = 0;
}

float EMData::get_amplitude_thres(float thres)
{

	if (thres < 0 || thres > 1){
	        LOGERR("threshold bust be between 0 and 1.");
	        throw InvalidValueException(thres, "thres: 0 <= thres <= 1");
	}
		
	EMData * amps = get_fft_amplitude();
	vector<float> ampvector = amps->get_data_as_vector();
	// yes I realize this may be slow if the map is big, but then again this function is only suited for tomo alignments, which if you have a big map will be VERY slow anyways!
	sort (ampvector.begin(), ampvector.end()); 
	int thresidx = int(thres * ampvector.size());
	float thresamp =  ampvector[thresidx];

	return thresamp;
}

vector<Vec3i > find_region(EMData* image,const vector<Vec3i >& coords, const float value, vector<Vec3i >& region)
{
	static vector<Vec3i> two_six_connected;
	if (two_six_connected.size() == 0) {
		for(int i = -1; i <= 1; ++i) {
			for(int j = -1; j <= 1; ++j) {
				for(int  k = -1; k <= 1; ++k) {
					if ( j != 0 || i != 0 || k != 0) {
						two_six_connected.push_back(Vec3i(i,j,k));
					}
				}
			}
		}
	}

	vector<Vec3i> ret;
	for(vector<Vec3i>::const_iterator it = two_six_connected.begin(); it != two_six_connected.end(); ++it ) {
		for(vector<Vec3i>::const_iterator it2 = coords.begin(); it2 != coords.end(); ++it2 ) {
			if  (image->get_value_at((*it2)[0],(*it2)[1],(*it2)[2]) != value) throw;
			Vec3i c = (*it)+(*it2);

			if ( c[0] < 0 || c[0] >= image->get_xsize()) continue;
			if ( c[1] < 0 || c[1] >= image->get_ysize()) continue;
			if ( c[2] < 0 || c[2] >= image->get_zsize()) continue;

			if( image->get_value_at(c[0],c[1],c[2]) == value ) {
				if (find(ret.begin(),ret.end(),c) == ret.end()) {
					if (find(region.begin(),region.end(),c) == region.end()) {
						region.push_back(c);
						ret.push_back(c);
					}
				}
			}
		}
	}
	return ret;
}

vector<Vec3i> EMData::mask_contig_region(const float& value, const Vec3i& seed) {
	Vec3i coord(seed[0],seed[1],seed[2]);
	vector<Vec3i> region;
	region.push_back(coord);
	vector<Vec3i> find_region_input = region;
	while (true) {
		vector<Vec3i> v = find_region(this,find_region_input, value, region);
		if (v.size() == 0 ) break;
		else find_region_input = v;
	}
	return region;
}
