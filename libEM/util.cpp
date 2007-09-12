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

#include "byteorder.h"
#include "emassert.h"
#include "emdata.h"
#include "util.h"
#include "marchingcubes.h"

#include <fcntl.h>
#include <iomanip>
#include <sstream>

#include <algorithm>
// using accumulate, inner_product, transform

#include <sys/types.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>

#ifndef WIN32
	#include <unistd.h>
	#include <sys/param.h>
#else
	#include <ctime>
	#include <io.h>
	#define access _access
	#define F_OK 00
#endif  //WIN32

using namespace std;
using namespace boost;
using namespace EMAN;

void Util::ap2ri(float *data, size_t n)
{
	Assert(n > 0);

	if (!data) {
		throw NullPointerException("pixel data array");
	}

	for (size_t i = 0; i < n; i += 2) {
		float f = data[i] * sin(data[i + 1]);
		data[i] = data[i] * cos(data[i + 1]);
		data[i + 1] = f;
	}
}

void Util::flip_complex_phase(float *data, size_t n)
{
	Assert(n > 0);

	if (!data) {
		throw NullPointerException("pixel data array");
	}

	for (size_t i = 0; i < n; i += 2) {
		data[i + 1] *= -1;
	}
}

int Util::file_lock_wait(FILE * file)
{
#ifdef WIN32
	return 1;
#else

	if (!file) {
		throw NullPointerException("Tried to lock NULL file");
	}

	int fdes = fileno(file);

	struct flock fl;
	fl.l_type = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start = 0;
	fl.l_len = 0;
#ifdef WIN32
	fl.l_pid = _getpid();
#else
	fl.l_pid = getpid();
#endif

#if defined __sgi
	fl.l_sysid = getpid();
#endif

	int err = 0;
	if (fcntl(fdes, F_SETLKW, &fl) == -1) {
		LOGERR("file locking error! NFS problem?");

		int i = 0;
		for (i = 0; i < 5; i++) {
			if (fcntl(fdes, F_SETLKW, &fl) != -1) {
				break;
			}
			else {
#ifdef WIN32
				Sleep(1000);
#else
				sleep(1);
#endif

			}
		}
		if (i == 5) {
			LOGERR("Fatal file locking error");
			err = 1;
		}
	}

	return err;
#endif
}



bool Util::check_file_by_magic(const void *first_block, const char *magic)
{
	if (!first_block || !magic) {
		throw NullPointerException("first_block/magic");
	}

	const char *buf = static_cast < const char *>(first_block);

	if (strncmp(buf, magic, strlen(magic)) == 0) {
		return true;
	}
	return false;
}

bool Util::is_file_exist(const string & filename)
{
	if (access(filename.c_str(), F_OK) == 0) {
		return true;
	}
	return false;
}


void Util::flip_image(float *data, size_t nx, size_t ny)
{
	if (!data) {
		throw NullPointerException("image data array");
	}
	Assert(nx > 0);
	Assert(ny > 0);

	float *buf = new float[nx];
	size_t row_size = nx * sizeof(float);

	for (size_t i = 0; i < ny / 2; i++) {
		memcpy(buf, &data[i * nx], row_size);
		memcpy(&data[i * nx], &data[(ny - 1 - i) * nx], row_size);
		memcpy(&data[(ny - 1 - i) * nx], buf, row_size);
	}

	if( buf )
	{
		delete[]buf;
		buf = 0;
	}
}

bool Util::sstrncmp(const char *s1, const char *s2)
{
	if (!s1 || !s2) {
		throw NullPointerException("Null string");
	}

	if (strncmp(s1, s2, strlen(s2)) == 0) {
		return true;
	}

	return false;
}

string Util::int2str(int n)
{
	char s[32] = {'\0'};
	sprintf(s, "%d", n);
	return string(s);
}

string Util::get_line_from_string(char **slines)
{
	if (!slines || !(*slines)) {
		throw NullPointerException("Null string");
	}

	string result = "";
	char *str = *slines;

	while (*str != '\n' && *str != '\0') {
		result.push_back(*str);
		str++;
	}
	if (*str != '\0') {
		str++;
	}
	*slines = str;

	return result;
}

vector<EMData *> Util::svdcmp(const vector<EMData *> &data,int nvec) {
	int nimg=data.size();
	if (nvec==0) nvec=nimg;
	vector<EMData *> ret(nvec);
	if (nimg==0) return ret;
	int pixels=data[0]->get_xsize()*data[0]->get_ysize()*data[0]->get_zsize();

	// Allocate the working space
	gsl_vector *work=gsl_vector_alloc(nimg);
	gsl_vector *S=gsl_vector_alloc(nimg);
	gsl_matrix *A=gsl_matrix_alloc(pixels,nimg);
	gsl_matrix *V=gsl_matrix_alloc(nimg,nimg);
	gsl_matrix *X=gsl_matrix_alloc(nimg,nimg);

	int im,x,y,z,i;
	for (im=0; im<nimg; im++) {
		for (z=0,i=0; z<data[0]->get_zsize(); z++) {
			for (y=0; y<data[0]->get_ysize(); y++) {
				for (x=0; x<data[0]->get_xsize(); x++,i++) {
					gsl_matrix_set(A,i,im,data[im]->get_value_at(x,y,z));
				}
			}
		}
	}

	// This calculates the SVD
	gsl_linalg_SV_decomp_mod (A,X, V, S, work);

	for (im=0; im<nvec; im++) {
		EMData *a=data[0]->copy_head();
		ret[im]=a;
		for (z=0,i=0; z<data[0]->get_zsize(); z++) {
			for (y=0; y<data[0]->get_ysize(); y++) {
				for (x=0; x<data[0]->get_xsize(); x++,i++) {
					a->set_value_at(x,y,z,gsl_matrix_get(A,i,im));
				}
			}
		}
	}
	return ret;
}

bool Util::get_str_float(const char *s, const char *float_var, float *p_val)
{
	if (!s || !float_var || !p_val) {
		throw NullPointerException("string float");
	}
	size_t n = strlen(float_var);
	if (strncmp(s, float_var, n) == 0) {
		*p_val = (float) atof(&s[n]);
		return true;
	}

	return false;
}

bool Util::get_str_float(const char *s, const char *float_var, float *p_v1, float *p_v2)
{
	if (!s || !float_var || !p_v1 || !p_v2) {
		throw NullPointerException("string float");
	}

	size_t n = strlen(float_var);
	if (strncmp(s, float_var, n) == 0) {
		sscanf(&s[n], "%f,%f", p_v1, p_v2);
		return true;
	}

	return false;
}

bool Util::get_str_float(const char *s, const char *float_var,
						 int *p_v0, float *p_v1, float *p_v2)
{
	if (!s || !float_var || !p_v0 || !p_v1 || !p_v2) {
		throw NullPointerException("string float");
	}

	size_t n = strlen(float_var);
	*p_v0 = 0;
	if (strncmp(s, float_var, n) == 0) {
		if (s[n] == '=') {
			*p_v0 = 2;
			sscanf(&s[n + 1], "%f,%f", p_v1, p_v2);
		}
		else {
			*p_v0 = 1;
		}
		return true;
	}
	return false;
}

bool Util::get_str_int(const char *s, const char *int_var, int *p_val)
{
	if (!s || !int_var || !p_val) {
		throw NullPointerException("string int");
	}

	size_t n = strlen(int_var);
	if (strncmp(s, int_var, n) == 0) {
		*p_val = atoi(&s[n]);
		return true;
	}
	return false;
}

bool Util::get_str_int(const char *s, const char *int_var, int *p_v1, int *p_v2)
{
	if (!s || !int_var || !p_v1 || !p_v2) {
		throw NullPointerException("string int");
	}

	size_t n = strlen(int_var);
	if (strncmp(s, int_var, n) == 0) {
		sscanf(&s[n], "%d,%d", p_v1, p_v2);
		return true;
	}
	return false;
}

bool Util::get_str_int(const char *s, const char *int_var, int *p_v0, int *p_v1, int *p_v2)
{
	if (!s || !int_var || !p_v0 || !p_v1 || !p_v2) {
		throw NullPointerException("string int");
	}

	size_t n = strlen(int_var);
	*p_v0 = 0;
	if (strncmp(s, int_var, n) == 0) {
		if (s[n] == '=') {
			*p_v0 = 2;
			sscanf(&s[n + 1], "%d,%d", p_v1, p_v2);
		}
		else {
			*p_v0 = 1;
		}
		return true;
	}
	return false;
}

string Util::change_filename_ext(const string & old_filename,
								 const string & ext)
{
	Assert(old_filename != "");
	if (ext == "") {
		return old_filename;
	}

	string filename = old_filename;
	size_t dot_pos = filename.rfind(".");
	if (dot_pos != string::npos) {
		filename = filename.substr(0, dot_pos+1);
	}
	else {
		filename = filename + ".";
	}
	filename = filename + ext;
	return filename;
}

string Util::remove_filename_ext(const string& filename)
{
    if (filename == "") {
        return "";
    }

	char *buf = new char[filename.size()+1];
	strcpy(buf, filename.c_str());
	char *old_ext = strrchr(buf, '.');
	if (old_ext) {
		buf[strlen(buf) - strlen(old_ext)] = '\0';
	}
	string result = string(buf);
	if( buf )
	{
		delete [] buf;
		buf = 0;
	}
	return result;
}

string Util::sbasename(const string & filename)
{
    if (filename == "") {
        return "";
    }

	char s = '/';
#ifdef WIN32
	s = '\\';
#endif
	char * c = strrchr(filename.c_str(), s);
    if (!c) {
        return filename;
    }
	else {
		c++;
	}
    return string(c);
}


string Util::get_filename_ext(const string& filename)
{
    if (filename == "") {
        return "";
    }

	string result = "";
	char *ext = strrchr(filename.c_str(), '.');
	if (ext) {
		ext++;
		result = string(ext);
	}
	return result;
}



void Util::calc_least_square_fit(size_t nitems, const float *data_x, const float *data_y,
								 float *slope, float *intercept, bool ignore_zero,float absmax)
{
	Assert(nitems > 0);

	if (!data_x || !data_y || !slope || !intercept) {
		throw NullPointerException("null float pointer");
	}
	double sum = 0;
	double sum_x = 0;
	double sum_y = 0;
	double sum_xx = 0;
	double sum_xy = 0;

	for (size_t i = 0; i < nitems; i++) {
		if ((!ignore_zero || (data_x[i] != 0 && data_y[i] != 0))&&(!absmax ||(data_y[i]<absmax && data_y[i]>-absmax))) {
			double y = data_y[i];
			double x = i;
			if (data_x) {
				x = data_x[i];
			}

			sum_x += x;
			sum_y += y;
			sum_xx += x * x;
			sum_xy += x * y;
			sum++;
		}
	}

	double div = sum * sum_xx - sum_x * sum_x;
	if (div == 0) {
		div = 0.0000001f;
	}

	*intercept = (float) ((sum_xx * sum_y - sum_x * sum_xy) / div);
	*slope = (float) ((sum * sum_xy - sum_x * sum_y) / div);
}

Vec3f Util::calc_bilinear_least_square(const vector<float> &p) {
unsigned int i;

// various sums used in the final solution
double Sx=0,Sy=0,Sxy=0,Sxx=0,Syy=0,Sz=0,Sxz=0,Syz=0,S=0;
for (i=0; i<p.size(); i+=3) {
	Sx+=p[i];
	Sy+=p[i+1];
	Sz+=p[i+2];
	Sxx+=p[i]*p[i];
	Syy+=p[i+1]*p[i+1];
	Sxy+=p[i]*p[i+1];
	S+=1.0;
	Sxz+=p[i]*p[i+2];
	Syz+=p[i+1]*p[i+2];
}
double d=S*Sxy*Sxy - 2*Sx*Sxy*Sy + Sxx*Sy*Sy  + Sx*Sx*Syy - S*Sxx*Syy;

Vec3f ret(0,0,0);

ret[0]=-((Sxy*Sxz*Sy - Sx*Sxz*Syy + Sx*Sxy*Syz - Sxx*Sy*Syz - Sxy*Sxy*Sz +Sxx*Syy*Sz)/d);
ret[1]=-((-Sxz*Sy*Sy  + S*Sxz*Syy - S*Sxy*Syz + Sx*Sy*Syz + Sxy*Sy*Sz -Sx*Syy*Sz) /d);
ret[2]=-((-S*Sxy*Sxz + Sx*Sxz*Sy - Sx*Sx*Syz + S*Sxx*Syz + Sx*Sxy*Sz -Sxx*Sy*Sz) /d);

return ret;
}

void Util::save_data(const vector < float >&x_array, const vector < float >&y_array,
					 const string & filename)
{
	Assert(x_array.size() > 0);
	Assert(y_array.size() > 0);
	Assert(filename != "");

	if (x_array.size() != y_array.size()) {
		LOGERR("array x and array y have different size: %d != %d\n",
			   x_array.size(), y_array.size());
		return;
	}

	FILE *out = fopen(filename.c_str(), "wb");
	if (!out) {
		throw FileAccessException(filename);
	}

	for (size_t i = 0; i < x_array.size(); i++) {
		fprintf(out, "%g\t%g\n", x_array[i], y_array[i]);
	}
	fclose(out);
}

void Util::save_data(float x0, float dx, const vector < float >&y_array,
					 const string & filename)
{
	Assert(dx != 0);
	Assert(y_array.size() > 0);
	Assert(filename != "");

	FILE *out = fopen(filename.c_str(), "wb");
	if (!out) {
		throw FileAccessException(filename);
	}

	for (size_t i = 0; i < y_array.size(); i++) {
		fprintf(out, "%g\t%g\n", x0 + dx * i, y_array[i]);
	}
	fclose(out);
}


void Util::save_data(float x0, float dx, float *y_array,
					 size_t array_size, const string & filename)
{
	Assert(dx > 0);
	Assert(array_size > 0);
	Assert(filename != "");

	if (!y_array) {
		throw NullPointerException("y array");
	}

	FILE *out = fopen(filename.c_str(), "wb");
	if (!out) {
		throw FileAccessException(filename);
	}

	for (size_t i = 0; i < array_size; i++) {
		fprintf(out, "%g\t%g\n", x0 + dx * i, y_array[i]);
	}
	fclose(out);
}


void Util::sort_mat(float *left, float *right, int *leftPerm, int *rightPerm)
// Adapted by PRB from a macro definition posted on SourceForge by evildave
{
	float *pivot ; int *pivotPerm;

	{
		float *pLeft  =   left; int *pLeftPerm  =  leftPerm;
		float *pRight =  right; int *pRightPerm = rightPerm;
		float scratch =  *left; int scratchPerm = *leftPerm;

		while (pLeft < pRight) {
			while ((*pRight > scratch) && (pLeft < pRight)) {
				pRight--; pRightPerm--;
			}
			if (pLeft != pRight) {
				*pLeft = *pRight; *pLeftPerm = *pRightPerm;
				pLeft++; pLeftPerm++;
			}
			while ((*pLeft < scratch) && (pLeft < pRight)) {
				pLeft++; pLeftPerm++;
			}
			if (pLeft != pRight) {
				*pRight = *pLeft; *pRightPerm = *pLeftPerm;
				pRight--; pRightPerm--;
			}
		}
		*pLeft = scratch; *pLeftPerm = scratchPerm;
		pivot = pLeft; pivotPerm= pLeftPerm;
	}
	if (left < pivot) {
		sort_mat(left, pivot - 1,leftPerm,pivotPerm-1);
	}
	if (right > pivot) {
		sort_mat(pivot + 1, right,pivotPerm+1,rightPerm);
	}
}

int Util::get_irand(int lo, int hi)
{

	int r = (int)((0.999999f + hi - lo) * rand() / (RAND_MAX + 1.0f) + lo);
	return r;
}

float Util::get_frand(int lo, int hi)
{
	return get_frand((float)lo, (float)hi);
}

/*the static variable does not work well in Python call for this function*/
float Util::get_frand(float lo, float hi)
{
/*	static bool inited = false;
	if (!inited) {
		std::cout << "(float)Not initiated..." << std::endl;
		srand(time(0));
		inited = true;
	}
*/
	float r = (hi - lo) * rand() / (RAND_MAX + 1.0f) + lo;
	return r;
}

/*the static variable does not work well in Python call for this function*/
float Util::get_frand(double lo, double hi)
{
/*	static bool inited = false;
	if (!inited) {
		std::cout << "(double)Not initiated..." << std::endl;
		srand(time(0));
		inited = true;
	}
*/
	double r = (hi - lo) * rand() / (RAND_MAX + 1.0) + lo;
	return (float)r;
}

float Util::get_gauss_rand(float mean, float sigma)
{
	float x = 0;
	float r = 0;
	bool valid = true;

	while (valid) {
		x = get_frand(-1.0, 1.0);
		float y = get_frand(-1.0, 1.0);
		r = x * x + y * y;

		if (r <= 1.0 && r > 0) {
			valid = false;
		}
	}

	float f = sqrt(-2.0f * log(r) / r);
	float result = x * f * sigma + mean;
	return result;
}

void Util::find_max(const float *data, size_t nitems, float *max_val, int *max_index)
{
	Assert(nitems > 0);

	if (!data || !max_val || !max_index) {
		throw NullPointerException("data/max_val/max_index");
	}
	float max = -FLT_MAX;
	int m = 0;

	for (size_t i = 0; i < nitems; i++) {
		if (data[i] > max) {
			max = data[i];
			m = (int)i;
		}
	}

	*max_val = (float)m;

	if (max_index) {
		*max_index = m;
	}
}

void Util::find_min_and_max(const float *data, size_t nitems,
							float *max_val, float *min_val,
							int *max_index, int *min_index)
{
	Assert(nitems > 0);

	if (!data || !max_val || !min_val || !max_index || !min_index) {
		throw NullPointerException("data/max_val/min_val/max_index/min_index");
	}
	float max = -FLT_MAX;
	float min = FLT_MAX;
	int max_i = 0;
	int min_i = 0;

	for (size_t i = 0; i < nitems; i++) {
		if (data[i] > max) {
			max = data[i];
			max_i = (int)i;
		}
		if (data[i] < min) {
			min = data[i];
			min_i = (int)i;
		}
	}

	*max_val = max;
	*min_val = min;

	if (min_index) {
		*min_index = min_i;
	}

	if (max_index) {
		*max_index = max_i;
	}

}

Dict Util::get_stats( const vector<float>& data )
{
	// Note that this is a heavy STL approach using generic algorithms - some memory could be saved
	// using plain c style code.
	
	if (data.size() == 0) EmptyContainerException("Error, attempting to call get stats on an empty container (vector<double>)");
	
	double sum = accumulate(data.begin(), data.end(), 0.0);
	
	double mean = sum / static_cast<double> (data.size());
	
	double std_dev = 0.0, skewness = 0.0, kurtosis = 0.0;
	
	if (data.size() > 1)
	{
		// read mm is "minus_mean"
		vector<double> data_mm(data.size());
		// read ts as "then squared"
		vector<double> data_mm_ts(data.size());
		
		// Subtract the mean from the data and store it in data_mm
		transform(data.begin(), data.end(), data_mm.begin(), std::bind2nd( std::minus<double>(), mean));
		
		// Get the square of the data minus the mean and store it in data_mm_ts
		transform(data_mm.begin(), data_mm.end(), data_mm.begin(), data_mm_ts.begin(), std::multiplies<double>());
		
		for ( unsigned int i = 0; i < data.size(); ++i )
		{
			cout << data[i] << " " << mean << " " << data_mm[i] << " " << data_mm_ts[i] << endl;
		}
		
		// Get the sum of the squares for the calculation of the standard deviation
		double square_sum = accumulate(data_mm_ts.begin(), data_mm_ts.end(), 0.0);
		
		std_dev = sqrtf(square_sum)/ static_cast<double>(data.size()-1);
		double std_dev_sq = std_dev * std_dev;
		
		double cubic_sum = inner_product(data_mm.begin(), data_mm.end(),data_mm_ts.begin(), 0.0);
		double quartic_sum = inner_product(data_mm_ts.begin(), data_mm_ts.end(),data_mm_ts.begin(), 0.0);
		
		cout << "Sums are " << square_sum << " " << cubic_sum << " " <<  quartic_sum << endl;
		
		// I got these definitions of skewness and kurtosis from
		// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
		skewness = cubic_sum / (static_cast<double>(data.size()-1) * std_dev_sq * std_dev );
		kurtosis = quartic_sum / (static_cast<double>(data.size()-1) * std_dev_sq * std_dev_sq );
		
	}
	
	Dict parms;
	parms["mean"] = mean;
	parms["std_dev"] = std_dev;
	parms["skewness"] = skewness;
	parms["kurtosis"] = kurtosis;
	
	return parms;
}

int Util::calc_best_fft_size(int low)
{
	Assert(low >= 0);

	//array containing valid sizes <1024 for speed
	static char *valid = NULL;

	if (!valid) {
		valid = (char *) calloc(4096, 1);

		for (float i2 = 1; i2 < 12.0; i2 += 1.0) {

			float f1 = pow((float) 2.0, i2);
			for (float i3 = 0; i3 < 8.0; i3 += 1.0) {

				float f2 = pow((float) 3.0, i3);
				for (float i5 = 0; i5 < 6.0; i5 += 1.0) {

					float f3 = pow((float) 5.0, i5);
					for (float i7 = 0; i7 < 5.0; i7 += 1.0) {

						float f = f1 * f2 * f3 * pow((float) 7.0, i7);
						if (f <= 4095.0) {
							int n = (int) f;
							valid[n] = 1;
						}
					}
				}
			}
		}
	}

	for (int i = low; i < 4096; i++) {
		if (valid[i]) {
			return i;
		}
	}

	LOGERR("Sorry, can only find good fft sizes up to 4096 right now.");

	return 1;
}

string Util::get_time_label()
{
	time_t t0 = time(0);
	struct tm *t = localtime(&t0);
	char label[32];
	sprintf(label, "%d/%02d/%04d %d:%02d",
			t->tm_mon + 1, t->tm_mday, t->tm_year + 1900, t->tm_hour, t->tm_min);
	return string(label);
}


void Util::set_log_level(int argc, char *argv[])
{
	if (argc > 1 && strncmp(argv[1], "-v", 2) == 0) {
		char level_str[32];
		strcpy(level_str, argv[1] + 2);
		Log::LogLevel log_level = (Log::LogLevel) atoi(level_str);
		Log::logger()->set_level(log_level);
	}
}

void Util::printMatI3D(MIArray3D& mat, const string str, ostream& out) {
	// Note: Don't need to check if 3-D because 3D is part of
	//       the MIArray3D typedef.
	out << "Printing 3D Integer data: " << str << std::endl;
	const multi_array_types::size_type* sizes = mat.shape();
	int nx = sizes[0], ny = sizes[1], nz = sizes[2];
	const multi_array_types::index* indices = mat.index_bases();
	int bx = indices[0], by = indices[1], bz = indices[2];
	for (int iz = bz; iz < nz+bz; iz++) {
		cout << "(z = " << iz << " slice)" << endl;
		for (int ix = bx; ix < nx+bx; ix++) {
			for (int iy = by; iy < ny+by; iy++) {
				cout << setiosflags(ios::fixed) << setw(5)
					 << mat[ix][iy][iz] << "  ";
			}
			cout << endl;
		}
	}
}
/*
Dict Util::get_isosurface(EMData * image, float surface_value, bool smooth)
{
	if (image->get_ndim() != 3) {
		throw ImageDimensionException("3D only");
	}
	
	MarchingCubes * mc = new MarchingCubes(image, smooth);
	mc->set_surface_value(surface_value);
	
	Dict d;
	if(smooth) {	
		d.put("points", *(mc->get_points()));
		d.put("faces", *(mc->get_faces()));
		d.put("normals", *(mc->get_normalsSm()));
	}
	else {
		d.put("points", *(mc->get_points()));
		d.put("faces", *(mc->get_faces()));
		d.put("normals", *(mc->get_normals()));
	}
	
	delete mc;
	mc = 0;
	
	return d;
}
*/
