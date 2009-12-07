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
#include "randnum.h"

#include <fcntl.h>
#include <iomanip>
#include <sstream>

#include <cstring>

#include <sys/types.h>
#include <gsl/gsl_linalg.h>
#include <algorithm> // using accumulate, inner_product, transform

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

void Util::rotate_phase_origin(float *data, size_t nx, size_t ny, size_t nz)
{
	if(ny==1 && nz==1) {	//1D, do nothing
		return;
	}
	else if(ny!=1 && nz==1) {	//2D, rotate vertically by ny/2
		size_t i, j, k, l;
		float re;
		l=ny/2*nx;
		for (i=0; i<ny/2; i++) {
			for (j=0; j<nx; j++) {
				k=j+i*nx;
				re=data[k];
				data[k]=data[k+l];
				data[k+l]=re;
			}
		}
	}
	else {	//3D, in the y,z plane, swaps quadrants I,III and II,IV, this is the 'rotation' in y and z
		size_t i, j, k, l, ii, jj;
		char * t=(char *)malloc(sizeof(float)*nx);

		k=nx*ny*(nz+1)/2;
		l=nx*ny*(nz-1)/2;
		jj=nx*sizeof(float);
		for (j=ii=0; j<nz/2; ++j) {
			for (i=0; i<ny; ++i,ii+=nx) {
				memcpy(t,data+ii,jj);
				if (i<ny/2) {
					memcpy(data+ii,data+ii+k,jj);
					memcpy(data+ii+k,t,jj);
				}
				else {
					memcpy(data+ii,data+ii+l,jj);
					memcpy(data+ii+l,t,jj);
				}
			}
		}
		free(t);
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

string Util::str_to_lower(const string& s) {
	string ret(s);
	std::transform(s.begin(),s.end(),ret.begin(), (int (*)(int) ) std::tolower);
	return ret;
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
					a->set_value_at(x,y,z,static_cast<float>(gsl_matrix_get(A,i,im)));
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
#ifdef _WIN32
	s = '\\';
#endif
	const char * c = strrchr(filename.c_str(), s);
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
	const char *ext = strrchr(filename.c_str(), '.');
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

ret[0]=static_cast<float>(-((Sxy*Sxz*Sy - Sx*Sxz*Syy + Sx*Sxy*Syz - Sxx*Sy*Syz - Sxy*Sxy*Sz +Sxx*Syy*Sz)/d));
ret[1]=static_cast<float>(-((-Sxz*Sy*Sy  + S*Sxz*Syy - S*Sxy*Syz + Sx*Sy*Syz + Sxy*Sy*Sz -Sx*Syy*Sz) /d));
ret[2]=static_cast<float>(-((-S*Sxy*Sxz + Sx*Sxz*Sy - Sx*Sx*Syz + S*Sxx*Syz + Sx*Sxy*Sz -Sxx*Sy*Sz) /d));

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

void Util::set_randnum_seed(unsigned long int seed)
{
	Randnum* randnum = Randnum::Instance();
	randnum->set_seed(seed);
}

unsigned long int Util::get_randnum_seed()
{
	Randnum* randnum = Randnum::Instance();
	return	randnum->get_seed();
}

int Util::get_irand(int lo, int hi)
{
	Randnum* randnum = Randnum::Instance();
	return randnum->get_irand(lo, hi);
}

float Util::get_frand(int lo, int hi)
{
	return get_frand((float)lo, (float)hi);
}

float Util::get_frand(float lo, float hi)
{
	Randnum* randnum = Randnum::Instance();
	return randnum->get_frand(lo, hi);
}

float Util::get_frand(double lo, double hi)
{
	Randnum* randnum = Randnum::Instance();
	return randnum->get_frand(lo, hi);
}

float Util::hypot_fast(int x, int y)
{
static float *mem = (float *)malloc(4*128*128);
static int dim = 0;
x=abs(x);
y=abs(y);

if (x>=dim || y>=dim) {
	if (x>4095 || y>4095) return hypot((float)x,(float)y);		// We won't cache anything bigger than 4096^2
	dim=dim==0?128:dim*2;
	mem=(float*)realloc(mem,4*dim*dim);
	for (int y=0; y<dim; y++) {
		for (int x=0; x<dim; x++) {
#ifdef	_WIN32
			mem[x+y*dim]=(float)_hypot((float)x,(float)y);
#else
			mem[x+y*dim]=hypot((float)x,(float)y);
#endif
		}
	}
}

return mem[x+y*dim];
}

short Util::hypot_fast_int(int x, int y)
{
static short *mem = (short *)malloc(2*128*128);
static int dim = 0;
x=abs(x);
y=abs(y);

if (x>=dim || y>=dim) {
	if (x>4095 || y>4095) return hypot((float)x,(float)y);		// We won't cache anything bigger than 4096^2
	dim=dim==0?128:dim*2;
	mem=(short*)realloc(mem,2*dim*dim);
	for (int y=0; y<dim; y++) {
		for (int x=0; x<dim; x++) {
#ifdef	_WIN32
			mem[x+y*dim]=(short)Util::round(_hypot((float)x,(float)y));
#else
			mem[x+y*dim]=(short)Util::round(hypot((float)x,(float)y));
#endif
		}
	}
}

return mem[x+y*dim];
}

// Uses a precached table to return a good approximate to exp(x)
// if outside the cached range, uses regular exp
float Util::fast_exp(const float &f) {
static float *mem = (float *)malloc(sizeof(float)*1000);
static bool needinit=true;

if (needinit) {
	needinit=false;
	for (int i=0; i<1000; i++) mem[i]=(float)exp(-i/50.0);
}
if (f>0 || f<-19.98) return exp(f);
int g=(int)(-f*50.0+0.5);

return mem[g];
}

// Uses a precached table to return a good approximate to acos(x)
// tolerates values outside the -1 - 1 domain by clamping to PI,0
float Util::fast_acos(const float &f) {
if (f>=1.0) return 0.0;
if (f<=-1.0) return M_PI;

static float *mem = (float *)malloc(sizeof(float)*2001);
static bool needinit=true;


if (needinit) {
	needinit=false;
	for (int i=0; i<=2000; i++) mem[i]=(float)acos(i/1000.0-1.0);
}
float f2=f*1000.0+1000.0;

int g=(int)(f2+.5);

return mem[g];

// This version interpolates, but is slower
/*int g=(int)f2;
f2-=g;
return mem[g+1]*f2+mem[g]*(1.0-f2);*/
}


float Util::get_gauss_rand(float mean, float sigma)
{
	Randnum* randnum = Randnum::Instance();
	return randnum->get_gauss_rand(mean, sigma);
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
	// using plain c style code, as in get_stats_cstyle below

	if (data.size() == 0) EmptyContainerException("Error, attempting to call get stats on an empty container (vector<double>)");

	double sum = accumulate(data.begin(), data.end(), 0.0);

	double mean = sum / static_cast<double> (data.size());

	double std_dev = 0.0, skewness = 0.0, kurtosis = 0.0;

	if (data.size() > 1)
	{
		// read mm is "minus_mean"
		vector<double> data_mm(data.size());
		// read ts as "then squared"
		vector<double> data_mm_sq(data.size());

		// Subtract the mean from the data and store it in data_mm
		transform(data.begin(), data.end(), data_mm.begin(), std::bind2nd(std::minus<double>(), mean));

		// Get the square of the data minus the mean and store it in data_mm_sq
		transform(data_mm.begin(), data_mm.end(), data_mm.begin(), data_mm_sq.begin(), std::multiplies<double>());

		// Get the sum of the squares for the calculation of the standard deviation
		double square_sum = accumulate(data_mm_sq.begin(), data_mm_sq.end(), 0.0);

		//Calculate teh standard deviation
		std_dev = sqrt(square_sum / static_cast<double>(data.size()-1));
		double std_dev_sq = std_dev * std_dev;

		// The numerator for the skewness fraction, as defined in http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
		double cubic_sum = inner_product(data_mm.begin(), data_mm.end(),data_mm_sq.begin(), 0.0);

		// The numerator for the kurtosis fraction, as defined in http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
		double quartic_sum = inner_product(data_mm_sq.begin(), data_mm_sq.end(),data_mm_sq.begin(), 0.0);

		// Finalize the calculation of the skewness and kurtosis, as defined in
		// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
		skewness = cubic_sum / ((data.size()-1) * std_dev_sq * std_dev );
		kurtosis = quartic_sum / ((data.size()-1) * std_dev_sq * std_dev_sq );

	}

	Dict parms;
	parms["mean"] = mean;
	parms["std_dev"] = std_dev;
	parms["skewness"] = skewness;
	parms["kurtosis"] = kurtosis;

	return parms;
}


Dict Util::get_stats_cstyle( const vector<float>& data )
{
	// Performs the same calculations as in get_stats, but uses a single pass, optimized c approach
	// Should perform better than get_stats

	if (data.size() == 0) EmptyContainerException("Error, attempting to call get stats on an empty container (vector<double>)");

	double square_sum = 0.0, sum = 0.0, cube_sum = 0.0, quart_sum = 0.0;
	for( vector<float>::const_iterator it = data.begin(); it != data.end(); ++it )
	{
		double val = *it;
		double square = val*val;
		quart_sum += square*square;
		cube_sum += square*val;
		square_sum += square;
		sum += val;
	}

	double mean = sum/(double)data.size();

	double std_dev = 0.0, skewness = 0.0, kurtosis = 0.0;

	if (data.size() > 1)
	{
		// The standard deviation is calculated here
		std_dev = sqrt( (square_sum - mean*sum)/(double)(data.size()-1));

		double square_mean = mean*mean;
		double cube_mean = mean*square_mean;

		double square_std_dev = std_dev*std_dev;

		// This is the numerator of the skewness fraction, if you expand the brackets, as defined in
		// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
		double cubic_sum = cube_sum - 3*square_sum*mean + 3*sum*square_mean - cube_mean*data.size();
		// Complete the skewness fraction
		skewness = cubic_sum/((data.size()-1)*square_std_dev*std_dev);

		// This is the numerator of the kurtosis fraction, if you expand the brackets, as defined in
		// http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
		double quartic_sum = quart_sum - 4*cube_sum*mean + 6*square_sum*square_mean - 4*sum*cube_mean  + square_mean*square_mean*data.size();
		// Complete the kurtosis fraction
		kurtosis = quartic_sum /( (data.size()-1)*square_std_dev*square_std_dev);
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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

void printmatrix( gsl_matrix* M, const int n, const int m, const string& message = "")
{
	cout << message << endl;
	for(int i = 0; i < n; ++i ){
		for (int j = 0; j < m; ++j ){
			cout << gsl_matrix_get(M,i,j) << "\t";
		}
		cout << endl;
	}
}

void printvector( gsl_vector* M, const int n, const string& message = "")
{
	cout << message << endl;
	for(int i = 0; i < n; ++i ){
		cout << gsl_vector_get(M,i) << "\t";
	}
	cout << endl;
}

float* Util::getBaldwinGridWeights( const int& freq_cutoff, const float& P, const float& r, const float& dfreq, const float& alpha, const float& beta)
{
	int i = 0;
	int discs = (int)(1+2*freq_cutoff/dfreq);

	float*  W = new float[discs];

	int fc = (int)(2*freq_cutoff + 1);
	gsl_matrix* M = gsl_matrix_calloc(fc,fc);

	gsl_vector* rhs = gsl_vector_calloc(fc);
	cout << i++ << endl;
	for(int k = -freq_cutoff; k <= freq_cutoff; ++k){
		for(int kp = -freq_cutoff; kp <= freq_cutoff; ++kp){
			int kdiff =abs( k-kp);
			int evenoddfac = ( kdiff % 2 == 0 ? 1 : -1);

			if (kdiff !=0){
				float val = sin(M_PI*(float)kdiff*r)/(sin(M_PI*(float)kdiff/(float)P))*(alpha+2.0f*beta*evenoddfac);
				gsl_matrix_set(M,int(k+freq_cutoff),int(kp+freq_cutoff),val);
			}
		}
		gsl_matrix_set(M,int(k+freq_cutoff),int(k+freq_cutoff),r*P* (alpha+2*beta));
		float val = alpha*sin(M_PI*k*r)/(sin(M_PI*(float)k/(float)P));
		if (k!=0) {
			gsl_vector_set(rhs,int(k+freq_cutoff),val);
		}
	}
	printmatrix(M,fc,fc,"M");

	gsl_vector_set(rhs,int(freq_cutoff),alpha*r*P);
	gsl_matrix* V = gsl_matrix_calloc(fc,fc);
	gsl_vector* S = gsl_vector_calloc(fc);
	gsl_vector* soln = gsl_vector_calloc(fc);
	gsl_linalg_SV_decomp(M,V,S,soln);

	gsl_linalg_SV_solve(M, V, S, rhs, soln); // soln now runs from -freq_cutoff to + freq_cutoff
	printvector(soln,fc,"soln");

	// we want to solve for W, which ranges from -freq_cutoff to +freq_cutoff in steps of dfreq                            2
	int Count=0;
	for(float q = (float)(-freq_cutoff); q <= (float)(freq_cutoff); q+= dfreq){
		float temp=0;
		for(int k = -freq_cutoff; k <= freq_cutoff; ++k){
			float dtemp;
			if (q!=k) {
				dtemp=(1/(float) P)* (float)gsl_vector_get(soln,int(k+freq_cutoff))  * sin(M_PI*(q-k))/sin(M_PI*(q-k)/((float) P));
			} else{
				dtemp = (1/(float) P)* (float)gsl_vector_get(soln,int(k+freq_cutoff))  * P;
			}
			temp +=dtemp;
		}
		W[Count]=temp;
		cout << W[Count] << " ";
		Count+=1;
	}
	cout << endl;
	return W;
}


void Util::equation_of_plane(const Vec3f& p1, const Vec3f& p2, const Vec3f& p3, float * plane )
{
	int x=0,y=1,z=2;
	plane[0] = p1[y]*(p2[z]-p3[z])+p2[y]*(p3[z]-p1[z])+p3[y]*(p1[z]-p2[z]);
	plane[1] = p1[z]*(p2[x]-p3[x])+p2[z]*(p3[x]-p1[x])+p3[z]*(p1[x]-p2[x]);
	plane[2] = p1[x]*(p2[y]-p3[y])+p2[x]*(p3[y]-p1[y])+p3[x]*(p1[y]-p2[y]);
	plane[3] = p1[x]*(p2[y]*p3[z]-p3[y]*p2[z])+p2[x]*(p3[y]*p1[z]-p1[y]*p3[z])+p3[x]*(p1[y]*p2[z]-p2[y]*p1[z]);
	plane[3] = -plane[3];
}


bool Util::point_is_in_triangle_2d(const Vec2f& p1, const Vec2f& p2, const Vec2f& p3,const Vec2f& point)
{

	Vec2f u = p2 - p1;
	Vec2f v = p3 - p1;
	Vec2f w = point - p1;

	float udotu = u.dot(u);
	float udotv = u.dot(v);
	float udotw = u.dot(w);
	float vdotv = v.dot(v);
	float vdotw = v.dot(w);

	float d = 1.0f/(udotv*udotv - udotu*vdotv);
	float s = udotv*vdotw - vdotv*udotw;
	s *= d;

	float t = udotv*udotw - udotu*vdotw;
	t *= d;

	// We've done a few multiplications, so detect when there are tiny residuals that may throw off the final
	// decision
	if (fabs(s) < Transform::ERR_LIMIT ) s = 0;
	if (fabs(t) < Transform::ERR_LIMIT ) t = 0;

	if ( fabs((fabs(s)-1.0)) < Transform::ERR_LIMIT ) s = 1;
	if ( fabs((fabs(t)-1.0)) < Transform::ERR_LIMIT ) t = 1;

// 	cout << "s and t are " << s << " " << t << endl;

	// The final decision, if this is true then we've hit the jackpot
	if ( s >= 0 && t >= 0 && (s+t) <= 1 ) return true;
	else return false;
}

bool Util::point_is_in_convex_polygon_2d(const Vec2f& p1, const Vec2f& p2, const Vec2f& p3, const Vec2f& p4,const Vec2f& actual_point)
{

	if (point_is_in_triangle_2d(p1,p2,p4,actual_point)) return true;
	else return point_is_in_triangle_2d(p3,p2,p4,actual_point);
}

/*
Dict Util::get_isosurface(EMData * image, float surface_palue, bool smooth)
{
	if (image->get_ndim() != 3) {
		throw ImageDimensionException("3D only");
	}

	MarchingCubes * mc = new MarchingCubes(image, smooth);
	mc->set_surface_palue(surface_palue);

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
