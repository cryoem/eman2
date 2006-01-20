/**
 * $Id$
 */
#include "byteorder.h"
#include "Assert.h"
#include "emdata.h"

#include <fcntl.h>
#include <iomanip>
#include <sstream>

#include <gsl/gsl_sf_bessel.h>
#include <sys/types.h>
#include <gsl/gsl_sf_bessel.h>
#include <algorithm>

#ifndef WIN32
	#include <unistd.h>
	#include <sys/param.h>
#else
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


void Util::spline_mat(float *x, float *y, int n,  float *xq, float *yq, int m) //PRB
{

	float x0= x[0];
	float x1= x[1];
	float x2= x[2];
	float y0= y[0];
	float y1= y[1];
	float y2= y[2];
	float yp1 =  (y1-y0)/(x1-x0) +  (y2-y0)/(x2-x0) - (y2-y1)/(x2-x1)  ;
	float xn  = x[n];
	float xnm1= x[n-1];
	float xnm2= x[n-2];
	float yn  = y[n];
	float ynm1= y[n-1];
	float ynm2= y[n-2];
	float ypn=  (yn-ynm1)/(xn-xnm1) +  (yn-ynm2)/(xn-xnm2) - (ynm1-ynm2)/(xnm1-xnm2) ;
	float *y2d = new float[n];
	Util::spline(x,y,n,yp1,ypn,y2d);
	Util::splint(x,y,y2d,n,xq,yq,m); //PRB
	delete [] y2d;
	return;
}


void Util::spline(float *x, float *y, int n, float yp1, float ypn, float *y2) //PRB
{
	int i,k;
	float p, qn, sig, un, *u;
	u=new float[n-1];

	if (yp1 > .99e30){
		y2[0]=u[0]=0.0;
	}else{
		y2[0]=-.5f;
		u[0] =(3.0f/ (x[1] -x[0]))*( (y[1]-y[0])/(x[1]-x[0]) -yp1);
	}

	for (i=1; i < n-1; i++) {
		sig= (x[i] - x[i-1])/(x[i+1] - x[i-1]);
		p = sig*y2[i-1] + 2.0f;
		y2[i]  = (sig-1.0f)/p;
		u[i] = (y[i+1] - y[i] )/(x[i+1]-x[i] ) -  (y[i] - y[i-1] )/(x[i] -x[i-1]);
		u[i] = (6.0f*u[i]/ (x[i+1]-x[i-1]) - sig*u[i-1])/p;
	}

	if (ypn>.99e30){
		qn=0; un=0;
	} else {
		qn= .5f;
		un= (3.0f/(x[n-1] -x[n-2])) * (ypn -  (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]= (un - qn*u[n-2])/(qn*y2[n-2]+1.0f);
	for (k=n-2; k>=0; k--){
		y2[k]=y2[k]*y2[k+1]+u[k];
	}
	delete [] u;
}

void Util::splint( float *xa, float *ya, float *y2a, int n,  float *xq, float *yq, int m) //PRB
{
	int klo, khi, k;
	float h, b, a;

//	klo=0; // can try to put here
	for (int j=0; j<m;j++){
		klo=0;
		khi=n-1;
		while (khi-klo >1) {
			k=(khi+klo) >>1;
			if  (xa[k]>xq[j]){ khi=k;}
			else { klo=k;}
		}
		h=xa[khi]- xa[klo];
		if (h==0.0) printf("Bad XA input to routine SPLINT \n");
		a =(xa[khi]-xq[j])/h;
		b=(xq[j]-xa[klo])/h;
		yq[j]=a*ya[klo] + b*ya[khi]
			+ ((a*a*a-a)*y2a[klo]
			     +(b*b*b-b)*y2a[khi]) *(h*h)/6.0f;
	}
//	printf("h=%f, a = %f, b=%f, ya[klo]=%f, ya[khi]=%f , yq=%f\n",h, a, b, ya[klo], ya[khi],yq[0]);
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

void Util::Radialize(int *PermMatTr, float *kValsSorted,   // PRB
               float *weightofkValsSorted, int Size, int *SizeReturned)
{
	int iMax = (int) floor( (Size-1.0)/2 +.01);
	int CountMax = (iMax+2)*(iMax+1)/2;
	int Count=-1;
	float *kVals     = new float[CountMax];
	float *weightMat = new float[CountMax];
	int *PermMat     = new   int[CountMax];
	SizeReturned[0] = CountMax;

//	printf("Aa \n");	fflush(stdout);
	for (int jkx=0; jkx< iMax+1; jkx++) {
		for (int jky=0; jky< jkx+1; jky++) {
			Count++;
			kVals[Count] = sqrtf((float) (jkx*jkx +jky*jky));
			weightMat[Count]=  1.0;
			if (jkx!=0)  { weightMat[Count] *=2;}
			if (jky!=0)  { weightMat[Count] *=2;}
			if (jkx!=jky){ weightMat[Count] *=2;}
			PermMat[Count]=Count+1;
	}}

	int lkVals = Count+1;
//	printf("Cc \n");fflush(stdout);

	sort_mat(&kVals[0],&kVals[Count],
	     &PermMat[0],  &PermMat[Count]);  //PermMat is
				//also returned as well as kValsSorted
	fflush(stdout);

	int newInd;

        for (int iP=0; iP < lkVals ; iP++ ) {
		newInd =  PermMat[iP];
		PermMatTr[newInd-1] = iP+1;
	}

//	printf("Ee \n"); fflush(stdout);

	int CountA=-1;
	int CountB=-1;

	while (CountB< (CountMax-1)) {
		CountA++;
		CountB++;
//		printf("CountA=%d ; CountB=%d \n", CountA,CountB);fflush(stdout);
		kValsSorted[CountA] = kVals[CountB] ;
		if (CountB<(CountMax-1) ) {
			while (fabs(kVals[CountB] -kVals[CountB+1])<.0000001  ) {
				SizeReturned[0]--;
				for (int iP=0; iP < lkVals; iP++){
//					printf("iP=%d \n", iP);fflush(stdout);
					if  (PermMatTr[iP]>CountA+1) {
						PermMatTr[iP]--;
		    			}
		 		}
				CountB++;
	    		}
		}
	}
	

	for (int CountD=0; CountD < CountMax; CountD++) {
	    newInd = PermMatTr[CountD];
	    weightofkValsSorted[newInd-1] += weightMat[CountD];
        }

}

float Util::get_frand(int lo, int hi)
{
	return get_frand((float)lo, (float)hi);
}

float Util::get_frand(float lo, float hi)
{
	static bool inited = false;
	if (!inited) {
		srand(time(0));
		inited = true;
	}

	float r = (hi - lo) * rand() / (RAND_MAX + 1.0f) + lo;
	return r;
}

float Util::get_frand(double lo, double hi)
{
	static bool inited = false;
	if (!inited) {
		srand(time(0));
		inited = true;
	}

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

vector<float>
Util::voea(float delta, float t1, float t2, float p1, float p2)
{
	vector<float> angles;
	float psi = 0.0;
	if ((0.0 == t1)&&(0.0 == t2)||(t1 >= t2)) {
		t1 = 0.0f;
		t2 = 90.0f;
	}
	if ((0.0 == p1)&&(0.0 == p2)||(p1 >= p2)) {
		p1 = 0.0f;
		p2 = 359.9f;
	}
	bool skip = ((t1 < 90.0)&&(90.0 == t2)&&(0.0 == p1)&&(p2 > 180.0));
	for (float theta = t1; theta <= t2; theta += delta) {
		float detphi;
		int lt;
		if ((0.0 == theta)||(180.0 == theta)) {
			detphi = 360.0f;
			lt = 1;
		} else {
			detphi = delta/sin(theta*static_cast<float>(dgr_to_rad));
			lt = int((p2 - p1)/detphi)-1;
			if (lt < 1) lt = 1;
			detphi = (p2 - p1)/lt;
		}
		for (int i = 0; i < lt; i++) {
			float phi = p1 + i*detphi;
			if (skip&&(90.0 == theta)&&(phi > 180.0)) continue;
			angles.push_back(phi);
			angles.push_back(theta);
			angles.push_back(psi);
		}
	}
	return angles;
}


float Util::triquad(double r, double s, double t, float f[]) {
	const float c2 = 1.0f / 2.0f;
	const float c4 = 1.0f / 4.0f;
	const float c8 = 1.0f / 8.0f;
	float rs = (float)(r*s);
	float st = (float)(s*t);
	float rt = (float)(r*t);
	float rst = (float)(r*st);
	float rsq = (float)(1 - r*r);
	float ssq = (float)(1 - s*s);
	float tsq = (float)(1 - t*t);
	float rm1 = (float)(1 - r);
	float sm1 = (float)(1 - s);
	float tm1 = (float)(1 - t);
	float rp1 = (float)(1 + r);
	float sp1 = (float)(1 + s);
	float tp1 = (float)(1 + t);

	return (float)(
		(-c8) * rst * rm1  * sm1  * tm1 * f[ 0] +
		( c4) * st	* rsq  * sm1  * tm1 * f[ 1] +
		( c8) * rst * rp1  * sm1  * tm1 * f[ 2] +
		( c4) * rt	* rm1  * ssq  * tm1 * f[ 3] +
		(-c2) * t	* rsq  * ssq  * tm1 * f[ 4] +
		(-c4) * rt	* rp1  * ssq  * tm1 * f[ 5] +
		( c8) * rst * rm1  * sp1  * tm1 * f[ 6] +
		(-c4) * st	* rsq  * sp1  * tm1 * f[ 7] +
		(-c8) * rst * rp1  * sp1  * tm1 * f[ 8] +

		( c4) * rs	* rm1  * sm1  * tsq * f[ 9] +
		(-c2) * s	* rsq  * sm1  * tsq * f[10] +
		(-c4) * rs	* rp1  * sm1  * tsq * f[11] +
		(-c2) * r	* rm1  * ssq  * tsq * f[12] +
					  rsq  * ssq  * tsq * f[13] +
		( c2) * r	* rp1  * ssq  * tsq * f[14] +
		(-c4) * rs	* rm1  * sp1  * tsq * f[15] +
		( c2) * s	* rsq  * sp1  * tsq * f[16] +
		( c4) * rs	* rp1  * sp1  * tsq * f[17] +

		( c8) * rst * rm1  * sm1  * tp1 * f[18] +
		(-c4) * st	* rsq  * sm1  * tp1 * f[19] +
		(-c8) * rst * rp1  * sm1  * tp1 * f[20] +
		(-c4) * rt	* rm1  * ssq  * tp1 * f[21] +
		( c2) * t	* rsq  * ssq  * tp1 * f[22] +
		( c4) * rt	* rp1  * ssq  * tp1 * f[23] +
		(-c8) * rst * rm1  * sp1  * tp1 * f[24] +
		( c4) * st	* rsq  * sp1  * tp1 * f[25] +
		( c8) * rst * rp1  * sp1  * tp1 * f[26]);
}

#define  fdata(i,j)      fdata  [((j)-1)*nxdata + (i)-1]
float Util::quadri(float xx, float yy, int nxdata, int nydata, float* fdata)
{
/*
c  purpose: quadratic interpolation 
c 
c  parameters:       xx,yy treated as circularly closed.
c                    fdata - image 1..nxdata, 1..nydata
c
c                    f3    fc       f0, f1, f2, f3 are the values
c                     +             at the grid points.  x is the
c                     + x           point at which the function
c              f2++++f0++++f1       is to be estimated. (it need
c                     +             not be in the first quadrant).
c                     +             fc - the outer corner point
c                    f4             nearest x.
c
c                                   f0 is the value of the fdata at
c                                   fdata(i,j), it is the interior mesh
c                                   point nearest  x.
c                                   the coordinates of f0 are (x0,y0),
c                                   the coordinates of f1 are (xb,y0),
c                                   the coordinates of f2 are (xa,y0),
c                                   the coordinates of f3 are (x0,yb),
c                                   the coordinates of f4 are (x0,ya),
c                                   the coordinates of fc are (xc,yc),
c
c                   o               hxa, hxb are the mesh spacings
c                   +               in the x-direction to the left
c                  hyb              and right of the center point.
c                   +
c            ++hxa++o++hxb++o       hyb, hya are the mesh spacings
c                   +               in the y-direction.
c                  hya
c                   +               hxc equals either  hxb  or  hxa
c                   o               depending on where the corner
c                                   point is located.
c
c                                   construct the interpolant
c                                   f = f0 + c1*(x-x0) +
c                                       c2*(x-x0)*(x-x1) +
c                                       c3*(y-y0) + c4*(y-y0)*(y-y1)
c                                       + c5*(x-x0)*(y-y0)
c
c
*/
    float x, y, dx0, dy0, f0, c1, c2, c3, c4, c5, dxb, dyb;
    float quadri;
    int   i, j, ip1, im1, jp1, jm1, ic, jc, hxc, hyc;
    
    x = xx;
    y = yy;

    // circular closure
    if (x < 1.0)               x = x+(1 - floor(x) / nxdata) * nxdata;
    if (x > (float)nxdata+0.5) x = fmod(x-1.0f,(float)nxdata) + 1.0f;
    if (y < 1.0)               y = y+(1 - floor(y) / nydata) * nydata;
    if (y > (float)nydata+0.5) y = fmod(y-1.0f,(float)nydata) + 1.0f;


    i   = (int) floor(x);
    j   = (int) floor(y);

    dx0 = x - i;
    dy0 = y - j;

    ip1 = i + 1;
    im1 = i - 1;
    jp1 = j + 1;
    jm1 = j - 1;

    if (ip1 > nxdata) ip1 = ip1 - nxdata;
    if (im1 < 1)      im1 = im1 + nxdata;
    if (jp1 > nydata) jp1 = jp1 - nydata;
    if (jm1 < 1)      jm1 = jm1 + nydata;

    f0  = fdata(i,j);
    c1  = fdata(ip1,j) - f0;
    c2  = (c1 - f0 + fdata(im1,j)) * 0.5;
    c3  = fdata(i,jp1) - f0;
    c4  = (c3 - f0 + fdata(i,jm1)) * 0.5;

    dxb = dx0 - 1;
    dyb = dy0 - 1;

    // hxc & hyc are either 1 or -1
    if (dx0 >= 0) {
       hxc = 1;
    }
    else {
       hxc = -1;
    }
    if (dy0 >= 0) {
       hyc = 1;
    }
    else {
       hyc = -1;
    }
 
    ic  = i + hxc;
    jc  = j + hyc;

    if (ic > nxdata) {
       ic = ic - nxdata;
    }
    else if (ic < 1) {
       ic = ic + nxdata;
    }

    if (jc > nydata) {
       jc = jc - nydata;
    }
    else if (jc < 1) {
       jc = jc + nydata;
    }

    c5  =  ( (fdata(ic,jc) - f0 - hxc * c1 - (hxc * (hxc - 1.0)) * c2 
            - hyc * c3 - (hyc * (hyc - 1.0)) * c4) * (hxc * hyc));

    quadri = f0 + dx0 * (c1 + dxb * c2 + dy0 * c5) + dy0 * (c3 + dyb * c4);

    return quadri; 
}
#undef fdata

Util::KaiserBessel::KaiserBessel(float alpha_, int K_, float r_, float v_,
		                         int N_, float vtable_, int ntable_) 
		: alpha(alpha_), v(v_), r(r_), N(N_), K(K_), vtable(vtable_), 
		  ntable(ntable_) {
	// Default values are alpha=1.25, K=6, r=0.5, v = K/2
	if (0.f == v) v = float(K)/2;
	if (0.f == vtable) vtable = v;
	fac = static_cast<float>(twopi)*alpha*r*v;
	alphar = alpha*r;
	vadjust = 1.0f*v;
	facadj = static_cast<float>(twopi)*alpha*r*vadjust;
	build_I0table();
}

float Util::KaiserBessel::i0win(float x) const {
	float val0 = float(gsl_sf_bessel_I0(facadj));
	float absx = fabs(x);
	if (absx > vadjust) return 0.f;
	float rt = sqrt(1.f - pow(absx/vadjust, 2));
	float res = float(gsl_sf_bessel_I0(facadj*rt))/val0;
	return res;
}

void Util::KaiserBessel::build_I0table() {
	i0table.resize(ntable+1); // i0table[0:ntable]
	int ltab = int(round(float(ntable)/1.25f));
	fltb = float(ltab)/(K/2);
	//float val0 = sqrt(facadj)*static_cast<float>(gsl_sf_bessel_I1(facadj));
	float val0 = gsl_sf_bessel_I0(facadj);
	for (int i=ltab+1; i <= ntable; i++) i0table[i] = 0.f;
	for (int i=0; i <= ltab; i++) {
		float s = float(i)/fltb/N;
		if (s < vadjust) {
			float rt = sqrt(1.f - pow(s/vadjust, 2));
			//i0table[i] = sqrt(facadj*rt)*static_cast<float>(gsl_sf_bessel_I1(facadj*rt))/val0;
			i0table[i] = gsl_sf_bessel_I0(facadj*rt)/val0;
		} else {
			i0table[i] = 0.f;
		}
	}
}

float Util::KaiserBessel::I0table_maxerror() {
	float maxdiff = 0.f;
	for (int i = 1; i <= ntable; i++) {
		float diff = fabs(i0table[i] - i0table[i-1]);
		if (diff > maxdiff) maxdiff = diff;
	}
	return maxdiff;
}

float Util::KaiserBessel::sinhwin(float x) const {
	float val0 = sinh(fac)/fac;
	float absx = fabs(x);
	if (0.0 == x) {
		float res = 1.0f;
		return res;
	} else if (absx == alphar) {
		return 1.0f/val0;
	} else if (absx < alphar) {
		float rt = sqrt(1.0f - pow((x/alphar), 2));
		float facrt = fac*rt;
		float res = (sinh(facrt)/facrt)/val0;
		return res;
	} else {
		float rt = sqrt(pow((x/alphar),2) - 1.f);
		float facrt = fac*rt;
		float res = (sin(facrt)/facrt)/val0;
		return res;
	}
}

float Util::FakeKaiserBessel::i0win(float x) const {
	float val0 = sqrt(facadj)*float(gsl_sf_bessel_I1(facadj));
	float absx = fabs(x);
	if (absx > vadjust) return 0.f;
	float rt = sqrt(1.f - pow(absx/vadjust, 2));
	float res = sqrt(facadj*rt)*float(gsl_sf_bessel_I1(facadj*rt))/val0;
	return res;
}

void Util::FakeKaiserBessel::build_I0table() {
	i0table.resize(ntable+1); // i0table[0:ntable]
	int ltab = int(round(float(ntable)/1.1f));
	fltb = float(ltab)/(K/2);
	float val0 = sqrt(facadj)*gsl_sf_bessel_I1(facadj);
	for (int i=ltab+1; i <= ntable; i++) i0table[i] = 0.f;
	for (int i=0; i <= ltab; i++) {
		float s = float(i)/fltb/N;
		if (s < vadjust) {
			float rt = sqrt(1.f - pow(s/vadjust, 2));
			i0table[i] = sqrt(facadj*rt)*gsl_sf_bessel_I1(facadj*rt)/val0;
		} else {
			i0table[i] = 0.f;
		}
	}
}

float Util::FakeKaiserBessel::sinhwin(float x) const {
	float val0 = sinh(fac)/fac;
	float absx = fabs(x);
	if (0.0 == x) {
		float res = 1.0f;
		return res;
	} else if (absx == alphar) {
		return 1.0f/val0;
	} else if (absx < alphar) {
		float rt = sqrt(1.0f - pow((x/alphar), 2));
		float facrt = fac*rt;
		float res = (sinh(facrt)/facrt)/val0;
		return res;
	} else {
		float rt = sqrt(pow((x/alphar),2) - 1.f);
		float facrt = fac*rt;
		float res = (sin(facrt)/facrt)/val0;
		return res;
	}
}

#if 0 // 1-st order KB window
float Util::FakeKaiserBessel::sinhwin(float x) const {
	//float val0 = sinh(fac)/fac;
	float prefix = 2*facadj*vadjust/float(gsl_sf_bessel_I1(facadj));
	float val0 = prefix*(cosh(facadj) - sinh(facadj)/facadj);
	float absx = fabs(x);
	if (0.0 == x) {
		//float res = 1.0f;
		float res = val0;
		return res;
	} else if (absx == alphar) {
		//return 1.0f/val0;
		return prefix;
	} else if (absx < alphar) {
		float rt = sqrt(1.0f - pow((x/alphar), 2));
		//float facrt = fac*rt;
		float facrt = facadj*rt;
		//float res = (sinh(facrt)/facrt)/val0;
		float res = prefix*(cosh(facrt) - sinh(facrt)/facrt);
		return res;
	} else {
		float rt = sqrt(pow((x/alphar),2) - 1.f);
		//float facrt = fac*rt;
		float facrt = facadj*rt;
		//float res = (sin(facrt)/facrt)/val0;
		float res = prefix*(sin(facrt)/facrt - cos(facrt));
		return res;
	}
}
#endif // 0


EMData* Util::Polar2D(EMData* image, vector<int> numr, string mode){
   int nsam = image->get_xsize();
   int nrow = image->get_ysize();
   int nring = numr.size()/3;
   int lcirc = numr[3*nring-2]+numr[3*nring-1]-1;
   EMData* out = new EMData();
   char cmode = (mode == "F" || mode == "f") ? 'f' : 'h';
   out->set_size(lcirc,1,1);
   alrq(image->get_data(), nsam, nrow, &numr[0], out->get_data(), lcirc, nring, cmode);
   return out;
}

#define  circ(i)         circ   [(i)-1]
#define  numr(i,j)       numr   [((j)-1)*3 + (i)-1]
#define  xim(i,j)        xim    [((j)-1)*nsam + (i)-1]
void Util::alrq(float *xim,  int nsam , int nrow , int *numr,
          float *circ, int lcirc, int nring, char mode)
{
/* 
c                                                                     
c  purpose:                                                          
c                                                                   
c  resmaple to polar coordinates
c                                                                  
*/
   //  dimension         xim(nsam,nrow),circ(lcirc)
   //  integer           numr(3,nring)

   double dfi, dpi;
   int    ns2, nr2, i, inr, l, nsim, kcirc, lt, j;
   float  yq, xold, yold, fi, x, y;

   ns2 = nsam/2+1;
   nr2 = nrow/2+1;
   dpi = 2.0*atan(1.0);

   for (i=1;i<=nring;i++) {
     // radius of the ring
     inr = numr(1,i);
     yq  = inr;
     l   = numr(3,i);
     if (mode == 'h' || mode == 'H') {
        lt = l/2;
     }
     else  {    //if (mode == 'f' || mode == 'F' )
        lt = l/4;
     }

     nsim           = lt-1;
     dfi            = dpi/(nsim+1);
     kcirc          = numr(2,i);
     xold           = 0.0;
     yold           = inr;
     circ(kcirc)    = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
     xold           = inr;
     yold           = 0.0;
     circ(lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

     if (mode == 'f' || mode == 'F') {
        xold              = 0.0;
        yold              = -inr;
        circ(lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        xold              = -inr;
        yold              = 0.0;
        circ(lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
     }

     for (j=1;j<=nsim;j++) {
        fi               = dfi*j;
        x                = sin(fi)*yq;
        y                = cos(fi)*yq;
        xold             = x;
        yold             = y;
        circ(j+kcirc)    = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        xold             =  y;
        yold             = -x;
        circ(j+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);

        if (mode == 'f' || mode == 'F')  {
           xold                = -x;
           yold                = -y;
           circ(j+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
           xold                = -y;
           yold                =  x;
           circ(j+lt+lt+lt+kcirc) = quadri(xold+ns2,yold+nr2,nsam,nrow,xim);
        };
     }
   }
 
}
