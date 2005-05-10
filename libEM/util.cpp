/**
 * $Id$
 */
#include "util.h"
#include "byteorder.h"
#include "log.h"
#include "exception.h"
#include "Assert.h"

#include <string.h>
#include <string>
#include <math.h>
#include <sys/types.h>
#ifndef WIN32
#include <unistd.h>
#include <sys/param.h>
#else
#include <io.h>
#define access _access
#define F_OK 00
#endif

#include <fcntl.h>
#include <float.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <sstream>

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

	delete[]buf;
	buf = 0;
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
	delete [] buf;
	buf = 0;
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
								 float *slope, float *intercept, bool ignore_zero)
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
		if (!ignore_zero || (data_x[i] != 0 && data_y[i] != 0)) {
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

namespace { // unnamed namespace to keep these functions private

	// Remove the '@' character if it is at the end of the string.
	// If there's a '@' elsewhere, then it's an error, and we return
	// an empty string.
	string rm_stack_char(string s) {
		string::size_type stackloc = s.find_first_of("@");
		if (stackloc != string::npos) {
			// there's a '@' somewhere
			if (stackloc == s.size()-1) {
				// it's at the end of the string
				return s.substr(0,s.size()-1);
			} else {
				// '@' not at the end, so it's an error
				throw InvalidStringException(s, "Invalid format: misplaced `@`.");
			}
		} else {
			// no '@' character
			return s;
		}
	}

	// Define the "location" of a field in the string.
	typedef struct { 
		int begin; 
		int end; 
		int count() { return end - begin + 1; }
	} fieldloc;

	// Find the locations of all {*...*} fields in the string.
	vector<fieldloc> findfields(string s) {
		vector<fieldloc> fields;
		string::size_type loc = 0;
		string::size_type begin;
		while ((begin = s.find_first_of("{", loc)) != string::npos) {
			int end = s.find_first_of("}", begin);
			fieldloc field; field.begin = begin; field.end = end;
			fields.push_back(field);
			// check to ensure that the field is valid
			string::size_type nasterisks = field.count() - 2;
			string asterisks = s.substr(field.begin+1, nasterisks);
			if (asterisks.find_first_not_of("*") != string::npos) {
				// there's something other than '*'s in the field
				throw InvalidStringException(s, "Malformed {*...*} field.");
			}
			loc = end;
		}
		return fields;
	}

} // un-named namespace

string EMAN::Util::parse_spider_fname(string s, vector<int> fieldvals) {
	// remove leading and trailing whitespace
	trim(s);
	// If there is a stack character, remove it.
	s = rm_stack_char(s);
	// where are the fields to replace?
	vector<fieldloc> fields = findfields(s);
	if (fields.size() != fieldvals.size()) {
		// wrong number of fields?!
		throw InvalidValueException(int(fields.size()), 
				"Number of fields provided "
				"differs from the number of {*...*} fields");
	}
	// create and return properly formatted string
	// foo{****} should return "foo0002" if 2 is the value to sub in.
	ostringstream news;
	news << setiosflags(ios::fixed) << setfill('0');
	for (string::size_type i = 0, loc = 0; i < fields.size(); i++) {
		news << s.substr(loc, fields[i].begin - loc)
			<< setw(fields[i].count() - 2)
			<< fieldvals[i];
		loc = fields[i].end + 1;
	}
	news << s.substr(fields[fields.size()-1].end+1);
	return news.str();
}

void Util::printMatI3D(MIArray3D& mat, int nx, int ny, int nz,
					   const string str, ostream& out) {
	out << "Printing 3D Integer data: " << str << std::endl;
	for (int iz = 1; iz <= nz; iz++) {
		cout << "(z = " << iz << " slice)" << endl;
		for (int ix = 0; ix <= nx; ix++) {
			for (int iy = 1; iy <= ny; iy++) {
				cout << setiosflags(ios::fixed) << setw(5) 
					 << mat[ix][iy][iz] << "  ";
			}
			cout << endl;
		}
	}
}

void Util::voea(string filename, float delta, float t1, float t2,
		float p1, float p2)
{
	const float QUADPI = 3.141592653589793238462643383279502884197;
	const float DGR_TO_RAD = QUADPI/180.;
	std::ofstream fout(filename.c_str());
	if (!fout.good()) {
		LOGERR("Cannot write angles output file");
		return;
	}
	float psi = 0.0;
	if ((0.0 == t1)&&(0.0 == t2)||(t1 >= t2)) {
		t1 = 0.0;
		t2 = 90.0;
	}
	if ((0.0 == p1)&&(0.0 == p2)||(p1 >= p2)) {
		p1 = 0.0;
		p2 = 359.9;
	}
	bool skip = ((t1 < 90.0)&&(90.0 == t2)&&(0.0 == p1)&&(p2 > 180.0));
	for (float theta = t1; theta <= t2; theta += delta) {
		float detphi;
		int lt;
		if ((0.0 == theta)||(180.0 == theta)) {
			detphi = 360.0;
			lt = 1;
		} else {
			detphi = delta/sin(theta*DGR_TO_RAD);
			lt = int((p2 - p1)/detphi)-1;
			if (lt < 1) lt = 1;
			detphi = (p2 - p1)/lt;
		}
		for (int i = 0; i < lt; i++) {
			float phi = p1 + i*detphi;
			if (skip&&(90.0 == theta)&&(phi > 180.0)) continue;
			fout << psi << '\t' << theta << '\t' << phi << std::endl;
		}
	}
	fout.close();
}
