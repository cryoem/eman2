/**
 * $Id$
 */
#include "util.h"
#include "assert.h"
#include "byteorder.h"
#include "log.h"

#include <string.h>
#include <math.h>
#include <sys/types.h>
#ifndef WIN32
#include <unistd.h>
#include <sys/param.h>
#endif

#include <fcntl.h>
#include <float.h>
#include <time.h>

using namespace EMAN;

void Util::ap2ri(float *data, size_t n)
{
    if (!data) {
	return;
    }

    for (size_t i = 0; i < n; i += 2) {
	float f = data[i] * sin(data[i + 1]);
	data[i] = data[i] * cos(data[i + 1]);
	data[i + 1] = f;
    }
}

void Util::flip_complex_phase(float *data, size_t n)
{
    if (!data) {
	return;
    }

    for (size_t i = 0; i < n; i += 2) {
	data[i+1] *= -1;
    }    
}

int Util::file_lock_wait(FILE * file)
{
#ifdef WIN32
    return 1;
#else

    if (!file) {
	Log::logger()->error("Tried to lock NULL file");
	return 1;
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

    if (fcntl(fdes, F_SETLKW, &fl) == -1) {
	Log::logger()->error("file locking error! NFS problem?");

	int i = 0;
	for (i = 0; i < 5; i++) {
	    if (fcntl(fdes, F_SETLKW, &fl) != -1) {
		break;
	    }
	    else {
#ifdef WIN32
		Sleep(0.001);
#else
		sleep(1);
#endif
		
	    }
	}
	if (i == 5) {
	    Log::logger()->error("Fatal file locking error");
	    return 1;
	}
    }

    return 0;
#endif
}

void Util::file_unlock(FILE *)
{
    
}

int Util::generate_machine_stamp()
{
    int stamp = 0;
    char *p = (char *) (&stamp);

    if (ByteOrder::is_host_big_endian()) {
	p[0] = 0x44;
	p[1] = 0x44;
	p[2] = 0;
	p[3] = 0;
    }
    else {
	p[0] = 0x11;
	p[1] = 0x11;
	p[2] = 0;
	p[3] = 0;
    }

    return stamp;

}

bool Util::check_file_by_magic(const void *first_block, const char *magic)
{
    if (!first_block || !magic) {
	return false;
    }

    const char *buf = static_cast<const char *>(first_block);

    if (strncmp(buf, magic, strlen(magic)) == 0) {
	return true;
    }
    return false;
}


void Util::flip_image(float *data, size_t nx, size_t ny)
{
    if (!data) {
	return;
    }
    
    float *buf = new float[nx];
    size_t row_size = nx * sizeof(float);

    for (size_t i = 0; i < ny / 2; i++) {
	memcpy(buf, &data[i * nx], row_size);
	memcpy(&data[i * nx], &data[(ny - 1 - i) * nx], row_size);
	memcpy(&data[(ny - 1 - i) * nx], buf, row_size);
    }

    delete [] buf;
    buf = 0;
}

bool Util::sstrncmp(const char *s1, const char *s2)
{
    if (strncmp(s1, s2, strlen(s2)) == 0) {
	return true;
    }
    return false;
}

bool Util::get_str_float(const char *s, const char *float_var, float *p_val)
{
    size_t n = strlen(float_var);    
    if (strncmp(s, float_var, n) == 0) {
	*p_val = atof(&s[n]);
	return true;
    }
    return false;
}

bool Util::get_str_float(const char *s, const char *float_var, float *p_v1, float *p_v2)
{
    size_t n = strlen(float_var);
    if (strncmp(s, float_var, n) == 0) {
	sscanf(&s[n], "%f,%f", p_v1, p_v2);
	return true;
    }
    return false;
}

bool Util::get_str_float(const char *s, const char *float_var, int *p_v0, float *p_v1, float *p_v2)
{
    size_t n = strlen(float_var);
    *p_v0 = 0;
    if (strncmp(s, float_var, n) == 0) {
	if (s[n] == '=') {
	    *p_v0 = 2;
	    sscanf(&s[n+1], "%f,%f", p_v1, p_v2);
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
    size_t n = strlen(int_var);    
    if (strncmp(s, int_var, n) == 0) {
	*p_val = atoi(&s[n]);
	return true;
    }
    return false;
}

bool Util::get_str_int(const char *s, const char *int_var, int *p_v1, int *p_v2)
{
    size_t n = strlen(int_var);    
    if (strncmp(s, int_var, n) == 0) {
	sscanf(&s[n], "%d,%d", p_v1, p_v2);
	return true;
    }    
    return false;
}

bool Util::get_str_int(const char *s, const char *int_var, int *p_v0, int *p_v1, int *p_v2)
{
    size_t n = strlen(int_var);
    *p_v0 = 0;
    if (strncmp(s, int_var, n) == 0) {
	if (s[n] == '=') {
	    *p_v0 = 2;
	    sscanf(&s[n+1], "%d,%d", p_v1, p_v2);
	}
	else {
	    *p_v0 = 1;
	}
	return true;
    }
    return false;
}

string Util::get_filename_by_ext(string old_filename, string ext)
{
    char buf[MAXPATHLEN];
    strcpy(buf, old_filename.c_str());
    char *old_ext = strrchr(buf, '.');
    buf[strlen(buf) - strlen(old_ext)] = '\0';
    strcat(buf, ext.c_str());
    return string(buf);
}

void Util::calc_least_square_fit(size_t nitems, const float *data_x, const float *data_y,
				 float *slope, float *intercept, bool ignore_zero)
{
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
	div = 0.0000001;
    }

    *intercept = (sum_xx * sum_y - sum_x * sum_xy) / div;
    *slope = (sum * sum_xy - sum_x * sum_y) / div;
}

void Util::save_data(const vector<float> & x_array, const vector<float> & y_array,
			     string filename)
{
    if (x_array.size() != y_array.size()) {
	Log::logger()->error("array x and array y have different size: %d != %d\n",
			     x_array.size(), y_array.size());
	return;
    }

    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) {
	Log::logger()->error("cannot open file to save data '%s'\n", filename.c_str());
	return;
    }

    for (size_t i = 0; i < x_array.size(); i++) {
	fprintf(out, "%g\t%g\n", x_array[i], y_array[i]);
    }
    fclose(out);
}

void Util::save_data(float x0, float dx, const vector<float> & y_array, string filename)
{
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) {
	Log::logger()->error("cannot open file to save data '%s'\n", filename.c_str());
	return;
    }

    for (size_t i = 0; i < y_array.size(); i++) {
	fprintf(out, "%g\t%g\n", x0 + dx * i, y_array[i]);
    }
    fclose(out);
}


void Util::save_data(float x0, float dx, float *y_array, size_t array_size, string filename)
{
    if (!y_array) {
	return;
    }
    
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) {
	Log::logger()->error("cannot open file to save data '%s'\n", filename.c_str());
	return;
    }

    for (size_t i = 0; i < array_size; i++) {
	fprintf(out, "%g\t%g\n", x0 + dx * i, y_array[i]);
    }
    fclose(out);
}

float Util::get_frand(float lo, float hi)
{
    static bool inited = false;
    if (!inited) {
	srand(time(0));
	inited = true;
    }

    float r = (hi - lo) * rand() / (RAND_MAX + 1.0) + lo;
    return r;
}

float Util::get_gauss_rand(float mean, float sigma)
{
    float x = 0;
    float r = 0;
    bool valid = true;
    
    while (valid) {
	x = get_frand(-1, 1);
	float y = get_frand(-1, 1);
	r = x * x + y * y;

	if (r <= 1.0 && r > 0) {
	    valid = false;
	}
    }

    float f = sqrt(-2.0 * log(r) / r);
    float result = x * f * sigma + mean;
    return result;
}

void Util::find_max(float *data, size_t nitems, float *max_val, int *max_index)
{
    if (!data) {
	return;
    }
    float max = -FLT_MAX;
    int m = 0;
    
    for (size_t i = 0; i < nitems; i++) {
	if (data[i] > max) {
	    max = data[i];
	    m = i;
	}
    }

    *max_val = m;
    
    if (max_index) {
	*max_index = m;
    }
}

void Util::find_min_and_max(float *data, size_t nitems, float *max_val, float *min_val,
			    int *max_index, int *min_index)
{
    if (!data) {
	return;
    }
    float max = -FLT_MAX;
    float min = FLT_MAX;
    int max_i = 0;
    int min_i = 0;

    for (size_t i = 0; i < nitems; i++) {
	if (data[i] > max) {
	    max = data[i];
	    max_i = i;
	}
	if (data[i] < min) {
	    min = data[i];
	    min_i = i;
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


// finds fft sizes with good primes. needs some work ...
int Util::calc_best_fft_size(int low)
{
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
			
			float f = f1 * f2  * f3 * pow((float) 7.0, i7);
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
    
    Log::logger()->error("Sorry, can only find good fft sizes up to 4096 right now.");

    return 1;
}

string Util::get_time_label()
{
    time_t t0 = time(0);
    struct tm *t = localtime(&t0);
    char label[32];
    sprintf(label, "%d/%02d/%04d %d:%02d",
	    t->tm_mon+1, t->tm_mday, t->tm_year+1900, t->tm_hour, t->tm_min);
    return string(label);
}

void Util::set_log_level(int argc, char* argv[])
{
    if (argc > 1 && strncmp(argv[1], "-v", 2) == 0) {
	char level_str[32];
	strcpy(level_str, argv[1]+2);
	Log::LogLevel log_level = (Log::LogLevel) atoi(level_str);
	Log::logger()->set_level(log_level);
    }
}
