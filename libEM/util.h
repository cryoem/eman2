#ifndef eman__util_h__
#define eman__util_h__ 1

#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>

using std::string;
using std::vector;

namespace EMAN
{
    class Util
    {
    public:
	static void ap2ri(float *data, size_t n);
	static void flip_complex_phase(float *data, size_t n);
	static int file_lock_wait(FILE * file);
	static int generate_machine_stamp();
	static bool check_file_by_magic(const void *first_block, const char *magic);
	static void flip_image(float *data, size_t nx, size_t ny);
	static bool is_sub_string(const char *s1, const char *s2);
	static string get_filename_by_ext(string old_filename, string ext);

	static void calc_least_square_fit(size_t nitems, const float *data_x, const float *data_y,
					  float *slope, float *intercept, bool ignore_zero);

	static void save_data_to_file(const vector<float> & x_array,
				      const vector<float> & y_array, string filename);
	static void save_data_to_file(float x0, float dx, const vector<float> & y_array,
				      string filename);
	static void save_data_to_file(float x0, float dx, float *y_array, size_t array_size,
				      string filename);

	static float get_frand(float low, float high);
	static float get_gaussian_rand(float mean, float sigma);

	static inline int round(float x)
	{
	    if (x < 0) {
		return (int) (x - 0.5);
	    }
	    return (int) (x + 0.5);
	}

	// p1=x0,y0, p2=x1,y0; p3=x1,y1; p4=x0,y1 
	static inline float bilinear_interpolate(float p1, float p2, float p3, float p4, float t,
						 float u)
	{
	    return (1 - t) * (1 - u) * p1 + t * (1 - u) * p2 + t * u * p3 + (1 - t) * u * p4;
	}

	// p1=x0,y0,z0; p2=x1,y0,z0; p3=x0,y1,z0, p4=x1,y1,z0
	// p5=x0,y0,z1; p6=x1,y0,z1; p7=x0,y1,z1, p8=x1,y1,z1
	static inline float trilinear_interpolate(float p1, float p2, float p3, float p4, float p5,
						  float p6, float p7, float p8, float t, float u,
						  float v)
	{
	    return ((1 - t) * (1 - u) * (1 - v) * p1 + t * (1 - u) * (1 - v) * p2
		    + (1 - t) * u * (1 - v) * p3 + t * u * (1 - v) * p4
		    + (1 - t) * (1 - u) * v * p5 + t * (1 - u) * v * p6
		    + (1 - t) * u * v * p7 + t * u * v * p8);
	}


	static void find_max(float *data, size_t nitems, float *max_val, int *max_index = 0);
	static void find_min_and_max(float *data, size_t nitems, float *max_val, float *min_val,
				     int *max_i = 0, int *min_i = 0);
	static int calc_best_fft_size(int low);

	static inline int square(int n)
	{
	    return n * n;
	}
	
	static inline float square(float x)
	{
	    return x * x;
	}
	static double square(double x)
	{
	    return x * x;
	}

	static double square_sum(double x, double y)
	{
	    return (x * x + y * y);
	}

	static inline double hypot3(double x, double y, double z)
	{
	    return sqrt(x * x + y * y + z * z);
	}

	static inline int fast_floor(float x)
	{
	    if (x < 0) {
		return ((int) x - 1);
	    }
	    return (int) x;
	}

	static inline float agauss(float a, float dx, float dy, float dz, float d)
	{
	    return (a * exp(-(dx * dx + dy * dy + dz * dz) / d));
	}

	
	static inline float min(float f1, float f2)
	{
	    return (f1 < f2 ? f1 : f2);
	}

	static inline float min(float f1, float f2, float f3)
	{
	    if (f1 <= f2 && f1 <= f3) {
		return f1;
	    }
	    if (f2 <= f1 && f2 <= f3) {
		return f2;
	    }
	    return f3;
	}

	static inline float min(float f1, float f2, float f3, float f4)
	{
	    float m = f1;
	    if (f2 < m) {
		m = f2;
	    }
	    if (f3 < m) {
		m = f3;
	    }
	    if (f4 < m) {
		m = f4;
	    }
	    return m;
	}

	static inline float max(float f1, float f2)
	{
	    return (f1 < f2 ? f2 : f1);
	}

	static inline float max(float f1, float f2, float f3)
	{
	    if (f1 >= f2 && f1 >= f3) {
		return f1;
	    }
	    if (f2 >= f1 && f2 >= f3) {
		return f2;
	    }
	    return f3;
	}

	static inline float max(float f1, float f2, float f3, float f4)
	{
	    float m = f1;
	    if (f2 > m) {
		m = f2;
	    }
	    if (f3 > m) {
		m = f3;
	    }
	    if (f4 > m) {
		m = f4;
	    }
	    return m;
	}

	static inline float angle_sub_2pi(float x, float y)
	{
	    float r = fmod(fabs(x - y), (float) (2 * M_PI));
	    if (r > M_PI) {
		r = 2.0 * M_PI - r;
	    }

	    return r;
	}

	static inline float angle_sub_pi(float x, float y)
	{
	    float r = fmod(fabs(x - y), (float) M_PI);
	    if (r > M_PI / 2.0) {
		r = M_PI - r;
	    }
	    return r;
	}

	static inline int goodf(float *f)
	{
	    // the first is abnormal zero the second is +-inf or NaN 
	    if ((((int *) f)[0] & 0x7f800000) == 0 || (((int *) f)[0] & 0x7f800000) == 255) {
		return 0;
	    }
	    return 1;
	}

	static string get_time_label();

	static void set_log_level(int argc, char* argv[]);

	static inline float eman_copysign(float a, float b)
	{
#ifndef _WIN32
	    return copysign(a, b);
#else
	    int flip = -1;
	    if ((a <= 0 && b <= 0) || ( a > 0 && b > 0)) {
		flip = 1;
	    }	    
	    return a*flip;
#endif
	}

	static inline double eman_erfc(double x)
	{
#ifndef _WIN32
	    return erfc(x);
#else
	    static double a[] = { -1.26551223, 1.00002368,
				   0.37409196, 0.09678418,
				  -0.18628806, 0.27886807,
				  -1.13520398, 1.48851587,
				  -0.82215223, 0.17087277 };
	    
	    double result = 1;
	    double z = fabs(x);
	    if (z > 0) {
		double t = 1 / (1 + 0.5*z);
		double f1 = t*(a[4]+t*(a[5]+t*(a[6]+t*(a[7]+t*(a[8]+t*a[9])))));
		result = t*exp((-z*z) + a[0] + t*(a[1]+t*(a[2]+t*(a[3]+f1))));
		
		if (x < 0) {
		    result = 2 - result;
		}
	    }
	    return result;
#endif
	}
	
    };
}

#endif
