#ifndef eman__util_h__
#define eman__util_h__ 1

#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;

namespace EMAN {
    class Util {
    public:
	static void ap2ri(float* data, int n);
	static int file_lock_wait(FILE *file);
	static int generate_machine_stamp();
	static bool check_file_by_magic(const void* first_block, const char* magic);
	static void flip_image(float* data, int nx, int ny);
	static bool is_sub_string(const char* s1, const char* s2);
	static string get_filename_by_ext(string old_filename, string ext);

	static void least_square_fit(int nitems, float* data_x, float* data_y, float* slope, float* intercept, bool ignore_zero);

	static void save_data_to_file(vector<float> x_array, vector<float> y_array, const char* filename);
	static void save_data_to_file(float x0, float dx, vector<float> y_array, const char* filename);


	static float get_frand(float lo, float hi);
	static float get_gaussian_rand(float avg, float std);
	
	
	static inline float bilinear_interpolate(float p1, float p2, float p3, float p4, float t, float u) {
	    return (1.0f-t)*(1.0f-u)*p1+t*(1.0f-u)*p2+t*u*p3+(1.0f-t)*u*p4;
	}

	static inline float trilinear_interpolate(float p1, float p2, float p3, float p4, float p5, float p6,
						  float p7, float p8, float t, float u, float v) {
	    return (1.0-t)*(1.0-u)*(1.0-v)*p1+t*(1.0-u)*(1.0-v)*p2+
		(1.0-t)*u*(1.0-v)*p3+t*u*(1.0-v)*p4+
		(1.0-t)*(1.0-u)*v*p5+t*(1.0-u)*v*p6+
		(1.0-t)*u*v*p7+t*u*v*p8;
	}
	
	static void find_max(float* data, int nitems, float* max_val, int* max_index = 0);
	static void find_min_and_max(float* data, int nitems, float* max_val, float* min_val,
				     int* max_i = 0, int* min_i = 0);
	static int calc_best_fft_size(int low);
    };
}

#endif
