/**
 * $Id$
 */
#ifndef eman__util_h__
#define eman__util_h__ 1

#include "emconstants.h"
#include "exception.h"
#include <vector>
#include <iostream>
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>

#ifdef WIN32
#include <windows.h>
#define M_PI 3.14159265358979323846f
//#define MAXPATHLEN (MAX_PATH * 4)
#endif

using std::string;
using std::vector;
using std::ostream;

namespace EMAN
{
	class EMData;
	class Dict;
	
	typedef boost::multi_array<int, 3> MIArray3D;

	/** @ingroup tested3c */
	/** Util is a collection of utility functions.
     */
	class Util
	{
	  public:
		/** convert complex data array from Amplitude/Phase format
		 * into Real/Imaginary format.
		 * @param data complex data array.
		 * @param n array size.
		 */
		static void ap2ri(float *data, size_t n);

		/** flip the phase of a complex data array.
		 * @param data complex data array.
		 * @param n array size.
		 */
		static void flip_complex_phase(float *data, size_t n);

		/** lock a file. If the lock fails, wait for 1 second; then
		 * try again. Repleat this wait-try for a maxinum of 5 times
		 * unless the lock succeeds.
		 *
		 * @param file The file to be locked.
		 * @return 0 if lock succeeds. 1 if fails.
		 */
		static int file_lock_wait(FILE * file);
		//static void file_unlock(FILE * file);

		/** check whether a file starts with certain magic string.
		 * @param first_block The first block of the file.
		 * @param magic The magic string to identify a file format.
		 * @return True if file matches magic. Otherwise, false.
		 */
		static bool check_file_by_magic(const void *first_block, const char *magic);

		/** check whether a file exists or not
		 * @return True if the file exists; False if not.
		 */
		static bool is_file_exist(const string & filename);

		/** Vertically flip the data of a 2D real image.
		 * @param data Data array of the 2D real image.
		 * @param nx Image Width.
		 * @param ny Image Height.
		 */
		static void flip_image(float *data, size_t nx, size_t ny);

		/** Safe string compare. It compares 's2' with the first N
		 * characters of 's1', where N is the length of 's2'.
		 * 
		 * @param s1 String 1. Its first strlen(s2) characters will be
		 *     used to do the comparison.
		 * @param s2 String 2. Its whole string will be used to do the comparison.
		 * @return True if the comparison is equal. False if not equal.
		 */
		static bool sstrncmp(const char *s1, const char *s2);

		/** Get a string format of an integer, e.g. 123 will be "123".
		 * @param[in] n The input integer.
		 * @return The string format of the given integer.
		 */
		static string int2str(int n);

		/** Extract a single line from a multi-line string. The line
		 * delimiter is '\n'. The multi-line string moves forward one
		 * line. If it is the last line, move to the end of the
		 * string.
		 *
		 * @param[in,out] str  A multiple-line string.
		 * @return A single line. 
		 */
		static string get_line_from_string(char **str);

		/* Extract the float value from a variable=value string
		 * with format like "XYZ=1.1", 
		 * where 'str' is "XYZ=1.1"; 'float_var' is "XYZ=";
		 * 'p_val' points to float number 1.1.
		 *
		 * @param[in] str A string like "XYZ=1.1";
		 * @param[in] float_var The variable name "XYZ=".
		 * @param[out] p_val The pointer to the float number.
		 * @return True if the extraction succeeds; False if
		 *         extraction fails.
		 */
		static bool get_str_float(const char *str, const char *float_var, float *p_val);

		
		/* Extract the float values from a variable=value1,value2 string
		 * with format like "XYZ=1.1,1.2", 
		 * where 'str' is "XYZ=1.1,1.2"; 'float_var' is "XYZ=";
		 * 'p_v1' points to  1.1; 'p_v2' points to  1.2.
		 *
		 * @param[in] str A string like "XYZ=1.1,1.2";
		 * @param[in] float_var The variable name "XYZ=".
		 * @param[out] p_v1 The pointer to the first float.
		 * @param[out] p_v2 The pointer to the second float.
		 * @return True if the extraction succeeds; False if
		 *         extraction fails.
		 */
		static bool get_str_float(const char *str, const char *float_var,
								  float *p_v1, float *p_v2);
	
		/* Extract number of values and the float values, if any,
		 * from a string whose format is
		 * either "variable=value1,value2 " or "variable".
		 *
		 * for example, if the string is "XYZ=1.1,1.2", then
		 * 'str' is "XYZ=1.1,1.2"; 'float_var' is "XYZ";
		 * 'p_nvalues' points to 2. 'p_v1' points to 1.1; 'p_v2' points to 1.2.
		 * If the string is "XYZ", then 'str' is "XYZ"; 'float_var' is "XYZ".
		 * 'p_nvalues' points to 0. 'p_v1' and 'p_v2' unchanged.
		 * 
		 * @param[in] str A string like "XYZ=1.1,1.2" or "XYZ".
		 * @param[in] float_var The variable name "XYZ".
		 * @param[out] p_nvalues Number of values in the string.
		 * @param[out] p_v1 The pointer to the first float, if any.
		 * @param[out] p_v2 The pointer to the second float, if any.
		 * @return True if the extraction succeeds; False if
		 *         extraction fails.
		 */
		static bool get_str_float(const char *str, const char *float_var, 
								  int *p_nvalues, float *p_v1, float *p_v2);

		/* Extract the int value from a variable=value string
		 * with format like "XYZ=1", 
		 * where 'str' is "XYZ=1"; 'int_var' is "XYZ=";
		 * 'p_val' points to float number 1.
		 *
		 * @param[in] str A string like "XYZ=1";
		 * @param[in] int_var The variable name "XYZ=".
		 * @param[out] p_val The pointer to the int number.
		 * @return True if the extraction succeeds; False if
		 *         extraction fails.
		 */
		static bool get_str_int(const char *str, const char *int_var, int *p_val);
		
		/* Extract the int value from a variable=value1,value2 string
		 * with format like "XYZ=1,2", 
		 * where 'str' is "XYZ=1,2"; 'int_var' is "XYZ=";
		 * 'p_val' points to float number 1.
		 *
		 * @param[in] str A string like "XYZ=1";
		 * @param[in] int_var The variable name "XYZ=".
		 * @param[out] p_v1 The pointer to the first int.
		 * @param[out] p_v2 The pointer to the second int.
		 * @return True if the extraction succeeds; False if
		 *         extraction fails.
		 */
		static bool get_str_int(const char *str, const char *int_var,
								int *p_v1, int *p_v2);
		
		/* Extract number of values and the int values, if any,
		 * from a string whose format is
		 * either "variable=value1,value2 " or "variable".
		 *
		 * for example, if the string is "XYZ=1,2", then
		 * 'str' is "XYZ=1,2"; 'int_var' is "XYZ";
		 * 'p_nvalues' points to 2. 'p_v1' points to 1; 'p_v2' points to 2.
		 * If the string is "XYZ", then 'str' is "XYZ"; 'int_var' is "XYZ".
		 * 'p_nvalues' points to 0. 'p_v1' and 'p_v2' unchanged.
		 * 
		 * @param[in] str A string like "XYZ=1,2" or "XYZ".
		 * @param[in] int_var The variable name "XYZ".
		 * @param[out] p_nvalues Number of values in the string.
		 * @param[out] p_v1 The pointer to the first int, if any.
		 * @param[out] p_v2 The pointer to the second int, if any.
		 * @return True if the extraction succeeds; False if
		 *         extraction fails.
		 */
		static bool get_str_int(const char *str, const char *int_var, 
								int *p_nvalues, int *p_v1, int *p_v2);

		/** Change a file's extension and return the new filename.
		 * If the given new extension is empty, the old filename is
		 * not changed. If the old filename has no extension, add the
		 * new extension to it.
		 *
		 * @param[in] old_filename Old filename.
		 * @param[in] new_ext  The new extension. It shouldn't have
		 * ".". e.g., for MRC file, it will be "mrc", not ".mrc".
		 * @return The new filename with the new extension.
		 */
		static string change_filename_ext(const string& old_filename,
										  const string & new_ext);

		/** Remove a filename's extension and return the new filename.
		 *
		 * @param[in] filename The old filename whose extension is
		 * going to be removed.
		 * @return The new filename without extension.
		 */
		static string remove_filename_ext(const string& filename);
		
		/** Get a filename's extension.
		 * @param[in] filename A given filename.
		 * @return The filename's extension, or empty string if the
		 * file has no extension.
		 */
		static string get_filename_ext(const string& filename);

		/** Get a filename's basename. For example, the basename of
		 * "hello.c" is still "hello.c"; The basename of
		 * "/tmp/abc/hello.c" is "hello.c".
		 * @param[in] filename The given filename, full path or relative path.
		 * @return The basename of the filename.
		 */
		static string sbasename(const string & filename);
		
		/** calculate the least square fit value.
		 *
		 * @param[in] nitems Number of items in array data_x and data_y.
		 * @param[in] data_x x data array.
		 * @param[in] data_y y data array. It should have the same
		 *        number of items to data_x.
		 * @param[out] p_slope pointer to the result slope.
		 * @param[out] p_intercept pointer to the result intercept.
		 * @param[in] ignore_zero If true, ignore data where either x
		 *        or y is 0. If false, includes all 0.
		 * @param[in] absmax Ignores values in y more than absmax from zero
		 */
		static void calc_least_square_fit(size_t nitems, const float *data_x,
										  const float *data_y, float *p_slope, 
										  float *p_intercept, bool ignore_zero,float absmax=0);

		/** Save (x y) data array into a file. Each line of the file
		 * have the format "x1TABy1", where x1, y1 are elements of x
		 * array and y array. The x, y arrays must have the same 
		 * number of items.
		 *
		 * @param[in] x_array The x array.
		 * @param[in] y_array The y array.
		 * @param[in] filename The output filename.
		 */
		 
		static void save_data(const vector < float >&x_array, 
							  const vector < float >&y_array,
							  const string & filename);

		/** Save x, y data into a file. Each line of the file have the
		 * format "x1TABy1", where x1 = x0 + dx*i; y1 = y_array[i].
		 *
		 * @param[in] x0 The starting point of x.
		 * @param[in] dx delta x. The increase step of x data.
		 * @param[in] y_array The y data array.
		 * @param[in] filename The output filename.
		 */
		static void save_data(float x0, float dx,
							  const vector < float >&y_array,
							  const string & filename);
		
		/** Save x, y data into a file. Each line of the file have the
		 * format "x1TABy1", where x1 = x0 + dx*i; y1 = y_array[i].
		 *
		 * @param[in] x0 The starting point of x.
		 * @param[in] dx delta x. The increase step of x data.
		 * @param[in] y_array The y data array.
		 * @param[in] array_size The y data array size.
		 * @param[in] filename The output filename.
		 */
		static void save_data(float x0, float dx, float *y_array, 
							  size_t array_size, const string & filename);

		/** Given a tabulated function y of x (n unordered points), and
		 * Given the values of the m values xq to be interpolated
		 * This routine returns the interpolated array yq, PRB
		 * This function is called by splint 
		 * @param[in] y of x is the tabulated function of length n
		 * @param[in] xq is the x values to be splined: has m points.
		 * @param[out]  yq are the splined values
		*/
		static void spline_mat(float *x, float *y, int n,  float *xq, float *yq, int m); //PRB

		/** Given a tabulated function y of x (unordered), and
		 * Given the values of the first derivatives at the end points
		 * This routine returns an array y2, that contains the second derivatives
		 * of the function at the tabulated points.   PRB
		 * This function is called by splint 
		 * @param[in] y of x is the tabulated function of length n
		 * @param[in] yp1, ypn : the derivatives of the first and last point.
		 * @param[out]  y2 is the value of the second derivatives
		*/
		static void spline(float *x, float *y, int n, float yp1, float ypn, float *y2);

		/** Given the arrays xa(ordered, ya of length n, which tabulate a function 
		 *  and given the array y2a which is the output of spline and an unordered array xq,
		 *  this routine returns a cubic-spline interpolated array yq.
		 * @param[in] y of x is the tabulated function of length n
		 * @param[in] y2a is returned from spline: second derivs
		 * @param[in] xq is the x values to be splined: has m points.
		 * @param[out]  yq are the splined values
		*/  // PRB
		static void splint( float *xa, float *ya, float *y2a, int n,  
		                                     float *xq, float *yq, int m);
		
		/** does a sort as in Matlab. Carries along the Permutation
		 * matrix
		 * @param[in] The array [left .. right] is sorted
		 * @param[in] The array [leftPerm rightPerm] is shuffled due to the sorting
		 * @param[out] Both arrays  are reshuffled.
		 */
		static void sort_mat(float *left, float *right, int *leftPerm, int
 				*rightPerm);

		/** list the sorted lengths of the  integer lattice 
		 * sites of a square sided image of  size Size. PRB
		 * @param[in] Size The length of the image 
		 * @param[out] PermMatTr  The matrix telling the ordering of the lattice 
		 *                        sites wrt the array
		 * @param[out] weightofkvalsSorted the number of sites at that distance
		 */
		 static void Radialize(int *PermMatTr,  float * kValsSorted,
		            float *weightofkvalsSorted, int Size, int *SizeReturned);
				     
		
		/** Get a float random number between low and high.
		 * @param[in] low The lower bound of the random number.
		 * @param[in] high The upper bound of the random number.
		 * @return The random number between low and high.
		 */
		static float get_frand(int low, int high);
		
		/** Get a float random number between low and high.
		 * @param[in] low The lower bound of the random number.
		 * @param[in] high The upper bound of the random number.
		 * @return The random number between low and high.
		 */
		static float get_frand(float low, float high);

		/** Get a float random number between low and high.
		 * @param[in] low The lower bound of the random number.
		 * @param[in] high The upper bound of the random number.
		 * @return The random number between low and high.
		 */
		static float get_frand(double low, double high);

		/** Get a Gaussian random number.
		 *
		 * @param[in] mean The gaussian mean
		 * @param[in] sigma The gaussian sigma
		 * @return the gaussian random number.
		 */
		static float get_gauss_rand(float mean, float sigma);


		/** Get ceiling round of a float number x.
		 * @param[in] x Given float number.
		 * @return Ceiling round of x.
		 */
		static inline int round(float x)
		{
			if (x < 0) {
				return (int) (x - 0.5f);
			}
			return (int) (x + 0.5f);
		}

		/** Get ceiling round of a float number x.
		 * @param[in] x Given float number.
		 * @return Ceiling round of x.
		 */
	   	static inline int round(double x)
		{
			if (x < 0) {
				return (int) (x - 0.5);
			}
			return (int) (x + 0.5);
		}

		/** Calculate bilinear interpolation.
		 * @param[in] p1 The first number. corresponding to (x0,y0).
		 * @param[in] p2 The second number. corresponding to (x1,y0).
		 * @param[in] p3 The third number. corresponding to (x1,y1).
		 * @param[in] p4 The fourth number. corresponding to (x0,y1).
		 * @param[in]  t t
		 * @param[in]  u u
		 * @return The bilinear interpolation value.
		 */
		static inline float bilinear_interpolate(float p1, float p2, float p3,
												 float p4, float t, float u)
		{
			return (1-t) * (1-u) * p1 + t * (1-u) * p2 + t * u * p3 + (1-t) * u * p4;
		}

		/** Calculate trilinear interpolation.
		 *
		 * @param[in] p1 The first number. corresponding to (x0,y0,z0).
		 * @param[in] p2 The second number. corresponding to (x1,y0,z0).
		 * @param[in] p3 The third number. corresponding to (x0,y1, z0).
		 * @param[in] p4 The fourth number. corresponding to (x1,y1,z0).
		 * @param[in] p5 The fifth number. corresponding to (x0,y0,z1).
		 * @param[in] p6 The sixth number. corresponding to (x1,y0,z1).
		 * @param[in] p7 The seventh number. corresponding to (x0,y1,z1).
		 * @param[in] p8 The eighth number. corresponding to (x1,y1,z1).
		 * @param[in] t t
		 * @param[in] u u
		 * @param[in] v v
		 * @return The trilinear interpolation value.
		 */
		static inline float trilinear_interpolate(float p1, float p2, float p3,
												  float p4, float p5, float p6, 
												  float p7, float p8, float t,
												  float u, float v)
		{
			return ((1 - t) * (1 - u) * (1 - v) * p1 + t * (1 - u) * (1 - v) * p2
					+ (1 - t) * u * (1 - v) * p3 + t * u * (1 - v) * p4
					+ (1 - t) * (1 - u) * v * p5 + t * (1 - u) * v * p6
					+ (1 - t) * u * v * p7 + t * u * v * p8);
		}

		/** Find the maximum value and (optional) its index in an array.
		 * @param[in] data data array.
		 * @param[in] nitems number of items in the data array.
		 * @param[out] p_max_val pointer to the maximum value.
		 * @param[out] p_max_index pointer to index of the maximum value.
		 */
		static void find_max(const float *data, size_t nitems,
							 float *p_max_val, int *p_max_index = 0);

		/** Find the maximum value and (optional) its index, minimum
		 * value and (optional) its index in an array.
		 *
		 * @param[in] data data array.
		 * @param[in] nitems number of items in the data array.
		 * @param[out] p_max_val  pointer to the maximum value.
		 * @param[out] p_min_val  pointer to the minimum value.
		 * @param[out] p_max_index pointer to index of the maximum value.
		 * @param[out] p_min_index pointer to index of the minimum value.
		 */
		static void find_min_and_max(const float *data, size_t nitems,
									 float *p_max_val, float *p_min_val,
									 int *p_max_index = 0, int *p_min_index = 0);

		/** Search the best FFT size with good primes. It supports
		 * FFT size up to 4096 now.
		 *
		 * @param[in] low low size the search starts with.
		 * @return The best FFT size.
		 */
		static int calc_best_fft_size(int low);

		/** Calculate a number's square.
		 * @param[in] n Given number.
		 * @return (n*n).
		 */
		static int square(int n)
		{
			return (n * n);
		}
		
		/** Calculate a number's square.
		 * @param[in] x Given number.
		 * @return (x*x).
		 */
		static float square(float x)
		{
			return (x * x);
		}

		/** Calculate a number's square.
		 * @param[in] x Given number.
		 * @return (x*x).
		 */
		static float square(double x)
		{
			return (float)(x * x);
		}

		/** Calcuate (x*x + y*y).
		 * @param[in] x The first number.
		 * @param[in] y The second number.
		 * @return (x*x + y*y).
		 */
		static float square_sum(float x, float y)
		{
			return (float)(x * x + y * y);
		}

		/** Euclidean distance function in 3D: f(x,y,z) = sqrt(x*x + y*y + z*z);
		 * @param[in] x The first number.
		 * @param[in] y The second number.
		 * @param[in] z The third number.
		 * @return sqrt(x*x + y*y + z*z);
		 */
		static inline float hypot3(int x, int y, int z)
		{
			return sqrtf((float)(x * x + y * y + z * z));
		}
		
		/** Euclidean distance function in 3D: f(x,y,z) = sqrt(x*x + y*y + z*z);
		 * @param[in] x The first number.
		 * @param[in] y The second number.
		 * @param[in] z The third number.
		 * @return sqrt(x*x + y*y + z*z);
		 */
		static inline float hypot3(float x, float y, float z)
		{
			return sqrtf(x * x + y * y + z * z);
		}
		
		/** Euclidean distance function in 3D: f(x,y,z) = sqrt(x*x + y*y + z*z);
		 * @param[in] x The first number.
		 * @param[in] y The second number.
		 * @param[in] z The third number.
		 * @return sqrt(x*x + y*y + z*z);
		 */
		static inline float hypot3(double x, double y, double z)
		{
			return (float) sqrt(x * x + y * y + z * z);
		}

		/** A fast way to calculate a floor, which is largest integral
		 * value not greater than argument.
		 *
		 * @param[in] x A float point number.
		 * @return floor of x.
		 */
		static inline int fast_floor(float x)
		{
			if (x < 0) {
				return ((int) x - 1);
			}
			return (int) x;
		}

		/** Calculate Gaussian value.
		 * @param[in] a 
		 * @param[in] dx 
		 * @param[in] dy 
		 * @param[in] dz 
		 * @param[in] d
		 * @return The Gaussian value.
		 */
		static inline float agauss(float a, float dx, float dy, float dz, float d)
		{
			return (a * exp(-(dx * dx + dy * dy + dz * dz) / d));
		}

		/** Get the minimum of 2 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @return The minimum of 2 numbers.
		 */
		static inline int get_min(int f1, int f2)
		{
			return (f1 < f2 ? f1 : f2);
		}
		
		/** Get the minimum of 3 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @param[in] f3 The third number.
		 * @return The minimum of 3 numbers.
		 */
		static inline int get_min(int f1, int f2, int f3)
		{
			if (f1 <= f2 && f1 <= f3) {
				return f1;
			}
			if (f2 <= f1 && f2 <= f3) {
				return f2;
			}
			return f3;
		}

		/** Get the minimum of 2 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @return The minimum of 2 numbers.
		 */
		static inline float get_min(float f1, float f2)
		{
			return (f1 < f2 ? f1 : f2);
		}
		
		/** Get the minimum of 3 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @param[in] f3 The third number.
		 * @return The minimum of 3 numbers.
		 */
		static inline float get_min(float f1, float f2, float f3)
		{
			if (f1 <= f2 && f1 <= f3) {
				return f1;
			}
			if (f2 <= f1 && f2 <= f3) {
				return f2;
			}
			return f3;
		}

		
		/** Get the minimum of 4 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @param[in] f3 The third number.
		 * @param[in] f4 The fourth number.
		 * @return The minimum of 4 numbers.
		 */
		static inline float get_min(float f1, float f2, float f3, float f4)
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
		
		/** Get the maximum of 2 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @return The maximum of 2 numbers.
		 */
		static inline float get_max(float f1, float f2)
		{
			return (f1 < f2 ? f2 : f1);
		}
		
		/** Get the maximum of 3 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @param[in] f3 The third number.
		 * @return The maximum of 3 numbers.
		 */
		static inline float get_max(float f1, float f2, float f3)
		{
			if (f1 >= f2 && f1 >= f3) {
				return f1;
			}
			if (f2 >= f1 && f2 >= f3) {
				return f2;
			}
			return f3;
		}
		
		/** Get the maximum of 4 numbers.
		 * @param[in] f1 The first number.
		 * @param[in] f2 The second number.
		 * @param[in] f3 The third number.
		 * @param[in] f4 The fourth number.
		 * @return The maximum of 4 numbers.
		 */
		static inline float get_max(float f1, float f2, float f3, float f4)
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

		/** Calculate the difference of 2 angles and makes the
		 * equivalent result to be less than Pi.
		 *
		 * @param[in] x The first angle.
		 * @param[in] y The second angle.
		 * @return The difference of 2 angles.
		 */
		static inline float angle_sub_2pi(float x, float y)
		{
			float r = fmod(fabs(x - y), (float) (2 * M_PI));
			if (r > M_PI) {
				r = (float) (2.0 * M_PI - r);
			}

			return r;
		}

		/** Calculate the difference of 2 angles and makes the
		 * equivalent result to be less than Pi/2.
		 *
		 * @param[in] x The first angle.
		 * @param[in] y The second angle.
		 * @return The difference of 2 angles.
		 */
		static inline float angle_sub_pi(float x, float y)
		{
			float r = fmod(fabs(x - y), (float) M_PI);
			if (r > M_PI / 2.0) {
				r = (float)(M_PI - r);
			}
			return r;
		}

		/** Check whether a number is a good float. A number is a good
		 * float if it is not abnormal zero, and not inf, -inf or NaN
		 *
		 * @param[in] p_f Pointer to the given float number.
		 * @return 1 if it is a good float; 0 if it's not good float.
		 */
		static inline int goodf(const float *p_f)
		{
			//This is the old way to judge a good float, which cause problem on 
			//Fedora Core 64 bit system. Now use isinff() and isnanf() on Linux.
			#if defined(_WIN32) || defined(__APPLE__)
			// the first is abnormal zero the second is +-inf or NaN 
			if ((((int *) p_f)[0] & 0x7f800000) == 0 ||
				(((int *) p_f)[0] & 0x7f800000) == 255) {
				return 0;
			}
			#else
			if(isinff(*p_f) || isnanf(*p_f)) return 0;
			#endif	//_WIN32 || __APPLE__
			
			return 1;
		}

		/** Get the current time in a string with format "mm/dd/yyyy hh:mm".
		 * @return The current time string.
		 */
		static string get_time_label();

		/** Set program logging level through command line option "-v N",
		 * where N is the level. 
		 *
		 * @param[in] argc Number of arguments.
		 * @param[in] argv Argument arrays.		 
		 */
		static void set_log_level(int argc, char *argv[]);
		
		/** copy sign of a number. return a value whose absolute value
		 * matches that of 'a', but whose sign matches that of 'b'.  If 'a'
		 * is a NaN, then a NaN with the sign of 'b' is returned.
		 *
		 * It is exactly copysign() on non-Windows system.		 
		 *
		 * @param[in] a The first number.
		 * @param[in] b The second number.
		 * @return Copy sign of a number.
		 */
		static inline float eman_copysign(float a, float b)
		{
#ifndef WIN32
			return copysign(a, b);
#else
			int flip = -1;
			if ((a <= 0 && b <= 0) || (a > 0 && b > 0)) {
				flip = 1;
			}
			return a * flip;
#endif
		}

		/** complementary error function. It is exactly erfc() on
		 * non-Windows system. On Windows, it tries to simulate erfc().
		 *
		 * The erf() function returns the error function of x; defined as
		 * erf(x) = 2/sqrt(pi)* integral from 0 to x of exp(-t*t) dt
		 *
		 * The erfc() function returns the complementary error function of x, that
		 * is 1.0 - erf(x).
		 *
		 * @param[in] x A float number.
		 * @return The complementary error function of x.
		 */
		static inline float eman_erfc(float x)
		{
#ifndef WIN32
			return (float)erfc(x);
#else
			static double a[] = { -1.26551223, 1.00002368,
								  0.37409196, 0.09678418,
								  -0.18628806, 0.27886807,
								  -1.13520398, 1.48851587,
								  -0.82215223, 0.17087277
			};

			double result = 1;
			double z = fabs(x);
			if (z > 0) {
				double t = 1 / (1 + 0.5 * z);
				double f1 = t * (a[4] + t * (a[5] + t * (a[6] +
							t * (a[7] + t * (a[8] + t * a[9])))));
				result = t * exp((-z * z) + a[0] + t * (a[1] + t * (a[2] + t * (a[3] + f1))));

				if (x < 0) {
					result = 2 - result;
				}
			}
			return (float)result;
#endif
		}

		/** Print a 3D integer matrix to a file stream (std out by default).
		 *
		 * @param[in] mat integer 3-d multi_array reference
		 * @param out Output stream; cout by default.
		 * @param[in] str Message string to be printed.
		 */
		static void printMatI3D(MIArray3D& mat, 
								const string str = string(""),
								ostream& out = std::cout);

		/** Compute a vector containing evenly-spaced euler angles
		 *
		 * The order of angles in the vector is phi, theta, psi.
		 * Note: We use the Spider construction for Euler angles here.
		 *
		 * @param[in] delta  Delta theta (spacing in theta).
		 * @param[in] t1  Starting (min) value of theta in degrees, default = 0.
		 * @param[in] t2  Ending (max) value of theta in degrees, default = 90.
		 * @param[in] p1  Starting (min) value of phi in degrees, default = 0.
		 * @param[in] p2  Ending (max) value of phi in degrees, default = 359.9.
		 * @return Vector of angles as a flat list of phi_0, theta_0, psi_0, ..., phi_N, theta_N, psi_N.
		 */
		static vector<float>
		voea(float delta, float t1=0, float t2=90, 
			 float p1=0, float p2=359.9);

		/** Tri-Quadratic interpolation.
		 *
		 *	@param[in] r x-coord value
		 *	@param[in] s y-coord value
		 *	@param[in] t z-coord value
		 *	@param[in] f 3x3x3 grid of measured values
		 *
		 *	@return Interpolated value
		 */
		static float triquad(double r, double s, double t, float f[]);

		/** Quadratic interpolation (2D).
		 *
		 *  Note:  This routine starts counting from 1, not 0!
		 *
		 *	This routine uses six image points for interpolation:
		 *
		 *@see M. Abramowitz & I.E. Stegun, Handbook of Mathematical
		 *     Functions (Dover, New York, 1964), Sec. 25.2.67.
		 *     http://www.math.sfu.ca/~cbm/aands/page_882.htm.
		 *
		 *@see http://www.cl.cam.ac.uk/users/nad/pubs/quad.pdf
		 *
		 *@verbatim
                f3    fc
                |
                | x
         f2-----f0----f1
                |
                |
                f4
		 *@endverbatim
		 *
		 *	f0 - f4 are image values near the interpolated point X.
		 *	f0 is the interior mesh point nearest x.
		 *
		 *	Coords:  
		 *@li        f0 = (x0, y0)
		 *@li        f1 = (xb, y0)
		 *@li        f2 = (xa, y0)
		 *@li        f3 = (x0, yb)
		 *@li        f4 = (x0, ya)
		 *@li        fc = (xc, yc)
		 *
		 *	Mesh spacings: 
		 *@li              hxa -- x- mesh spacing to the left of f0
		 *@li              hxb -- x- mesh spacing to the right of f0
		 *@li              hyb -- y- mesh spacing above f0
		 *@li              hya -- y- mesh spacing below f0
		 *
		 *	Interpolant:
		 *	  f = f0 + c1*(x-x0) + c2*(x-x0)*(x-x1)
		 *			 + c3*(y-y0) + c4*(y-y0)*(y-y1)
		 *			 + c5*(x-x0)*(y-y0)
		 *
		 *	@param[in] x x-coord value
		 *	@param[in] y y-coord value
		 *	@param[in] image Image object (pointer)
		 *
		 *	@return Interpolated value
		 */
		static float quadri(float x, float y, int nx, int ny, float* image);

		/** 1-D Kaiser-Bessel window function class.
		 *  (It's a class so that the windowing parameters may be
		 *   instantiated and held in the instance object.)
		 *
		 *  The I0 version can be tabulated and interpolated upon
		 *  demand, but the max error needs to be checked.  The 
		 *  "vtable" parameter corresponds to the maximum value of x
		 *  for which the I0 window is non-zero.  Setting "vtable"
		 *  different from "v" corresponds to a change in units of x.
		 *  In practice, it is often handy to replace x in some sort
		 *  of absolute units with x described in terms of grid
		 *  intervals.
		
		 *
		 *  The get_kbsinh_win and get_kbi0_win functions return
		 *  single-argument function objects, which is what a 
		 *  generic routine is likely to want.
		 *
		 *  @see P. A. Penczek, R. Renka, and H. Schomberg,
		 *      J. Opt. Soc. Am. _21_, 449 (2004)
		 * 
		 */
		class KaiserBessel 
		{
			protected:
				float alpha, v, r; /** Kaiser-Bessel parameters */
				int N; /** size in Ix-space */
				int K; /** I0 window size */
				float vtable; /** table I0 non-zero domain maximum */
				int ntable;
				vector<float> i0table;
				float dtable; /** table spacing */
				float alphar; /** alpha*r */
				float fac; /** 2*pi*alpha*r*v */
				float vadjust; 
				float facadj; /** 2*pi*alpha*r*vadjust */
				virtual void build_I0table(); /** Tabulate I0 window for speed */
				float fltb;
			public:
				KaiserBessel(float alpha_, int K, float r_,
						     float v_, int N_, float vtable_=0.f, 
							 int ntable_ = 5999);
				virtual ~KaiserBessel() {};
				/** Compute the maximum error in the table */
				float I0table_maxerror();
				vector<float> dump_table() {
					return i0table;
				}
				/** Kaiser-Bessel Sinh window function */
				virtual float sinhwin(float x) const;
				/** Kaiser-Bessel I0 window function */
				virtual float i0win(float x) const;
				/** Kaiser-Bessel I0 window function (uses table lookup) */
				float i0win_tab(float x) const {
					float absx = fabs(x);
					int loc = int(round(absx*fltb));
					return i0table[loc];
#if 0 // old version
					if (absx > vtable) return 0.f;
					float loc = absx/dtable;
					return i0table[int(loc + 0.5f)];
#endif // 0
				}
				/** Return the size of the I0 window */
				int get_window_size() const { return K; }
				/** Sinh window function object */
				class kbsinh_win {
					KaiserBessel& kb;
					public:
					kbsinh_win(KaiserBessel& kb_) : kb(kb_) {}
					float operator()(float x) const {
						return kb.sinhwin(x);
					}
					int get_window_size() const {return kb.get_window_size();}
				};
				/** Sinh window function object factory */
				kbsinh_win get_kbsinh_win() { 
					return kbsinh_win(*this);
				}
				/** I0 window function object */
				class kbi0_win {
					KaiserBessel& kb;
					public:
					kbi0_win(KaiserBessel& kb_) : kb(kb_) {}
					float operator()(float x) const {
						return kb.i0win(x);
					}
					int get_window_size() const {return kb.get_window_size();}
				};
				/** I0 window function object factory */
				kbi0_win get_kbi0_win() { 
					return kbi0_win(*this);
				}
		};

		class FakeKaiserBessel : public KaiserBessel {
			public:
				FakeKaiserBessel(float alpha, int K, float r_,
						         float v_, int N_, float vtable_=0.f,
								 int ntable_ = 5999) 
		        : KaiserBessel(alpha, K, r_, v_, N_, vtable_, ntable_) {
					build_I0table();
				}
				float sinhwin(float x) const;
				float i0win(float x) const;
				void build_I0table();
		};

		/** Gaussian function class.
		 *
		 *  Usage:
		 *
		 *     Gaussian gauss(sigma);
		 *     float g = gauss(x);
		 */
		class Gaussian {
			float sigma;
			float rttwopisigma;
			float twosigma2;
			public:
			Gaussian(float sigma_ = 1.0) : sigma(sigma_) {
				rttwopisigma = sqrtf(static_cast<float>(twopi)*sigma);
				twosigma2 = 2*sigma*sigma;
			}
			inline float operator()(float x) const {
				return exp(-x*x/(twosigma2))/rttwopisigma;
			}
		};

		/** Sign function
		 *
		 *	@param[in] val value who's sign is to be checked
		 *
		 *	@return +1 if the value is positive, -1 if the value is negative.
		 */
		template<class T> static inline T sgn(T& val) {
			return (val > 0) ? T(+1) : T(-1);
		}
		static void alrq(float *xim,  int nsam , int nrow , int *numr,
                                 float *circ, int lcirc, int nring, char mode);
                static EMData* Polar2D(EMData* image, vector<int> numr, string mode);
                static void alrq_ms(float *xim, int    nsam, int  nrow, float cns2, float cnr2,
                                    int  *numr, float *circ, int lcirc, int  nring, char  mode);
                static EMData* Polar2Dm(EMData* image, float cns2, float cnr2, vector<int> numr, string mode);
                static void  fftr_q(float  *xcmplx, int nv);
                static void  fftr_d(double *xcmplx, int nv);
                static void  fftc_q(float  *br, float  *bi, int ln, int ks);
                static void  fftc_d(double *br, double *bi, int ln, int ks);
                static void  Frngs(EMData* circ, vector<int> numr);
                static void  frngs(float *circ, int *numr, int nring);
		static boost::tuple<double, float, int>
			Crosrng_e(EMData* circ1, EMData* circ2, vector<int> numr, int neg);
                static void  crosrng_e(float *circ1, float *circ2, int lcirc,
                                       int    nring, int   maxrin, int *numr,
                                       double *qn, float *tot, int neg);
//                static void  Crosrng_ms(EMData* circ1, EMData* circ2, vector<int> numr,
//                                        double *qn, float *tot, double   *qm, double *tmt);
		static boost::tuple<double, float, double, float>
                         Crosrng_ms(EMData* circ1, EMData* circ2, vector<int> numr);
                static void  crosrng_ms(float *circ1, float *circ2, int  lcirc, int  nring,
                                        int   maxrin, int   *numr , double *qn, float *tot,
                                        double   *qm, float *tmt);
                static void  prb1d(double *b, int npoint, float *pos);
	};
}

#endif

/* vim: set ts=4 noet nospell: */
