#ifndef __util_h__
#define __util_h__

#include <stdio.h>
#include <string>

using std::string;

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
    };
}

#endif
