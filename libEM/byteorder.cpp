/**
 * $Id$
 */
#include "byteorder.h"
#include "util.h"

using namespace EMAN;

bool ByteOrder::is_host_endian_checked = false;
bool ByteOrder::host_big_endian = false;

bool ByteOrder::is_host_big_endian()
{
	if (!is_host_endian_checked) {
		int one = 1;
		char *p_one = (char *) (&one);

		if (p_one[0] == 1 && p_one[1] == 0 && p_one[2] == 0 && p_one[3] == 0) {
			host_big_endian = false;
		}
		else {
			host_big_endian = true;
		}

		is_host_endian_checked = true;
	}

	return host_big_endian;
}

bool ByteOrder::is_float_big_endian(float f)
{
	bool is_big = false;
	
	if (Util::goodf(&f) && f > 0 && f < 65535.0 && f == floor(f)) {
		is_big = is_host_big_endian();
	}
	else {
		is_big = !is_host_big_endian();
	}

	return is_big;
}
