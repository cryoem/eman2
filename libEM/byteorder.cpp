/**
 * $Id$
 */
#include "byteorder.h"

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
