/**
 * $Id$
 */
#include "byteorder.h"

using namespace EMAN;

bool ByteOrder::is_machine_endian_checked = false;
bool ByteOrder::machine_big_endian = false;

bool ByteOrder::is_machine_big_endian()
{
    if (!is_machine_endian_checked) {
	int one = 1;
	char* p_one = (char*) (&one);
	
	if(p_one[0] == 1 && p_one[1] == 0 && p_one[2] == 0 && p_one[3] == 0) {
	    machine_big_endian = false;
	}
	else {
	    machine_big_endian = true;
	}

	is_machine_endian_checked = true;
    }

    return machine_big_endian;
}
