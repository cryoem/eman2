/**
 * $Id$
 */
 
/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */
 
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
