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

#include <cstdio>
#include "geometry.h"

using namespace EMAN;

IntPoint EMAN::operator -( const IntPoint& p)
{
	return IntPoint(-p[0],-p[1],-p[2]);
}

bool EMAN::operator<(const Pixel& p1, const Pixel& p2)
{
	if (p1.value < p2.value) {
		return true;
	}
	return false;
}

bool EMAN::operator==(const Pixel& p1, const Pixel& p2)
{
	if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z && p1.value == p2.value) {
		return true;
	}
	return false;
}

bool EMAN::operator!=(const Pixel& p1, const Pixel& p2)
{
	return !(p1 == p2);
}


bool Region::inside_region() const
{
	if (size[0] >= 0 && size[1] >= 0 && size[2] >= 0) {
		return true;
	}

	return false;
}

bool Region::inside_region(const FloatPoint & p) const
{
	if (p.get_ndim() == 1) {
		return inside_region(p[0]);
	}


	if (p.get_ndim() == 2) {
		return inside_region(p[0], p[1]);
	}

	return inside_region(p[0], p[1], p[2]);
}

bool Region::inside_region(float x) const
{
	if (size[0] >= 0 && origin[0] <= x &&
		(origin[0] + size[0]) > x ) {
		return true;
	}
	return false;
}


bool Region::inside_region(float x, float y) const
{
	if (size[0] >= 0 && size[1] >= 0 &&
		origin[0] <= x && origin[1] <= y &&
		(origin[0] + size[0]) > x && (origin[1] + size[1]) > y) {
		return true;
	}
	return false;
}

#include <iostream>
using std::cout;
using std::endl;

bool Region::inside_region(float x, float y, float z) const
{
	if (size[0] >= 0 && size[1] >= 0 && size[2] >= 0 &&
		origin[0] <= x && origin[1] <= y && origin[2] <= z &&
		(origin[0] + size[0]) > x &&
		(origin[1] + size[1]) > y && (origin[2] + size[2]) > z) {
		return true;
	}
	return false;
}


bool Region::is_region_in_box(const FloatSize & box) const
{
	if (size[0] >= 0 && size[1] >= 0 && size[2] >= 0 &&
		origin[0] >= 0 && origin[1] >= 0 && origin[2] >= 0 &&
		(origin[0] + size[0]) <= box[0] &&
		(origin[1] + size[1]) <= box[1] &&
		(origin[2] + size[2]) <= box[2]) {
		return true;
	}

	return false;
}



string Region::get_string() const
{
	char str[1028];
	int ndim = origin.get_ndim();

	if (ndim == 2) {
		sprintf(str, "(%2.1f, %2.1f; %2.1f, %2.1f)",
				origin[0], origin[1], size[0], size[1]);
	}
	else if (ndim == 3) {
		sprintf(str, "(%2.1f, %2.1f, %2.1f; %2.1f, %2.1f, %2.1f)",
				origin[0], origin[1], origin[2], size[0], size[1], size[2]);
	}

	return string(str);
}
