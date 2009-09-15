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

#include "xydata.h"
#include <algorithm>
#ifndef WIN32
#include <sys/param.h>
#else
#include <windows.h>
#define MAXPATHLEN (MAX_PATH*4)
#endif
#include <cfloat>
#include <cmath>
#include <cstdio>
#include "log.h"
#include "exception.h"

using namespace EMAN;

XYData::XYData()
	: ymin(FLT_MAX), ymax(-FLT_MAX), mean_x_spacing(1.0)
{
}

void XYData::update()
{
	std::sort(data.begin(), data.end());

	ymin = FLT_MAX;
	ymax = -FLT_MAX;

	typedef vector < Pair >::const_iterator Ptype;
	for (Ptype p = data.begin(); p != data.end(); p++) {
		if (p->y > ymax) {
			ymax = p->y;
		}
		if (p->y < ymin) {
			ymin = p->y;
		}
	}

	size_t n = data.size();
	mean_x_spacing = (data[n - 1].x - data[0].x) / (float) n;
}


int XYData::read_file(const string & filename)
{
	FILE *in = fopen(filename.c_str(), "rb");
	if (!in) {
		LOGERR("cannot open xydata file '%s'", filename.c_str());
		return 1;
	}

	char buf[MAXPATHLEN];
	char tmp_str[MAXPATHLEN];

	while (fgets(buf, MAXPATHLEN, in)) {
		if (buf[0] != '#') {
			float x = 0;
			float y = 0;

			if (sscanf(buf, " %f%[^.0-9-]%f", &x, tmp_str, &y) != 3) {
				break;
			}
			data.push_back(Pair(x, y));
		}
	}

	fclose(in);
	in = 0;

	update();

	return 0;
}

int XYData::write_file(const string & filename) const
{
	FILE *out = fopen(filename.c_str(), "wb");
	if (!out) {
		LOGERR("cannot open xydata file '%s' to write", filename.c_str());
		return 1;
	}

	typedef vector < Pair >::const_iterator Ptype;
	for (Ptype p = data.begin(); p != data.end(); p++) {
		fprintf(out, "%1.6g\t%1.6g\n", p->x, p->y);
	}

	fclose(out);
	out = 0;

	return 0;
}


float XYData::calc_correlation(XYData * xy, float minx, float maxx) const
{
	size_t n = data.size();
	float x0 = data[0].x;
	float xn = data[n - 1].x;

	if (maxx <= minx || minx >= xn || maxx <= x0) {
		LOGERR("incorrect minx, maxx=%f,%f for this XYData range [%f,%f]",
			   minx, maxx, x0, xn);
		return 0;
	}

	float scc = 0;
	float norm1 = 0;
	float norm2 = 0;

	xy->update();
	for (size_t i = 0; i < n; i++) {
		float x = data[i].x;
		if (x >= minx && x <= maxx && xy->is_validx(x)) {
			float selfy = data[i].y;
			float xyy = xy->get_yatx(x);

			scc += selfy * xyy;
			norm1 += selfy * selfy;
			norm2 += xyy * xyy;
		}
	}

	float result = scc / sqrt(norm1 * norm2);
	return result;
}


float XYData::get_yatx(float x)
{
// You MUST be kidding. Do you know how expensive update() is !?!?
// the correct answer is to call update after you change the data, just like with EMData objects. --steve
//	update();	//update to set the mean_x_spacing value

	if (data.size()==0 || mean_x_spacing==0) return 0.0;

	int s = (int) floor((x - data[0].x) / mean_x_spacing);
	int nx = (int) data.size();

	if (data[s].x > x) {
		while (data[s].x > x && s >= 0) {
			s--;
		}
	}
	else if (data[s + 1].x < x) {
		while (data[s + 1].x < x && s < (nx - 1)) {
			s++;
		}
	}

	if (s >= nx || s < 0) {
		return 0.0;
	}

	float f = (x - data[s].x) / (data[s + 1].x - data[s].x);
	float y = data[s].y * (1 - f) + data[s + 1].y * f;
	return y;
}

void XYData::set_xy_list(const vector<float>& xlist, const vector<float>& ylist)
{
	if(xlist.size() != ylist.size()) {
		throw InvalidParameterException("xlist and ylist size does not match!");
	}

	for(unsigned int i=0; i<xlist.size(); ++i) {
		data.push_back(Pair(xlist[i], ylist[i]));
	}
}

void XYData::set_size(size_t n)
{
	data.resize(n, Pair(0.0f, 0.0f));
}

vector<float> XYData::get_xlist() const
{
	vector<float> xlist;
	vector<Pair>::const_iterator cit;
	for(cit=data.begin(); cit!=data.end(); ++cit) {
		xlist.push_back( (*cit).x);
	}

	return xlist;
}

vector<float> XYData::get_ylist() const
{
	vector<float> ylist;
	vector<Pair>::const_iterator cit;
	for(cit=data.begin(); cit!=data.end(); ++cit) {
		ylist.push_back( (*cit).y);
	}

	return ylist;
}
