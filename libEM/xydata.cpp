/**
 * $Id$
 */
#include "xydata.h"
#include <algorithm>
#ifndef WIN32
#include <sys/param.h>
#else
#define MAXPATHLEN 1024
#endif
#include <float.h>
#include <math.h>
#include "log.h"

using namespace EMAN;

XYData::XYData()
{
	ymin = FLT_MAX;
	ymax = -FLT_MAX;
	mean_x_spacing = 0;
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


int XYData::read_file(string filename)
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

int XYData::write_file(string filename) const
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


float XYData::get_yatx(float x) const
{
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
		return 0;
	}

	float f = (x - data[s].x) / (data[s + 1].x - data[s].x);
	float y = data[s].y * (1 - f) + data[s + 1].y * f;
	return y;
}
