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
 
#include "ctf.h"
#include "emdata.h"
#include "xydata.h"
#include "Assert.h"

using namespace EMAN;

SimpleCtf::SimpleCtf()
{
	defocus = 0;
	bfactor = 0;
	amplitude = 0;
	ampcont = 0;
	noise1 = 0;
	noise2 = 0;
	noise3 = 0;
	noise4 = 0;
	voltage = 0;
	cs = 0;
	apix = 0;
}


SimpleCtf::~SimpleCtf()
{
}


int SimpleCtf::from_string(const string & ctf)
{
	Assert(ctf != "");
	int i = sscanf(ctf.c_str(), "%f %f %f %f %f %f %f %f %f %f %f",
				   &defocus, &bfactor, &amplitude, &ampcont, &noise1,
				   &noise2, &noise3, &noise4, &voltage, &cs, &apix);
	if (i != 11) {
		return 1;
	}
	return 0;
}

void SimpleCtf::from_dict(const Dict & dict)
{
	defocus = dict["defocus"];
	bfactor = dict["bfactor"];
	amplitude = dict["amplitude"];
	ampcont = dict["ampcont"];
	noise1 = dict["noise1"];
	noise2 = dict["noise2"];
	noise3 = dict["noise3"];
	noise4 = dict["noise4"];
	voltage = dict["voltage"];
	cs = dict["cs"];
	apix = dict["apix"];
}

Dict SimpleCtf::to_dict() const
{
	Dict dict;
	dict["defocus"] = defocus;
	dict["bfactor"] = bfactor;
	dict["amplitude"] = amplitude;
	dict["ampcont"] = ampcont;
	dict["noise1"] = noise1;
	dict["noise2"] = noise2;
	dict["noise3"] = noise3;
	dict["noise4"] = noise4;
	dict["voltage"] = voltage;
	dict["cs"] = cs;
	dict["apix"] = apix;

	return dict;
}

void SimpleCtf::from_vector(const vector<float>& vctf)
{
	defocus = vctf[0];
	bfactor = vctf[1];
	amplitude = vctf[2];
	ampcont = vctf[3];
	noise1 = vctf[4];
	noise2 = vctf[5];
	noise3 = vctf[6];
	noise4 = vctf[7];
	voltage = vctf[8];
	cs = vctf[9];
	apix = vctf[10];
}

vector<float> SimpleCtf::to_vector() const
{
	vector<float> vctf;
	
	vctf.push_back(defocus);
	vctf.push_back(bfactor);
	vctf.push_back(amplitude);
	vctf.push_back(ampcont);
	vctf.push_back(noise1);
	vctf.push_back(noise2);
	vctf.push_back(noise3);
	vctf.push_back(noise4);
	vctf.push_back(voltage);
	vctf.push_back(cs);
	vctf.push_back(apix);
	
	return vctf;
}
		

string SimpleCtf::to_string() const
{
	char ctf[1024];
	sprintf(ctf, "%1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g",
			defocus, bfactor, amplitude, ampcont, noise1, noise2, noise3, noise4, voltage, cs,
			apix);

	return string(ctf);
}

void SimpleCtf::copy_from(const Ctf * new_ctf)
{
	if (new_ctf) {
		SimpleCtf *c = (SimpleCtf *) (new_ctf);
		defocus = c->defocus;
		bfactor = c->bfactor;
		amplitude = c->amplitude;
		ampcont = c->ampcont;
		noise1 = c->noise1;
		noise2 = c->noise2;
		noise3 = c->noise3;
		noise4 = c->noise4;
		voltage = c->voltage;
		cs = c->cs;
		apix = c->apix;
	}
}


vector < float >SimpleCtf::compute_1d(int size, CtfType type, XYData * sf)
{
	Assert(size > 0);
	
	float tmp_f1 = CTFOS * sqrt((float) 2) * size / 2;
	int np = (int) ceil(tmp_f1) + 2;
	vector < float >r;

	r.resize(np);

	float ds = 1 / (apix * size * CTFOS);
	float s = 0;
	float g1 = calc_g1();
	float g2 = calc_g2();
	float amp1 = calc_amp1();

	switch (type) {
	case CTF_AMP_S:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf1(amp1, gamma, s);
			s += ds;
		}
		break;

	case CTF_NOISE_S:
		for (int i = 0; i < np; i++) {
			r[i] = calc_noise(s);
			s += ds;
		}
		break;

	case CTF_ABS_AMP_S:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = fabs(calc_ctf1(amp1, gamma, s));
			s += ds;
		}
		break;

	case CTF_RELATIVE_SNR:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf(amp1, gamma, s);
			s += ds;
		}
		break;

	case CTF_ABS_SNR:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf(amp1, gamma, s);
			if (s) {
				r[i] *= pow(10.0f, sf->get_yatx(s));
			}
			s += ds;
		}

		break;

	case CTF_SNR_WIENER:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf(amp1, gamma, s);
			if (s) {
				r[i] *= pow(10.0f, sf->get_yatx(s));
			}

			r[i] = 1.0f / (1.0f + 1.0f / r[i]);
			s += ds;
		}
		break;

	case CTF_WIENER_CTF_CORRECTION1:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf(amp1, gamma, s);
			if (s) {
				r[i] *= pow(10.0f, sf->get_yatx(s));
			}


			float v = fabs(calc_ctf1(amp1, gamma, s));
			if (r[i] == 0 || v == 0) {
				r[i] = 0;
			}
			else {
				r[i] = (1.0f / (1.0f + 1.0f / r[i])) / v;
			}
			s += ds;
		}
		break;

	case CTF_WIENER_CTF_CORRECTION2:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf(amp1, gamma, s);
			if (s) {
				r[i] *= pow(10.0f, sf->get_yatx(s));
			}
			float v = calc_ctf1(amp1, gamma, s);
			if (r[i] == 0 || v == 0) {
				r[i] = 0;
			}
			else {
				r[i] = (1.0f / (1.0f + 1.0f / r[i])) / v;
			}
			s += ds;
		}
		break;

	case CTF_TOTAL_CURVE:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			if (sf->is_validx(s)) {
				r[i] = calc_ctf1(amp1, gamma, s);
				r[i] = r[i] * r[i] * pow(10.0f, sf->get_yatx(s)) + calc_noise(s);
			}
			else {
				r[i] = 0;
			}
			s += ds;
		}
		break;
	default:
		break;
	}

	return r;
}


void SimpleCtf::compute_2d_real(EMData *, CtfType, XYData *)
{


}



void SimpleCtf::compute_2d_complex(EMData * image, CtfType type, XYData * sf)
{
	if (!image) {
		LOGERR("image is null. cannot computer 2D complex CTF");
		return;
	}

	if (image->is_complex() == false) {
		LOGERR("compute_2d_complex can only work on complex images");
		return;
	}

	int nx = image->get_xsize();
	int ny = image->get_ysize();

	if (nx != ny + 2) {
		LOGERR("compute_2d_complex only works on (nx, nx-2) images");
		return;
	}

	float ds = 1.0f / (apix * ny);
	image->to_one();

	float *d = image->get_data();
	float g1 = calc_g1();
	float g2 = calc_g2();

	if (type == CTF_NOISE) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {
				float s = (float) hypot(x, y - ny / 2.0f) * ds;
				d[x * 2 + ynx] = calc_noise(s);
				d[x * 2 + ynx + 1] = d[x * 2 + ynx];
			}
		}
	}
	else if (type == CTF_BFACTOR) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float v = calc_amplitude(gamma);
				d[x * 2 + ynx] *= v;
				d[x * 2 + ynx + 1] = d[x * 2 + ynx];
			}
		}
	}
	else if (type == CTF_AMP) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {
				float s = (float)hypot((float) x, (float) y - ny / 2) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float v = fabs(calc_amplitude(gamma));
				d[x * 2 + ynx] *= v;
				d[x * 2 + ynx + 1] = d[x * 2 + ynx];
			}
		}
	}
	else if (type == CTF_SIGN) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;
			for (int x = 0; x < nx / 2; x++) {
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float v = calc_amplitude(gamma);
				v = v > 0 ? 1.0f : -1.0f;
				d[x * 2 + ynx] *= v;
				d[x * 2 + ynx + 1] = d[x * 2 + ynx];
			}
		}

	}
	else if (type == CTF_BFACTOR) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
				float v = exp(-(bfactor * s * s));
				d[x * 2 + ynx] *= v;
				d[x * 2 + ynx + 1] = d[x * 2 + ynx];
			}
		}
	}
	else if (type == CTF_WIENER_FILTER) {
		float amp1 = calc_amp1();

		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {

				float s = (float)hypot(x, y - ny / 2.0f) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float f = calc_ctf1(amp1, gamma, s);
				float noise = calc_noise(s);
				f = f * f / noise;

				if (s) {
					f *= pow(10.0f, sf->get_yatx(s));
				}
				f = 1.0f / (1.0f + 1.0f / f);
				d[x * 2 + ynx] *= f;
				d[x * 2 + ynx + 1] = d[x * 2 + ynx];
			}
		}
	}

	image->done_data();
}



bool SimpleCtf::equal(const Ctf * ctf1) const
{
	if (ctf1) {
		SimpleCtf *c = (SimpleCtf *) (ctf1);

		if (defocus == c->defocus &&
			bfactor == c->bfactor &&
			amplitude == c->amplitude &&
			ampcont == c->ampcont &&
			noise1 == c->noise1 &&
			noise2 == c->noise2 &&
			noise3 == c->noise3 &&
			noise4 == c->noise4 && voltage == c->voltage && cs == c->cs && apix == c->apix) {
			return true;
		}
	}
	return false;
}
