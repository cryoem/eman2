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
#include "emassert.h"

using namespace EMAN;

EMAN1Ctf::EMAN1Ctf()
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


EMAN1Ctf::~EMAN1Ctf()
{
}


int EMAN1Ctf::from_string(const string & ctf)
{
	Assert(ctf != "");
	char type;
	int i = sscanf(ctf.c_str(), "%c%f %f %f %f %f %f %f %f %f %f %f",
				   &type,&defocus, &bfactor, &amplitude, &ampcont, &noise1,
				   &noise2, &noise3, &noise4, &voltage, &cs, &apix);
	if (i != 11) {
		return 1;
	}
	return 0;
}

void EMAN1Ctf::from_dict(const Dict & dict)
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

Dict EMAN1Ctf::to_dict() const
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

void EMAN1Ctf::from_vector(const vector<float>& vctf)
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

vector<float> EMAN1Ctf::to_vector() const
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


string EMAN1Ctf::to_string() const
{
	char ctf[1024];
	sprintf(ctf, "O%1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g %1.3g",
			defocus, bfactor, amplitude, ampcont, noise1, noise2, noise3, noise4, voltage, cs,
			apix);

	return string(ctf);
}

void EMAN1Ctf::copy_from(const Ctf * new_ctf)
{
	if (new_ctf) {
		const EMAN1Ctf *c = static_cast<const EMAN1Ctf *>(new_ctf);
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


vector < float >EMAN1Ctf::compute_1d(int size, float ds, CtfType type, XYData * sf)
{
	Assert(size > 0);

	float tmp_f1 = CTFOS * sqrt((float) 2) * size / 2;
	int np = (int) ceil(tmp_f1) + 2;
	vector < float >r;

	r.resize(np);

//	float ds = 1 / (apix * size * CTFOS);
	float s = 0;
	float g1 = calc_g1();
	float g2 = calc_g2();
	float amp1 = calc_amp1();

	switch (type) {
	case CTF_AMP:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf1(amp1, gamma, s);
			s += ds;
		}
		break;

	case CTF_SIGN:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf1(amp1, gamma, s)>0?1.0f:-1.0f;
			s += ds;
		}
		break;

	case CTF_BACKGROUND:
		for (int i = 0; i < np; i++) {
			r[i] = calc_noise(s);
			s += ds;
		}
		break;

	case CTF_SNR:
	case CTF_SNR_SMOOTH:
// 		if (!sf) {
// 			LOGERR("CTF computation error, no SF found\n");
// 			return r;
// 		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_snr(amp1, gamma, s);
			if (s && sf) {
				r[i] *= sf->get_yatx(s);
			}
			s += ds;
		}

		break;

	case CTF_WIENER_FILTER:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_snr(amp1, gamma, s);
			if (s && sf) {
				r[i] *= sf->get_yatx(s);
			}

			r[i] = 1.0f / (1.0f + 1.0f / r[i]);
			s += ds;
		}
		break;

	case CTF_TOTAL:
		if (!sf) {
			LOGERR("CTF computation error, no SF found\n");
			return r;
		}

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			if (sf) {
				r[i] = calc_ctf1(amp1, gamma, s);
				r[i] = r[i] * r[i] * sf->get_yatx(s) + calc_noise(s);
			}
			else {
				r[i] = calc_ctf1(amp1, gamma, s);
				r[i] = r[i] * r[i] + calc_noise(s);
			}
			s += ds;
		}
		break;
	default:
		break;
	}

	return r;
}


void EMAN1Ctf::compute_2d_real(EMData *, CtfType, XYData *)
{


}



void EMAN1Ctf::compute_2d_complex(EMData * image, CtfType type, XYData * sf)
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

	if (type == CTF_BACKGROUND) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {
#ifdef	_WIN32
				float s = (float) _hypot(x, y - ny / 2.0f) * ds;
#else
				float s = (float) hypot(x, y - ny / 2.0f) * ds;
#endif
				d[x * 2 + ynx] = calc_noise(s);
				d[x * 2 + ynx + 1] = 0;			// The phase is somewhat arbitrary
			}
		}
	}
	else if (type == CTF_AMP) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {
#ifdef	_WIN32
				float s = (float)_hypot((float) x, (float) y - ny / 2) * ds;
#else
				float s = (float)hypot((float) x, (float) y - ny / 2) * ds;
#endif	//_WIN32
				float gamma = calc_gamma(g1, g2, s);
				float v = fabs(calc_amplitude(gamma));
				d[x * 2 + ynx] = v;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}
	else if (type == CTF_SIGN) {
		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;
			for (int x = 0; x < nx / 2; x++) {
#ifdef	_WIN32
				float s = (float)_hypot(x, y - ny / 2.0f) * ds;
#else
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
#endif
				float gamma = calc_gamma(g1, g2, s);
				float v = calc_amplitude(gamma);
				d[x * 2 + ynx] = v > 0 ? 1.0f : -1.0f;
				d[x * 2 + ynx + 1] = 0;
			}
		}

	}
	else if (type == CTF_SNR || type == CTF_SNR_SMOOTH) {
		float amp1 = calc_amp1();

		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {

#ifdef	_WIN32
				float s = (float)_hypot(x, y - ny / 2.0f) * ds;
#else
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
#endif
				float gamma = calc_gamma(g1, g2, s);
				float f = calc_ctf1(amp1, gamma, s);
				float noise = calc_noise(s);
				f = f * f / noise;

				if (s && sf) {
					f *= sf->get_yatx(s);
				}
				d[x * 2 + ynx] *= f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}
	else if (type == CTF_WIENER_FILTER) {
		float amp1 = calc_amp1();

		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {

#ifdef	_WIN32
				float s = (float)_hypot(x, y - ny / 2.0f) * ds;
#else
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
#endif
				float gamma = calc_gamma(g1, g2, s);
				float f = calc_ctf1(amp1, gamma, s);
				float noise = calc_noise(s);
				f = f * f / noise;

				if (s) {
					f *= sf->get_yatx(s);
				}
				f = 1.0f / (1.0f + 1.0f / f);
				d[x * 2 + ynx] *= f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}
	else if (type == CTF_TOTAL) {
		float amp1 = calc_amp1();

		for (int y = 0; y < ny; y++) {
			int ynx = y * nx;

			for (int x = 0; x < nx / 2; x++) {

#ifdef	_WIN32
				float s = (float)_hypot(x, y - ny / 2.0f) * ds;
#else
				float s = (float)hypot(x, y - ny / 2.0f) * ds;
#endif
				float gamma = calc_gamma(g1, g2, s);
				float f = calc_ctf1(amp1, gamma, s);
				float noise = calc_noise(s);
				f = f * f;

				if (sf && s) {
					f *= sf->get_yatx(s);
				}
				f+=noise;

				d[x * 2 + ynx] *= f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}

	image->update();
}



bool EMAN1Ctf::equal(const Ctf * ctf1) const
{
	if (ctf1) {
		const EMAN1Ctf *c = static_cast<const EMAN1Ctf *>(ctf1);
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

/*************************************
EMAN2Ctf
*************************************/

EMAN2Ctf::EMAN2Ctf()
{
	defocus = 0;
	dfdiff = 0;
	dfang = 0;
	bfactor = 0;
	ampcont = 0;
	voltage = 0;
	cs = 0;
	apix = 1.0;
	dsbg=-1;
	background.clear();
	snr.clear();
}


EMAN2Ctf::~EMAN2Ctf()
{
}


int EMAN2Ctf::from_string(const string & ctf)
{
	Assert(ctf != "");
	char type=' ';
	int pos,i,j;
	int bglen=0,snrlen=0;
	float v;
	const char *s=ctf.c_str();

	sscanf(s, "%c%f %f %f %f %f %f %f %f %f %d%n",
				   &type,&defocus, &dfdiff,&dfang,&bfactor,&ampcont,&voltage, &cs, &apix,&dsbg,&bglen,&pos);
	if (type!='E') throw InvalidValueException(type,"Trying to initialize Ctf object with bad string");


	background.resize(bglen);
	for (i=0; i<bglen; i++) {
		if (sscanf(s+pos,",%f%n",&v,&j)<1) return(1);
		background[i]=v;
		pos+=j;
	}

	sscanf(s+pos," %d%n",&snrlen,&j);
	pos+=j;
	snr.resize(snrlen);
	for (i=0; i<snrlen; i++) {
		if (sscanf(s+pos,",%f%n",&v,&j)<1) return(1);
		snr[i]=v;
		pos+=j;
	}

	return 0;

}

string EMAN2Ctf::to_string() const
{
	char ctf[256];
	sprintf(ctf, "E%1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %d",
			defocus, dfdiff, dfang, bfactor, ampcont, voltage, cs, apix, dsbg,(int)background.size());

	string ret=ctf;
	for (int i=0; i<(int)background.size(); i++) {
		sprintf(ctf,",%1.3g",background[i]);
		ret+=ctf;
	}

	sprintf(ctf, " %d",(int)snr.size());
	ret+=ctf;
	for (int i=0; i<(int)snr.size(); i++) {
		sprintf(ctf,",%1.3g",snr[i]);
		ret+=ctf;
	}


	return ret;
}

void EMAN2Ctf::from_dict(const Dict & dict)
{
	defocus = (float)dict["defocus"];
	dfdiff = (float)dict["dfdiff"];
	dfang = (float)dict["dfang"];
	bfactor = (float)dict["bfactor"];
	ampcont = (float)dict["ampcont"];
	voltage = (float)dict["voltage"];
	cs = (float)dict["cs"];
	apix = (float)dict["apix"];
	dsbg = (float)dict["dsbg"];
	background = dict["background"];
	snr = dict["snr"];
}

Dict EMAN2Ctf::to_dict() const
{
	Dict dict;
	dict["defocus"] = defocus;
	dict["dfdiff"] = dfdiff;
	dict["dfang"] = dfang;
	dict["bfactor"] = bfactor;
	dict["ampcont"] = ampcont;
	dict["voltage"] = voltage;
	dict["cs"] = cs;
	dict["apix"] = apix;
	dict["dsbg"] = dsbg;
	dict["background"] = background;
	dict["snr"] = snr;

	return dict;
}

void EMAN2Ctf::from_vector(const vector<float>& vctf)
{
	int i;
	defocus = vctf[0];
	dfdiff = vctf[1];
	dfang = vctf[2];
	bfactor = vctf[3];
	ampcont = vctf[4];
	voltage = vctf[5];
	cs = vctf[6];
	apix = vctf[7];
	dsbg = vctf[8];
	background.resize((int)vctf[9]);
	for (i=0; i<(int)vctf[9]; i++) background[i]=vctf[i+10];
	snr.resize((int)vctf[i+10]);
	for (int j=0; j<(int)vctf[i+10]; j++) snr[j]=vctf[i+j+11];
}

vector<float> EMAN2Ctf::to_vector() const
{
	vector<float> vctf;

	vctf.push_back(defocus);
	vctf.push_back(dfdiff);
	vctf.push_back(dfang);
	vctf.push_back(bfactor);
	vctf.push_back(ampcont);
	vctf.push_back(voltage);
	vctf.push_back(cs);
	vctf.push_back(apix);
	vctf.push_back(dsbg);
	vctf.push_back((float)background.size());
	for (unsigned int i=0; i<background.size(); i++) vctf.push_back(background[i]);
	vctf.push_back((float)snr.size());
	for (unsigned int j=0; j<snr.size(); j++) vctf.push_back(snr[j]);

	return vctf;
}



void EMAN2Ctf::copy_from(const Ctf * new_ctf)
{
	if (new_ctf) {
		const EMAN2Ctf *c = static_cast<const EMAN2Ctf *>(new_ctf);
		defocus = c->defocus;
		dfdiff = c->dfdiff;
		dfang = c->dfang;
		bfactor = c->bfactor;
		ampcont = c->ampcont;
		voltage = c->voltage;
		cs = c->cs;
		apix = c->apix;
		dsbg = c->dsbg;
		background = c->background;
		snr = c->snr;
	}
}

inline int max_int(int a,int b) { return a>b?a:b; }
inline int min_int(int a,int b) { return a<b?a:b; }

vector < float >EMAN2Ctf::compute_1d(int size,float ds, CtfType type, XYData * sf)
{
	Assert(size > 0);

//	float tmp_f1 =  sqrt((float) 2) * size / 2;
//	int np = (int) ceil(tmp_f1) + 2;
	int np=size/2;
	vector < float >r;

	r.resize(np);

//	float ds = 1 / (apix * size);
	float s = 0;
	float g1 = calc_g1();
	float g2 = calc_g2();
	float amp1 = calc_amp1();

	switch (type) {
	case CTF_AMP:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf1(amp1, gamma, s);
			s += ds;
		}
		break;

	case CTF_SIGN:
		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			r[i] = calc_ctf1(amp1, gamma, s)>=0?1.0f:-1.0f;
			s += ds;
		}
		break;

	case CTF_BACKGROUND:
		for (int i = 0; i < np; i++) {
			float f = s/dsbg;
			int j = (int)floor(f);
			f-=j;
			if (j>(int)background.size()-2) r[i]=background.back();
			else r[i]=background[j]*(1.0f-f)+background[j+1]*f;
			s+=ds;
		}
		break;

	case CTF_SNR:
		for (int i = 0; i < np; i++) {
			float f = s/dsbg;
			int j = (int)floor(f);
			f-=j;
			if (j>(int)snr.size()-2) r[i]=snr.back();
			else r[i]=snr[j]*(1.0f-f)+snr[j+1]*f;
//			printf("%d\t%f\n",j,snr[j]);
			s+=ds;
		}
		break;
	case CTF_SNR_SMOOTH:
		// This apparently complicated routine tries to make a nice smooth and accurate SNR curve. It does this
		// by fitting local regions of the SNR vs the theoretical SNR (theoretical CTF^2/measured background),
		// then taking the slope of the result times the theoretical SNR to produce a local SNR estimate

		{ // <- is to permit new temporary value allocation
			vector < float >tsnr;	// theoretical SNR
			tsnr.resize(np);
			vector < float >dsnr;	// data SNR
			dsnr.resize(np);
			
			float s0=s;
			
			for (int i = 0; i < np; i++) {
				float gamma = calc_gamma(g1, g2, s);
				tsnr[i] = calc_ctf1(amp1, gamma, s);		// ctf amp

				// background value
				float f = s/dsbg;
				int j = (int)floor(f);
				f-=j;
				float bg;
				if (j>(int)background.size()-2) bg=background.back();
				else bg=background[j]*(1.0f-f)+background[j+1]*f;
				if (bg <=0) bg=.001f;

				tsnr[i] = tsnr[i]*tsnr[i]/bg;		// This is now a SNR curve
				if (sf && s) {
					tsnr[i] *= sf->get_yatx(s);
				}

				
				// This is the SNR computed from the data without fitting
				if (j>(int)snr.size()-2) dsnr[i]=snr.back();
				else dsnr[i]=snr[j]*(1.0f-f)+snr[j+1]*f;
				
				s+=ds;
			}

			int npsm=np/25;			// 1/2 number of points to smooth over, 25 is arbitrary
			if (npsm<2) npsm=2;
			
			s=s0;
			for (int i = 1; i < np; i++) {
				// simple linear regression embedded here
				double sum = 0;
				double sum_x = 0;
				double sum_y = 0;
				double sum_xx = 0;
				double sum_xy = 0;

				for (int k=max_int(i-npsm,1); k<=min_int(i+npsm,np-1); k++) {
					double y = dsnr[k];
					double x = tsnr[k];

					sum_x += x;
					sum_y += y;
					sum_xx += x * x;
					sum_xy += x * y;
					sum++;
				}

				double div = sum * sum_xx - sum_x * sum_x;
// 				if (div == 0) {
// 					div = 0.0000001f;
// 				}

	//			*intercept = (float) ((sum_xx * sum_y - sum_x * sum_xy) / div);
	//			*slope = (float) ((sum * sum_xy - sum_x * sum_y) / div);

				if (div!=0.0) r[i]=(float) ((sum * sum_xy - sum_x * sum_y) / div)*tsnr[i];
				else r[i]=0.0;
				if (r[i]<0) r[i]=0;
				
				s+=ds;
			}
			r[0]=0;
		}
		break;

	case CTF_WIENER_FILTER:
// 		if (!sf) {
// 			LOGERR("CTF computation error, no SF found\n");
// 			return r;
// 		}

		for (int i = 0; i < np; i++) {
			float f = s/dsbg;
			int j = (int)floor(f);
			float bg;
			f-=j;
			if (j>(int)snr.size()-2) {
/*				r[i]=snr.back();
				bg=background.back();*/
				r[i]=0;
			}
			else {
				r[i]=snr[j]*(1.0f-f)+snr[j+1]*f;
				bg=background[j]*(1.0f-f)+background[j+1]*f;
			}
			if (r[i]<0) r[i]=0;
			r[i]=r[i]/(r[i]+1.0f);		// Full Wiener filter assuming perfect image with noise
//			r[i]=sqrt(r[i]/bg)/(r[i]+1.0);	// Wiener filter with 1/CTF term (sort of) to restore image then filter
			s+=ds;
		}
		r[0]=0;
		break;

	case CTF_TOTAL:

		for (int i = 0; i < np; i++) {
			float gamma = calc_gamma(g1, g2, s);
			if (sf) {
				r[i] = calc_ctf1(amp1, gamma, s);
				r[i] = r[i] * r[i] * sf->get_yatx(s) + calc_noise(s);
			}
			else {
				r[i] = calc_ctf1(amp1, gamma, s);
				r[i] = r[i] * r[i] + calc_noise(s);
			}
			s += ds;
		}
		break;
	default:
		break;
	}

	return r;
}


void EMAN2Ctf::compute_2d_real(EMData *, CtfType, XYData *)
{


}



void EMAN2Ctf::compute_2d_complex(EMData * image, CtfType type, XYData * sf)
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

	if ((ny%2==1 && nx!=ny+1) || (ny%2==0 && nx != ny + 2)) {
		printf("nx=%d  ny=%d\n",nx,ny);
		LOGERR("compute_2d_complex only works on (nx, nx-2) images");
		return;
	}

	float ds = 1.0f / (apix * ny);
	image->to_one();

	float *d = image->get_data();
	float g1 = calc_g1();
	float g2 = calc_g2();
	float amp1 = calc_amp1();

	if (type == CTF_BACKGROUND) {
		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;

			for (int x = 0; x < nx / 2; x++) {
				float s = (float) Util::hypot_fast(x, y ) * ds;
				d[x * 2 + ynx] = calc_noise(s);
				d[x * 2 + ynx + 1] = 0;			// The phase is somewhat arbitrary
			}
		}
	}
	else if (type == CTF_AMP) {
		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;

			for (int x = 0; x < nx / 2; x++) {
				float s = (float)Util::hypot_fast(x,y ) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float v = calc_ctf1(amp1, gamma, s);
//				float v = calc_amplitude(gamma);
				d[x * 2 + ynx] = v;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}
	else if (type == CTF_SIGN) {
		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;
			for (int x = 0; x < nx / 2; x++) {
				float s = (float)Util::hypot_fast(x,y ) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float v = calc_amplitude(gamma);
				d[x * 2 + ynx] = v >= 0 ? 1.0f : -1.0f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}
	else if (type == CTF_SNR) {

		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;

			for (int x = 0; x < nx / 2; x++) {

				float s = (float)Util::hypot_fast(x,y ) * ds;
				float f = s/dsbg;
				int j = (int)floor(f);
				f-=j;
				if (j>(int)snr.size()-2) d[x*2+ynx]=snr.back();
				else d[x*2+ynx]=snr[j]*(1.0f-f)+snr[j+1]*f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
		d[0]=0;
	}
	else if (type == CTF_SNR_SMOOTH) {
		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;

			for (int x = 0; x < nx / 2; x++) {

				float s = (float)Util::hypot_fast(x,y ) * ds;
				float f = s/dsbg;
				int j = (int)floor(f);
				f-=j;
				if (j>(int)snr.size()-2) d[x*2+ynx]=snr.back();
				else d[x*2+ynx]=snr[j]*(1.0f-f)+snr[j+1]*f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
		d[0]=0;
	}
	else if (type == CTF_WIENER_FILTER) {
		if (dsbg==0) printf("Warning, DSBG set to 0\n");
		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;

			for (int x = 0; x < nx / 2; x++) {

				float s = (float)Util::hypot_fast(x,y ) * ds;
				float f = s/dsbg;
				int j = (int)floor(f);
				float bg,snrf;
				f-=j;
				if (j>(int)snr.size()-2) {
/*					bg=background.back();
					d[x*2+ynx]=snr.back()/(snr.back()+1.0);*/
					d[x*2+ynx]=0;
				}
				else {
					bg=background[j]*(1.0f-f)+background[j+1]*f;
					snrf=snr[j]*(1.0f-f)+snr[j+1]*f;

//					printf("%d\t%f\n",x,sqrt(snrf/bg)/(snrf+1.0));
					if (snrf<0) snrf=0.0;
//					d[x*2+ynx]=sqrt(snrf/bg)/(snrf+1.0);	// Note that this is a Wiener filter with a 1/CTF term to compensate for the filtration already applied, but needs to be multiplied by the structure factor
					d[x*2+ynx]=snrf/(snrf+1);	// This is just the simple Wiener filter

				}
				d[x * 2 + ynx + 1] = 0;
			}
		}
		d[0]=0;
	}
	else if (type == CTF_TOTAL) {
		float amp1 = calc_amp1();

		for (int y = -ny/2; y < ny/2; y++) {
			int y2=(y+ny)%ny;
			int ynx = y2 * nx;

			for (int x = 0; x < nx / 2; x++) {

				float s = (float)Util::hypot_fast(x,y ) * ds;
				float gamma = calc_gamma(g1, g2, s);
				float f = calc_ctf1(amp1, gamma, s);
				float noise = 0;
				f = f * f;

				if (sf && s) {
					f *= sf->get_yatx(s);
				}
				f+=noise;

				d[x * 2 + ynx] *= f;
				d[x * 2 + ynx + 1] = 0;
			}
		}
	}
	else printf("Unknown CTF image mode\n");

	image->update();
}



bool EMAN2Ctf::equal(const Ctf * ctf1) const
{
	if (ctf1) {
		const EMAN2Ctf *c = static_cast<const EMAN2Ctf *>(ctf1);
		if (defocus != c->defocus ||
			dfdiff != c->dfdiff ||
			dfang != c->dfang ||
			bfactor != c->bfactor ||
			ampcont != c->ampcont ||
			voltage != c->voltage ||
			cs != c->cs ||
			apix != c->apix ||
			dsbg != c->dsbg ||
			background.size() != c->background.size() ||
			snr.size() != c->snr.size()
			) return false;

		for (unsigned int i=0; i<background.size(); i++) {
			if (background[i]!=c->background[i]) return false;
		}
		for (unsigned int i=0; i<snr.size(); i++) {
			if (snr[i]!=c->snr[i]) return false;
		}
		return true;
	}
	return false;
}
