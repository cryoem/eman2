/**
 * $Id$
 */
#include "ctf.h"
#include "log.h"
#include "emdata.h"
#include "xydata.h"


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


int SimpleCtf::from_string(string ctf)
{
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
	SimpleCtf *c = (SimpleCtf*)(new_ctf);
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


vector<float> SimpleCtf::compute_1d(EMData * image, CtfType type, XYData * sf)
{
    if (!image) {
	Log::logger()->error("image is null. cannot computer 1D CTF");
	return vector<float>();
    }
    
    int ny = image->get_ysize();
    float tmp_f1 = CTFOS * sqrt((float)2) * ny / 2;
    int np = (int) ceil(tmp_f1) + 2;
    vector<float> r;

    if (image->has_ctff()) {
	if (type == CTF_AMP_S ||
	    type == CTF_NOISE_S ||
	    type == CTF_SNR_WIENER ||
	    type == CTF_WIENER_CTF_CORRECTION1 ||
	    type == CTF_WIENER_CTF_CORRECTION2 || type == CTF_TOTAL_CURVE) {
	    Log::logger()->error("%d: Cannot generate this CTF info from a file.", type);
	    return r;
	}

	char imagefile[80];
	float pix = 0;
	int cn = 0;

	if (sscanf(image->get_name().c_str(), "!$ %d,%f,%s", &cn, &pix, imagefile) != 3) {
	    Log::logger()->error("Reported HASCTFF, but filename not present in name.");
	    return r;
	}

	EMData tmp_img;

	if (tmp_img.read_image(imagefile, cn)) {
	    Log::logger()->error("Cannot read HASCTFF file: %s,%d: ", imagefile, cn);
	    return r;
	}

	float scale = tmp_img.get_xsize() / (CTFOS * ny / 2.0);
	float *sr = tmp_img.get_data();
	r.resize(np);

	if (type == CTF_RELATIVE_SNR) {
	    sr += tmp_img.get_xsize();
	}
	else if (type == CTF_ABS_SNR) {
	    sr += 2 * tmp_img.get_xsize();
	}

	for (int i = 0; i < CTFOS * ny / 2; i++) {
	    float iscale = i * scale;
	    int fl = (int) floor(iscale);
	    r[i] = sr[fl + 1] * (iscale - fl) + sr[fl] * (1 + fl - iscale);
	}

	tmp_img.done_data();
    }
    else {
	if (image->get_ctf() == 0) {
	    Log::logger()->error("ctfCurve with no CTF parameters\n");
	    return r;
	}
	
	r.resize(np);
	
	float ds = 1 / (apix * ny * CTFOS);
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
		Log::logger()->error("CTF computation error, no SF found\n");
		return r;
	    }

	    for (int i = 0; i < np; i++) {
		float gamma = calc_gamma(g1, g2, s);
		r[i] = calc_ctf(amp1, gamma, s);
		if (s) {
		    r[i] *= pow((float) 10.0, sf->get_yatx(s));
		}
		s += ds;
	    }

	    break;

	case CTF_SNR_WIENER:
	    if (!sf) {
		Log::logger()->error("CTF computation error, no SF found\n");
		return r;
	    }

	    for (int i = 0; i < np; i++) {
		float gamma = calc_gamma(g1, g2, s);
		r[i] = calc_ctf(amp1, gamma, s);
		if (s) {
		    r[i] *= pow((float) 10.0, sf->get_yatx(s));
		}

		r[i] = 1.0 / (1.0 + 1.0 / r[i]);
		s += ds;
	    }
	    break;

	case CTF_WIENER_CTF_CORRECTION1:
	    if (!sf) {
		Log::logger()->error("CTF computation error, no SF found\n");
		return r;
	    }

	    for (int i = 0; i < np; i++) {
		float gamma = calc_gamma(g1, g2, s);
		r[i] = calc_ctf(amp1, gamma, s);
		if (s) {
		    r[i] *= pow((float) 10.0, sf->get_yatx(s));
		}


		float v = fabs(calc_ctf1(amp1, gamma, s));
		if (r[i] == 0 || v == 0) {
		    r[i] = 0;
		}
		else {
		    r[i] = (1.0 / (1.0 + 1.0 / r[i])) / v;
		}
		s += ds;
	    }
	    break;

	case CTF_WIENER_CTF_CORRECTION2:
	    if (!sf) {
		Log::logger()->error("CTF computation error, no SF found\n");
		return r;
	    }

	    for (int i = 0; i < np; i++) {
		float gamma = calc_gamma(g1, g2, s);
		r[i] = calc_ctf(amp1, gamma, s);
		if (s) {
		    r[i] *= pow((float) 10.0, sf->get_yatx(s));
		}
		float v = calc_ctf1(amp1, gamma, s);
		if (r[i] == 0 || v == 0) {
		    r[i] = 0;
		}
		else {
		    r[i] = (1.0 / (1.0 + 1.0 / r[i])) / v;
		}
		s += ds;
	    }
	    break;

	case CTF_TOTAL_CURVE:
	    if (!sf) {
		Log::logger()->error("CTF computation error, no SF found\n");
		return r;
	    }

	    for (int i = 0; i < np; i++) {
		float gamma = calc_gamma(g1, g2, s);
		if (sf->is_validx(s)) {
		    r[i] = calc_ctf1(amp1, gamma, s);
		    r[i] = r[i] * r[i] * pow((float)10.0, sf->get_yatx(s)) + calc_noise(s);
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
    }

    return r;
}


void SimpleCtf::compute_2d_real(EMData * image, CtfType type, XYData * sf)
{


}



void SimpleCtf::compute_2d_complex(EMData * image, CtfType type, XYData * sf)
{
    if (!image) {
	Log::logger()->error("image is null. cannot computer 2D complex CTF");
	return;
    }
    
    if (image->get_ctf() == 0) {
	Log::logger()->error("image has no CTF parameters");
	return;
    }

    if (image->is_complex() == false) {
	Log::logger()->error("compute_2d_complex can only work on complex images");
	return;
    }

    int nx = image->get_xsize();
    int ny = image->get_ysize();
    
    if (nx!=ny+2) {
	Log::logger()->error("compute_2d_complex only works on (nx, nx-2) images");
	return;
    }
    
    float ds=1.0/(apix*ny);
    image->to_one();

    float *d=image->get_data();
    float g1 = calc_g1();
    float g2 = calc_g2();
	
    if (type==CTF_NOISE) {
	for (int y=0; y<ny; y++) {
	    int ynx = y * nx;
	    
	    for (int x=0; x<nx/2; x++) {
		float s=hypot(x,y-ny/2.0)*ds;
		d[x*2+ynx]=calc_noise(s);
		d[x*2+ynx+1]=d[x*2+ynx];
	    }
	}
    }
    else if (type==CTF_BFACTOR) {
	for (int y=0; y<ny; y++) {
	    int ynx = y * nx;
	    
	    for (int x=0; x<nx/2; x++) {
		float s=hypot(x,y-ny/2.0)*ds;
		float gamma=calc_gamma(g1, g2, s);
		float v=calc_amplitude(gamma);
		d[x*2+ynx]*=v;
		d[x*2+ynx+1]=d[x*2+ynx];
	    }
	}
    }
    else if (type == CTF_AMP) {
	for (int y=0; y<ny; y++) {
	    int ynx = y * nx;
		
	    for (int x=0; x<nx/2; x++) {
		float s=hypot((float)x,(float)y-ny/2)*ds;
		float gamma=calc_gamma(g1, g2, s);
		float v=fabs(calc_amplitude(gamma));
		d[x*2+ynx]*=v;
		d[x*2+ynx+1]=d[x*2+ynx];
	    }
	}
    }
    else if (type == CTF_SIGN) {	    
	for (int y=0; y<ny; y++) {
	    int ynx = y * nx;
	    for (int x=0; x<nx/2; x++) {
		float s=hypot(x,y-ny/2.0)*ds;
		float gamma=calc_gamma(g1, g2, s);
		float v = calc_amplitude(gamma);
		v = v>0?1.0:-1.0;
		d[x*2+ynx]*=v;
		d[x*2+ynx+1]=d[x*2+ynx];
	    }
	}
	    
    }    
    else if (type==CTF_BFACTOR) {
	for (int y=0; y<ny; y++) {
	    int ynx = y * nx;
	    
	    for (int x=0; x<nx/2; x++) {
		float s=hypot(x,y-ny/2.0)*ds;
		float v=exp(-(bfactor*s*s));
		d[x*2+ynx]*=v;
		d[x*2+ynx+1]=d[x*2+ynx];
	    }
	}
    }	
    else if (type==CTF_WIENER_FILTER) {
	float amp1 = calc_amp1();
	
	for (int y=0; y<ny; y++) {
	    int ynx = y * nx;
	    
	    for (int x=0; x<nx/2; x++) {
		
		float s=hypot(x,y-ny/2.0)*ds;
		float gamma=calc_gamma(g1, g2, s);
		float f=calc_ctf1(amp1, gamma, s);
		float noise = calc_noise(s);
		f=f*f/noise;

		if (s) {
		    f*=pow((float)10.0, sf->get_yatx(s));
		}
		f=1.0/(1.0+1.0/f);
		d[x*2+ynx]*=f;
		d[x*2+ynx+1]=d[x*2+ynx];
	    }
	}
    }

    image->done_data();
}



bool SimpleCtf::equal(const Ctf *ctf1) const
{
    if (ctf1) {
	SimpleCtf *c = (SimpleCtf*)(ctf1);
	
	if (defocus == c->defocus &&
	    bfactor == c->bfactor &&
	    amplitude == c->amplitude &&
	    ampcont == c->ampcont &&
	    noise1 == c->noise1 &&
	    noise2 == c->noise2 &&
	    noise3 == c->noise3 &&
	    noise4 == c->noise4 &&
	    voltage == c->voltage &&
	    cs == c->cs &&
	    apix == c->apix) {
	    return true;
	}
    }
    return false;
}

   
