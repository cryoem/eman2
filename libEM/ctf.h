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

#ifndef eman__ctf__h__
#define eman__ctf__h__ 1

#include <cmath>

#include "emobject.h"

#ifdef WIN32
	#ifndef M_PI
		#define M_PI 3.14159265358979323846f
	#endif	//M_PI
#endif	//WIN32

using std::string;
using std::map;

namespace EMAN
{
	class EMData;
	class XYData;

	/** Ctf is the base class for all CTF model.
     *
     * Contrast transfer function (CTF) is the  function that
     * describes the transfer of information from the object
     * to the contrast observed in the image for electron microscopy.
     */
	class Ctf
	{
	  public:
		// NOTE: ctf is positive for the first peak, instead of negative
		enum CtfType
		{
			CTF_AMP,			// ctf ampltidue only
			CTF_SIGN,			// ctf sign (+-1)
			CTF_BACKGROUND,		// Background, no ctf oscillation
			CTF_SNR,			// Signal to noise ratio
			CTF_SNR_SMOOTH,		// Signal to noise ratio, smoothed, algorithm may vary, but this should be more suitable for weighting
			CTF_WIENER_FILTER,	// Weiner Filter = 1/(1+1/snr)
			CTF_TOTAL			// AMP*AMP+NOISE
		};
	  public:
		virtual ~ Ctf()
		{
		}

		float defocus;			// 	Defocus in microns, note that postitive is now underfocus, whereas values in EMAN1 are negative overfocus
		float bfactor;			// 	B-factor expressed using the x-ray convention (e^-B/4 s^2 in amplitude space) EMAN1 used E^-B s^2
		float voltage;			//  in kV
		float cs;				//  in mm
		float apix;				//

		virtual int from_string(const string & ctf) = 0;	// The first letter of the string indicates the subclass type
		virtual string to_string() const = 0;

		virtual void from_dict(const Dict & dict) = 0;
		virtual Dict to_dict() const = 0;

		virtual void from_vector(const vector<float>& vctf) = 0;
		virtual vector<float> to_vector() const = 0;

		virtual vector < float >compute_1d(int size,float ds, CtfType t, XYData * struct_factor = 0) = 0;
		virtual void compute_2d_real(EMData * img, CtfType t, XYData * struct_factor = 0) = 0;
		virtual void compute_2d_complex(EMData * img, CtfType t, XYData * struct_factor = 0) = 0;

		virtual void copy_from(const Ctf * new_ctf) = 0;
		virtual bool equal(const Ctf * ctf1) const = 0;

	  public:
		enum
		{ CTFOS = 5 };

	};

	/** EMAN1Ctf is the CTF model used in EMAN1.
     */
	class EMAN1Ctf:public Ctf
	{
	  public:
//		float defocus;			// 0	Defocus in microns, note that postitive is now underfocus, whereas values in EMAN1 are negative overfocus
//		float bfactor;			// 1	B-factor expressed using the x-ray convention (e^-B/4 s^2 in amplitude space) EMAN1 used E^-B s^2
		float amplitude;		// 2
		float ampcont;			// 3
		float noise1;			// 4
		float noise2;			// 5
		float noise3;			// 6
		float noise4;			// 7
//		float voltage;			// 8
//		float cs;				// 9
//		float apix;				// 10

	  public:
		EMAN1Ctf();
		EMAN1Ctf(const vector<float>& vf) {from_vector(vf);}	//for unpickling
		~EMAN1Ctf();

		vector < float >compute_1d(int size,float ds, CtfType type, XYData * struct_factor = 0);
		void compute_2d_real(EMData * image, CtfType type, XYData * struct_factor = 0);
		void compute_2d_complex(EMData * image, CtfType type, XYData * struct_factor = 0);

		int from_string(const string & ctf);
		string to_string() const;

		void from_dict(const Dict & dict);
		Dict to_dict() const;

		void from_vector(const vector<float>& vctf);
		vector<float> to_vector() const;

		void copy_from(const Ctf * new_ctf);
		bool equal(const Ctf * ctf1) const;

		float get_defocus() const
		{
			return defocus;
		}
		float get_bfactor() const
		{
			return bfactor;
		}

	  private:
		inline float calc_amp1()
		{
			return (sqrt(1 - ampcont * ampcont/10000.0f));
		}

		inline float calc_lambda()
		{
			float lambda = 12.2639f / sqrt(voltage * 1000.0f + 0.97845f * voltage * voltage);
			return lambda;
		}

		inline float calc_g1()
		{
			float lambda = calc_lambda();
			float g1 = 2.5e6f * cs * lambda * lambda * lambda;
			return g1;
		}

		inline float calc_g2()
		{
			float lambda = calc_lambda();
			float g2 = 5000.0f * -defocus * lambda;
			return g2;
		}

		inline float calc_gamma(float g1, float g2, float s)
		{
			float s2 = s * s;
			float gamma = (float) (-2 * M_PI * (g1 * s2 * s2 + g2 * s2));
			return gamma;
		}

		inline float calc_ctf1(float g, float gamma, float s)
		{
			float r = amplitude * exp(-(bfactor/4.0f * s * s)) * (g * sin(gamma) + ampcont/100.0f * cos(gamma));
			return r;
		}

		inline float calc_amplitude(float gamma)
		{
			float t1 = sqrt(1.0f - ampcont * ampcont/10000.0f) * sin(gamma);
			float v = amplitude * (t1 + ampcont/100.0f * cos(gamma));
			return v;
		}

		inline float calc_noise(float s)
		{
			float ns = (float) M_PI / 2 * noise4 * s;
			float ns2 = ns * ns;
			float n = noise3 * exp(-ns2 - s * noise2 - sqrt(fabs(s)) * noise1);
			return n;
		}

		inline float calc_snr(float g1, float gamma, float s)
		{
			float ctf1 = calc_ctf1(g1, gamma, s);
			float ctf2 = ctf1 * ctf1 / calc_noise(s);
			return ctf2;
		}

	};

	/** EMAN2Ctf is the default CTF model used in EMAN2
     */
	class EMAN2Ctf:public Ctf
	{
	  public:
//		float defocus;		// defocus in microns, positive underfocus
		float dfdiff;		// defocus difference for astigmatism, defocus is the major elliptical axis
		float dfang;		// angle of the major elliptical axis in degrees measured counterclockwise from x
//		float bfactor;		// B-factor in 1/A^2 expressed using the x-ray convention (e^-B/4 s^2 in amplitude space) EMAN1 used E^-B s^2
		float ampcont;		// amplitude contrast as a percentage ie- this should be 10 for 10% amp contrast
//		float voltage;		// microscope voltage in kV
//		float cs;			// Cs in mm
//		float apix;			// A/pix value used when generating 2D results
		float dsbg;			// ds value for background and SNR
		vector<float> background;	// background intensity, 1 value per radial pixel (NX/2, corners ignored)
		vector<float> snr;			// SNR, 1 value per radial pixel (NX/2, corners assumed 0)

		vector<float> get_snr(){ return snr;}
		void set_snr(const vector<float>& vf) {snr = vf;}
		vector<float> get_background(){ return background;}
		void set_background(const vector<float>& vf) {background = vf;}

	  public:
		EMAN2Ctf();
		EMAN2Ctf(const vector<float>& vf) {from_vector(vf);}	//for unpickling
		~EMAN2Ctf();

		vector < float >compute_1d(int size,float ds, CtfType type, XYData * struct_factor = 0);
		void compute_2d_real(EMData * image, CtfType type, XYData * struct_factor = 0);
		void compute_2d_complex(EMData * image, CtfType type, XYData * struct_factor = 0);

		int from_string(const string & ctf);
		string to_string() const;

		void from_dict(const Dict & dict);
		Dict to_dict() const;

		void from_vector(const vector<float>& vctf);
		vector<float> to_vector() const;

		void copy_from(const Ctf * new_ctf);
		bool equal(const Ctf * ctf1) const;

	  private:
		inline float calc_amp1()
		{
			return (sqrt(1 - ampcont * ampcont/10000.0f));
		}

		inline float calc_lambda()
		{
			float lambda = 12.2639f / sqrt(voltage * 1000.0f + 0.97845f * voltage * voltage);
			return lambda;
		}

		inline float calc_g1()
		{
			float lambda = calc_lambda();
			float g1 = 2.5e6f * cs * lambda * lambda * lambda;
			return g1;
		}

		inline float calc_g2()
		{
			float lambda = calc_lambda();
			float g2 = 5000.0f * -defocus * lambda;
			return g2;
		}

		inline float calc_gamma(float g1, float g2, float s)
		{
			float s2 = s * s;
			float gamma = (float) (-2 * M_PI * (g1 * s2 * s2 + g2 * s2));
			return gamma;
		}

		inline float calc_ctf1(float g, float gamma, float s)
		{
			float r = exp(-(bfactor/4.0f * s * s)) * (g * sin(gamma) + ampcont/100.0f * cos(gamma));
			return r;
		}

		inline float calc_amplitude(float gamma)
		{
			float t1 = sqrt(1.0f - ampcont * ampcont/10000.0f) * sin(gamma);
			float v = (t1 + ampcont/100.0f * cos(gamma));
			return v;
		}

		inline float calc_noise(float s)
		{
			int si=(int)(s/dsbg);
			if (si>(int)background.size()||si<0) return background.back();
			return background[si];
		}

	};

}



#endif	//eman__ctf__h__
