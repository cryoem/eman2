#ifndef eman__ctf__h__
#define eman__ctf__h__ 1

#include <string>
#include <map>
#include <math.h>
#include "emobject.h"

using std::string;
using std::map;

namespace EMAN {

    class Ctf {
    public:
	virtual ~Ctf() {}
	virtual bool cmp() const = 0;
	//virtual CtfMapType get_maptype() = 0;
    
	virtual int from_string(string ctf) = 0;
	virtual string to_string() const = 0;

	virtual int from_dict(map<string, EMObject>& dict) = 0;
	virtual int to_dict(map<string, EMObject>& dict) const = 0;

    public:
	enum CtfMapType {
	    CTF_MAP_CTF,	// the true ctf with B decay and with positive and negative peaks, 
	    // NOTE: ctf is positive for the first peak, instead of negative
	    CTF_MAP_CTF_NO_B,	// the true ctf without B decay and with positive and negative peaks,
	    // NOTE: ctf is positive for the first peak, instead of negative
	    CTF_MAP_AMP,	// ctf ampltidue only = fabs(CTF_MAP_CTF)
	    CTF_MAP_AMP_NO_B,	// ctf ampltidue only = fabs(CTF_MAP_CTF_NO_B)
	    CTF_MAP_SIGN,	// ctf sign (+-1)	= sign(CTF_MAP_CTF)
	    CTF_MAP_B_FACTOR,	// B factor decay only, no ctf oscillation
	    CTF_MAP_BACKGROUND,	// Background, no ctf oscillation
	    CTF_MAP_SNR,	// Signal to noise ratio
	    CTF_MAP_SNR_SIGN,	// Signal to noise ratio with sign = CTF_MAP_SNR*CTF_MAP_SIGN
	    CTF_MAP_WIENER_FILTER,	// Weiner Filter = 1/(1+1/snr)
	    CTF_MAP_WIENER_CTF_CORRECTION	// ctf correction with Weiner Filter = 1/(ctf*exp(-b*s^2)*(1+1/snr))
	};
    
	enum CtfCurveType {
	    CTF_CURVE_AMP_S = 0,
	    CTF_CURVE_NOISE_S = 1,
	    CTF_CURVE_ABS_AMP_S = 2,
	    CTF_CURVE_RELATIVE_SNR = 3,
	    CTF_CURVE_ABS_SNR = 4,
	    CTF_CURVE_SNR_WIENER = 5,
	    CTF_CURVE_WIENER_CTF_CORRECTION1 = 6,
	    CTF_CURVE_WIENER_CTF_CORRECTION2 = 7,
	    CTF_CURVE_TOTAL_CURVE = 8
	};
    };
    

    class EMData;
    class XYData;


    class SimpleCtf : public Ctf {
    public:
	float defocus; // 0
	float bfactor; // 1
	float amplitude; // 2
	float ampcont; // 3
	float noise1; // 4
	float noise2; // 5
	float noise3; // 6
	float noise4; // 7
	float voltage; // 8
	float cs; // 9
	float apix; // 10

	Ctf::CtfMapType ctfmaptype;
	float astig_amp;
	float astig_ang;
	float drift_amp;
	float drift_ang;
    
    public:
	SimpleCtf();
	~SimpleCtf();

	bool cmp() const;
	Ctf::CtfMapType get_maptype();
	void compute_map(Ctf::CtfMapType maptype, EMData* power_spectrum = 0);
	void compute_curve(Ctf::CtfCurveType type, XYData* xy);
    
	bool is_changed() const;
    
	//void compute_bg(EMData* power_spectrum);
	bool is_set_properly();

	int from_string(string ctf);
	string to_string() const;
	
	int from_dict(map<string, EMObject>& dict);
	int to_dict(map<string, EMObject>& dict) const;

	static float* calc_ctf_curve(EMData* image, CtfCurveType type, XYData *sf);

    private:
	static inline float calc_gamma(float g1, float g2, float s)
	{
	    float s2 = s * s;
	    float gamma = -2 * M_PI * (g1 * s2 * s2 + g2 * s2);
	    return gamma;
	}

	static inline float calc_ctf1(SimpleCtf* ctf, float g, float gamma, float s)
	{
	    float r = ctf->amplitude * exp(-(ctf->bfactor*s*s)) *
		(g * sin(gamma) + ctf->ampcont * cos(gamma));
	    return r;
	}

	static inline float calc_noise(SimpleCtf* ctf, float n4, float s)
	{
	    float ns = n4 * s;
	    float ns2 = ns * ns;
	    float n = ctf->noise3 * exp(-ns2 - s * ctf->noise2 - sqrt(fabs(s)) * ctf->noise1);
	    return n;
	}
	
	static inline float calc_ctf(SimpleCtf* ctf, float g1, float n4, float gamma, float s)
	{
	    float ctf1 = calc_ctf1(ctf, g1, gamma, s);
	    float ctf2 = ctf1 * ctf1 / calc_noise(ctf, n4, s);
	    return ctf2;
	}

    };

}



#endif


