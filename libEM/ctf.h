#ifndef eman__ctf__h__
#define eman__ctf__h__ 1

#include <string>
#include <map>
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
	    CTF_CURVE_NOISE_S,
	    CTF_CURVE_ABS_AMP_S,
	    CTF_CURVE_RELATIVE_SNR,
	    CTF_CURVE_ABS_SNR,
	    CTF_CURVE_SNR_WIENER_FILTER,
	    CTF_CURVE_WIENER_CTF_CORRECTION1,
	    CTF_CURVE_WIENER_CTF_CORRECTION2,
	    CTF_CURVE_TOTAL_CURVE
	};
    };
    

    class EMData;
    class XYData;


    class SimpleCtf : public Ctf {
    public:
	float defocus;
	float bfactor;
	float amplitude;
	float ampcont;
	float noise1;
	float noise2;
	float noise3;
	float noise4;
	float voltage;
	float cs;
	float apix;

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
    };

}



#endif


