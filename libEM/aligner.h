#ifndef eman__aligner_h__
#define eman__aligner_h__ 1


#include "emobject.h"
#include "transform.h"

namespace EMAN
{
    class EMData;
    class Cmp;

    class Aligner
    {
    public:
	virtual ~Aligner() { }

	virtual EMData *align(EMData * this_img, string cmp_name = "") const = 0;

	virtual Dict get_params() const
	{
	    return params;
	}
	
	virtual void set_params(const Dict & new_params)
	{
	    params = new_params;
	}

	virtual TypeDict get_param_types() const = 0;
	virtual string get_name() const = 0;

    protected:
	mutable Dict params;
    };

    class TranslateAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "(TranslateAligner";
	}
	
	static Aligner *NEW()
	{
	    return new TranslateAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("intonly", EMObject::INT);
	    d.put("maxshift", EMObject::INT);
	    return d;
	}
    };

    class Translate3DAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "Translate3DAligner";
	}
	
	static Aligner *NEW()
	{
	    return new Translate3DAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("useparent", EMObject::INT);
	    return d;
	}

    };

    class RotateAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotateAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    return d;
	}
    };


    class RotatePrecenterAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotatePrecenterAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotatePrecenterAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    return d;
	}
    };

    class RotateCHAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotateCHAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateCHAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("irad", EMObject::INT);
	    d.put("orad", EMObject::INT);
	    return d;
	}
    };

    class RotateTranslateAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotateTranslateAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateTranslateAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("usedot", EMObject::INT);
	    d.put("maxshift", EMObject::INT);
	    return d;
	}
    };

    class RotateTranslateBestAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotateTranslateBestAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateTranslateBestAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("maxshift", EMObject::INT);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
    };


    class RotateTranslateRadonAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotateTranslateRadonAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateTranslateRadonAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("maxshift", EMObject::INT);
	    d.put("radonwith", EMObject::EMDATA);
	    d.put("radonthis", EMObject::EMDATA);
	    return d;
	}
    };


    class RotateFlipAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RotateFlipAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateFlipAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("flip", EMObject::EMDATA);
	    d.put("imask", EMObject::INT);
	    return d;
	}
    };

    class RotateTranslateFlipAligner: public Aligner
    {
    public:	
	EMData *align(EMData * this_img, string cmp_name = "") const;

	string get_name() const
	{
	    return "RotateTranslateFlipAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RotateTranslateFlipAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("flip", EMObject::EMDATA);
	    d.put("usedot", EMObject::INT);
	    d.put("maxshift", EMObject::INT);
	    return d;
	}
    };

    class RTFSlowAligner: public Aligner
    {
    public:	
	EMData *align(EMData * this_img, string cmp_name = "") const;

	string get_name() const
	{
	    return "RTFSlowAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RTFSlowAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("flip", EMObject::EMDATA);
	    d.put("maxshift", EMObject::INT);
	    return d;
	}
    };

    class RTFSlowestAligner: public Aligner
    {
    public:	
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RTFSlowestAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RTFSlowestAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("flip", EMObject::EMDATA);
	    d.put("maxshift", EMObject::INT);
	    return d;
	}
    };


    class RTFBestAligner: public Aligner
    {
    public:
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RTFBestAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RTFBestAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("flip", EMObject::EMDATA);
	    d.put("maxshift", EMObject::INT);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
    };


    class RTFRadonAligner: public Aligner
    {
    public:	
	EMData *align(EMData * this_img, string cmp_name = "") const;

	string get_name() const
	{
	    return "RTFRadonAligner";
	}
	
	static Aligner *NEW()
	{
	    return new RTFRadonAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("maxshift", EMObject::INT);
	    d.put("thisf", EMObject::EMDATA);
	    d.put("radonwith", EMObject::EMDATA);
	    d.put("radonthis", EMObject::EMDATA);
	    d.put("radonthisf", EMObject::EMDATA);
	    return d;
	}
    };


    class RefineAligner: public Aligner
    {
    public:	
	EMData *align(EMData * this_img, string cmp_name = "") const;
	
	string get_name() const
	{
	    return "RefineAligner";
	}
	static Aligner *NEW()
	{
	    return new RefineAligner();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("mode", EMObject::INT);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
    };
    
}

#endif
