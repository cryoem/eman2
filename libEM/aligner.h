#ifndef eman__aligner_h__
#define eman__aligner_h__ 1


#include "emobject.h"
#include "transform.h"

namespace EMAN {

    class EMData;
    class Cmp;
    
    class Aligner {
    public:
	virtual ~Aligner() {}

	virtual EMData* align(EMData* this_img, string cmp_name = "") const = 0;
	
	virtual Dict get_params() const { return params; }
	virtual void set_params(const Dict& new_params) { params = new_params; }
    
	virtual string get_name() const = 0;
    
    protected:
	mutable Dict params;	
    };

    /*
     * params: EMData* with, int intonly, int useparent, int maxshift
     */
    class TranslateAligner : public Aligner {
    public:
	string get_name() const { return "(TranslateAligner"; }
	static Aligner* NEW() { return new TranslateAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };

    class RotateAligner : public Aligner {
    public:
	string get_name() const { return "RotateAligner"; }
	static Aligner* NEW() { return new RotateAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    
    class RotatePrecenterAligner : public Aligner {
    public:
	string get_name() const { return "RotatePrecenterAligner"; }
	static Aligner* NEW() { return new RotatePrecenterAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };

    class RefineAligner : public Aligner {
    public:
	string get_name() const { return "RefineAligner"; }
	static Aligner* NEW() { return new RefineAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RotateTranslateAligner : public Aligner {
    public:
	string get_name() const { return "RotateTranslateAligner"; }
	static Aligner* NEW() { return new RotateTranslateAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RotateFlipAligner : public Aligner {
    public:
	string get_name() const { return "RotateFlipAligner"; }
	static Aligner* NEW() { return new RotateFlipAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RotateTranslateFlipAligner : public Aligner {
    public:
	string get_name() const { return "RotateTranslateFlipAligner"; }
	static Aligner* NEW() { return new RotateTranslateFlipAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RTFSlowestAligner : public Aligner {
    public:
	string get_name() const { return "RTFSlowestAligner"; }
	static Aligner* NEW() { return new RTFSlowestAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RTFSlowAligner : public Aligner {
    public:
	string get_name() const { return "RTFSlowAligner"; }
	static Aligner* NEW() { return new RTFSlowAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RTFRadonAligner : public Aligner {
    public:
	string get_name() const { return "RTFRadonAligner"; }
	static Aligner* NEW() { return new RTFRadonAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class Translate3DAligner : public Aligner {
    public:
	string get_name() const { return "Translate3DAligner"; }
	static Aligner* NEW() { return new Translate3DAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    class RotateCHAligner : public Aligner {
    public:
	string get_name() const { return "RotateCHAligner"; }
	static Aligner* NEW() { return new RotateCHAligner(); }
	
	EMData* align(EMData* this_img, string cmp_name = "") const;
    };
    
    
    typedef Aligner* (*AlignerType)();
    
    class AlignerFactory {
    public:
	static AlignerFactory* instance();
	
	void add(AlignerType aligner);
	Aligner* get(string align_name);
	Aligner* get(string align_name, const Dict& params);
	vector<string> get_list();
	
    private:
	AlignerFactory();
	AlignerFactory(const AlignerFactory& a);
	~AlignerFactory();

	void force_add(AlignerType aligner);
	
	static AlignerFactory* my_instance;
	map<string, AlignerType> my_dict;
    };


}

#endif
