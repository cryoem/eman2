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

	virtual int align(EMData* em, Vec3f& result, string comp_name = "") const = 0;
	
	virtual Dict get_params() const { return params; }
	virtual void set_params(const Dict& new_params) { params = new_params; }
    
	virtual string get_name() const = 0;
    
    protected:
	mutable Dict params;
    };

    class TransAligner : public Aligner {
    public:
	string get_name() const { return "(TransAligner"; }
	static Aligner* NEW() { return new TransAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };

    class RotAligner : public Aligner {
    public:
	string get_name() const { return "RotAligner"; }
	static Aligner* NEW() { return new RotAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    
    class RotPreCenteredAligner : public Aligner {
    public:
	string get_name() const { return "RotPreCenteredAligner"; }
	static Aligner* NEW() { return new RotPreCenteredAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };

    class RefineAligner : public Aligner {
    public:
	string get_name() const { return "RefineAligner"; }
	static Aligner* NEW() { return new RefineAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    class RTAligner : public Aligner {
    public:
	string get_name() const { return "RTAligner"; }
	static Aligner* NEW() { return new RTAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    class RFAligner : public Aligner {
    public:
	string get_name() const { return "RFAligner"; }
	static Aligner* NEW() { return new RFAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
        
    class RTFAligner : public Aligner {
    public:
	string get_name() const { return "RTFAligner"; }
	static Aligner* NEW() { return new RTFAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    class RTFSlowestAligner : public Aligner {
    public:
	string get_name() const { return "RTFSlowestAligner"; }
	static Aligner* NEW() { return new RTFSlowestAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    class RTFSlowAligner : public Aligner {
    public:
	string get_name() const { return "RTFSlowAligner"; }
	static Aligner* NEW() { return new RTFSlowAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    class RTFRadonAligner : public Aligner {
    public:
	string get_name() const { return "RTFRadonAligner"; }
	static Aligner* NEW() { return new RTFRadonAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    
        
    class Trans3DAligner : public Aligner {
    public:
	string get_name() const { return "Trans3DAligner"; }
	static Aligner* NEW() { return new Trans3DAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
    };
    
    class RotCHAligner : public Aligner {
    public:
	string get_name() const { return "RotCHAligner"; }
	static Aligner* NEW() { return new RotCHAligner(); }
	
	int align(EMData* em, Vec3f& result, string comp_name = "") const;
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
