#ifndef eman_cmp__h__
#define eman_cmp__h__ 1


#include "emobject.h"

namespace EMAN
{
    
    class EMData;
    class Transform;
    
    // todo: 
    // define a function to return native value; another function to
    // return normlized value from 0-1, where bigger means better.

    class Cmp
    {
    public:
	virtual ~Cmp() { }
	virtual float cmp(EMData * em, Transform * transform = 0) const = 0;
	virtual TypeDict get_param_types() const = 0;
	virtual string get_name() const = 0;
	
	virtual Dict get_params() const
	{
	    return params;
	}

	virtual void set_params(const Dict & new_params)
	{
	    params = new_params;
	}
	
    protected:
	mutable Dict params;
    };

    class DotCmp : public Cmp
    {
    public:
	float cmp(EMData * em, Transform * transform = 0) const;

	string get_name() const
	{
	    return "DotCmp";
	}

	static Cmp *NEW() 
	{
	    return new DotCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("evenonly", EMObject::INT);
	    return d;
	}

    };

    class LinearCmp : public Cmp
    {
    public:
	float cmp(EMData * em, Transform * transform = 0) const;

	string get_name() const
	{
	    return "LinearCmp";
	}

	static Cmp *NEW()
	{
	    return new LinearCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("invert", EMObject::INT);
	    d.put("keepzero", EMObject::INT);
	    return d;
	}
    };

    class PhaseCmp : public Cmp
    {
    public:
	float cmp(EMData * em, Transform * transform = 0) const;

	string get_name() const
	{
	    return "PhaseCmp";
	}

	static Cmp *NEW()
	{
	    return new PhaseCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);

	    return d;
	}
    };

    class FSCCmp : public Cmp
    {
    public:
	float cmp(EMData * em, Transform * transform = 0) const;

	string get_name() const
	{
	    return "FSCCmp";
	}

	static Cmp *NEW()
	{
	    return new FSCCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
    };


    typedef Cmp *(*CmpType) ();

    class CmpFactory
    {
    public:
	static CmpFactory *instance();

	void add(CmpType cmp);
	Cmp *get(string comp_name);
	Cmp *get(string comp_name, const Dict & params);

	vector<string> get_list();

    private:
	CmpFactory();
	CmpFactory(const CmpFactory &) { }
	~CmpFactory();

	void force_add(CmpType cmp);

	static CmpFactory *my_instance;
	map<string, CmpType> my_dict;
    };

}


#endif
