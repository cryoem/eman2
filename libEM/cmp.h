/**
 * $Id$
 */
#ifndef eman_cmp__h__
#define eman_cmp__h__ 1


#include "emobject.h"

namespace EMAN
{
    
    class EMData;
    class Transform;
    
    /** Image Comparison methods. The bigger the result, the better.
     */
    
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
	    return "Dot";
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

    class VarianceCmp : public Cmp
    {
    public:
	float cmp(EMData * em, Transform * transform = 0) const;

	string get_name() const
	{
	    return "Variance";
	}

	static Cmp *NEW()
	{
	    return new VarianceCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
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
	    return "Phase";
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

    class FRCCmp : public Cmp
    {
    public:
	float cmp(EMData * em, Transform * transform = 0) const;

	string get_name() const
	{
	    return "FRC";
	}

	static Cmp *NEW()
	{
	    return new FRCCmp();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("with", EMObject::EMDATA);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}
    };

    template<> Factory<Cmp>::Factory();
}


#endif
