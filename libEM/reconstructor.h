/**
 * $Id$
 */
#ifndef eman_reconstructor_h__
#define eman_reconstructor_h__ 1

#include "emobject.h"
#include <vector>
#include <map>
#include <string>
#include <math.h>

using std::vector;
using std::map;
using std::string;

namespace EMAN
{

    /**
     * these magic numbers and macros need to be removed later
     */    
#define I5G	(10.4 / (M_PI*M_PI))	// used for 5x5x5 interpolation
#define I4G	(8.8 / (M_PI*M_PI))	// used for 4 interpolation
#define I3G	(6.4 / (M_PI*M_PI))	// used for 3 and 5x5x5 interpolation
#define I2G	(4.0 / (M_PI*M_PI))	// 2 interpolation


    class EMData;
    class Rotation;


    /** Reconstructor is the base class for all reconstructors.
     * Each subclass Reconstructor defines a way to do 3D recontruction.
     * Each specific Reconstructor type has a unique ID name. This name
     * is used to create a Reconstructor instance or do a reconstruction.
     *
     * Typical usages of Reconstructors are as follows:
     *
     * 1. How to get all the Reconstructor types:
     *
     *    vector<string> all_reconstructors = Factory<Reconstructor>.instance()->get_list();
     *
     * 2. How to use a Reconstructor
     *
     *    Reconstructor* r = Factory<Reconstructor>.instance()->get("Fourier");
     *    r->setup;
     *    r->insert->slice(slice1, euler1);
     *    insert more
     *    EMData* result = r->finish();
     *
     * 3. How to define a new Reconstructor type
     *
     *    A new XYZReconstructor class should implement the following functions:
     * 
     *        int setup();
     *        int insert_slice(EMData * slice, const Rotation & euler);
     *        EMData * finish();
     *        string get_name() const { return "XYZ"; }
     *        static Reconstructor *NEW() { return new XYZReconstructor(); }
     *        TypeDict get_param_types() const;
     */
    class Reconstructor
    {
    public:
	virtual ~Reconstructor()
	{
	}
	virtual int setup() = 0;
	virtual int insert_slice(EMData * slice, const Rotation & euler) = 0;
	virtual EMData *finish() = 0;

	virtual string get_name() const = 0;

	virtual Dict get_params() const
	{
	    return params;
	}

	virtual void set_params(const Dict & new_params)
	{
	    params = new_params;
	}
	
	virtual TypeDict get_param_types() const = 0;

      protected:
	mutable Dict params;
    };

    /** Fourier space 3D reconstruction
     */
    class FourierReconstructor : public Reconstructor
    {
      public:
	FourierReconstructor();
	~FourierReconstructor();

	int setup();
	int insert_slice(EMData * slice, const Rotation & euler);
	EMData *finish();

	string get_name() const
	{
	    return "Fourier";
	}
	static Reconstructor *NEW()
	{
	    return new FourierReconstructor();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
	    d.put("mode", EMObject::INT);
	    d.put("weight", EMObject::FLOAT);
	    d.put("dlog", EMObject::INT);
	    return d;
	}

      private:
	EMData * image;
	int nx;
	int ny;
	int nz;
    };

    /** Fourier space 3D reconstruction with slices already Wiener filtered.
     */
    class WienerFourierReconstructor : public Reconstructor
    {
      public:
	WienerFourierReconstructor();
	~WienerFourierReconstructor();

	int setup();
	int insert_slice(EMData * slice, const Rotation & euler);
	EMData *finish();

	string get_name() const
	{
	    return "WienerFourier";
	}
	static Reconstructor *NEW()
	{
	    return new WienerFourierReconstructor();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
	    d.put("mode", EMObject::INT);
	    d.put("padratio", EMObject::FLOAT);
	    d.put("snr", EMObject::FLOATARRAY);
	    return d;
	}

      private:
	EMData * image;
	int nx;
	int ny;
	int nz;
    };

    /** Real space 3D reconstruction using back projection.
     * 
     * Back-projection is a method of 3D reconstruction from 2D
     * projections. It is based on superposing 3D functions
     * ("back-projection bodies") obtained by translating the
     * 2D projections along the directions of projection. 
     */
    class BackProjectionReconstructor : public Reconstructor
    {
      public:
	BackProjectionReconstructor();
	~BackProjectionReconstructor();

	int setup();
	int insert_slice(EMData * slice, const Rotation & euler);
	EMData *finish();

	string get_name() const
	{
	    return "BackProjection";
	}
	static Reconstructor *NEW()
	{
	    return new BackProjectionReconstructor();
	}

	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
	    d.put("weight", EMObject::FLOAT);
	    return d;
	}

      private:
	EMData * image;
	int nx;
	int ny;
	int nz;
    };

    template<> Factory<Reconstructor>::Factory();
}

#endif
