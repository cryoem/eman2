#ifndef eman_reconstructor_template_h__
#define eman_reconstructor_template_h__ 1

#include "reconstructor.h"

namespace EMAN {
    
    /** XYZReconstructor is a reconstructor template for defining new
     * reconstructors. Please add your own code at the proper place.
     *
     * 1) Replace all 'XYZ' with your new reconstructor name.
     * 2) Define the reconstructor parameter names and types in get_param_types().
     * 3) Implement the reconstructor in setup(), insert_slice(), and finish();
     */
    class XYZReconstructor : public Reconstructor
    {
    public:
	XYZReconstructor();
	~XYZReconstructor();

	/** initialize the reconstructor
	 */
	int setup();
	
	/** insert each image slice to the reconstructor. You may call
	 * this function multiple times.
	 */
	int insert_slice(EMData * slice, const Rotation & euler);

	/** finish reconstruction and return the complete model.
	 */
	EMData *finish();

	string get_name() const
	{
	    return "XYZ";
	}
	static Reconstructor *NEW()
	{
	    return new XYZReconstructor();
	}
	
	/** Add your reconstructor parameter names and types in
	 * get_param_types(). For available parameter types, please
	 * refer class EMObject.
	 * 
	 * As an example, XYZReconstructor has 3 parameters:
	 *    int size;
	 *    float patratio;
	 *    vector<float> snr;
	 */
	TypeDict get_param_types() const
	{
	    TypeDict d;
	    d.put("size", EMObject::INT);
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

    /** Add your new reconstructor to ReconstructorFactoryExt().
     */
    class ReconstructorFactoryExt {
    public:
	ReconstructorFactoryExt() {
	    Factory<Reconstructor>::add(&XYZReconstructor::NEW);
	}
    };

    static ReconstructorFactoryExt rf_ext;
}


#endif
