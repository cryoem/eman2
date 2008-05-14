
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <typeconverter.h>
#include <xydata.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyTypeConverter2)
{
    class_< EMAN::EMNumPy >("EMNumPy", init<  >())
        .def(init< const EMAN::EMNumPy& >())
        .def("em2numpy", &EMAN::EMNumPy::em2numpy)
        .def("numpy2em", &EMAN::EMNumPy::numpy2em)
        .staticmethod("em2numpy")
        .staticmethod("numpy2em")
    ;


	import_array();
        python::numeric::array::set_module_and_type("numpy", "ndarray");


	EMAN::vector_to_python<int>();
	EMAN::vector_to_python<long>();
	EMAN::vector_to_python<float>();
	EMAN::vector_to_python<double>();
	EMAN::vector_to_python<EMAN::Transform3D>();
	EMAN::vector_to_python<std::string>();
	EMAN::vector_to_python<EMAN::EMData*>();
	EMAN::vector_to_python<EMAN::Pixel>();
	EMAN::vector_to_python<EMAN::EMObject>();
	EMAN::vector_to_python<EMAN::Vec3f>();
	EMAN::vector_to_python<EMAN::IntPoint>();
	EMAN::vector_to_python< std::vector<EMAN::Vec3f> >();
	EMAN::vector_from_python<int>();
	EMAN::vector_from_python<long>();
	EMAN::vector_from_python<float>();
	EMAN::vector_from_python<EMAN::Transform3D>();
	EMAN::vector_from_python<std::string>();
	EMAN::vector_from_python<EMAN::EMData*>();
	EMAN::vector_from_python<EMAN::Pixel>();
	EMAN::vector_from_python<EMAN::EMObject>();
	EMAN::vector_from_python<EMAN::Vec3f>();
	EMAN::vector_from_python<std::vector<float> >();
	EMAN::map_to_python_2<unsigned int, unsigned int>();
	EMAN::map_to_python<int>();
	EMAN::map_to_python<long>();
	EMAN::map_to_python<float>();
	EMAN::map_to_python<double>();
	EMAN::map_to_python<std::string>();
	EMAN::map_to_python<vector<string> >();

	EMAN::map_from_python<int>();
	EMAN::map_from_python<long>();
	EMAN::map_from_python<float>();
	EMAN::map_from_python<std::string>();
	EMAN::map_from_python<vector<string> >();

    EMAN::EMObject_to_python();
	EMAN::Dict_to_python();
	EMAN::Dict_from_python();

	EMAN::tuple3_to_python<EMAN::IntPoint>();
	EMAN::tuple3_to_python<EMAN::FloatPoint>();

	EMAN::tuple3_to_python<EMAN::IntSize>();
	EMAN::tuple3_to_python<EMAN::FloatSize>();

	EMAN::tuple3_from_python<EMAN::IntPoint, int>();
	EMAN::tuple3_from_python<EMAN::FloatPoint, float>();

	EMAN::tuple3_from_python<EMAN::IntSize, int>();
	EMAN::tuple3_from_python<EMAN::FloatSize, float>();

	EMAN::tuple3_from_python<EMAN::Vec3i, int>();
	EMAN::tuple3_from_python<EMAN::Vec3f, float>();

	EMAN::emobject_array_from_python();
	EMAN::emobject_emdata_from_python();
	EMAN::emobject_string_from_python();
	EMAN::emobject_xydata_from_python();
	EMAN::emobject_transform3d_from_python();
	EMAN::emobject_null_from_python();

	implicitly_convertible<int, EMAN::EMObject>();
	implicitly_convertible<float, EMAN::EMObject>();
	implicitly_convertible<double, EMAN::EMObject>();
	implicitly_convertible<const char*, EMAN::EMObject>();
	implicitly_convertible<void*, EMAN::EMObject>();
	implicitly_convertible<EMAN::EMData*, EMAN::EMObject>();
	implicitly_convertible<EMAN::XYData*, EMAN::EMObject>();
	implicitly_convertible<EMAN::Transform3D*, EMAN::EMObject>();

	EMAN::MArrayND_to_python<2>();
	EMAN::MArrayND_to_python<3>();
	EMAN::MCArrayND_to_python<2>();
	EMAN::MCArrayND_to_python<3>();

}

