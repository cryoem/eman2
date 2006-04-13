
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <emdata.h>
#include <polardata.h>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpyPolarData2)
{
    class_< EMAN::UnevenMatrix, boost::noncopyable >("UnevenMatrix", init<  >())
        .def("get_xsize", &EMAN::UnevenMatrix::get_xsize)
        .def("get_xmin", &EMAN::UnevenMatrix::get_xmin)
        .def("get_xmax", &EMAN::UnevenMatrix::get_xmax)
        .def("get_size", &EMAN::UnevenMatrix::get_size)
    ;

    class_< EMAN::PolarData, bases< EMAN::UnevenMatrix > , boost::noncopyable >("PolarData", init<  >())
        .def(init< EMAN::EMData*, int, int, std::string >())
    ;

}

