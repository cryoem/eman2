
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
    class_< EMAN::UnevenMatrix >("UnevenMatrix", init<  >())
        .def(init< const EMAN::UnevenMatrix& >())
        .def("get_xsize", &EMAN::UnevenMatrix::get_xsize)
        .def("get_xmin", &EMAN::UnevenMatrix::get_xmin)
        .def("get_xmax", &EMAN::UnevenMatrix::get_xmax)
        //.def("print_UnevenMatrix", &EMAN::UnevenMatrix::print_UnevenMatrix)
    ;

    class_< EMAN::PolarData, bases< EMAN::UnevenMatrix >  >("PolarData", init<  >())
        .def(init< const EMAN::PolarData& >())
        .def(init< EMAN::EMData*, int, int, string >())
        //.def("print_polar", &EMAN::PolarData::print_polar)
    ;

}

