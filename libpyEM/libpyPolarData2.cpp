
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
    class_< EMAN::PolarData >("PolarData", init<  >())
        .def(init< const EMAN::PolarData& >())
        .def(init< EMAN::EMData*, int, int >())
        .def("print_polar", &EMAN::PolarData::print_polar)
    ;

}

