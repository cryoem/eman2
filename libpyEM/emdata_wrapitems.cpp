#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include "emdata.h"

using namespace boost::python;
using namespace EMAN;

object emdata_getitem(object self, object key) {
    EMData& s = extract<EMData&>(self);
    // Check if the key is a single int
    extract<int> x(key);
    if (x.check()) {
        // it is
        int i = x();
        if (s.is_complex())
            return object(s.cmplx(i));
        return object(s(i));
    }
    // Check if the key is a tuple 
    // (which it should be, if more than one index is passed)
    extract<tuple> t(key);
    if (t.check()) {
        int size = len(key);
        if (3 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            int iz = extract<int>(key[2]);
            if (s.is_complex())
                return object(s.cmplx(ix,iy,iz));
            return object(s(ix,iy,iz));
        } else if (2 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            if (s.is_complex())
                return object(s.cmplx(ix,iy));
            return object(s(ix,iy));
        } else if (1 == size) {
            int ix = extract<int>(key[0]);
            if (s.is_complex())
                return object(s.cmplx(ix));
            return object(s(ix));
        } else {
            throw ImageDimensionException("Need 1, 2, or 3 indices.");
        }
    }
    // not a tuple or an int, so bail
    throw TypeException("Expected 1, 2, or 3 *integer* indices", "");
}

void emdata_setitem(object self, object key, object val) {
    EMData& s = extract<EMData&>(self);
    extract<int> x(key);
    if (x.check()) {
        int i = x();
        if (s.is_complex()) 
            s.cmplx(i) = extract<std::complex<float> >(val);
        else
            s(i) = extract<float>(val);
        return;
    }
    extract<tuple> t(key);
    if (t.check()) {
        int size = len(key);
        if (3 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            int iz = extract<int>(key[2]);
            if (s.is_complex())
                s.cmplx(ix,iy,iz) = extract<std::complex<float> >(val);
            else
                s(ix,iy,iz) = extract<float>(val);
            return;
        } else if (2 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            if (s.is_complex())
                s.cmplx(ix,iy) = extract<std::complex<float> >(val);
            else
                s(ix,iy) = extract<float>(val);
            return;
        } else if (1 == size) {
            int ix = extract<int>(key[0]);
            if (s.is_complex())
                s.cmplx(ix) = extract<std::complex<float> >(val);
            else
                s(ix) = extract<float>(val);
            return;
        } else {
            throw ImageDimensionException("Need 1, 2, or 3 indices.");
        }
    }
    throw TypeException("Expected 1, 2, or 3 integer indices", "");
}


