/*
 * Author: Grant Goodyear, 11/15/2003 (sludtke@bcm.edu)
 * Copyright (c) 2000-2006 Baylor College of Medicine
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * */

#ifdef _WIN32
	#pragma warning(disable:4819)
#endif	//_WIN32

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
            return object(s.sget_value_at(ix,iy,iz));
        } else if (2 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            if (s.is_complex())
                return object(s.cmplx(ix,iy));
            return object(s.sget_value_at(ix,iy));
        } else if (1 == size) {
            int ix = extract<int>(key[0]);
            if (s.is_complex())
                return object(s.cmplx(ix));
            return object(s.sget_value_at(ix));
        } else {
            throw ImageDimensionException("Need 1, 2, or 3 indices.");
        }
    }
	// If we got a string, then we'll check the metadata dictionary
	extract<string> k(key);
	if (k.check()) return object(s.get_attr(k));

	extract<std::wstring> k2(key);
	if (k2.check()) return object(s.get_attr(extract<string>(str(key).encode())));
	
    // not a tuple or an int, so bail
    throw TypeException("Need x,y,z or attribute key", "");
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
//		printf("good2\n");
        if (3 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            int iz = extract<int>(key[2]);
            if (s.is_complex())
                s.set_complex_at(ix,iy,iz,extract<std::complex<float> >(val));
            else
                s.set_value_at(ix,iy,iz,extract<float>(val));
            return;
        } else if (2 == size) {
            int ix = extract<int>(key[0]);
            int iy = extract<int>(key[1]);
            if (s.is_complex())
                s.set_complex_at(ix,iy,extract<std::complex<float> >(val));
            else
                s.set_value_at(ix,iy,extract<float>(val));
            return;
        } else if (1 == size) {
            int ix = extract<int>(key[0]);
            if (s.is_complex())
                s.set_complex_at(ix,0,extract<std::complex<float> >(val));
            else
                s.set_value_at(ix,extract<float>(val));
            return;
        } else {
            throw ImageDimensionException("Need 1, 2, or 3 indices.");
        }
    }
	extract<std::string> k(key);
	extract<EMAN::EMObject> v(val);
	if (k.check() && v.check()) {
		s.set_attr(extract<string>(key),v);
		return;
	}
	
	extract<std::wstring> k2(key);
	if (k2.check() && v.check()) {
		s.set_attr(extract<string>(str(key).encode()),v);
		return;
	}

	
    throw TypeException("Need x,y,z or attribute key", "");
}

