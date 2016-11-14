/**
 * $Header$
 */

/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
 */

#ifndef EMARRAY_H_
#define EMARRAY_H_

#include "exception.h"
#include <vector>

using std::vector;

namespace EMAN {

    /** EMArray -- 1-, 2-, or 3-D array of types T.
     *
     * Sometimes you need an array of something.
     *
     * Data is ordered with the first element increasing the fastest.
     */
    template<class T> 
    class EMArray {
        public:
            /** Constructor */
            explicit EMArray(size_t xsize, size_t ysize=1, size_t zsize=1) 
                : nx(xsize),ny(ysize),nz(zsize),
                  xoff(0),yoff(0),zoff(0) {
                size = (size_t)nx*ny*nz;
                array = new T[size];
            }
            /** Destructor */
            ~EMArray() { delete [] array; }
            /** Set the array offsets.  By default an EMArray ranges
             *  from 0..nx-1,0..ny-1,0..nz-1.  Setting the offsets to (1,1,1)
             *  changes the ranges to 1..nx,1..ny,1..nz, for example.
             */
            void set_array_offsets(const size_t xoff_=0, const size_t yoff_=0,
                                   const size_t zoff_=0) {
                xoff = xoff_; yoff = yoff_; zoff = zoff_;
            }
            /** Set the array offsets using an offset vector.  This
             *  method exists mainly for restoring the original offsets.
             */
            void set_array_offsets(const vector<int> offsets) {
                set_array_offsets(offsets[0],offsets[1],offsets[2]);
            }
            /** Get an integer array containing the array offsets. */
            vector<int> get_array_offsets() {
                vector<int> offsets;
                offsets.push_back(xoff);
                offsets.push_back(yoff);
                offsets.push_back(zoff);
                return offsets;
            }
            /** Index into the array using Fortran-style indexing:
             *
             *    EMArray<double> arr(2,2) // 2-D 2x2 array of doubles
             *    arr(0,0) = 1.f;
             */
            T& operator()(const size_t ix, const size_t iy=0, 
                          const size_t iz=0) {
                long pos = (ix-xoff) + ((iy-yoff)+(iz-zoff)*ny)*nx;
#ifdef BOUNDS_CHECKING
                if (pos < 0 || pos >= long(size))
                    throw OutofRangeException(0, size-1, pos, "emarray");
#endif // BOUNDS_CHECKING
                return *(array + pos);
            }

            int get_xsize() const { return nx; }
	    int get_ysize() const { return ny; }
	    int get_zsize() const { return nz; }

        private:
            const size_t nx, ny, nz;
            size_t size;
            size_t xoff, yoff, zoff;
            T* array;
            // disable copy and assignment
            EMArray(const EMArray&); 
            EMArray& operator=(const EMArray&);
    };
}
#endif // EMARRAY_H_
