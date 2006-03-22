/**
 * $Header$
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
            EMArray(size_t xsize, size_t ysize=1, size_t zsize=1) 
                : nx(xsize),ny(ysize),nz(zsize),
                  xoff(0),yoff(0),zoff(0) {
                size = nx*ny*nz;
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
