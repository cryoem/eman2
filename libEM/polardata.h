/**
* $Id$
*/

/*
 * Author: Grant Tang, 02/06/2006 (gtang@bcm.edu)
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

#ifndef eman__polardata_h__
#define eman__polardata_h__

#include <map>
#include <vector>
#include "log.h"
#include "exception.h"

using std::map;
using std::vector; 

namespace EMAN
{
	static const int MAXFFT=32768;
	class EMData;
	
	/** a general data structure for a matrix with variable x dim size for different y*/
	class UnevenMatrix
	{
	public:
		UnevenMatrix() : data(0) {
			printf("Welcome to UnevenMatrix\n");
		}
		virtual ~UnevenMatrix() {
			if(data) {
				delete data;
				data = 0;
			}
			printf("Destructor of UnevenMatrix...\n");
		}
		
		/** get the x dim size for a given y 
		 *  @param y int the y value for which we need the x dim size 
		 *  @return in tthe x dim size */
		inline int get_xsize(int y) {
			int xsize = desc_data[y].x1 - desc_data[y].x0;
			if( xsize <= 0 ) {
				throw InvalidValueException(xsize, "xsize <= 0");
			}
			else {
				return xsize;
			}
		}
		
		/** get the minimal x dim value for a given y 
		 *  @param y int the y value for which we need the corresponding minimal x dim value 
		 *  @return int the minimal x dim value for a given y */
		inline int get_xmin(int y) {
			return desc_data[y].x0;
		}
		
		/** get the maximal x dim value for a given y, note: x1 is one out of the max 
		 *  @param y int the y value for which we need the corresponding maximal x dim value 
		 *  @return int the maximal x dim value for a given y */
		inline int get_xmax(int y) {
			return desc_data[y].x1 - 1;
		}
		
		/**	get the total size of the data block
		 *  @return int the number of float stored in data */
		inline int get_size() {
			return tot_size;
		}

#ifdef DEBUG
		void print_UnevenMatrix() {
			printf("print function in UnevenMatrix\n");
		} 
#endif 	//DEBUG

	protected:
		/** allocation memory for data array
		 * @return int 
		 * @exception InvalidValueException if the desc_data map size is zero*/
		void alloc_data();
		
		/** struct to define x dimension size for each y, x0 is inclusive, 
		 * x1 is one after the maximum, [x0, x1), so the corresponding x dim size is (x1-x0) */
		struct Xdim {
			int x0;
			int x1;
			Xdim() {}
			Xdim(int i0, int i1) : x0(i0), x1(i1) {
				if( x1 <= x0 ) {
					LOGERR("x dimension error ... x0(%d) > x1(%d)", x0, x1);
				}
			}
		};
		
		/** store all data in one dimension float array for cache efficiency, 
		 * we calculate the offset for x, y dimension */
		float * data;
		
		/** describe for each y, the x dimension's size and range */
		map< int, Xdim >	desc_data;
		
		/** the total size of the data*/
		int tot_size;
	};
	
	/** a specialized image class for storing the results of a transform 
	 *  from EMData to polar coordinates, currently support 2D only.
	 *  data on x dimension may be variable size, which is defined in 
	 *  map< int, Xdim > desc_data*/
	class PolarData : public UnevenMatrix
	{
	public:
		PolarData() {printf("Welcome to PolarData class... \n");}
		
		/** Construct a PolarData object from a EMData
		 * @param image   the EMData object to be converted
		 * @param xcen the x dimension of the center
		 * @param ycen the y dimension of the center
		 * @param mode
		 */
		PolarData(EMData * image, int xcen, int ycen, string mode);
		
		virtual ~PolarData() {
			printf("Destructor of PolarData...\n");
		}
		
#ifdef DEBUG		
		/** a testing function */
		void print_polar() {
			printf("PolarData class is a specialized image class for storing the \n");
			printf("results of a transform from EMData to polar coordinates, currently \n");
			printf("support 2D only. Data on x dimension may be variable size, which \n");
			printf("is defined in map< int, Xdim > desc_data\n");
		}
		
		int test_init_desc_data();
#endif	//DEBUG

	private:
		/** the ring weights for each radius r */
		map< int, float >	weight;
		
		/** calculate the number of element for each ring
		 * @param first_ring the ring number for the first ring
		 * @param last_ring the ring number for the last ring
		 * @param skip step of ring 
		 * @param mode half mode('H'/'h') or full mode('F'/'f')
		 * @return the vector for number of elements of each ring */
		vector<int> Numrinit(int first_ring, int last_ring, int skip, string mode);
	
		/**Returns the smallet power by which 2 has to be raised to obtain an integer kess equal n
	 	* @param n int
	 	* @return int */
		int log2(int n);
		
		/** calculate ring weights for rotational alignment 
		 * @param numr number of element for each ring
		 * @param mode half mode('H'/'h') or full mode('F'/'f')
		 * @return vector<float> the weight for each ring */
		vector<float> ringwe( vector<int> numr, string mode );
	};

}


#endif	//eman__polardata_h__
