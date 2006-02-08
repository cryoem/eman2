/**
 * $Id$
 */
 
#ifndef eman__polardata_h__
#define eman__polardata_h__

#include <map>

using std::map;

namespace EMAN
{
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
		 *  @param y the y value for which we need the x dim size 
		 *  @return the x dim size */
		inline int get_xsize(int y) {
			return desc_data[y].x1 - desc_data[y].x0;
		}
		
		/** get the minimal x dim value for a given y 
		 *  @param y the y value for which we need the corresponding minimal x dim value 
		 *  @return the minimal x dim value for a given y */
		inline int get_xmin(int y) {
			return desc_data[y].x0;
		}
		
		/** get the maximal x dim value for a given y, note: x1 is one out of the max 
		 *  @param y the y value for which we need the corresponding maximal x dim value 
		 *  @return the maximal x dim value for a given y */
		inline int get_xmax(int y) {
			return desc_data[y].x1 - 1;
		}
		
		void print_UnevenMatrix() {
			printf("print function in UnevenMatrix\n");
		} 
		
	protected:
		/** struct to define x dimension size for each y, x0 is inclusive, 
		 * x1 is one after the maximum, [x0, x1), so the corresponding x dim size is (x1-x0) */
		struct Xdim {
			int x0;
			int x1;
		};
		
		/** store all data in one dimension float array for cache efficiency, 
		 * we calculate the offset for x, y dimension */
		float * data;
		
		/** describe for each y, the x dimension's size and range */
		map< int, Xdim >	desc_data;
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
		 * @param em   the EMData object to be converted
		 * @param xcen the x dimension of the center
		 * @param ycen the y dimension of the center
		 */
		PolarData(EMData * em, int xcen, int ycen);
		
		virtual ~PolarData() {
			printf("Destructor of PolarData...\n");
		}
		
		/** a testing function */
		void print_polar() {
			printf("PolarData class is a specialized image class for storing the \n");
			printf("results of a transform from EMData to polar coordinates, currently \n");
			printf("support 2D only. Data on x dimension may be variable size, which \n");
			printf("is defined in map< int, Xdim > desc_data\n");
		}
	
	private:
		/** the weight for each y dimension */
		map< int, float >	weight;
	};

}

#endif	//eman__polardata_h__
