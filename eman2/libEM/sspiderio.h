/**
 * $Id$
 */

/*
 * Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

#ifndef eman__single_spider_io_h__
#define eman__single_spider_io_h__ 1

#include "spiderio.h"


namespace EMAN
{
	/** Single Spider Image I/O class. For Spider and Single Spider
	 * image format, please refer spiderio.h.
	 * @see spiderio.h
	 */	
	class SingleSpiderIO:public SpiderIO
	{
	  public:
		/** SingleSpiderIO constructor.
		 * @param filename The filename of a single spider image file.
		 * @param rw_mode Read/Write file open mode.
		 */
		explicit SingleSpiderIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SingleSpiderIO();
		
		/** Write a single SPIDER format header to an image.
		 *
		 * @param dict A keyed-dictionary storing the header information.
		 * @param image_index The index of the image to write.
		 * @param area The region to write data to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian Whether to use the host machine
		 *        endian to write out or not. If false, use the
		 *        endian opposite to the host machine's endian.
		 * @return 0 if OK; 1 if error.
		 */
		int write_header(const Dict & dict, int image_index = 0, const Region* area = 0,
						 EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
						 bool use_host_endian = true);
		
		/** Write data to an image.
		 *
		 * @param data An array storing the data.
		 * @param image_index The index of the image to write.
		 * @param area The region to write data to.
		 * @param filestoragetype The image data type used in the output file.
		 * @param use_host_endian Whether to use the host machine
		 *        endian to write out or not. If false, use the
		 *        endian opposite to the host machine's endian.
		 * @return 0 if OK; 1 if error.
		 */
		int write_data(float *data, int image_index = 0, const Region* area = 0,
					   EMUtil::EMDataType filestoragetype = EMUtil::EM_FLOAT,
					   bool use_host_endian = true);

		static bool is_valid(const void *first_block);
		
		bool is_single_image_format() const
		{
			return true;
		}
		
	  protected:
		  bool is_valid_spider(const void *first_block);
	};

}

#endif
