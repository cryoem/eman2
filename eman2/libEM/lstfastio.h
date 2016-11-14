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

#ifndef eman__lstiofast_h__
#define eman__lstiofast_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** A LSX file is a high performance ASCII file that contains a list of image
	 * numbers and file names. Each line of a LSX file has the following format:
	 * after the #LSX magic number are 2 lines, the second line containing
	 * # comment
	 * # line-length
	 * Denoting the length of each line in the file. If writing additional
	 * lines to the file, this length must not be exceeded.
	 *
	 * Additional '#' comment lines are NOT permitted within the file
	 *
	 * lines have the form:
	 * reference_image_index  reference-image-filename comments
	 */
		
		
	class LstFastIO:public ImageIO
	{
	  public:
		explicit LstFastIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~LstFastIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

		bool is_single_image_format() const
		{
			return false;
		}
		int get_nimg();
	  private:
		string filename;
		IOMode rw_mode;
		FILE *lst_file;

		bool is_big_endian;
		bool initialized;
		int nimg;
		unsigned int line_length;
		unsigned int head_length;

		ImageIO *imageio;
		string ref_filename;

		int last_lst_index;
		int last_ref_index;

		int calc_ref_image_index(int image_index);
		static const char *MAGIC;
	};

}


#endif	//eman__lstio_h__
