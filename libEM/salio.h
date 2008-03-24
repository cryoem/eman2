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

#ifndef eman__salio_h__
#define eman__salio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/** A SAL image is an image from Perkin Elmer PDS Microdensitometer.
	 * A SAL image consists of 2 files: 1 header file "X.hdr" and a
	 * data file "X.img". Header file is in ASCII format. Data file is
	 * in binary format.
	 *
	 * Each pair of hdr/img SAL files contains 1 2D image.
	*/

	class SalIO:public ImageIO
	{
	  public:
		explicit SalIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SalIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);

	  private:
		static const char *HDR_EXT;
		static const char *IMG_EXT;
		static const char *MAGIC;

		enum ScanAxis
		{
			X_SCAN_AXIS,
			Y_SCAN_AXIS
		};

		enum ScanMode
		{
			NON_RASTER_SCAN,
			RASTER_SCAN
		};

	  private:
		string filename;
		IOMode rw_mode;
		FILE *sal_file;
		bool initialized;
		int nx;
		int ny;
		int record_length;
		ScanMode scan_mode;
		float pixel;
	};

}


#endif	//eman__salio_h__
