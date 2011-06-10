/**
 * $Id$
 */

/*
 * Author: Grant Tang, 06/07/2011 (gtang@bcm.edu)
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

#ifndef eman__omapio_h__
#define eman__omapio_h__ 1

#include "imageio.h"

namespace EMAN
{
	/**DSN6 MAP is composed of a series of records which are all 512 bytes long.
	 * The first is a header which contains all the information required to
	 * intercept the rest of the file. The subsequent bricks contain blocks of
	 * electron density. Each density sample is one byte in size. The order of
	 * the samples within a brick is "x fast, y medium, and z slow".
	 * http://www.uoxray.uoregon.edu/tnt/manual/node104.html
	 *
	 * - BRICK (BINARY; input only) - treated as DSN6
	 * - BRIX (BINARY; input and output) - similar to DSN6, but has a human-readable first record that contains the cell constants etc.
	 * - OMAP (input and output) - treated as DSN6
	 *  */
	class OmapIO:public ImageIO
	{
	public:
		explicit OmapIO(const string & omapname, IOMode rw_mode = READ_ONLY);
		~OmapIO();

		DEFINE_IMAGEIO_FUNC;

		static bool is_valid(const void *first_block, off_t file_size = 0);

	private:
		/**The data in the header is composed of 256 short integers.*/
		struct OmapHeader
		{
			short xstart;		//x start
			short ystart;		//y start
			short zstart;		//z start
			short nx;			//x extent
			short ny;			//y extent
			short nz;			//z extent
			short apix_x;		//x sampling rate
			short apix_y;		//y sampling rate
			short apix_z;		//z sampling rate
			short header10;		//header18*A Cell Edge
			short header11;		//header18*B Cell Edge
			short header12;		//header18*C Cell Edge
			short alpha;		//alpha
			short beta;			//beta
			short gamma;		//gamma
			short header16;		//header19*(253-3)/(Rhomax-Rhomin)
			short header17;		//(3*Rhomax-253*Rhomin)/(Rhomax-Rhomin)
			short scale;		//Cell Constant Scaling Factor
			short header19;		//constant 100

			short unused[237];	//unused space
		};

		string filename;
		IOMode rw_mode;
		FILE *omapfile;
		OmapHeader omaph;

		bool is_big_endian;
		bool initialized;
		bool is_new_file;
	};

}

#endif	//eman__omapio_h__
