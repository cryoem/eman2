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

#ifndef eman__serio_h__
#define eman__serio_h__ 1

#include "imageio.h"

#ifdef _WIN32
typedef unsigned long long uint64_t;
#else
#include "stdint.h"
#endif

namespace EMAN
{
	/** SER (Series File Format) is a file format created by Dr. Chris Boothroyd.
	 * The file format is used by TIA (Tecnai imaging and analysis), which is the
	 * program used on FEI Tecnai and Titan microscopes for acquiring and displaying
	 * scanned images and spectra.
	 *
	 * Each .ser file stores a number of 1D or 2D images.
	 * We do not support complex SER image for now.
	 *
	 * http://www.microscopy.cen.dtu.dk/~cbb/info/TIAformat/index.html
	 * */
	class SerIO : public ImageIO
	{
	public:
		explicit SerIO(const string & filename, IOMode rw_mode = READ_ONLY);
		~SerIO();

		DEFINE_IMAGEIO_FUNC;
		static bool is_valid(const void *first_block);
		int get_nimg();

	private:
		enum SerImgMode {
			oneD = 0x4120,
			twoD = 0x4122
		};

		enum SerTagType {
			timeOnly = 0x4152,
			posTime  = 0x4122
		};

		enum SerDataMode {
			SER_UCHAR = 1,
			SER_USHORT,
			SER_UINT,
			SER_CHAR,
			SER_SHORT,
			SER_INT,
			SER_FLOAT,
			SER_DOUBLE,
			SER_COMPLEX8,
			SER_COMPLEX16,
			UNKNOWN
		};

		//This header struct has an alignment issue. It's size is 32 instead of 30.
		//So we need read the item one by one.
		struct SerHeader {
			short 	ByteOrder;
			short 	SeriesID;
			short 	SeriesVersion;
			int		DataTypeID;
			int		TagTypeID;
			int 	TotalNumberElements;
			int		ValidNumberElements;
			int		OffsetArrayOffset;		// could also be 64 bits in new version, we don't actually use this header
			int		NumberDimensions;
		};

		string 	filename;
		IOMode 	rw_mode;
		FILE *	serfile;
		bool 	initialized;
		bool 	is_new_file;

		SerHeader serh;

		uint64_t *  data_offset_array;
		uint64_t *  tag_offset_array;

		int nimg;	//total image number in this file
		int nx;
		int ny;
		int nz;
		int datatypeid;	//1D or 2D image
		int datamode;	//data value mode (int, float, etc.)

		/**helper function to read attributes in dimension array*/
		void read_dim_arr(Dict & dict, int idx);

		/**helper function to read header attributes in data element */
		void read_data_element(Dict & dict);

		/**helper function to read header attributes in data tag*/
		void read_data_tag(Dict & dict);
	};

}

#endif	//eman__serio_h__
