/**
* $Id$
*/

/*
 * Author: Muthu Alagappan, 07/14/2009, (m.alagappan901@gmail.com)
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

#ifndef eman_pdbreader_h_
#define eman_pdbreader_h_

#include "emdata.h"
#include "transform.h"
#include "pointarray.h"

#if defined NFFT || NFFT2
extern "C"
{
	#include "nfft.h"
	#include "utils.h"
}
#endif

#include <sys/stat.h>
#include <vector>


namespace EMAN
{
	/** PointArray defines a double array of points with values in a 3D space. */
	class PDBReader
	{
		public:
		enum Density2PointsArrayAlgorithm
		{
			PEAKS_SUB, PEAKS_DIV, KMEANS
		};

		public:
		PDBReader();
		explicit PDBReader(int nn);
		~PDBReader();
		void zero();
		PDBReader *copy();
		PDBReader & operator=(PDBReader & pa);
		size_t get_number_points() const;
		void set_number_points(size_t nn);

		/** Reads and parses all information from file
		*
		* @param The .pdb file that you want to read from
		*/
		bool read_from_pdb(const char *file);

		/** Saves all atom information into a pdb in the official format
		*
		* @param The file that you want the pdb info written to
		*/
		void save_to_pdb(const char *file) const;

		/** Returns the double array of points
		*
		* @return A double array of points
		*/
		double *get_points_array();

		/** Allows the user to set the double array of points
		*
		* @param A double array of points
		*/
		void set_points_array(double *p);

		/** Returns all x,y,z triplets packed into a vector<float>
		*
		* @return All points packed into a vector<float>
		*/
		vector<float> get_points();

		/** Does Transform*v as opposed to v*Transform (as in the transform function)
		 * @param transform an EMAN2 Transform object 
		*/
		void right_transform(const Transform& transform);
		PointArray* makePointArray (const PDBReader& p);

		vector<float> get_x();
		vector<float> get_y();
		vector<float> get_z();
		vector<string> get_atomName();
		vector<string> get_resName();
		vector<int> get_resNum();
		

		private:
		double *points;
		vector<int> pointInfo;
		vector<string> pWords;
		vector<string> atomName;
		vector<string> residueName;
		vector<string> chainId;
		vector<string> elementSym;
		vector<string> tail;
		vector<string> head;
		vector<string> lines;
		size_t n;
		int ter_stop;

		vector<float> x;
		vector<float> y;
		vector<float> z;
		vector<int> resNum;
	};
}

#endif
