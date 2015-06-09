/*
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
 * it under the terms of the GNU General Public License as published by.edge
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
 */

#ifndef eman_tomoseg_h_
#define eman_tomoseg_h_

#include "emdata.h"

namespace EMAN
{
	class TomoObject{
	public:
		TomoObject(vector<Vec3i> allpt, float md=1, int slice=0);
		void write_imod(FILE *fp);
		int get_size();
		
		/* Starting from connecting two terminis and iteratively bend the polyline by breaking at the 
		 * most distant point from the polyline.*/
		int bend();
		
		/* Enlong the line segments and connect the adjecent objects. */
		void enlong(EMData *bwmap,EMData *skelmap);
		
	public:
		vector<int> ptid;	// ID of the cornor points in the object
		vector<int> segid;	// ID of the segments in the labeled map that are in this object.
		vector<Vec3i> allpoints;// Cooridinates of all points in this object
		float maxdist;		// Max distance to bend the polyline.
		int nowslice;
	};
	class TomoSeg
	{
	public:
		TomoSeg(){verb=false;};
		void set_verb(){verb=!verb;};
		
		/* Read the black-white map. Each segment in the map should be labeled. */
		int read_bwmap(EMData *map);
		
		/* Read the skeleton map. The map should be labeled the same as the bwmap */
		int read_skelmap(EMData *map);
		
		/*  Generate objects from numo largest segments in the skeleton map. */
		int generate_objects(int numo, float maxdist, int nowslice);
		
		/* Check the value of 8 neighbors. Return the number of white neighbors. */
		int check_neighbors(int x,int y);
		
		/* Write objects to imod format */
		void write_imod(const char *file);
		
	private:
		EMData *skelmap;
		EMData *bwmap;
		bool verb;
		vector<TomoObject> objs;
	};
}

#endif