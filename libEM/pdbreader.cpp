/**
 * $Id$
 */

/*
 * Author: Muthu Alagappan, 08/11/2004 (m.alagappan901@gmail.com)
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

#include "pdbreader.h"
#include <vector>
#include <cstring>
#include <string>
#include <stdio.h>
#include <iostream>

using namespace EMAN;


PDBReader::PDBReader()
{
	points = 0;
	n = 0;
}

PDBReader::PDBReader( int nn)
{
	n = nn;
	points = (double *) calloc(4 * n, sizeof(double));
	ter_stop = 0;
}

PDBReader::~PDBReader()
{
	if( points )
	{
		free(points);
		points = 0;
	}
}

void PDBReader::zero()
{
	memset((void *) points, 0, 4 * n * sizeof(double));
}


//copies all relevant information that exists into a new PDBReader object.
PDBReader *PDBReader::copy()
{
	PDBReader *pa2 = new PDBReader();
	pa2->set_number_points(get_number_points());
	double *pa2data = pa2->get_points_array();
	memcpy(pa2data, get_points_array(), sizeof(double) * 4 * get_number_points());
	pa2->pWords = pWords;
	pa2->atomName = atomName;
	pa2->residueName = residueName;
	pa2->chainId = chainId;
	pa2->elementSym = elementSym;
	pa2->tail = tail;
	pa2->head = head;
	pa2->pointInfo = pointInfo;
	pa2->lines = lines;
	return pa2;
}


PDBReader & PDBReader::operator=(PDBReader & pa)
{
	if (this != &pa) {
		set_number_points(pa.get_number_points());
		memcpy(get_points_array(), pa.get_points_array(), sizeof(double) * 4 * get_number_points());
	}
	return *this;
}

size_t PDBReader::get_number_points() const
{
	return n;
}

//reallocates the number of points
void PDBReader::set_number_points(size_t nn)
{
	if (n != nn) {
		n = nn;
		points = (double *) realloc(points, 4 * n * sizeof(double));
	}
}

double *PDBReader::get_points_array()
{
	return points;
}

void PDBReader::set_points_array(double *p)
{
	points = p;
}



// adds purely the x values to a vector
vector<float> PDBReader::get_x() {
	if (count_stop == 0) {count_stop = atomName.size();}
	for (int i=0; i<count_stop; i++) {
        	x.push_back((float)points[4*i]);
       	 	y.push_back((float)points[4*i + 1]);
        	z.push_back((float)points[4*i + 2]);
		resNum.push_back(pointInfo[2*i+1]);
    	}

	return x;
}

// adds purely the y values to a vector
vector<float> PDBReader::get_y() {
	//cout << y.size() << endl;
	return y;
}

// adds purely the z values to a vector
vector<float> PDBReader::get_z() {
	return z;
}

vector<string> PDBReader::get_atomName() {
    	return atomName;
}

vector<string> PDBReader::get_resName() {
    	return residueName;
}

vector<int> PDBReader::get_resNum() {
	return resNum;
}


//Accurately parses a pdb file for all information
bool PDBReader::read_from_pdb(const char *file)
{

	pWords.clear();
	atomName.clear();
	residueName.clear();
	chainId.clear();
	elementSym.clear();
	tail.clear();
	head.clear();
	lines.clear();
	pointInfo.clear();
	ter_stop = 0;
	count_stop = 0;

	struct stat filestat;
	stat(file, &filestat);
	set_number_points(( int)(filestat.st_size / 80 + 1)); //80 bytes per line
#ifdef DEBUG
	printf("PDBReader::read_pdb(): try %4lu atoms first\n", get_number_points());
#endif

	FILE *fp = fopen(file, "r");
	if(!fp) {
		fprintf(stderr,"ERROR in PDBReader::read_pdb(): cannot open file %s\n",file);
		throw;
	}
	char s[200];
	size_t count = 0;

	while ((fgets(s, 200, fp) != NULL)) {
		lines.push_back(s);
		if ((strncmp(s, "ENDMDL", 6) == 0) && (ter_stop==0)){
			ter_stop =1;
			count_stop = count;
		}
		if (strncmp(s, "END", 6) == 0){
			break;
		}
		//if ((strncmp(s, "TER", 3) == 0) && (ter_stop ==0)){
			//if (count !=0) {ter_stop = 1;}
			//count_stop = count;
		//}
		if (strncmp(s, "ATOM", 4) != 0) {
			if (count == 0) {head.push_back(s);}
			else {tail.push_back(s);}
		}
		else {
			pWords.push_back(s);
			atomName.push_back(pWords[count].substr(12,4));
			residueName.push_back(pWords[count].substr(17,3));
			chainId.push_back(pWords[count].substr(21,1));
			elementSym.push_back(pWords[count].substr(76,2));

			float x, y, z, tf;
			int an, sn;
			sscanf(&s[6], "%d", &an);
			sscanf(&s[23], "%d %f %f %f %*f %f", &sn, &x, &y, &z, &tf);

			if (count + 1 > get_number_points()) {
				set_number_points(2 * (count + 1)); 
			}   //makes sure point array is big enough

			points[4 * count] = x;
			points[4 * count + 1] = y;
			points[4 * count + 2] = z;
			points[4 * count + 3] = tf;
			pointInfo.push_back(an);
			pointInfo.push_back(sn);
			count++;
		}
	}

	fclose(fp);
	set_number_points(count);
	return true;
}


//Accurately writes a pdb file once the information is alread stored.
/*
void PDBReader::save_to_pdb(const char *file) const {
	FILE *fp = fopen(file, "w");
	for (int i = 0; i <head.size(); i++) {
		char head2 [200];
		strcpy(head2, head[i].c_str());
		fprintf (fp, head2);
	}
	int m = 0;
	for ( size_t i = 0; i < get_number_points(); i++) {
		string curr = pWords[m];
		string mid, final;
		mid = curr.substr(12, 10);
		final = curr.substr(76,2);
		char mid2 [12];
		strcpy(mid2, mid.c_str());
		char final2 [4];
		strcpy(final2, final.c_str());
		m++;
		fprintf(fp, "ATOM  %5d %10s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\n", pointInfo[2*i], mid2, pointInfo[2*i +1], points[4 * i], points[4 * i + 1], points[4 * i + 2], points[4*i + 3], final2);
	}
	for (int i = 0; i <tail.size(); i++) {
		char tail2 [200];
		strcpy(tail2, tail[i].c_str());
		fprintf (fp, tail2);
	}
	fprintf (fp, "END");
	fclose(fp);

}
*/


//Accurately writes a pdb file once the information is alread stored.
void PDBReader::save_to_pdb(const char *file) const {
	FILE *fp = fopen(file, "w");
	int m = 0;
	for (size_t i =0; i< lines.size(); i++) {
		char liner [200];
		strcpy (liner, lines[i].c_str());
		if (strncmp(liner, "ATOM", 4) != 0) {
			fprintf (fp, liner);
		}
		else {
			string curr = pWords[m];
			string mid, final;
			mid = curr.substr(12, 10);
			final = curr.substr(76,2);
			char mid2 [12];
			strcpy(mid2, mid.c_str());
			char final2 [4];
			strcpy(final2, final.c_str());
			fprintf(fp, "ATOM  %5d %10s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\n", pointInfo[2*m], mid2, pointInfo[2*m +1], points[4 * m], points[4 * m + 1], points[4 * m + 2], points[4*m + 3], final2);
			m++;
		}
	}	
	fclose(fp);

}


//returns a vector of the x,y,z values 
vector<float> PDBReader::get_points() {
vector<float> ret;
for (unsigned int i=0; i<n; i++) {
	ret.push_back((float)points[i*4]);
	ret.push_back((float)points[i*4+1]);
	ret.push_back((float)points[i*4+2]);
}

return ret;
}

//performs a right transform of all x,y,z points given a Transform object
void PDBReader::right_transform(const Transform& transform) {
	for ( unsigned int i = 0; i < 4 * n; i += 4) {
		Transform s = transform.transpose();
		Vec3f v((float)points[i],(float)points[i+1],(float)points[i+2]);
		v= s*v;
		points[i]  =v[0];
		points[i+1]=v[1];
		points[i+2]=v[2];
	}
}

PointArray* PDBReader::makePointArray (const PDBReader& p) {
	PointArray* pArray  = new PointArray;
	p.save_to_pdb("thisFile3.txt");
	pArray->read_from_pdb("thisFile3.txt");
	remove("thisFile3.txt");
	
	return pArray;
}


