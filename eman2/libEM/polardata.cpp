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

#include <new>
#include <cmath>
#include "emdata.h"
#include "polardata.h"
#include "exception"
#include "util.h"

#ifdef DEBUG
	#include <iostream>
	using std::cout;
	using std::endl;
#endif

#ifdef WIN32
	#ifndef M_PI
		#define M_PI 3.14159265358979323846f
	#endif	//M_PI
#endif	//WIN32

using namespace EMAN;
using std::bad_alloc; 


void UnevenMatrix::alloc_data()
{
	if( desc_data.size() == 0 ) {
		LOGERR("No data space need to be allocated for UnevenMatrix, check you desc_data...");
		throw InvalidValueException((float)desc_data.size(), "desc_data size == 0");
	}
	
	int size = 0;	//total number of float need to be stored in data block
	map< int, Xdim >::const_iterator iter;
	int y = 0;	
	for( iter = desc_data.begin(); iter != desc_data.end(); ++iter ) {
		y = (*iter).first;
		size += get_xsize(y);
	}
	
	this->tot_size = size;
	try {
		this->data = new float[size];
	}
	catch( bad_alloc exception ) {
		LOGERR("memory allocation for UnevenMatrix failed");
	}
	catch(...) {
		LOGERR("Unknown error in memory allocation for UnevenMatrix");
	}
}

PolarData::PolarData(EMData * , int , int , string )
{
//	int nsam = image->get_xsize();
//	int nrow = image->get_ysize();
	
	
//	int nring = numr.size()/3; 
		
}
/*
EMData* PolarData::calc_ccf(EMData * em)
{
	em = 0;
}
*/
/*  I had to block this section, it is in conflict with my python code PAP
vector<int> PolarData::Numrinit(int first_ring, int last_ring, int skip, string mode)
{
	float dpi;
	if(mode == "f" || mode == "F") {
		dpi = 2 * M_PI;
	}
	else {
		dpi = M_PI;
	}
	
	vector<int>	numr;
	int lcirc = 1;
	int jp = 0; 
	int ip = 0;
	for(int k=first_ring; k<=last_ring; ++k) {
		numr.push_back(k);
		jp = int(dpi * k + 0.99999999);
		ip = (int)pow(2, (double)log2(jp));
		if ( k+skip <= last_ring && jp > ip+std::floor((double)ip/2) ) {
#ifdef _WIN32
			ip = _cpp_min(MAXFFT, 2*ip);
#else
			ip = std::min<int>(MAXFFT, 2*ip);
#endif	//_WIN32
		}
		if ( k+skip > last_ring && jp > ip+std::floor((double)ip/5) ) {
#ifdef _WIN32
			ip = _cpp_min(MAXFFT, 2*ip);	
#else
			ip = std::min<int>(MAXFFT, 2*ip);
#endif
		}
		
		numr.push_back(lcirc);
		numr.push_back(ip);
		lcirc += ip;
	}
	
	--lcirc;
	return numr;
}

int PolarData::log2(int n)
{
	int m = 1;
    int k = -1;
    int i;
    while ( m <= n) {
       i = m;
       ++k;
       m = 2*i;
    }
    
    return k;
}

vector<float> PolarData::ringwe( vector<int> numr, string mode )
{
	float dpi;
	if(mode == "f" || mode == "F") {
		dpi = 2 * M_PI;
	}
	else {
		dpi = M_PI;
	}
	
	vector<float>	wr;
	int nring = numr.size()/3;
	wr.resize(nring);
	float maxrin = (float)numr[numr.size()-1];
	for(int i=0; i<nring; ++i) {
		wr[i] = (numr[i*3]*dpi)*maxrin / (float)(numr[2+i*3]); 
	}
	
	return wr;
}
*/
#ifdef DEBUG
int PolarData::test_init_desc_data()
{
	desc_data[0] = Xdim(0, 10);
	desc_data[3] = Xdim(20, 30);
	desc_data[6] = Xdim(3, 5);
	desc_data[10] = Xdim(2, 12);
	//desc_data[10] = Xdim(22, 12);	//error test
	
	cout << "allocation of data success" << endl;
	
	//int j = 10;
	//cout << "log2(10) = " << log2(j) << endl;
	
	return 0;
}
#endif	//DEBUG
