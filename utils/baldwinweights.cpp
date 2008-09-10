/*
 * Author: Baldwin and Woolford 2008
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

#include "emdata.h"
#include "util.h"
#include "transform.h"

using namespace EMAN;

#include <iostream>
using std::cout;
using std::endl;
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

void test_shared_pointer() {
	shared_ptr<Transform> p( new Transform);
	cout << p.use_count() << endl;
	Transform* t = p.get();
	cout << p.use_count() << endl;
	shared_ptr<Transform> p2(p);
	cout << p.use_count() << endl;
	cout << p2.use_count() << endl;
	shared_ptr<Transform> p3 = p;
	cout << p.use_count() << endl;
	cout << p3.use_count() << endl;
}

int main(int argc, char *argv[])
{
// 	int nx = 64;
// 	int P = (int)((1.0+0.25)*nx+1);
// 	float r = (float)(nx+1)/(float)P;
// 	int mFreqCutoff = 2;
// 	float mDFreq = 0.2;
	
// 	float* W = Util::getBaldwinGridWeights(mFreqCutoff, P, r, mDFreq,0.5,0.2);
// 	cout << "Test 2" << endl;
// 	W = Util::getBaldwinGridWeights(3, 35, 0.9, 1,0.5,0.2);
	
	test_shared_pointer();
// 	delete [] W;
    return 0;
}




