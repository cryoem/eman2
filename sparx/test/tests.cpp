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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "emdata.h"
#include <iostream>

using namespace EMAN;
using namespace std;
int main(int argc, char** argv) {
float x,y,w,z,s1,s2,ss;
ss=0.f;
for (int k=0; k<200; k++) {
s1=0.f;s2=0.f;
for (int j=0; j<30000; j++) {
x=1.002f;y=1.002;w=1.002;z=1.002;
  for (int i=0; i<1000; i++) { x=x*z*1.001f;y=y*w*1.02f;}
//cout <<"   " << x <<"   "<<y <<"   "<< w <<"   "<<z << endl;
//cout << "   "<<j<<"   "<<"  x  " << x <<endl;
s1+=x+y;
s2=z+w;
}
ss=ss+s1+s2;
}
cout <<"   " << x <<"   "<<y <<"   "<< w <<"   "<<z <<"   "<<ss << endl;
}
 
