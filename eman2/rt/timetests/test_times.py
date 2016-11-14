#!/usr/bin/env python

#
# Author: David Woolford, 09/12/2007 (woolford@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# this was added by d.woolford. It is envisaged that a variety of timing tests will be added to this script to generate useful
# data for benchmarking purposes.

from EMAN2 import *
from math import *
from random import *
from time import *

def main():
	precision_test()

def timetest():
	n = 10000000
	it = 10
	a = [ random() for i in range(n) ]
	
	time1 = clock()
	for j in range(it):
		b = Util.get_stats_cstyle(a)
	time2 = clock()
	dt1 = time2 - time1
	
	print "It took %f seconds to gets the stats (in c style) for %d random numbers, %d times" %(dt1,n,it)
	
	time1 = clock()
	for j in range(it):
		b = Util.get_stats(a)
	time2 = clock()
	dt1 = time2 - time1
	
	print "It took %f seconds to gets the stats (in c++ style) for %d random numbers, %d times" %(dt1,n,it)

def precision_test():
        """test RotateTranslateAligner ....................."""
		
	n = 32
	err_dx = 0
	err_dy = 0
	err_az = 0
	
	for i in range(0,1):
		for dx in range(-n/4+1,n/4-1):
			for dy in range(-n/4+1,n/4-1):
				for az in range(0,181,5):
					t3d = Transform({'type':'eman', 'az':az, 'alt':alt, 'phi':phi})
					t3d.set_trans(dx,dy)
					e3 = EMData()
					e3.set_size(n+i,n+i,1)
					e3.process_inplace('testimage.squarecube', {'axis':'x', 'edge_length':1, 'fill':1} )
					if ( i == 0 ) : e3.write_image("testc.img")
					
					
					e4 = e3.copy()
					
					e3.transform(t3d)
					
					e5 = e4.align('rotate_translate', e3, {'maxshift':n/2})
					
					soln_dx = e5.get_attr("align.dx")
					soln_dy = e5.get_attr("align.dy")
					soln_az = e5.get_attr("align.az")
					
					if (abs(soln_dx-dx) > err_dx ): err_dx = abs(soln_dx-dx)
					if (abs(soln_dy-dy) > err_dy ): err_dy = abs(soln_dy-dy)
					if (abs(soln_az-az)%180 > err_az ): err_az = abs(soln_az-az)%180
	
	
	print "RotateTranslateAlign precision is %f %f %f)" %(err_dx,err_dy,err_az)

if __name__ == "__main__":
    main()
