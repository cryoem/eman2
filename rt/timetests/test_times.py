#!/usr/bin/env python

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
					t3d = Transform3D(EULER_EMAN,az,0,0)
					t3d.set_posttrans(dx,dy)
					e3 = EMData()
					e3.set_size(n+i,n+i,1)
					e3.process_inplace('testimage.squarecube', {'axis':'x', 'edge_length':1, 'fill':1} )
					if ( i == 0 ) : e3.write_image("testc.img")
					
					
					e4 = e3.copy()
					
					e3.rotate_translate(t3d)
					
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