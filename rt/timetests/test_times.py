#!/usr/bin/env python

# this was added by d.woolford. It is envisaged that a variety of timing tests will be added to this script to generate useful
# data for benchmarking purposes.

from EMAN2 import *
from math import *
from random import *
from time import *

def main():
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
	
if __name__ == "__main__":
    main()