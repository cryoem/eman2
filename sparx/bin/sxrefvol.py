#!/usr/bin/env python

import os
import sys
import global_def
from global_def import *
from optparse import OptionParser
from sys import argv, exit

def main():
	from sys import argv
	progname = os.path.basename(argv[0])
	usage = progname + " rawvol_1 rawvol_2 ... resolution_1 resolution_2 ... refvol.hdf mask"
	
	if( len(argv)%2==0 ):
	    print "Usage: ", usage
	    exit(-1)
 
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

        nvol = (len(argv)-3)/2

	from applications import refvol
	vollist = argv[1:1+nvol]
	fsclist = argv[1+nvol:1+2*nvol]
	output  = argv[1+2*nvol]
	mask    = argv[2+2*nvol]

        refvol( vollist, fsclist, output, mask )

if __name__ == "__main__":
	main()
