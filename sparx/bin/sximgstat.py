#!/usr/bin/env python
import global_def
from   global_def     import *

def main():
	import os
	import sys
	from optparse    import OptionParser

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + "  stack1 <stack2> <mask> --ccc --fsc file --inf --rad=r"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option( "--ccc", action="store_true", default=False, help="print cross corelation coefficient" )
	parser.add_option( "--fsc", type="string",       default="",    help="calculate resolution curve" )
	parser.add_option( "--inf", action="store_true", default=False, help="print basic infomation of the img" )
	parser.add_option( "--rad", type="int",          default=-1,    help="radius of operation" )

        (options,args) = parser.parse_args( arglist[1:] )
     
	if len(args)<1 or len(args)>3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
		sys.exit(-1)


	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	from applications import imgstat
	global_def.BATCH = True
	imgstat( args, options.ccc, options.fsc, options.inf, options.rad )

if __name__=="__main__":
	main()
