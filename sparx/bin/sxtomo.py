#!/usr/bin/env python




import global_def
from   global_def import *
import sys
from   optparse import OptionParser
import os
from utilities import get_image
def main():

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + " prj_stack vol_stack fsc_curve <mask> --box "
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--box",  type="string", help="box size")
	(options,args) = parser.parse_args(arglist[1:])     

        from applications import tomo
        tomo( options.box )


if __name__ == "__main__":
	main()
