#!/bin/env python

# $Id$

from EMAN2 import *
from optparse import OptionParser
import sys
import os.path

def main():
    progname = os.path.basename(sys.argv[0])
    usage = progname + " options inputfile outputfile"
    parser = OptionParser(usage)

    parser.add_option("--shrink", type="int", action="append", help=".")
    parser.add_option("--scale", type="float", action="append", help="")
    parser.add_option("--clip", type="float", nargs=3, action="append", help="")
    parser.add_option("--clipc", type="float", nargs=3, action="append", help="")
    parser.add_option("--fftclip", type="float", nargs=3, action="append", help="")
    
    # todo: rot, trans, apix
 
    
    
if __name__ == "__main__":
    main()
    

