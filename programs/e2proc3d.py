#!/bin/env python

# $Id$

# todo: verify the filters who have the same names in proc3d
#       and proc2d have the same implementation

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
    
    parser.add_option("--apix", type="float", help="")
    parser.add_option("--mult", type="float", help="")
    parser.add_option("--add", type="float", help="")
    
    # next: signoise
    
    
    
if __name__ == "__main__":
    main()
    

