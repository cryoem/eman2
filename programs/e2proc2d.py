#!/bin/env python

# $Id$

from EMAN2 import *
from optparse import OptionParser
import sys
import os.path
import math

def main():
    usage = os.path.basename(sys.argv[0]) + " options inputfile outputfile"

    parser = OptionParser(usage)

    parser.add_option("--apix", type="float", help="")
    parser.add_option("--average", action="store_true", help="")
    parser.add_option("--calcsf", type="int", help="")

    parser.add_option("--clip", type="float", nargs=2, help="2D clip")
    parser.add_option("--ctfsplit", action="store_true", help="")
    parser.add_option("--exclude", type="string", help="")

    parser.add_option("--fileavg", type="string", help="")
    parser.add_option("--filter", type="string", action="append", help="filter name")
    parser.add_option("--first", type="int", help="first image")

    parser.add_option("--inplace", action="store_true", help="inplace ")
    parser.add_option("--interlv", type="string", help="")
    parser.add_option("--last", type="int", help="last image")

    parser.add_option("--list", type="string", help="")
    parser.add_option("--meanshrink", type="int", help="")
    parser.add_option("--mraprep",  action="store_true", help="")

    parser.add_option("--norefs", action="store_true", help="")
    parser.add_option("--outtype", type="string", help="")
    parser.add_option("--phot",  action="store_true", help="")

    parser.add_option("--plt", type="string", help="")
    parser.add_option("--radon",  action="store_true", help="")
    parser.add_option("--randomize", type="float", nargs=3, help="")

    parser.add_option("--rfp",  action="store_true", help="")
    parser.add_option("--rot", type="float", help="rotation angle")
    parser.add_option("--scale", type="float", help="scale")

    parser.add_option("--selfcl", type="int", nargs=2, help="")
    parser.add_option("--setsfpairs",  action="store_true", help="")
    parser.add_option("--sfout", type="string", help="")

    parser.add_option("--shrink", type="int", help="shrink factor")
    parser.add_option("--split", type="int", help="")
    parser.add_option("--sym", type="string", help="")

    parser.add_option("--trans", type="float", nargs=2, help="translation")
    parser.add_option("--verbose", type="int", help="verbose level [1-5]")


    (options, args) = parser.parse_args()

    inputfile = args[0]
    outputfile = args[1]

    csym = 0
    if options.sym:
        if options.sym[0] != 'c':
            print "Error: Only C symmetry currently supported!"
            sys.exit(1)
        else:
            csym = int(options.sym[1:])
        
    if options.randomize:
        options.randomize = (options.randomize[0], options.randomize[1] * math.pi / 180.0, 
                             options.randomize[2])

    scl = 90
    if options.selfcl:
        scl = options.selfcl[0] / 2

    n0 = 0
    n1 = -1

    if options.first:
        n0 = options.first
    if options.last:
        n1 = options.last

    MAXMICROCTF = 1000
    defocus_val = [0] * MAXMICROCTF
    bfactor_val = [0] * MAXMICROCTF

    
    d = EMData()
    nimg = EMUtil.get_image_count(inputfile)
    if nimg <= n1 or n1 < 0:
        n1 = nimg - 1

    ld = EMData()
    print "nimg = ", nimg
    
    for i in range(n0, n1+1):
        d.read_image(inputfile, i)
        nx = d.get_xsize()

        if options.ctfsplit and (i == n0 or (not EMUtil.is_same_ctf(d, ld))):
            ctf = d.get_ctf()

            for j in range(1, options.ctfsplit):
                if defocus_val[j] == ctf.get_defocus() and bfactor_val[j] == ctf.get_bfactor():
                    break
            if options.ctfsplit <= j:
                options.ctfsplit = j + 1
                print "New CTF at " + i
                
            defocus_val[j] = ctf.get_defocus()
	    bfactor_val[j] = ctf.get_bfactor()
        


    
if __name__ == "__main__":
    main()
    
