#!/bin/env python

# $Id$

# usage: proc2d inputfile outputfile [first=<n>] [last=<n>] [inplace]  [rot=<angle>] [trans=<dx>,<dy>]   [scale=<sca>] [clip=<x,y>] [shrink=<n>] [meanshrink=<n>] [apix=<A/pix>]   [norefs] [average]  [split=n] [ctfsplit] [interlv=<file 2 interleave>] [sym=<Cn>] [plt[=<plt file>]] [setsfpairs]  [randomize=<dr>,<da>,<flip>]  [selfcl[=<steps>][,<mode>]] [radon] [rfp] [mraprep] [phot] rfilt=<filtername><:param1=value1><...> [mrc]   [pif] [hdf] [png] [em] [spider] [spider-single] [pgm[=<low>,<hi>]] 

#from EMAN2 import *
from optparse import OptionParser

def filter_callback (option, opt, value, parser):
    print "opt=", opt
    print "value=", value
    
def main():
    parser = OptionParser()
    parser.add_option("--verbose", type="int", help="verbose level [1-5]")
    
    parser.add_option("--first", type="int", help="first image")
    parser.add_option("--last", type="int", help="last image")

    parser.add_option("--apix", type="float", help="")
    parser.add_option("--inplace", action="store_true", help="inplace ")
    parser.add_option("--rot", type="float", help="rotation angle")
    
    parser.add_option("--trans", type="float", nargs=2, help="translation")
    parser.add_option("--scale", type="float", help="scale")
    parser.add_option("--clip", type="float", nargs=2, help="2D clip")

    parser.add_option("--shrink", type="int", help="shrink factor")
    parser.add_option("--meanshrink", type="int", help="")

    parser.add_option("--list", type="string", help="")
    parser.add_option("--exclude", type="string", help="")
    parser.add_option("--filefilt", type="string", help="")

    parser.add_option("--calcsf", type="int", help="")
    parser.add_option("--sfout", type="string", help="")
    
    parser.add_option("--fileavg", type="string", help="")
    
    parser.add_option("--norefs", action="store_true", help="")
    parser.add_option("--average", action="store_true", help="")
    parser.add_option("--ctfsplit", action="store_true", help="")

    parser.add_option("--snrfilt", string="string", help="")
    parser.add_option("--wiener", string="string", help="")
    
    parser.add_option("--split", type="int", help="")
    parser.add_option("--interlv", type="string", help="")
    parser.add_option("--sym", type="string", help="")

    parser.add_option("--plt", type="string", help="")
    parser.add_option("--randomize", type="float", nargs=3, help="")
    parser.add_option("--selfcl", type="int", nargs=2, help="")
    
    parser.add_option("--setsfpairs",  action="store_true", help="")
    parser.add_option("--radon",  action="store_true", help="")
    parser.add_option("--rfp",  action="store_true", help="")
    parser.add_option("--mraprep",  action="store_true", help="")
    parser.add_option("--phot",  action="store_true", help="")
 
    parser.add_option("--filter", type="string", action="callback", help="filter name", callback=filter_callback)
    
    (options, args) = parser.parse_args()
    
    
    print "verbose=", options.verbose
    
    
if __name__ == "__main__":
    main()
    
