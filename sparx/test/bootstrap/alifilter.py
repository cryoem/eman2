#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit
from string import split, atof

def read_fsc( fscfile ) :

    f = open( fscfile, 'r' )

    fsc = []
    fsc.append( [] )
    fsc.append( [] )
    fsc.append( [] )

    line = f.readline()
    while len(line) > 0 :
 
        items = split( line )

        if len(items) == 3 :
            fsc[0].append( atof(items[0]) )
            fsc[1].append( atof(items[1]) )
            fsc[2].append( atof(items[2]) )
        line = f.readline()
     
    return fsc


if len(argv) != 4:
    print "ali_filter vol_in fsc vol_out"
    exit(-1)

vol = getImage( argv[1] )
fsc = read_fsc( argv[2] )


if(fsc[1][0] < 0.5) : fsc[1][0] = 1.0
filt = filt_from_fsc2(fsc, 0.05)

refvol = filt_table(vol, filt)

cs = refvol.phase_cog()
refvol = fshift(refvol, -cs[0], -cs[1] -cs[2])

refvol = refine_with_mask(refvol)
dropImage( refvol, argv[3], 's' )



