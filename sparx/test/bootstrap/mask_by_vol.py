
from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atof

if len(argv) != 5 :
    print "Usage: mask_by_vol input mask_vol thresh output"
    exit(-1)


input = getImage( argv[1] )
mask  = getImage( argv[2] )
th    = atof( argv[3] )
output = argv[4]


mask.process_inplace( "eman1.threshold.binary", {"value":th} )

input *= mask

dropImage( input, output )



