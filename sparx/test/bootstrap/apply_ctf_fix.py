#!/usr/bin/env python

from sparx import *
from EMAN2 import *
from sys import argv
from sys import exit
from sys import stdout
from random import random

if( len( argv ) != 3 ) :
    print "usage: apply_ctf image_file_in image_file_out"
    exit(-1)

i_file = argv[1]
o_file = argv[2]

Cs = 2.0 # in mm
amp_contrast = 0.1
voltage = 300
pixel = 4.88

nimage = EMUtil.get_image_count(i_file)

vol = EMData()

for i in xrange(nimage) :

    if( i%50 == 0 ) : 
        stdout.write( "%4d " % i )
	stdout.flush()

    vol.read_image( i_file, i )

    r = random()

    if r < 0.3333:
        defocus = 22500
    elif r < 0.6666:
        defocus = 30000
    else:
        defocus = 40000

    ctfvol = filt_ctf(vol,defocus,Cs,voltage,pixel,amp_contrast)

    ctfvol.set_attr_dict( {"ctf_applied":0.0,"defocus":defocus,"amp_contrast":amp_contrast} )
    ctfvol.set_attr_dict( {"voltage":voltage, "Cs":Cs, "Pixel_size":pixel} )

    ctfvol.write_image( o_file, i )

    stdout.write( "." )
    stdout.flush()

    if( i%50 == 49 ) :
        stdout.write( "\n" )
	stdout.flush()

stdout.write( "\n" )
stdout.flush( )


