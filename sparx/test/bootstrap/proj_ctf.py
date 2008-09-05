#!/usr/bin/env python


from sparx import *
from EMAN2 import *
from sys import argv
from sys import exit
from sys import stdout
from string import atof
from string import split

if( len(argv) != 4 ) :
    print "Usage: proj vol_input angle_input proj_output"
    exit(-1)


vol_in = argv[1]
angle_in = open( argv[2], "r" )
proj_out = argv[3]

# ignore the first line
line = angle_in.readline()

vol = getImage(vol_in)
volft,kb=prep_vol(vol)

Cs = 2.0
pixel = 4.88
voltage = 100.0
amp_contrast = 0.1
defocus = 13000 

iang = 0

while angle_in :

    if( iang %50 ==0 ) :
        stdout.write( "%4d " % iang )
	stdout.flush()


    items = split( angle_in.readline() )
    if( len(items) < 5 ) :
        break

    psi = atof( items[2] )
    theta = atof( items[3] )
    phi = atof( items[4] )

#    if iang < 1000:
#        defocus = 21000
#    elif iang < 2000:
#        defocus = 25000
#    else:
#        defocus = 29000
#    ctfvol = filt_ctf(vol,defocus,Cs,voltage,pixel,amp_contrast)
#    volft,kb = prep_vol(ctfvol)

    proj=prgs(volft, kb, [ phi,theta,psi, 0.0, 0.0])


    #proj = filt_ctf(proj, defocus, Cs, voltage, pixel, amp_contrast)
    #proj = filt_ctf(proj, defocus, Cs, voltage, pixel, amp_contrast)


    proj.set_attr_dict( {'phi':phi,'theta':theta,'psi':psi,'s2x':0.0,'s2y':0.0,'mirror':0.0} )
    proj.set_attr_dict( {'active':1} )

    proj.set_attr_dict( {'ctf_applied':1.0, 'defocus':defocus, 'amp_contrast':amp_contrast} )
    proj.set_attr_dict( {'voltage':voltage,'Cs':Cs, 'Pixel_size':pixel } )


    proj.write_image( proj_out, iang )

    stdout.write( "." )
    stdout.flush()

    
    if( iang % 50 == 49 ) :
        stdout.write( "\n" )
 
    iang = iang + 1

stdout.write( "\n" )

