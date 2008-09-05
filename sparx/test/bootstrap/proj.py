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

curtvol = getImage( vol_in )
volft,kb=prep_vol(curtvol)

iproj = 0

while angle_in :

    if( iproj %50 ==0 ) :
        stdout.write("%4d " % iproj)
	stdout.flush()

    items = split( angle_in.readline() )
    if( len(items) < 5 ) :
        break

    psi = atof( items[2] )
    theta = atof( items[3] )
    phi = atof( items[4] )

    proj=prgs(volft, kb, [ phi,theta,psi, 0.0, 0.0])
    proj.set_attr_dict( {'phi':phi,'theta':theta,'psi':psi,'sx':0.0,'sy':0.0,'mirror':0.0} )
    proj.set_attr_dict( {'active':1} )

#    Cs = curt.get_attr( "Cs" )
#    pixel = curt.get_attr( "pixel" )
#    defocus = curt.get_attr( "defocus" )
#    voltage = curt.get_attr( "voltage" )
#    amp_contrast = curt.get_attr( "amp_contrast" )
#    proj.set_attr_dict( {'ctf_applied':0.0, 'defocus':defocus, 'amp_contrast':amp_contrast} )
#    proj.set_attr_dict( {'voltage':voltage,'Cs':Cs, 'pixel':pixel } )

    proj.write_image( proj_out, iproj)
    stdout.write( "." )
    stdout.flush()

    if( iproj % 50 == 49 ) :
        stdout.write( "\n" )

    iproj = iproj + 1

stdout.write( "\n" )

