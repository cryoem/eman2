#!/usr/bin/env python

from sparx import *
from EMAN2 import *
from sys   import exit
from sys   import argv
from string import split
from string import atof
from string import atoi
from random import normalvariate
from math   import exp

if( len( argv ) != 4 ) :
    print "Usage: make_special_model angledoc output_volume output_proj"
    exit(-1)

angledoc = argv[1]
volume_file = argv[2]
proj_file = argv[3]

size = 64
sigma_x = 8
sigma_y = 12
sigma_z = 10
maxium  = 10.0
base = model_gauss( sigma_x, size, size, size, sigma_y,sigma_z )
[mean,sigma,min,max]=Util.infomask(base,None, True)
base = base * (maxium/max)

threshold = exp(-1.5)

#point 1
xcenter = 45
ycenter = 39
zcenter = 43
xsigma  = 6
ysigma  = 5
zsigma  = 4
maxium  = 3.0

print "make point at %4d, %4d, %4d" % ( xcenter, ycenter, zcenter )
p1 = model_gauss( xsigma, size,size,size, ysigma, zsigma, xcenter,ycenter,zcenter )
[mean,sigma,min,max]=Util.infomask(p1,None, True)
p1 = p1/max
p1.process_inplace( "threshold.belowtozero", {"minval": threshold} )
p1 = p1*maxium

#point 2
xcenter = 37
ycenter = 27
zcenter = 19
xsigma  = 6
ysigma  = 8
zsigma  = 7
maxium  = 2.0

print "make point at %4d, %4d, %4d" % ( xcenter, ycenter, zcenter )
p2 = model_gauss( xsigma, size,size,size, ysigma, zsigma, xcenter,ycenter,zcenter )
[mean,sigma,min,max]=Util.infomask(p2,None, True)
p2 = p2/max
p2.process_inplace( "threshold.belowtozero", {"minval": threshold} )
p2 = p2*maxium

#point 3
xcenter = 21
ycenter = 25
zcenter = 43
xsigma  = 9
ysigma  = 8
zsigma  = 6
maxium  = -2.0

print "make point at %4d, %4d, %4d" % ( xcenter, ycenter, zcenter )
p3 = model_gauss( xsigma, size,size,size, ysigma, zsigma, xcenter,ycenter,zcenter )
[mean,sigma,min,max]=Util.infomask(p3,None, True)
p3 = p3/max
p3.process_inplace( "threshold.belowtozero", {"minval": threshold} )
p3 = p3*maxium

#point 4
xcenter = 19
ycenter = 35
zcenter = 27
xsigma  = 6
ysigma  = 4
zsigma  = 5
maxium  = 1.2

print "make point at %4d, %4d, %4d" % ( xcenter, ycenter, zcenter )
p4 = model_gauss( xsigma, size,size,size, ysigma, zsigma, xcenter,ycenter,zcenter )
[mean,sigma,min,max]=Util.infomask(p4,None,True)
p4 = p4/max
p4.process_inplace( "threshold.belowtozero", {"minval": threshold} )
p4 = p4*maxium

angles = open( angledoc, "r" )
line = angles.readline() # ignore the first line

iproj = 0
line = angles.readline()
while( line ) :
    print "proj # %4d" % iproj
    whole = base

    r1 = normalvariate(0.0, 1.0)
    r2 = normalvariate(0.0, 1.0)
    r3 = normalvariate(0.0, 1.0)
    r4 = normalvariate(0.0, 1.0)

    fr1 = r1
    fr2 = 1.2 * r2
    fr3 = 2.5 * (r1+r3)
    fr4 = (r1+2.0*r4)

    sp1 = p1 + (fr1/9.0*0.577) * p1*p1*p1
    sp2 = p2 + (fr2/4.0*0.555) * p2*p2*p2
    sp3 = p3 + (fr3/16.*0.5) * p3*p3*p3
    sp4 = p4 + (fr4/1.2*0.566) * p4*p4*p4

    whole = base + sp1 + sp2 + sp3 + sp4
    whole.write_image( volume_file, iproj )
        
    items = split( line )
    psi = atof( items[2] )
    theta = atof( items[3] )
    phi = atof( items[4] )
    volft,kb=prep_vol(whole)
    proj=prgs(volft, kb, [ phi,theta,psi, 0.0, 0.0])
    proj.set_attr_dict( {'phi':phi,'theta':theta,'psi':psi,'s2x':0.0,'s2y':0.0,'mirror':0.0} )
    proj.set_attr_dict( {'active':1, 'ctf_applied':0.0} )
    proj.write_image( proj_file, iproj )

    line = angles.readline()
    iproj = iproj + 1



