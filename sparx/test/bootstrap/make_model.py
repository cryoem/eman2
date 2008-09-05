from sparx import *
from EMAN2 import *
from sys   import exit
from sys   import argv
from string import split
from string import atof
from string import atoi
from random import normalvariate
from math   import exp

if( len( argv ) != 5 ) :
    print "Usage: model input angledoc output_volume output_proj"
    exit(-1)

input = open( argv[1], "r" )
angledoc = argv[2]
volume_file = argv[3]
proj_file = argv[4]

size = atoi( input.readline() )

items = split( input.readline() )
sigma_x = atoi( items[0] )
sigma_y = atoi( items[1] )
sigma_z = atoi( items[2] )
maxium  = atof( items[3] )
base = model_gauss( sigma_x, size,size,size, sigma_y,sigma_z )
[mean,sigma,min,max]=Util.infomask(base,None)
base = base * (maxium/max)

threshold = exp(-1.5)

line = input.readline()

points = []
varias = []

model = base
while line :
    items = split( line )
    xcenter = atoi( items[0] )
    ycenter = atoi( items[1] )
    zcenter = atoi( items[2] )

    xsigma  = atoi( items[3] )
    ysigma  = atoi( items[4] )
    zsigma  = atoi( items[5] )

    maxium  = atof( items[6] )
    varian  = atof( items[7] )

    print "make point at %4d, %4d, %4d" % ( xcenter, ycenter, zcenter )

    p1 = model_gauss( xsigma, size,size,size, ysigma, zsigma, xcenter,ycenter,zcenter )
    [mean,sigma,min,max]=Util.infomask(p1,None)
    p1 = p1/max
    p1.process_inplace( "eman1.threshold.belowtozero", {"minval": threshold} )
    p1 = p1*maxium

    points.append( p1 )
    varias.append( varian )
    model = model + p1

    line = input.readline()


dropImage( model, "model002.spi" )


angles = open( angledoc, "r" )
line = angles.readline() # ignore the first line

iproj = 0
line = angles.readline()
while( line ) :
    print "proj # %4d" % iproj
    whole = base
    for i in xrange( len(points) ) :
        ran = normalvariate(0.0, 1.0)
        p  = points[i] + points[i]*points[i]*points[i] *( ran*varias[i] )
	whole = whole + p
	whole.write_image( volume_file, iproj )
        
    items = split( line )
    psi = atof( items[2] )
    theta = atof( items[3] )
    phi = atof( items[4] )
    volft,kb=prep_vol(whole)
    proj=prgs(volft, kb, [ phi,theta,psi, 0.0, 0.0])
    proj.set_attr_dict( {'phi':phi,'theta':theta,'psi':psi,'sx':0.0,'sy':0.0,'mirror':0.0} )
    proj.set_attr_dict( {'active':1, 'ctf_applied':0.0} )
    proj.write_image( proj_file, iproj )

    line = angles.readline()
    iproj = iproj + 1



