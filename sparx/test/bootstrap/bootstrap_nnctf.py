#!/usr/bin/env python

from EMAN2  import *
from sparx  import *
from random import randint
from string import atoi
from string import atof

import sys
import string
import pypar

size = pypar.size()
myid = pypar.rank()

if myid == 0 :
    argv = sys.argv
    for i in xrange(1,size) :
       pypar.send( argv, destination=i)
else:
    argv = pypar.receive( source=0 )



proj_stack = None
volume_prefix = None
status_prefix = None
niter = None
media = None
npad = 4
ctf = 'NONE'
snr = None
sign = 1

iarg = 1
while iarg < len(argv) :
    tag = argv[iarg]
    if( tag == "-h" ) :
        print "usage: bootstrap_nn [-h] -i proj_file -vp volume_prefix -sp status_prefix -n niter -media memory|filename -npad npad -snr signal_noise_raito -sign sign"
        exit(-1)
    elif( tag == "-i" ) :
        iarg=iarg+1
        proj_stack = argv[iarg]
    elif( tag == "-vp" ) :
        iarg=iarg+1
        volume_prefix = argv[iarg]
    elif( tag == "-sp" ) :
        iarg=iarg+1
        status_prefix = argv[iarg]
    elif( tag == "-n" ) :
        iarg=iarg+1
        niter = atoi( argv[iarg] )
    elif( tag == "-media" ) :
        iarg=iarg+1
        media = argv[iarg]
    elif( tag == "-npad" ) :
        iarg=iarg+1
        npad = atoi( argv[iarg] )
    elif( tag == "-snr" ) :
        iarg=iarg+1
        snr = atof( argv[iarg] )
    elif( tag == "-sign" ) :
        iarg=iarg+1
	sign = atoi( argv[iarg] )
    else :
       # print "unknown option ", tag, "use -h for help"
       pass

    iarg=iarg+1

if( proj_stack == None ) :
    print "no input projection specified"
    exit(-1)

if( volume_prefix == None ) :
    print "volume file name not specified"
    exit(-1)

if( status_prefix == None ) :
    print "status file name not specified"
    exit(-1)

if( media == None ) :
    print "media type not specified"
    exit(-1)

if( snr == None ) :
    print "signal noise raito not specified for ctf"
    exit(-1)

myvolume_file = "%s%4d.hdf" % ( volume_prefix, myid )
mystatus_file = "%s%4d.inf" % ( status_prefix, myid )

myvolume_file = string.replace( myvolume_file, ' ', '0' )
mystatus_file = string.replace( mystatus_file, ' ', '0' )

mystatus = open( mystatus_file, 'w' )

myiter = niter / size

if( media == "memory" ) :
    mymedia = media
else :
    mymedia = "%s%4d.hdf" % (media,myid)
    mymedia = string.replace(mymedia, ' ', '0')

nproj = EMUtil.get_image_count(proj_stack)
mystatus.write( "# of projs: %d\n" % nproj )

myseed = 1000+myid

list_proj = range(nproj)
recons3d_bootstrap( proj_stack, myvolume_file, list_proj, myiter, mymedia, npad, "c1", mystatus, "yes", snr, sign, myseed )



