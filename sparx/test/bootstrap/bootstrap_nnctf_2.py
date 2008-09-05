#!/usr/bin/env python

from EMAN2  import *
from sparx  import *
from random import randint
from string import atoi
from string import atof
from mpi import *

import sys
import string

sys.argv = mpi_init(len(sys.argv), sys.argv)
size = mpi_comm_size(MPI_COMM_WORLD)
myid = mpi_comm_rank(MPI_COMM_WORLD)


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
while iarg < len(sys.argv) :
    tag = sys.argv[iarg]
    if( tag == "-h" ) :
        print "usage: bootstrap_nnctf [-h] -i proj_file -vp volume_prefix -sp status_prefix -n niter -media memory|filename -npad npad -snr signal_noise_raito -sign sign"
        exit(-1)
    elif( tag == "-i" ) :
        iarg=iarg+1
        proj_stack = sys.argv[iarg]
    elif( tag == "-vp" ) :
        iarg=iarg+1
        volume_prefix = sys.argv[iarg]
    elif( tag == "-sp" ) :
        iarg=iarg+1
        status_prefix = sys.argv[iarg]
    elif( tag == "-n" ) :
        iarg=iarg+1
        niter = atoi( sys.argv[iarg] )
    elif( tag == "-media" ) :
        iarg=iarg+1
        media = sys.argv[iarg]
    elif( tag == "-npad" ) :
        iarg=iarg+1
        npad = atoi( sys.argv[iarg] )
    elif( tag == "-snr" ) :
        iarg=iarg+1
        snr = atof( sys.argv[iarg] )
    elif( tag == "-sign" ) :
        iarg=iarg+1
	sign = atoi( sys.argv[iarg] )
    else :
        #print "unknown option ", tag, "use -h for help"
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
    #mymedia = "%s%4d.hdf" % ( media, myid )
    #mymedia = string.replace( mymedia, ' ', '0' )
    mymedia = media

nproj = EMUtil.get_image_count(proj_stack)
mystatus.write( "# of projs: %d\n" % nproj )

myseed = 1000 + myid

list_proj = range(nproj)
bootstrap_nnctf( proj_stack, myvolume_file, list_proj, myiter, mymedia, npad, "c1", mystatus, snr, sign, myseed)



