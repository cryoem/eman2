#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from mpi import *
from string import replace
from random import random
from math import sqrt,pow

import sys

sys.argv = mpi_init( len(sys.argv), sys.argv )
nproc = mpi_comm_size( MPI_COMM_WORLD )
myid  = mpi_comm_rank( MPI_COMM_WORLD )

prj_stack = sys.argv[1]
outdir = sys.argv[2]
out_stack = "%s%4d.hdf" % (sys.argv[3], myid)
out_stack = replace( out_stack, ' ', '0' )

    
img1st = getImage( prj_stack )
nx = img1st.get_xsize()
ny = img1st.get_ysize()

img_number     = EMUtil.get_image_count( prj_stack )#len(img_list)
img_per_node   = img_number/nproc
img_node_start = img_per_node*myid

if myid == nproc-1: 
    img_node_end = img_number
else:               
    img_node_end = img_node_start + img_per_node

if myid==0:
    os.system( "rm -rf " + outdir )
    os.system( "mkdir " + outdir )
mpi_barrier( MPI_COMM_WORLD )

infofile = "progress%4d.txt" % (myid+1)
infofile = replace(infofile, ' ', '0')
finfo = open( infofile, 'w' )

imgdata = []
for i in xrange(img_node_start, img_node_end):
    img = get_im(prj_stack, i)
    imgdata.append(img)
finfo.write( ' all imgs loaded\n' )
finfo.flush( )

mask2d = model_circle( 25, nx, ny )
mask3d = None
niter = 10
snr = 50.0
sym = "c1"
sign = 1

fftvol_all,weight_all = prepare_recons_ctf_fftvol( imgdata, snr, sym, myid, 0, range(len(imgdata)) )
bcast_EMData_to_all( fftvol_all, myid )
bcast_EMData_to_all( weight_all, myid )

vol_all = recons_ctf_from_fftvol( nx, fftvol_all.copy(), weight_all.copy(), snr, sym, weighting=0 )
[vavg, vstd, vmin, vmax] = Util.infomask( vol_all, None, True )
finfo.write( "Info of reconstructed volume: %10.3e %10.3e %10.3e %10.3e\n" % (vavg, vstd, vmin, vmax) )
finfo.flush()

for i in xrange( len(imgdata) ) :

        curt_img = imgdata[i]

        fftvol_out = fftvol_all.copy()
        weight_out = weight_all.copy()
        params =  {"size":nx, "npad":4, "snr":snr, "sign":1, "symmetry":sym, "fftvol":fftvol_out, "weight":weight_out, "weighting":0}
        r = Reconstructors.get( "nn4_ctf", params )
        r.setup()

        curt_img.set_attr( "remove", 1 )
        Ttype = Transform3D.EulerType.SPIDER
        phi = curt_img.get_attr( 'phi' )
        theta = curt_img.get_attr( 'theta' )
        psi = curt_img.get_attr( 'psi' )    
        r.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
 
   
        vol_out = recons_ctf_from_fftvol(nx, fftvol_out.copy(), weight_out.copy(), snr, sym, weighting=0)

        cccoa = ccc( vol_out, vol_all )
        dotoo = vol_out.dot( vol_out )
        dotoa = vol_out.dot( vol_all )
        a = sqrt( dotoa/dotoo )


        add = (1.0-a)*100000.0
        if add > 0 :
            scale = 1.0 + add
        else :
            scale = 1.0 / (1.0 -add) 

        finfo.write( 'imgid, ccc, a, scale, round 1: %6d %10.8f %10.8f %10.8f\n' % ( img_node_start+i, cccoa, a, scale ) )
        finfo.flush()

        curt_img *= scale

        curt_img.set_attr( 'remove', 0 )
        r.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
        fftvol_all = fftvol_out
        weight_all = weight_out
        vol_all = recons_ctf_from_fftvol( nx, fftvol_all.copy(), weight_all.copy(), snr, sym, weighting=0)

        cccoa = ccc( vol_out, vol_all )
        dotoa = vol_out.dot( vol_all )
        a = sqrt( dotoa/dotoo )

        finfo.write( 'imgid, ccc, a,        round 2: %6d %10.8f %10.8f\n' % ( img_node_start+i, cccoa, a ) )
        finfo.flush()

        curt_img.write_image( out_stack, i  )

