#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from mpi import *
from string import replace
from random import random
from math import sqrt,pow

import sys

def decimate(img, ndeci):
    btwl = filt_btwl( img, 0.23, 0.27 )
    deci = Util.decimate( btwl, ndeci, ndeci, 1 )

    phi = img.get_attr( 'phi' )
    theta = img.get_attr( 'theta' )
    psi = img.get_attr( 'psi' )
    s2x = img.get_attr( 's2x' ) / ndeci
    s2y = img.get_attr( 's2y' ) / ndeci

    Cs = img.get_attr( 'Cs' )
    ctfed = img.get_attr( 'ctf_applied' )
    pixel = img.get_attr( 'Pixel_size' ) * ndeci
    defocus = img.get_attr( 'defocus' )
    voltage = img.get_attr( 'voltage' )
    amp_contrast = img.get_attr( 'amp_contrast' )

    deci.set_attr_dict( {'phi':phi, 'theta':theta, 'psi':psi, 's2x':s2x, 's2y':s2y, 'mirror':0.0} )
    deci.set_attr_dict( {'active':1, 'ctf_applied':ctfed} )
    deci.set_attr_dict( {'defocus':defocus, 'amp_contrast':amp_contrast, 'voltage':voltage, 'Cs':Cs} )
    deci.set_attr_dict( {'Pixel_size':pixel, 'B_factor':0.0} )

    return deci    



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

img_number     = EMUtil.get_image_count( prj_stack )
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

finfo.write( ' image range: %6d - %6d\n' %(img_node_start, img_node_end) )
finfo.flush()

imgdata = []
for i in xrange(img_node_start, img_node_end):
    img = get_im(prj_stack, i)
    imgdata.append(img)

totscale = [1.0] * len(imgdata)

finfo.write( ' all imgs loaded\n' )
finfo.flush( )

mask2d = model_circle( 25, nx, ny )
mask3d = None
niter = 10
snr = 1.0
sym = "c1"
sign = 1
niter = 10

for iter in xrange(niter):
    finfo.write( "Iteration #%2d\n" % iter )
    finfo.flush()

    odd_start = img_node_start%2
    eve_start = 1 - odd_start

    fftvol_odd,weight_odd = prepare_recons_ctf_fftvol( imgdata, snr, sym, myid, 0, range(odd_start, len(imgdata), 2) )
    fftvol_eve,weight_eve = prepare_recons_ctf_fftvol( imgdata, snr, sym, myid, 0, range(eve_start, len(imgdata), 2) )
    fftvol_all = fftvol_odd.copy()
    weight_all = weight_odd.copy()
    fftvol_all += fftvol_eve
    weight_all += weight_eve


    finfo.write( "bcasting fftvol\n" )
    finfo.flush()

    bcast_EMData_to_all( fftvol_odd, myid )
    bcast_EMData_to_all( weight_odd, myid )
    bcast_EMData_to_all( fftvol_eve, myid )
    bcast_EMData_to_all( weight_eve, myid )
    bcast_EMData_to_all( fftvol_all, myid )
    bcast_EMData_to_all( weight_all, myid )


    finfo.write( "bcasting fftvol done\n" )
    finfo.flush()


    vol_odd = recons_ctf_from_fftvol( nx, fftvol_odd.copy(), weight_odd.copy(), snr, sym, weighting=0 ) 
    [vavg, vstd, vmin, vmax] = Util.infomask( vol_odd, None, True )
    finfo.write( "Info of odd volume: %10.3e %10.3e %10.3e %10.3e\n" % (vavg, vstd, vmin, vmax) )
    finfo.flush()

    vol_eve = recons_ctf_from_fftvol( nx, fftvol_eve.copy(), weight_eve.copy(), snr, sym, weighting=0 )
    [vavg, vstd, vmin, vmax] = Util.infomask( vol_eve, None, True )
    finfo.write( "Info of eve volume: %10.3e %10.3e %10.3e %10.3e\n" % (vavg, vstd, vmin, vmax) )
    finfo.flush()

    vol_all = recons_ctf_from_fftvol( nx, fftvol_all.copy(), weight_all.copy(), snr, sym, weighting=0 )
    [vavg, vstd, vmin, vmax] = Util.infomask( vol_eve, None, True )
    finfo.write( "Info of all volume: %10.3e %10.3e %10.3e %10.3e\n" % (vavg, vstd, vmin, vmax) )
    finfo.flush()


    fscfile = outdir + ( "/fsc%2d_%6d.dat" % (iter,img_node_start) )
    fscfile = replace(fscfile, ' ', '0')

    fscc = fsc( vol_odd, vol_eve, filename=fscfile )


    newscale = [] 

    for i in xrange( img_node_start, img_node_end ) :
        curt_img = imgdata[i-img_node_start]

        fftvol_all_new = fftvol_all.copy()
        weight_all_new = weight_all.copy()
        allparams =  {"size":nx, "npad":4, "snr":snr, "sign":1, "symmetry":sym, "fftvol":fftvol_all_new, "weight":weight_all_new}
        rall = Reconstructors.get( "nn4_ctf", allparams )
        rall.setup()

        if i%2==0 :
            fftvol_hlf_new = fftvol_odd.copy()
            weight_hlf_new = weight_odd.copy()
        else:
            fftvol_hlf_new = fftvol_eve.copy()
            weight_hlf_new = weight_eve.copy()

        hlfparams =  {"size":nx, "npad":4, "snr":snr, "sign":1, "symmetry":sym, "fftvol":fftvol_hlf_new, "weight":weight_hlf_new}
        rhlf = Reconstructors.get( "nn4_ctf", hlfparams )
        rhlf.setup()

        Ttype = Transform3D.EulerType.SPIDER
        phi = curt_img.get_attr( 'phi' )
        theta = curt_img.get_attr( 'theta' )
        psi = curt_img.get_attr( 'psi' )    

        curt_img.set_attr( "remove", 1 )
        rall.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
        rhlf.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
 
        vol_all_new = recons_ctf_from_fftvol(nx, fftvol_all_new.copy(), weight_all_new.copy(), snr, sym, weighting=0 )
        dot_new_new = vol_all_new.dot( vol_all_new )
        dot_new_old = vol_all_new.dot( vol_all )
        a = sqrt( dot_new_new/dot_new_old )

        add = (1.0-a)*50000.0
        if add > 0:
            scale = 1.0 + add
        else:
            scale = 1.0/ (1.0-add)

        finfo.write( 'imgid, a, scale, round 1: %6d %10.8f %10.8f\n' % ( i, a, scale ) )
        finfo.flush()

        curt_img *= scale
        curt_img.set_attr( "remove", 0 )
        rall.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
        rhlf.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
  
        vol_all_new = recons_ctf_from_fftvol(nx, fftvol_all_new.copy(), weight_all_new.copy(), snr, sym, weighting=0)
        dot_new_new = vol_all_new.dot( vol_all_new )
        dot_new_old = vol_all_new.dot( vol_all )
        a = sqrt( dot_new_new/dot_new_old )


        vol_hlf_new = recons_ctf_from_fftvol(nx, fftvol_hlf_new.copy(), weight_hlf_new.copy(), snr, sym, weighting=0)
        if i%1000==0:
            fscfile = outdir + ( "/fsc%2d_%6d.dat" % (iter, i+1) )
            fscfile = replace(fscfile, ' ', '0')
        else:
            fscfile = None

        if i%2==0:
            fscc_new = fsc( vol_hlf_new, vol_eve, filename=fscfile )
        else:
            fscc_new = fsc( vol_hlf_new, vol_odd, filename=fscfile )

        finfo.write( 'imgid, a, scale, fsc    : %6d %10.8f %10.8f %10.8f %10.8f\n' % ( i, a, scale, sum(fscc[1]), sum(fscc_new[1]) ) )
        finfo.flush()

        if sum(fscc_new[1]) < sum(fscc[1]):
            curt_img /= scale
            scale = 1.0
            finfo.write( 'img %d remain unchanged \n' % i) 
            finfo.flush()
        else:
            if i%2==0:
                fftvol_odd = fftvol_hlf_new
                weight_odd = weight_hlf_new
                vol_odd    = vol_hlf_new
            else:
                fftvol_eve = fftvol_hlf_new
                weight_eve = weight_hlf_new
                vol_eve    = vol_hlf_new
 
            finfo.write( 'img %d new a: %10.8f\n' %(i, a)  )
 
        
            fftvol_all = fftvol_all_new
            weight_all = weight_all_new
            vol_all    = vol_all_new
            fscc       = fscc_new
        
        newscale.append(scale)

    assert len(totscale)==len(newscale)
    for ii in xrange( len(newscale) ):
        totscale[ii] *= newscale[ii]
        curt_img.write_image( out_stack, ii )

    newscalefile = outdir + ( "/newscale%2d_%4d.txt" % (iter,myid) )
    newscalefile = replace( newscalefile, ' ', '0' )
    dropSpiderDoc( newscalefile, totscale )

