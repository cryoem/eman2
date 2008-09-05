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

finfo.write( ' image range: %6d - %6d\n' %(img_node_start, img_node_end) )
finfo.flush()

imgdata = []
for i in xrange(img_node_start, img_node_end):
    img = get_im(prj_stack, i)
    dec = decimate( img, 2 )
    imgdata.append(img)

finfo.write( ' all imgs loaded\n' )
finfo.flush( )

mask2d = model_circle( 25, nx, ny )
mask3d = None
niter = 10
snr = 1.0
sym = "c1"
sign = 1



fftvol_odd,weight_odd = prepare_recons_ctf_fftvol( imgdata, snr, sym, myid, 0, range(0, len(imgdata), 2), finfo )
fftvol_eve,weight_eve = prepare_recons_ctf_fftvol( imgdata, snr, sym, myid, 0, range(1, len(imgdata), 2), finfo )

finfo.write( "bcasting fftvol\n" )
finfo.flush()

bcast_EMData_to_all( fftvol_odd, myid )
bcast_EMData_to_all( weight_odd, myid )
bcast_EMData_to_all( fftvol_eve, myid )
bcast_EMData_to_all( weight_eve, myid )

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



fscfile = outdir + ( "/fsc%6d.dat" % img_node_start )
fscfile = replace(fscfile, ' ', '0')

fscc = fsc( vol_odd, vol_eve, filename=fscfile )


newscale = [] 

for i in xrange( img_node_start, img_node_end ) :
    curt_img = imgdata[i-img_node_start]

    if i%2==0 :
        fftvol_half_new = fftvol_odd.copy()
        weight_half_new = weight_odd.copy()
    else:
        fftvol_half_new = fftvol_eve.copy()
        weight_half_new = weight_eve.copy()

    params =  {"size":nx, "npad":4, "snr":snr, "sign":1, "symmetry":sym, "fftvol":fftvol_half_new, "weight":weight_half_new, "weighting":0}
    rhalf = Reconstructors.get( "nn4_ctf", params )
    rhalf.setup()


    Ttype = Transform3D.EulerType.SPIDER
    phi = curt_img.get_attr( 'phi' )
    theta = curt_img.get_attr( 'theta' )
    psi = curt_img.get_attr( 'psi' )    

    curt_img.set_attr( "remove", 1 )
    rhalf.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
 
        
    curt_img *= 0.95
    curt_img.set_attr( "remove", 0 )
    rhalf.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
  
    vol_half_new = recons_ctf_from_fftvol(nx, fftvol_half_new.copy(), weight_half_new.copy(), snr, sym, weighting=0)

    if i%100==0:
        fscfile = outdir + ( "/fsc%6d.dat" % (i+1) )
        fscfile = replace(fscfile, ' ', '0')
    else:
        fscfile = None

    if i%2==0:
        fscc_new = fsc( vol_half_new, vol_eve, filename=fscfile )
    else:
        fscc_new = fsc( vol_half_new, vol_odd, filename=fscfile )

    scale = 0.95      
    finfo.write( 'round 1 imgid scale, fsc_before, fsc_after: %6d %11.8f %11.8f %11.8f\n' % (i, scale, sum(fscc[1]), sum(fscc_new[1])) )
    finfo.flush()
       
        
    if sum(fscc_new[1]) < sum(fscc[1]):
        curt_img.set_attr( "remove", 1 )
        rhalf.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
           
        curt_img *= 1.1 
        curt_img.set_attr( "remove", 0 )
        rhalf.insert_slice( curt_img, Transform3D(Ttype, phi, theta, psi) )
  
        vol_half_new = recons_ctf_from_fftvol(nx, fftvol_half_new.copy(), weight_half_new.copy(), snr, sym, weighting=0)

        if i%2==0:
            fscc_new = fsc( vol_half_new, vol_eve, filename=fscfile )
        else:
            fscc_new = fsc( vol_half_new, vol_odd, filename=fscfile )

        scale = 1.05
        finfo.write( 'round 2 imgid scale, fsc_before, fsc_after: %6d %11.8f %11.8f %11.8f\n' % (i, scale, sum(fscc[1]), sum(fscc_new[1])) )
        finfo.flush()
 
    if sum(fscc_new[1]) < sum(fscc[1]):
        scale = 1.0
        finfo.write( 'img %d already normalized \n' % i) 
        finfo.flush()
    else:
        if i%2==0:
            fftvol_odd = fftvol_half_new
            weight_odd = weight_half_new
            vol_odd    = vol_half_new
        else:
            fftvol_eve = fftvol_half_new
            weight_eve = weight_half_new
            vol_eve    = vol_half_new

        fscc       = fscc_new
        
    newscale.append(scale)


newscalefile = outdir + ( "/newscale%4d.txt" % myid )
newscalefile = replace( newscalefile, ' ', '0' )
dropSpiderDoc( newscalefile, newscale )

