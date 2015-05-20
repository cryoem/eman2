#!/usr/bin/env python
#
# Author: 
# Copyright (c) 2012 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#


from __future__ import print_function
from EMAN2 import *
from sparx import *
from logger import Logger, BaseLogger_Files
import global_def

from mpi   import  *
from math  import  *



import os
import sys
import subprocess
import time
import string
from   sys import exit
from   time import localtime, strftime

def subdict(d,u):
	# substitute values in dictionary d by those given by dictionary u
	for q in u:  d[q] = u[q]

def cmdexecute(cmd):
	from   time import localtime, strftime
	import subprocess
	outcome = subprocess.call(cmd, shell=True)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(outcome == 1):
		print(  line,"ERROR!!   Command failed:  ", cmd)
		exit()
	else:  print(line,"Executed successfully: ",cmd)

def fuselowf(vs, fq):
	n = len(vs)
	for i in xrange(n): fftip(vs[i])
	a = vs[0].copy()
	for i in xrange(1,n):
		Util.add_img(a, vs[i])
	Util.mul_scalar(a, 1.0/float(n))
	a = filt_tophatl(a, fq)
	for i in xrange(n):
		vs[i] = fft(Util.addn_img(a, filt_tophath(vs[i], fq)))
	return

global passlastring, mempernode
def hlfmem(x):
	return (len(even_angles(x))*4.*passlastring**2*4./ mempernode - 0.5)**2


def comparetwoalis(params1, params2, thresherr=1.0, radius = 1.0):
	#  Find errors per image
	nn = len(params1)
	perr = 0
	for k in xrange(nn):
		if(max_3D_pixel_error(params1[k], params2[k], r=radius) < thresherr):
			perr += 1
	return perr/float(nn)*100.0


def checkstep(item, keepchecking, myid, main_node):
	if(myid == main_node):
		if keepchecking:
			if(os.path.exists(item)):
				doit = 0
			else:
				doit = 1
				keepchecking = False
		else:
			doit = 1
	else:
		doit = 1
	doit = bcast_number_to_all(doit, source_node = main_node)
	return doit, keepchecking


def getindexdata(stack, partids, partstack, myid, nproc):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	lpartids  = map(int, read_text_file(partids) )
	ndata = len(lpartids)
	partstack = read_text_row(partstack)
	if( ndata < nproc):
		if(myid<ndata):
			image_start = myid
			image_end   = myid+1
		else:
			image_start = 0
			image_end   = 1			
	else:
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
	lpartids  = lpartids[image_start:image_end]
	partstack = partstack[image_start:image_end]
	data = EMData.read_images(stack, lpartids)
	for i in xrange(len(partstack)):  set_params_proj(data[i], partstack[i])
	return data


def getalldata(stack, myid, nproc):
	if(myid == 0):  ndata = EMUtil.get_image_count(stack)
	else:           ndata = 0
	ndata = bcast_number_to_all(ndata)	
	if( ndata < nproc):
		if(myid<ndata):
			image_start = myid
			image_end   = myid+1
		else:
			image_start = 0
			image_end   = 1			
	else:
		image_start, image_end = MPI_start_end(ndata, nproc, myid)
	data = EMData.read_images(stack, range(image_start, image_end))
	return data



class ali3d_options:
	ir     = 1
	rs     = 1
	ou     = -1
	xr     = "-1"
	yr     = "-1"
	ts     = "1"
	an     = "-1"
	sym    = "d2"
	delta  = "2"
	npad   = 2
	center = 0
	CTF    = True
	ref_a  = "S"
	snr    = 1.0
	mask3D = "startm.hdf"
	fl     = 0.4
	aa     = 0.1
	initfl = 0.4
	pwreference = "rotpw3i3.txt"
#################################



def run3Dalignment(paramsdict, partids, partstack, outputdir, procid, myid, main_node, nproc):
	#  Reads from paramsdict["stack"] particles partids set parameters in partstack
	#    and do refinement as specified in paramsdict
	#
	#  Will create outputdir
	#  Will write to outputdir output parameters: params-chunk0.txt and params-chunk1.txt
	if(myid == main_node):
		#  Create output directory
		log = Logger(BaseLogger_Files())
		log.prefix = os.path.join(outputdir)
		#cmd = "mkdir "+log.prefix
		#cmdexecute(cmd)
		log.prefix += "/"
	else:  log = None
	mpi_barrier(MPI_COMM_WORLD)

	ali3d_options.delta  = paramsdict["delta"]
	ali3d_options.ts     = paramsdict["ts"]
	ali3d_options.xr     = paramsdict["xr"]
	#  low pass filter is applied to shrank data, so it has to be adjusted
	ali3d_options.fl     = paramsdict["lowpass"]/paramsdict["shrink"]
	ali3d_options.initfl = paramsdict["initialfl"]/paramsdict["shrink"]
	ali3d_options.aa     = paramsdict["falloff"]
	ali3d_options.maxit  = paramsdict["maxit"]
	ali3d_options.mask3D = paramsdict["mask3D"]
	ali3d_options.an	 = paramsdict["an"]
	ali3d_options.ou     = paramsdict["radius"]  #  This is changed in ali3d_base, but the shrank value is needed in vol recons, fixt it!
	shrinkage            = paramsdict["shrink"]

	projdata = getindexdata(paramsdict["stack"], partids, partstack, myid, nproc)
	onx = projdata[0].get_xsize()
	last_ring = ali3d_options.ou
	if last_ring < 0:	last_ring = int(onx/2) - 2
	mask2D  = model_circle(last_ring,onx,onx) - model_circle(ali3d_options.ir,onx,onx)
	if(shrinkage < 1.0):
		# get the new size
		masks2D = resample(mask2D, shrinkage)
		nx = masks2D.get_xsize()
		masks2D  = model_circle(int(last_ring*shrinkage+0.5),nx,nx) - model_circle(max(int(ali3d_options.ir*shrinkage+0.5),1),nx,nx)
	nima = len(projdata)
	oldshifts = [0.0,0.0]*nima
	for im in xrange(nima):
		#data[im].set_attr('ID', list_of_particles[im])
		ctf_applied = projdata[im].get_attr_default('ctf_applied', 0)
		phi,theta,psi,sx,sy = get_params_proj(projdata[im])
		projdata[im] = fshift(projdata[im], sx, sy)
		set_params_proj(projdata[im],[phi,theta,psi,0.0,0.0])
		#  For local SHC set anchor
		#if(nsoft == 1 and an[0] > -1):
		#	set_params_proj(data[im],[phi,tetha,psi,0.0,0.0], "xform.anchor")
		oldshifts[im] = [sx,sy]
		if ali3d_options.CTF :
			ctf_params = projdata[im].get_attr("ctf")
			if ctf_applied == 0:
				st = Util.infomask(projdata[im], mask2D, False)
				projdata[im] -= st[0]
				projdata[im] = filt_ctf(projdata[im], ctf_params)
				projdata[im].set_attr('ctf_applied', 1)
		if(shrinkage < 1.0):
			#phi,theta,psi,sx,sy = get_params_proj(projdata[im])
			projdata[im] = resample(projdata[im], shrinkage)
			st = Util.infomask(projdata[im], None, True)
			projdata[im] -= st[0]
			st = Util.infomask(projdata[im], masks2D, True)
			projdata[im] /= st[1]
			#sx *= shrinkage
			#sy *= shrinkage
			#set_params_proj(projdata[im], [phi,theta,psi,sx,sy])
			if ali3d_options.CTF :
				ctf_params.apix /= shrinkage
				projdata[im].set_attr('ctf', ctf_params)
		else:
			st = Util.infomask(projdata[im], None, True)
			projdata[im] -= st[0]
			st = Util.infomask(projdata[im], mask2D, True)
			projdata[im] /= st[1]
	del mask2D
	if(shrinkage < 1.0): del masks2D

	"""
	if(paramsdict["delpreviousmax"]):
		for i in xrange(len(projdata)):
			try:  projdata[i].del_attr("previousmax")
			except:  pass
	"""
	if(myid == main_node):
		print_dict(paramsdict,"3D alignment parameters")
		print("                    =>  actual lowpass      :  ",ali3d_options.fl)
		print("                    =>  actual init lowpass :  ",ali3d_options.initfl)
		if(len(ali3d_options.pwreference)>0): \
		print("                    =>  PW adjustment       :  ",ali3d_options.pwreference)
		print("                    =>  partids             :  ",partids)
		print("                    =>  partstack           :  ",partstack)
		
	if(ali3d_options.fl > 0.46):  ERROR("Low pass filter in 3D alignment > 0.46 on the scale of shrank data","sxcenter_projections",1,myid) 

	#  Run alignment command, it returns params per CPU
	params = center_projections_3D(projdata, paramsdict["refvol"], \
									ali3d_options, onx, shrinkage, \
									mpi_comm = MPI_COMM_WORLD,  myid = myid, main_node = main_node, log = log )
	del log, projdata

	params = wrap_mpi_gatherv(params, main_node, MPI_COMM_WORLD)

	#  store params
	if(myid == main_node):
		for im in xrange(nima):
			params[im][0] = params[im][0]/shrinkage +oldshifts[im][0]
			params[im][1] = params[im][1]/shrinkage +oldshifts[im][1]
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"Executed successfully: ","3D alignment","  number of images:%7d"%len(params))
		write_text_row(params, os.path.join(outputdir,"params.txt") )

def print_dict(dict,theme):
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	print(line,theme)
	spaces = "                    "
	for q in dict:  print("                    => ",q+spaces[len(q):],":  ",dict[q])


def main():

	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	import user_functions
	from applications import MPI_start_end
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc, do_volume
	from logger import Logger, BaseLogger_Files
	import sys
	import os
	import time
	import socket

	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  [output_directory]  initial_volume  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --CTF  --fl --aa --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",      		type= "int",   default= 1,			help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",      		type= "int",   default= -1,			help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",      		type= "int",   default= 1,			help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",      		type="string", default= "-1",		help="range for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",      		type="string", default= "-1",		help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",      		type="string", default= "1",		help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",   		type="string", default= "-1",		help="angular step of reference projections during initialization step (default automatically selected based on radius of the structure.)")
	parser.add_option("--an",      		type="string", default= "-1",		help="angular neighborhood for local searches (phi and theta) (Default exhaustive searches)")
	parser.add_option("--CTF",     		action="store_true", default=False,	help="Use CTF (Default no CTF correction)")
	parser.add_option("--shrink",     	type="float",  default= 1.0,		help="Reduce data size by shrink factor (default 1.0)")
	parser.add_option("--snr",     		type="float",  default= 1.0,		help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",   		type="string", default= "S",		help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",     		type="string", default= "c1",		help="symmetry of the refined structure")
	parser.add_option("--npad",    		type="int",    default= 2,			help="padding size for 3D reconstruction (default=2)")

	#options introduced for the do_volume function
	parser.add_option("--fl",			type="float",	default=0.12,		help="cut-off frequency of hyperbolic tangent low-pass Fourier filte (default 0.12)")
	parser.add_option("--aa",			type="float",	default=0.1,		help="fall-off of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--pwreference",	type="string",	default="",			help="text file with a reference power spectrum (default no power spectrum adjustment)")
	parser.add_option("--mask3D",		type="string",	default=None,		help="3D mask file (default a sphere  WHAT RADIUS??)")


	(options, args) = parser.parse_args(sys.argv[1:])

	#print( "  args  ",args)
	if( len(args) == 3):
		volinit = args[2]
		masterdir = args[1]
	elif(len(args) == 2):
		volinit = args[1]
		masterdir = ""
	else:
		print( "usage: " + usage)
		print( "Please run '" + progname + " -h' for detailed options")
		return 1

	stack = args[0]

	#  INPUT PARAMETERS
	radi  = options.ou
	global_def.BATCH = True
	ali3d_options.ir     = options.ir
	ali3d_options.rs     = options.rs
	ali3d_options.ou     = options.ou
	ali3d_options.xr     = options.xr
	ali3d_options.yr     = options.yr
	ali3d_options.ts     = options.ts
	ali3d_options.an     = "-1"
	ali3d_options.sym    = options.sym
	ali3d_options.delta  = options.delta
	ali3d_options.npad   = options.npad
	ali3d_options.CTF    = options.CTF
	ali3d_options.ref_a  = options.ref_a
	ali3d_options.snr    = options.snr
	ali3d_options.mask3D = options.mask3D
	ali3d_options.pwreference = ""  #   It will have to be turned on after exhaustive done by setting to options.pwreference
	ali3d_options.fl     = 0.4
	ali3d_options.initfl = 0.4
	ali3d_options.aa     = 0.1

	mpi_init(0, [])



	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	# Get the pixel size, if none set to 1.0, and the original image size
	if(myid == main_node):
		total_stack = EMUtil.get_image_count(stack)
		a = get_im(stack)
		nxinit = a.get_xsize()
		if ali3d_options.CTF:
			i = a.get_attr('ctf')
			pixel_size = i.apix
			fq = pixel_size/fq
		else:
			pixel_size = 1.0
			#  No pixel size, fusing computed as 5 Fourier pixels
			fq = 5.0/nxinit
		del a
	else:
		total_stack = 0
		nxinit = 0
		pixel_size = 1.0
	total_stack = bcast_number_to_all(total_stack, source_node = main_node)
	pixel_size  = bcast_number_to_all(pixel_size, source_node = main_node)
	nxinit      = bcast_number_to_all(nxinit, source_node = main_node)

	if(radi < 1):  radi = nxinit//2-2
	elif((2*radi+2)>nxinit):  ERROR("Particle radius set too large!","sxcenter_projections",1,myid)
	ali3d_options.ou = radi

	shrink = options.shrink
	nxshrink = int(nxinit*shrink+0.5)
	angular_neighborhood = "-1"

	#  MASTER DIRECTORY
	if(myid == main_node):
		print( "   masterdir   ",masterdir)
		if( masterdir == ""):
			timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
			masterdir = "master"+timestring
		li = len(masterdir)
		cmd = "{} {}".format("mkdir", masterdir)
		cmdexecute(cmd)
	else:
		li = 0

	li = mpi_bcast(li,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]

	if( li > 0 ):
		masterdir = mpi_bcast(masterdir,li,MPI_CHAR,main_node,MPI_COMM_WORLD)
		masterdir = string.join(masterdir,"")


	nnxo        = nxinit

	#  INITIALIZATION

	initdir = masterdir

	#  This is initial setting, has to be initialized here, we do not want it to run too long.
	#    INITIALIZATION THAT FOLLOWS WILL HAVE TO BE CHANGED SO THE USER CAN PROVIDE INITIAL GUESS OF RESOLUTION
	#  If we new the initial resolution, it could be done more densely
	if(options.xr == "-1"):  xr = "%d"%((nnxo - (2*radi-1))//2)
	else:  xr = options.xr
	if(options.yr == "-1"):  yr = xr
	else:  yr = options.yr

	delta = float(options.delta)
	if(delta <= 0.0):  delta = "%f"%round(degrees(atan(1.0/float(radi))), 2)
	else:    delta = "%f"%delta

	paramsdict = {	"stack":stack,"delta":delta, "ts":"1.0", "xr":xr, "an":angular_neighborhood, \
					"center":"0", "maxit":1, "local":False,\
					"lowpass":options.fl, "initialfl":0.4, "falloff":options.aa, "radius":radi, \
					"nsoft":0, "delpreviousmax":True, "shrink":options.shrink, "saturatecrit":1.0, "pixercutoff":2.0,\
					"refvol":volinit, "mask3D":options.mask3D}

	partids = os.path.join(masterdir, "ids.txt")
	partstack = os.path.join(masterdir, "paramszero.txt")


	if( myid == main_node ):
		write_text_file(range(total_stack), partids)
		write_text_row([[0.0,0.0,0.0,0.0,0.0] for i in xrange(total_stack) ], partstack)

	run3Dalignment(paramsdict, partids, partstack, initdir, 0, myid, main_node, nproc)

	mpi_barrier(MPI_COMM_WORLD)

	mpi_finalize()


if __name__=="__main__":
	main()

