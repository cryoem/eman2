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



def cmdexecute(cmd):
	from   time import localtime, strftime
	import subprocess
	outcome = subprocess.call(cmd, shell=True)
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	if(outcome == 1):
		print(  line,"ERROR!!   Command failed:  ", cmd)
		exit()
	else:  print(line,"Executed successfully: ",cmd)

def volshiftali(vv, mask3d=None):
	nv = len(vv)
	ni = vv[0].get_xsize()
	if(mask3d == None):
		mask3d = model_circle(ni//2-2, ni, ni, ni)

	for i in xrange(nv):
		Util.mul_img(vv[i], mask3d)
		fftip(vv[i])

	del mask3d
	ps = [[0.,0.,0.] for i in xrange(nv)]
	changed = True
	while(changed):
		a = EMData(ni, ni, ni, False)
		for i in vv:
			Util.add_img(a,i)

		changed = False
		for i in xrange(nv):
			pp = peak_search(ccf(vv[i],a), 1)
			for j in xrange(3):
				if(pp[0][5+j] != ps[i][j]):
					ps[i][j] = pp[0][5+j]
					changed = True

	return ps

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

def doXfiles(path, source = "chunk", inparams = "params", params = "params", dest = "X"):
	#  will produce X*.txt and paramsX*.txt
	#  Generate six Xfiles from four chunks and export parameters.  This is hardwired as it is always done in the same way
	#  AB
	#     indices
	write_text_file( \
		map(int, read_text_file(os.path.join(path,source+"0.txt")))+map(int, read_text_file(os.path.join(path,source+"1.txt"))),   \
		os.path.join(path,dest+"0.txt"))
	#  params
	write_text_row( \
		read_text_row(os.path.join(path,inparams+"00.txt"))+read_text_row(os.path.join(path,inparams+"10.txt")), \
		os.path.join(path,params+dest+"0.txt"))
	#  AC
	write_text_file( \
		map(int, read_text_file(os.path.join(path,source+"0.txt")))+map(int, read_text_file(os.path.join(path,source+"2.txt"))),   \
		os.path.join(path,dest+"1.txt"))
	write_text_row( \
		read_text_row(os.path.join(path,inparams+"01.txt"))+read_text_row(os.path.join(path,inparams+"20.txt")), \
		os.path.join(path,params+dest+"1.txt"))
	#  AD
	write_text_file( \
		map(int, read_text_file(os.path.join(path,source+"0.txt")))+map(int, read_text_file(os.path.join(path,source+"3.txt"))),   \
		os.path.join(path,dest+"2.txt"))
	write_text_row( \
		read_text_row(os.path.join(path,inparams+"02.txt"))+read_text_row(os.path.join(path,inparams+"30.txt")), \
		os.path.join(path,params+dest+"2.txt"))
	#  BC
	write_text_file( \
		map(int, read_text_file(os.path.join(path,source+"1.txt")))+map(int, read_text_file(os.path.join(path,source+"2.txt"))),   \
		os.path.join(path,dest+"3.txt"))
	write_text_row( \
		read_text_row(os.path.join(path,inparams+"11.txt"))+read_text_row(os.path.join(path,inparams+"21.txt")), \
		os.path.join(path,params+dest+"3.txt"))
	#  BD
	write_text_file( \
		map(int, read_text_file(os.path.join(path,source+"1.txt")))+map(int, read_text_file(os.path.join(path,source+"3.txt"))),   \
		os.path.join(path,dest+"4.txt"))
	write_text_row( \
		read_text_row(os.path.join(path,inparams+"12.txt"))+read_text_row(os.path.join(path,inparams+"31.txt")), \
		os.path.join(path,params+dest+"4.txt"))
	#  CD
	write_text_file( \
		map(int, read_text_file(os.path.join(path,source+"2.txt")))+map(int, read_text_file(os.path.join(path,source+"3.txt"))),   \
		os.path.join(path,dest+"5.txt"))
	write_text_row( \
		read_text_row(os.path.join(path,inparams+"22.txt"))+read_text_row(os.path.join(path,inparams+"32.txt")), \
		os.path.join(path,params+dest+"5.txt"))
	return

def	mergeparfiles(i1,i2,io,p1,p2,po):
	#  1 - rescued
	#  2 - good old
	l1 = map(int, read_text_file( i1 ))
	l2 = map(int, read_text_file( i2 ))
	if(l1[0] == 0):
			write_text_file( l2, io)
			for ll in xrange(3):
				p = read_text_row(p2+"%01d.txt"%ll)
				write_text_row(2, po+"%01d.txt"%ll)
	else:
		t = l1 +l2
		for i in xrange(len(t)):
			t[i] = [t[i],i]
		t.sort()
		write_text_file( [t[i][0] for i in xrange(len(t))], io)
		for ll in xrange(3):
			p = read_text_row(p1+"%01d.txt"%ll) + read_text_row(p2+"%01d.txt"%ll)
			write_text_row([p[t[i][1]] for i in xrange(len(t))], po+"%01d.txt"%ll)
	return


def getindexdata(stack, partids, partstack, myid, nproc):
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
	

def compute_fscs(stack, outputdir, chunkname, newgoodname, fscoutputdir, doit, keepchecking, nproc, myid, main_node):
	#  Compute reconstructions per group from good particles only to get FSC curves
	#  We will compute two FSC curves - from not averaged parameters and from averaged parameters
	#     So, we have to build two sets:
	#    not averaged  (A2+C3) versus (B0+D5)
	#          averaged  (A0+C1) versus (B3+D4)
	#    This requires pulling good subsets given by goodX*;  I am not sure why good, sxconsistency above produced newgood text files.
	#                                                                 Otherwise, I am not sure what newbad will contain.
	# Input that should vary:  
	#    "bdb:"+os.path.join(outputdir,"chunk%01d%01d"%(procid,i))
	#    os.path.join(outputdir,"newgood%01d.txt"%procid)
	#  Output should be in a separate directory "fscoutputdir"

	if(myid == main_node):
		if keepchecking:
			if(os.path.exists(fscoutputdir)):
				doit = 0
				print("Directory  ",fscoutputdir,"  exists!")
			else:
				doit = 1
				keepchecking = False
		else:
			doit = 1
		if doit:
			cmd = "{} {}".format("mkdir", fscoutputdir)
			cmdexecute(cmd)
	mpi_barrier(MPI_COMM_WORLD)
	
	#  not averaged
	doit, keepchecking = checkstep(os.path.join(fscoutputdir,"volfscn0.hdf"), keepchecking, myid, main_node)
	if doit:
		if(myid == main_node):
			#  A2+C3
			#     indices
			write_text_file( \
				map(int, read_text_file(os.path.join(outputdir,chunkname+"0.txt")))+map(int, read_text_file(os.path.join(outputdir,chunkname+"2.txt"))),   \
				os.path.join(fscoutputdir,"chunkfn0.txt"))
			#  params
			write_text_row( \
				read_text_row(os.path.join(outputdir,newgoodname+"02.txt"))+read_text_row(os.path.join(outputdir,newgoodname+"22.txt")), \
				os.path.join(fscoutputdir,"params-chunkfn0.txt"))

		mpi_barrier(MPI_COMM_WORLD)

		projdata = getindexdata(stack, os.path.join(fscoutputdir,"chunkfn0.txt"), os.path.join(fscoutputdir,"params-chunkfn0.txt"), myid, nproc)
		if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		del projdata
		if(myid == main_node):
			vol.write_image(os.path.join(fscoutputdir,"volfscn0.hdf"))
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"Executed successfully: 3D reconstruction of", os.path.join(fscoutputdir,"volfscn0.hdf"))
		del vol


	doit, keepchecking = checkstep(os.path.join(fscoutputdir,"volfscn1.hdf"), keepchecking, myid, main_node)
	if doit:
		if(myid == main_node):
			#  B0+D5
			#     indices
			write_text_file( \
				map(int, read_text_file(os.path.join(outputdir,chunkname+"1.txt")))+map(int, read_text_file(os.path.join(outputdir,chunkname+"3.txt"))),   \
				os.path.join(fscoutputdir,"chunkfn1.txt"))
			#  params
			write_text_row( \
				read_text_row(os.path.join(outputdir,newgoodname+"10.txt"))+read_text_row(os.path.join(outputdir,newgoodname+"32.txt")), \
				os.path.join(fscoutputdir,"params-chunkfn1.txt"))

		mpi_barrier(MPI_COMM_WORLD)

		projdata = getindexdata(stack, os.path.join(fscoutputdir,"chunkfn1.txt"), os.path.join(fscoutputdir,"params-chunkfn1.txt"), myid, nproc)
		if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		del projdata
		if(myid == main_node):
			vol.write_image(os.path.join(fscoutputdir,"volfscn1.hdf"))
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"Executed successfully: 3D reconstruction of", os.path.join(fscoutputdir,"volfscn1.hdf"))
		del vol

	#      averaged
	doit, keepchecking = checkstep(os.path.join(fscoutputdir,"volfsca0.hdf"), keepchecking, myid, main_node)
	if doit:
		if(myid == main_node):
			#  A0+C1
			#     indices
			write_text_file( \
				map(int, read_text_file(os.path.join(outputdir,chunkname+"0.txt")))+map(int, read_text_file(os.path.join(outputdir,chunkname+"2.txt"))),   \
				os.path.join(fscoutputdir,"chunkfa0.txt"))
			#  params
			write_text_row( \
				read_text_row(os.path.join(outputdir,newgoodname+"00.txt"))+read_text_row(os.path.join(outputdir,newgoodname+"20.txt")), \
				os.path.join(fscoutputdir,"params-chunkfa0.txt"))
		mpi_barrier(MPI_COMM_WORLD)

		projdata = getindexdata(stack, os.path.join(fscoutputdir,"chunkfa0.txt"), os.path.join(fscoutputdir,"params-chunkfa0.txt"), myid, nproc)
		if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		del projdata
		if(myid == main_node):
			vol.write_image(os.path.join(fscoutputdir,"volfsca0.hdf"))
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"Executed successfully: 3D reconstruction of", os.path.join(fscoutputdir,"volfsca0.hdf"))
		del vol


	doit, keepchecking = checkstep(os.path.join(fscoutputdir,"volfsca1.hdf"), keepchecking, myid, main_node)
	if doit:
		if(myid == main_node):
			#  B3+D4
			write_text_file( \
				map(int, read_text_file(os.path.join(outputdir,chunkname+"1.txt")))+map(int, read_text_file(os.path.join(outputdir,chunkname+"3.txt"))),   \
				os.path.join(fscoutputdir,"chunkfa1.txt"))
			#  params
			write_text_row( \
				read_text_row(os.path.join(outputdir,newgoodname+"11.txt"))+read_text_row(os.path.join(outputdir,newgoodname+"31.txt")), \
				os.path.join(fscoutputdir,"params-chunkfa1.txt"))
		mpi_barrier(MPI_COMM_WORLD)

		projdata = getindexdata(stack, os.path.join(fscoutputdir,"chunkfa1.txt"), os.path.join(fscoutputdir,"params-chunkfa1.txt"), myid, nproc)
		if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
		del projdata
		if(myid == main_node):
			vol.write_image(os.path.join(fscoutputdir,"volfsca1.hdf"))
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"Executed successfully: 3D reconstruction of", os.path.join(fscoutputdir,"volfsca1.hdf"))
		del vol


 
	#  Get updated FSC curves
	if(myid == main_node):
		if(ali3d_options.mask3D is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
		else:
			mask = get_im(ali3d_options.mask3D)
		if keepchecking:
			if(os.path.exists(os.path.join(fscoutputdir,"fscn.txt"))):
				doit = 0
			else:
				doit = 1
				keepchecking = False
		else:  doit = 1
		if  doit:  fsc(get_im(os.path.join(fscoutputdir,"volfscn0.hdf"))*mask,\
				get_im(os.path.join(fscoutputdir,"volfscn1.hdf"))*mask,\
				1.0,os.path.join(fscoutputdir,"fscn.txt") )
		if keepchecking:
			if(os.path.exists(os.path.join(fscoutputdir,"fsca.txt"))):
				doit = 0
			else:
				doit = 1
				keepchecking = False
		else:  doit = 1
		if  doit:  fsc(get_im(os.path.join(fscoutputdir,"volfsca0.hdf"))*mask,\
				get_im(os.path.join(fscoutputdir,"volfsca1.hdf"))*mask,\
				1.0,os.path.join(fscoutputdir,"fsca.txt") )

		nfsc = read_text_file(os.path.join(fscoutputdir,"fscn.txt") ,-1)
		currentres = 0.5
		ns = len(nfsc[1])
		for i in xrange(1,ns-1):
			if ( (2*nfsc[1][i]/(1.0+nfsc[1][i]) ) < 0.5):
				currentres = nfsc[0][i-1]
				break
		print("  Current resolution ",i,currentres)
	else:
		currentres = 0.0
	currentres = bcast_number_to_all(currentres, source_node = main_node)
	if(currentres < 0.0):
		if(myid == main_node):
			print("  Something wrong with the resolution, cannot continue")
		mpi_finalize()
		exit()

	mpi_barrier(MPI_COMM_WORLD)
	return  currentres, doit, keepchecking



def compute_resolution(vol, radi, nnxo, fscoutputdir):
	# this function is single processor
	#  Get updated FSC curves
	if(ali3d_options.mask3D is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
	else:                              mask = get_im(ali3d_options.mask3D)
	nfsc = fsc(vol[0]*mask,vol[1]*mask, 1.0,os.path.join(fscoutputdir,"fsc.txt") )
	currentres = 0.5
	ns = len(nfsc[1])
	for i in xrange(1,ns-1):
		if ( (2*nfsc[1][i]/(1.0+nfsc[1][i]) ) < 0.5):
			currentres = nfsc[0][i-1]
			break
	#print("  Current resolution ",i,currentres)
	if(currentres < 0.0):
		print("  Something wrong with the resolution, cannot continue")
		mpi_finalize()
		exit()

	return  currentres



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
	pwreference = "rotpw3i3.txt"
#################################



def metamove(paramsdict, partids, partstack, outputdir, procid, myid, main_node, nproc):
	if(myid == main_node):
		log = Logger(BaseLogger_Files())
		log.prefix = os.path.join(outputdir)
		cmd = "mkdir "+log.prefix
		cmdexecute(cmd)
		log.prefix += "/"
	else:  log = None
	mpi_barrier(MPI_COMM_WORLD)

	ali3d_options.delta = paramsdict["delta"]
	ali3d_options.center = paramsdict["center"]
	ali3d_options.ts = paramsdict["ts"]
	ali3d_options.xr = paramsdict["xr"]
	ali3d_options.fl = paramsdict["currentres"]
	ali3d_options.aa = paramsdict["aa"]
	ali3d_options.maxit = paramsdict["maxit"]
	ali3d_options.mask3D = paramsdict["mask3D"]
	projdata = getindexdata(paramsdict["stack"], partids, partstack, myid, nproc)
	if(paramsdict["delpreviousmax"]):
		for i in xrange(len(projdata)):
			try:  projdata[i].del_attr("previousmax")
			except:  pass
	ali3d_options.ou = paramsdict["radius"]  #  This is changed in ali3d_base, but the shrank value is needed in vol recons, fixt it!
	if(myid == main_node):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"METAMOVE parameters")
		spaces = "                 "
		for q in paramsdict:  print("                    => ",q+spaces[len(q):],":  ",paramsdict[q])
		print("                    =>    partids:              ",partids)
		print("                    =>    partstack:            ",partstack)

	#  Run alignment command
	params = ali3d_base(projdata, get_im(paramsdict["refvol"]), \
				ali3d_options, paramsdict["shrink"], mpi_comm = MPI_COMM_WORLD, log = log, \
				nsoft = paramsdict["nsoft"], saturatecrit = paramsdict["saturatecrit"] )
	del log, projdata
	#  store params
	if(myid == main_node):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"Executed successfully: ","ali3d_base_MPI %d"%paramsdict["nsoft"])
		write_text_row(params, os.path.join(outputdir,"params-chunk%01d.txt"%procid) )



def main():

	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	from logger import Logger, BaseLogger_Files
	import sys
	import os
	import time
	import socket
	import user_functions
	from applications import MPI_start_end
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc, do_volume

	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  [output_directory]  initial_volume  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1  --fl --aa --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",       type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "-1",               help="range for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",                help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",                help="angular step of reference projections (default 2)")
	#parser.add_option("--an",       type="string", default= "-1",              help="angular neighborhood for local searches (phi and theta)")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit",    type="int",  default= 400,                  help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--outlier_percentile",     type="float",    default= 95, help="percentile above which outliers are removed every iteration")
	parser.add_option("--iteration_start",     type="int",    default= 0,       help="starting iteration for rviper, 0 means go to the most recent one (default).")
	parser.add_option("--CTF",      action="store_true", default=False,         help="NOT IMPLEMENTED Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")

	#options introduced for the do_volume function
	parser.add_option("--fl",      type="float",  default=0.12,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filte (default 0.12)")
	parser.add_option("--aa",      type="float",  default=0.1,    help="fall-off of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--pwreference",      type="string",  default="",    help="text file with a reference power spectrum (default no power spectrum adjustment)")
	parser.add_option("--mask3D",      type="string",  default=None,    help="3D mask file (default a sphere  WHAT RADIUS??)")
			


	(options, args) = parser.parse_args(sys.argv[1:])

	#print( args)
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

	orgstack = args[0]
	#print(  orgstack,masterdir,volinit )

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
	ali3d_options.center = options.center
	ali3d_options.CTF    = options.CTF
	ali3d_options.ref_a  = options.ref_a
	ali3d_options.snr    = options.snr
	ali3d_options.mask3D = options.mask3D
	ali3d_options.pwreference = options.pwreference
	ali3d_options.fl     = 0.4
	ali3d_options.aa     = 0.1

	if( ali3d_options.xr == "-1" ):  ali3d_options.xr = "2"
	"""
	print( options)

	print( 'ali3d_options',  ali3d_options.ir    ,\
	ali3d_options.rs        ,\
	ali3d_options.ou        ,\
	ali3d_options.xr        ,\
	ali3d_options.yr        ,\
	ali3d_options.ts        ,\
	ali3d_options.an        ,\
	ali3d_options.sym       ,\
	ali3d_options.delta     ,\
	ali3d_options.npad      ,\
	ali3d_options.center    ,\
	ali3d_options.CTF       ,\
	ali3d_options.ref_a     ,\
	ali3d_options.snr       ,\
	ali3d_options.mask3D    ,\
	ali3d_options.fl        ,\
	ali3d_options.aa    \
	)

		#exit()
"""



	mpi_init(0, [])



	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	#mpi_finalize()
	#exit()

	nxinit = -1  #int(280*0.3*2)
	nsoft = 0

	mempernode = 4.0e9


	#  PARAMETERS OF THE PROCEDURE 
	#  threshold error
	thresherr = 0
	fq = 0.11 # low-freq limit to which fuse ref volumes.  Should it be estimated somehow?

	# Get the pixel size, if none set to 1.0, and the original image size
	if(myid == main_node):
		a = get_im(orgstack)
		nnxo = a.get_xsize()
		if ali3d_options.CTF:
			i = a.get_attr('ctf')
			pixel_size = i.apix
		else:
			pixel_size = 1.0
		del a
	else:
		nnxo = 0
		pixel_size = 1.0
	pixel_size = bcast_number_to_all(pixel_size, source_node = main_node)
	nnxo = bcast_number_to_all(nnxo, source_node = main_node)

	if(radi < 1):  radi = nnxo//2-2
	elif((2*radi+2)>nnxo):  ERROR("HERE","particle radius set too large!",1)
	ali3d_options.ou = radi
	if(nxinit < 0):  nxinit = min(32, nnxo)

	nxshrink = nxinit
	shrink = float(nxshrink)/float(nnxo)
	#  Randomize initial parameters (??)
	#  cmd = "{} {} {}".format("sxheader.py", stack, "--rand_alpha  --params=xform.projection")
	#  cmdexecute(cmd)

	#  MASTER DIRECTORY
	if(myid == main_node):
		if( masterdir == ""):
			timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
			masterdir = "master"+timestring
			li = len(masterdir)
			cmd = "{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)
			keepchecking = 0
		else:
			li = 0
			keepchecking = 1
	else:
		li = 0
		keepchecking = 1

	li = mpi_bcast(li,1,MPI_INT,main_node,MPI_COMM_WORLD)[0]

	if( li > 0 ):
		masterdir = mpi_bcast(masterdir,li,MPI_CHAR,main_node,MPI_COMM_WORLD)
		masterdir = string.join(masterdir,"")

	#  create a vstack from input stack to the local stack in masterdir
	#  Stack name set to default
	stack = "bdb:"+masterdir+"/rdata"
	# Initialization of stacks
	if(myid == main_node):
		if keepchecking:
			if(os.path.exists(os.path.join(masterdir,"EMAN2DB/rdata.bdb"))):  doit = False
			else:  doit = True
		else:  doit = True
		if  doit:
			if(orgstack[:4] == "bdb:"):	cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+stack)
			else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, stack)
			cmdexecute(cmd)
			cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", stack)
			cmdexecute(cmd)
			keepchecking = False
		total_stack = EMUtil.get_image_count(stack)
		junk = get_im(stack)
		nnxo = junk.get_xsize()
		del junk
	else:
		total_stack = 0
		nnxo = 0

	total_stack = bcast_number_to_all(total_stack, source_node = main_node)
	nnxo        = bcast_number_to_all(nnxo, source_node = main_node)

	#  INITIALIZATION

	#  Run exhaustive projection matching to get initial orientation parameters
	#  Estimate initial resolution
	initdir = os.path.join(masterdir,"main000")
		

	xr = min(8,(nnxo - (2*radi+1))//2)
	if(xr > 3):  ts = "2"
	else:  ts = "1"

	paramsdict = {	"stack":stack,"delta":"2.0", "ts":ts, "xr":"%f"%xr, "an":"-1", "center":options.center, "maxit":1, \
					"currentres":0.4, "aa":0.1, "radius":radi, "nsoft":0, "delpreviousmax":True, "shrink":1.0, "saturatecrit":1.0, \
					"refvol":volinit, "mask3D":options.mask3D}
	doit, keepchecking = checkstep(initdir, keepchecking, myid, main_node)
	if  doit:

		partids = os.path.join(masterdir, "ids.txt")
		partstack = os.path.join(masterdir, "paramszero.txt")
		
		if( myid == main_node ):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"INITIALIZATION")
			write_text_file(range(total_stack), partids)
			write_text_row([[0.0,0.0,0.0,0.0,0.0] for i in xrange(total_stack) ], partstack)

		metamove(paramsdict, partids, partstack, initdir, 0, myid, main_node, nproc)

		#  store params
		partids = [None]*2
		for procid in xrange(2):  partids[procid] = os.path.join(initdir,"chunk%01d.txt"%procid)
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(initdir,"params-chunk%01d.txt"%procid)
		from random import shuffle
		if(myid == main_node):
			print(line,"Executed successfully: ","initialization ali3d_base_MPI  %d"%nsoft)
			#  split randomly
			params = read_text_row(os.path.join(initdir,"params-chunk0.txt"))
			ll = range(total_stack)
			shuffle(ll)
			l1 = ll[:total_stack//2]
			l2 = ll[total_stack//2:]
			del ll
			l1.sort()
			l2.sort()
			write_text_file(l1,partids[0])
			write_text_file(l2,partids[1])
			write_text_row([params[i] for i in l1], partstack[0]) 
			write_text_row([params[i] for i in l2], partstack[1])
			del params, l1, l2
		mpi_barrier(MPI_COMM_WORLD)
		#  Now parallel
		vol = [None]*2
		for procid in xrange(2):
			projdata = getindexdata(stack, partids[procid], partstack[procid], myid, nproc)
			if ali3d_options.CTF:  vol[procid] = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
			else:                  vol[procid] = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
			del projdata
			if( myid == main_node):
				vol[procid].write_image(os.path.join(initdir,"vol%01d.hdf"%procid) )
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				print(  line,"Generated inivol #%01d "%procid)


		if(myid == main_node):
			currentres = compute_resolution(vol, radi, nnxo, initdir)		
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(  line,"Initial resolution %6.4f"%currentres)
			write_text_file([currentres],os.path.join(initdir,"current_resolution.txt"))
		else:  currentres = 0.0
		currentres = bcast_number_to_all(currentres, source_node = main_node)
	else:
		if(myid == main_node): currentres = read_text_file(os.path.join(initdir,"current_resolution.txt"))[0]		
		else:  currentres = 0.0
		currentres = bcast_number_to_all(currentres, source_node = main_node)

	# set for the first iteration
	nxshrink = min(max(32, int((currentres+paramsdict["aa"]/2.)*2*nnxo + 0.5)), nnxo)
	shrink = float(nxshrink)/nnxo
	tracker = {"previous-resolution":currentres, "movedup":False,"eliminated-outliers":False,\
				"previous-nx":nxshrink, "previous-shrink":shrink, "extention":0, "bestsolution":0}
	
	previousoutputdir = initdir
	#  MAIN ITERATION
	mainiteration = 0
	keepgoing = 1
	while(keepgoing):
		mainiteration += 1


		#  prepare output directory
		mainoutputdir = os.path.join(masterdir,"main%03d"%mainiteration)

		if(myid == main_node):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"MAIN ITERATION #",mainiteration, shrink, nxshrink)
			if keepchecking:
				if(os.path.exists(mainoutputdir)):
					doit = 0
					print("Directory  ",mainoutputdir,"  exists!")
				else:
					doit = 1
					keepchecking = False
			else:
				doit = 1

			if doit:
				cmd = "{} {}".format("mkdir", mainoutputdir)
				cmdexecute(cmd)

		# prepare names of input file names
		partids = [None]*2
		for procid in xrange(2):  partids[procid] = os.path.join(previousoutputdir,"chunk%01d.txt"%procid)
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(previousoutputdir,"params-chunk%01d.txt"%procid)

		mpi_barrier(MPI_COMM_WORLD)


		#mpi_finalize()
		#exit()

		#print("RACING  A ",myid)
		outvol = [os.path.join(previousoutputdir,"vol%01d.hdf"%procid) for procid in xrange(2)]
		for procid in xrange(2):
			doit, keepchecking = checkstep(outvol[procid], keepchecking, myid, main_node)

			if  doit:
				from multi_shc import do_volume
				projdata = getindexdata(stack, partids[procid], partstack[procid], myid, nproc)
				if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
				else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
				del projdata
				if( myid == main_node):
					vol.write_image(outvol[procid])
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(  line,"Generated inivol #%01d "%procid)
				del vol

		if(myid == main_node):
			if keepchecking:
				procid = 1
				if(os.path.join(mainoutputdir,"fusevol%01d.hdf"%procid)):
					doit = 0
				else:
					doit = 1
					keepchecking = False
			else:
				doit = 1
			if doit:
				vol = [get_im(outvol[procid]) for procid in xrange(2) ]
				fq = 0.11 # which part to fuse
				fuselowf(vol, fq)
				for procid in xrange(2):  vol[procid].write_image(os.path.join(mainoutputdir,"fusevol%01d.hdf"%procid) )
				del vol
		else:  doit = 0
		mpi_barrier(MPI_COMM_WORLD)
		doit = bcast_number_to_all(doit, source_node = main_node)

		#  Refine two groups at a current resolution
		lastring = int(shrink*radi + 0.5)
		if(lastring < 2):
			print(  line,"ERROR!!   lastring too small  ", radi, shrink, lastring)
			break

		#  REFINEMENT
		for procid in xrange(2):
			coutdir = os.path.join(mainoutputdir,"loga%01d"%procid)
			doit, keepchecking = checkstep(coutdir  , keepchecking, myid, main_node)

			paramsdict = {	"stack":stack,"delta":"%f"%round(degrees(atan(1.0/lastring)), 2) , "ts":"1", "xr":"2", "an":"-1", "center":options.center, "maxit":1500,  \
							"currentres":currentres, "aa":0.1, "radius":radi, "nsoft":1, "saturatecrit":0.75, "delpreviousmax":True, "shrink":shrink, \
							"refvol":os.path.join(mainoutputdir,"fusevol%01d.hdf"%procid),"mask3D":options.mask3D }

			if  doit:

				metamove(paramsdict, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)

		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(mainoutputdir, "loga%01d"%procid, "params-chunk%01d.txt"%procid)

		for procid in xrange(2):
			outvol = os.path.join(mainoutputdir,"loga%01d"%procid,"shcvol%01d.hdf"%procid)
			doit, keepchecking = checkstep(outvol, keepchecking, myid, main_node)

			if  doit:
				from multi_shc import do_volume
				projdata = getindexdata(stack, partids[procid], partstack[procid], myid, nproc)
				if ali3d_options.CTF:  vol = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
				else:                  vol = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
				del projdata
				if( myid == main_node):
					vol.write_image(outvol)
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(  line,"Generated shcvol #%01d "%procid)
				del vol

		if(myid == main_node):
			if keepchecking:
				procid = 1
				if(os.path.join(mainoutputdir,"loga%01d"%procid,"fusevol%01d.hdf"%procid) ):
					doit = 0
				else:
					doit = 1
					keepchecking = False
			else:
				doit = 1
			if doit:
				vol = []
				for procid in xrange(2):  vol.append(get_im(os.path.join(mainoutputdir,"loga%01d"%procid,"shcvol%01d.hdf"%procid) ))
				fq = 0.11 # which part to fuse
				fuselowf(vol, fq)
				for procid in xrange(2):  vol[procid].write_image( os.path.join(mainoutputdir,"loga%01d"%procid,"fusevol%01d.hdf"%procid) )
				del vol
		else:  doit = 0
		mpi_barrier(MPI_COMM_WORLD)
		doit = bcast_number_to_all(doit, source_node = main_node)

			
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(mainoutputdir,"loga%01d"%procid,"params-chunk%01d.txt"%procid)

		for procid in xrange(2):
			coutdir = os.path.join(mainoutputdir,"logb%01d"%procid)
			doit, keepchecking = checkstep(coutdir, keepchecking, myid, main_node)

			#  Run exhaustive to finish up matching
			paramsdict = {	"stack":stack,"delta":"%f"%round(degrees(atan(1.0/lastring)), 2) , "ts":"1", "xr":"2", "an":"-1", "center":options.center, "maxit":10, \
							"currentres":currentres, "aa":0.1, "radius":radi, "nsoft":0, "saturatecrit":0.95, "delpreviousmax":True, "shrink":shrink, \
							"refvol":os.path.join(mainoutputdir,"loga%01d"%procid,"fusevol%01d.hdf"%procid), "mask3D":options.mask3D }

			if  doit:

				metamove(paramsdict, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)

		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(mainoutputdir,"logb%01d"%procid,"params-chunk%01d.txt"%procid)
		vol = [None]*2

		for procid in xrange(2):
			doit, keepchecking = checkstep(os.path.join(mainoutputdir,"vol%01d.hdf"%procid), keepchecking, myid, main_node)

			#  sxrecons3d.py  (full size)
			if  doit:
				projdata = getindexdata(stack, partids[procid], partstack[procid], myid, nproc)
				if ali3d_options.CTF:  vol[procid] = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
				else:                  vol[procid] = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
				del projdata
				if( myid == main_node):
					vol[procid].write_image(os.path.join(mainoutputdir,"vol%01d.hdf"%procid))
					line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
					print(  line,"Generated vol #%01d "%procid)
			else:
				if(myid == main_node):
					vol[procid] = get_im(os.path.join(mainoutputdir,"vol%01d.hdf"%procid))
				else:
					#  The other nodes do not need volumes
					vol[procid] = model_blank(nnxo,nnxo,nnxo)


		doit, keepchecking = checkstep(os.path.join(mainoutputdir,"current_resolution.txt"), keepchecking, myid, main_node)
		newres = 0.0
		if  doit:
			if(myid == main_node):
				newres = compute_resolution(vol, radi, nnxo, mainoutputdir)		
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				print(  line,"Current resolution %6.4f"%newres)
				write_text_file([newres],os.path.join(mainoutputdir,"current_resolution.txt"))
		else:
			if(myid == main_node): newres = read_text_file( os.path.join(mainoutputdir,"current_resolution.txt") )[0]		
		newres = bcast_number_to_all(newres, source_node = main_node)

		doit, keepchecking = checkstep(os.path.join(mainoutputdir,"volf.hdf"), keepchecking, myid, main_node)
		if  doit:
			volf = (vol[0]+vol[1])*0.5
			ali3d_options.fl = newres
			volf = do_volume(volf,ali3d_options, mainiteration, mpi_comm = MPI_COMM_WORLD)
			if(myid == main_node): volf.write_image(os.path.join(mainoutputdir,"volf.hdf"))
			del volf, vol
		

		mpi_barrier(MPI_COMM_WORLD)

		#print("RACING  X ",myid)
		if(newres == currentres):

			for procid in xrange(2):
				coutdir = os.path.join(mainoutputdir,"logc%01d"%procid)
				doit, keepchecking = checkstep(coutdir, keepchecking, myid, main_node)

				if  doit:
					#
					paramsdict = {	"stack":stack,"delta":"%f"%round(degrees(atan(1.0/lastring)), 2) , "ts":"1", "xr":"2", "an":"-1", "center":options.center, "maxit":1,  \
									"currentres":newres, "aa":0.1, "radius":radi, "nsoft":0, "saturatecrit":0.95, "delpreviousmax":True, "shrink":shrink, \
									"refvol":os.path.join(mainoutputdir,"vol%01d.hdf"%(1-procid)), "mask3D":options.mask3D }

					metamove(paramsdict, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)

			# identify bad apples
			doit, keepchecking = checkstep(os.path.join(mainoutputdir,"badapples.txt"), keepchecking, myid, main_node)
			if  doit:
				if(myid == main_node):
					from utilities import get_symt
					from pixel_error import max_3D_pixel_error
					ts = get_symt(ali3d_options.sym)
					badapples = []
					deltaerror = 2.0
					for procid in xrange(2):
						bad = []
						ids  = map(int,read_text_file( partids[procid] ))
						oldp = read_text_row(partstack[procid])
						newp = read_text_row(os.path.join(mainoutputdir,"logc%01d"%procid,"params-chunk%01d.txt"%procid))
						for i in xrange(len(ids)):

							t1 = Transform({"type":"spider","phi":oldp[i][0],"theta":oldp[i][1],"psi":oldp[i][2]})
							t1.set_trans(Vec2f(-oldp[i][3]*shrink, -oldp[i][4]*shrink))
							t2 = Transform({"type":"spider","phi":newp[i][0],"theta":newp[i][1],"psi":newp[i][2]})
							t2.set_trans(Vec2f(-newp[i][3]*shrink, -newp[i][4]*shrink))
							if(len(ts) > 1):
								# only do it if it is not c1
								pixel_error = +1.0e23
								for kts in ts:
									ut = t2*kts
									# we do not care which position minimizes the error
									pixel_error = min(max_3D_pixel_error(t1, ut, lastring), pixel_error)
							else:
								pixel_error = max_3D_pixel_error(t1, t2, lastring)

							if(pixel_error > deltaerror):
								bad.append(i)
						if(len(bad)>0):
							badapples += [ids[bad[i]] for i in xrange(len(bad))]
							for i in xrange(len(bad)-1,-1,-1):
								del oldp[bad[i]],ids[bad[i]]
						write_text_file(ids,os.path.join(mainoutputdir,"chunk%01d.txt"%procid))
						write_text_row(oldp,os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))
					if(len(badapples)>0):
						badapples.sort()
						write_text_file(badapples,os.path.join(mainoutputdir,"badapples.txt"))
						eli = 100*float(len(badapples))/float(len(ids))
						line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
						print(line,"Elimination of outliers: %5.1f percent"%eli )
					else:  eli = 0.0
					del badapples, oldp,ids,bad,newp,ts
				else:  eli =0.0
				eli = bcast_number_to_all(eli, source_node = main_node)	
				if(eli > 0.0):  tracker["eliminated-outliers"] = True
				else:  tracker["eliminated-outliers"] = False
		else:
			tracker["eliminated-outliers"] = False
			#mpi_finalize()
			#exit()


		keepgoing = 0
		if( newres > currentres or tracker["eliminated-outliers"]):
			if(myid == main_node):  print("  Resolution improved, full steam ahead!")

			if( newres > currentres ):  tracker["movedup"] = True
			else:   tracker["movedup"] = False
			shrink = min(2*newres + paramsdict["aa"], 1.0)
			tracker["extension"] = 4
			nxshrink = min(int(nnxo*shrink + 0.5) + tracker["extension"],nnxo)
			tracker["previous-resolution"] = newres
			currentres = newres
			tracker["bestsolution"] = mainiteration
			tracker["eliminated-outliers"] = False
			keepgoing = 1
		
		elif(newres < currentres):
			if(not tracker["movedup"] and tracker["extension"] < 2):
				keepgoing = 0
				if(myid == main_node):  print("  Cannot improve resolution, the best result is in the directory main%03d"%tracker["bestsolution"])
			else:
				if(not tracker["movedup"] and tracker["extension"] > 1 and mainiteration > 1):
					if(myid == main_node):  print("  Resolution decreased.  Will decrease target resolution and will fall back on the best so far:  main%03d"%tracker["bestsolution"])
					bestoutputdir = os.path.join(masterdir,"main%03d"%tracker["bestsolution"])
				elif( tracker["movedup"] and tracker["extension"] > 1 and mainiteration > 1):
					if(myid == main_node):  print("  Resolution decreased.  Will decrease target resolution and will try starting from previous stage:  main%03d"%(mainiteration - 1))
					bestoutputdir = os.path.join(masterdir,"main%03d"%(mainiteration-1))
				else:  # missing something here?
					if(myid == main_node):  print(" Should not be here, ERROR 175!")
					break
					mpi_finale()
					exit()
					
				cmd = "{} {} {}".format("cp -p ",os.path.join(bestoutputdir,"chunk%01d.txt"%procid) , os.path.join(mainoutputdir,"chunk%01d.txt"%procid))
				cmdexecute(cmd)
				cmd = "{} {} {}".format("cp -p ",os.path.join(bestoutputdir,"params-chunk%01d.txt"%procid), os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))
				cmdexecute(cmd)
				if(myid == main_node):  currentres = read_text_file( os.path.join(bestoutputdir,"current_resolution.txt") )[0]	
				currentres = bcast_number_to_all(currentres, source_node = main_node)	

				shrink = min(2*currentres + paramsdict["aa"], 1.0)
				tracker["extension"] -= 1
				nxshrink = min(int(nnxo*shrink + 0.5) + tracker["extension"],nnxo)
				tracker["previous-resolution"] = newres
				currentres = newres
				tracker["eliminated-outliers"] = False
				tracker["movedup"] = False
				keepgoing = 1
			

		elif(newres == currentres):
			if( tracker["extension"] > 0 ):
				if(myid == main_node):  print("The resolution did not improve. This is look ahead move.  Let's try to relax slightly and hope for the best")
				tracker["extension"] -= 1

				tracker["movedup"] = False
			
				shrink = min(2*newres + paramsdict["aa"], 1.0)
				nxshrink = min(int(nnxo*shrink + 0.5) + tracker["extension"],nnxo)

				shrink = min(2*currentres + paramsdict["aa"], 1.0)
				nxshrink = min(int(nnxo*shrink + 0.5) + tracker["extension"],nnxo)
				tracker["previous-resolution"] = newres
				currentres = newres
				tracker["eliminated-outliers"] = False
				tracker["movedup"] = False
				keepgoing = 1
			else:
				if(myid == main_node):  print("The resolution did not improve.")

			

		if( keepgoing == 1 ):
			if(myid == main_node):
				print("  New shrink and image dimension :",shrink,nxshrink)
				#  Will continue, so update the params files
				for procid in xrange(2):
					if(not os.path.exists(os.path.join(mainoutputdir,"chunk%01d.txt"%procid))): 
						cmd = "{} {} {}".format("cp -p ", partids[procid] , os.path.join(mainoutputdir,"chunk%01d.txt"%procid))
						cmdexecute(cmd)
					if(not os.path.exists(os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))):
						cmd = "{} {} {}".format("cp -p ", partstack[procid], os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))
						cmdexecute(cmd)
			
			previousoutputdir = mainoutputdir
			tracker["previous-shrink"]     = shrink
			tracker["previous-nx"]         = nxshrink

		else:
			if(myid == main_node):
				print("  Terminating, the best solution is in the directory main%03d"%tracker["bestsolution"])
		mpi_barrier(MPI_COMM_WORLD)

	mpi_finalize()





if __name__=="__main__":
	main()

