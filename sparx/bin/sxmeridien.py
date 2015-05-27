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

def get_pixercutoff(radius, delta = 2.0, dsx = 0.5):
	#  Estimate tolerable error based on current delta and shrink.
	#  Current radius (radi*shrink)
	#  delta - current angular step
	#  dsx   - expected pixel error (generally, for polar searches it is 0.5, for gridding 0.1.
	t1 = Transform({"type":"spider","phi":0.0,"theta":0.0,"psi":0.0})
	t1.set_trans(Vec2f(0.0, 0.0))
	t2 = Transform({"type":"spider","phi":0.0,"theta":delta,"psi":delta})
	t2.set_trans(Vec2f(dsx, dsx))
	return max_3D_pixel_error(t1, t2, radius)


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


def get_resolution(vol, radi, nnxo, fscoutputdir):
	# this function is single processor
	#  Get updated FSC curves, user can also provide a mask using radi variable
	import types
	if(type(radi) == int):
		if(ali3d_options.mask3D is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
		else:                              mask = get_im(ali3d_options.mask3D)
	else:  mask = radi
	nfsc = fsc(vol[0]*mask,vol[1]*mask, 1.0,os.path.join(fscoutputdir,"fsc.txt") )
	currentres = -1.0
	ns = len(nfsc[1])
	#  This is actual resolution, as computed by 2*f/(1+f)
	for i in xrange(1,ns-1):
		if ( nfsc[1][i] < 0.333333333333333333333333):
			currentres = nfsc[0][i-1]
			break
	if(currentres < 0.0):
		print("  Something wrong with the resolution, cannot continue")
		mpi_finalize()
		exit()
	"""
	lowpass = 0.5
	ns = len(nfsc[1])
	#  This is resolution used to filter half-volumes
	for i in xrange(1,ns-1):
		if ( nfsc[1][i] < 0.5 ):
			lowpass = nfsc[0][i-1]
			break
	"""
	lowpass, falloff = fit_tanh1(nfsc, 0.01)

	return  round(lowpass,4), round(falloff,4), round(currentres,2)

def compute_resolution(stack, outputdir, partids, partstack, radi, nnxo, CTF, myid, main_node, nproc):
	vol = [None]*2
	fsc = [None]*2
	if(ali3d_options.mask3D is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
	else:                              mask = get_im(ali3d_options.mask3D)

	for procid in xrange(2):
		#  sxrecons3d.py  (full size)
		if(len(stack) > 2):
			projdata = getindexdata(stack, partids[procid], partstack[procid], myid, nproc)
			#if CTF:  vol[procid] = recons3d_4nn_ctf_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)
			#else:    vol[procid] = recons3d_4nn_MPI(myid, projdata, symmetry=ali3d_options.sym, npad = 2)

			if CTF:
				from reconstruction import rec3D_MPI
				vol[procid],fsc[procid] = rec3D_MPI(projdata, snr = 1.0, symmetry = ali3d_options.sym, \
					mask3D = mask, fsc_curve = os.path.join(outputdir,"within-fsc%01d.txt"%procid), \
					myid = myid, main_node = main_node, odd_start = 1, eve_start = 0, finfo = None, npad = 2)
			else :
				from reconstruction import rec3D_MPI_noCTF
				vol[procid],fsc[procid] = rec3D_MPI_noCTF(projdata, symmetry = ali3d_options.sym, \
					mask3D = mask, fsc_curve = os.path.join(outputdir,"within-fsc%01d.txt"%procid), \
					myid = myid, main_node = main_node, odd_start = 1, eve_start = 0, finfo = None, npad = 2)


			del projdata
		else:
			#  There are only two entries, these have to be volumes
			vol = stack
		if( myid == main_node):
			vol[procid].write_image(os.path.join(outputdir,"vol%01d.hdf"%procid))
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(  line,"Generated vol #%01d "%procid)
	del fsc
	lowpass    = 0.0
	falloff    = 0.0
	currentres = 0.0

	if(myid == main_node):
		lowpass, falloff, currentres = get_resolution(vol, mask, nnxo, outputdir)
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(  line,"Current resolution %6.2f, low-pass filter cut-off %6.4f and fall-off %6.4f"%(currentres,lowpass,falloff))
		write_text_row([lowpass, falloff, currentres],os.path.join(outputdir,"current_resolution.txt"))
	#  Returns: low-pass filter cutoff;  low-pass filter falloff;  current resolution
	currentres = bcast_number_to_all(currentres, source_node = main_node)
	lowpass    = bcast_number_to_all(lowpass, source_node = main_node)
	falloff    = bcast_number_to_all(falloff, source_node = main_node)
	return round(lowpass,4), round(falloff,4), round(currentres,2)

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



def metamove(paramsdict, partids, partstack, outputdir, procid, myid, main_node, nproc):
	#  Reads from paramsdict["stack"] particles partids set parameters in partstack
	#    and do refinement as specified in paramsdict
	#
	#  Will create outputdir
	#  Will write to outputdir output parameters: params-chunk0.txt and params-chunk1.txt
	if(myid == main_node):
		#  Create output directory
		log = Logger(BaseLogger_Files())
		log.prefix = os.path.join(outputdir)
		cmd = "mkdir "+log.prefix
		cmdexecute(cmd)
		log.prefix += "/"
	else:  log = None
	mpi_barrier(MPI_COMM_WORLD)

	ali3d_options.delta  = paramsdict["delta"]
	ali3d_options.center = paramsdict["center"]
	ali3d_options.ts     = paramsdict["ts"]
	ali3d_options.xr     = paramsdict["xr"]
	#  low pass filter is applied to shrank data, so it has to be adjusted
	ali3d_options.fl     = paramsdict["lowpass"]/paramsdict["shrink"]
	ali3d_options.initfl = paramsdict["initialfl"]/paramsdict["shrink"]
	ali3d_options.aa     = paramsdict["falloff"]
	ali3d_options.maxit  = paramsdict["maxit"]
	ali3d_options.mask3D = paramsdict["mask3D"]
	ali3d_options.an	 = paramsdict["an"]
	projdata = getindexdata(paramsdict["stack"], partids, partstack, myid, nproc)
	if(paramsdict["delpreviousmax"]):
		for i in xrange(len(projdata)):
			try:  projdata[i].del_attr("previousmax")
			except:  pass
	ali3d_options.ou = paramsdict["radius"]  #  This is changed in ali3d_base, but the shrank value is needed in vol recons, fixt it!
	if(myid == main_node):
		print_dict(paramsdict,"METAMOVE parameters")
		print("                    =>  actual lowpass      :  ",ali3d_options.fl)
		print("                    =>  actual init lowpass :  ",ali3d_options.initfl)
		if(len(ali3d_options.pwreference)>0): \
		print("                    =>  PW adjustment       :  ",ali3d_options.pwreference)
		print("                    =>  partids             :  ",partids)
		print("                    =>  partstack           :  ",partstack)
		
	if(ali3d_options.fl > 0.46):  ERROR("Low pass filter in metamove > 0.46 on the scale of shrank data","sxmeridien",1,myid) 

	#  Run alignment command
	if(paramsdict["local"]): params = local_ali3d_base_MPI(projdata, get_im(paramsdict["refvol"]), \
				ali3d_options, paramsdict["shrink"], mpi_comm = MPI_COMM_WORLD, log = log, \
		    	chunk = 0.25, saturatecrit = paramsdict["saturatecrit"], pixercutoff =  paramsdict["pixercutoff"])
	else: params = ali3d_base(projdata, get_im(paramsdict["refvol"]), \
				ali3d_options, paramsdict["shrink"], mpi_comm = MPI_COMM_WORLD, log = log, \
				nsoft = paramsdict["nsoft"], saturatecrit = paramsdict["saturatecrit"],  pixercutoff =  paramsdict["pixercutoff"] )
	del log, projdata
	#  store params
	if(myid == main_node):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"Executed successfully: ","ali3d_base_MPI %d"%paramsdict["nsoft"],"  number of images:%7d"%len(params))
		write_text_row(params, os.path.join(outputdir,"params-chunk%01d.txt"%procid) )

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
	usage = progname + " stack  [output_directory]  initial_volume  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --fl --aa --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",      		type= "int",   default= 1,			help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",      		type= "int",   default= -1,			help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",      		type= "int",   default= 1,			help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",      		type="string", default= "-1",		help="range for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",      		type="string", default= "-1",		help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",      		type="string", default= "1",		help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",   		type="string", default= "-1",		help="angular step of reference projections during initialization step (default automatically selected based on radius of the structure.)")
	parser.add_option("--an",      		type="string", default= "-1",		help="angular neighborhood for local searches (phi and theta) (Default exhaustive searches)")
	parser.add_option("--center",  		type="float",  default= -1,			help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit",   		type="int",  	default= 400,		help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--outlier_percentile",type="float",    default= 95,	help="percentile above which outliers are removed every iteration")
	parser.add_option("--iteration_start",type="int",    default= 0,		help="starting iteration for rviper, 0 means go to the most recent one (default).")
	parser.add_option("--CTF",     		action="store_true", default=False,	help="Use CTF (Default no CTF correction)")
	parser.add_option("--snr",     		type="float",  default= 1.0,		help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",   		type="string", default= "S",		help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",     		type="string", default= "c1",		help="symmetry of the refined structure")
	parser.add_option("--npad",    		type="int",    default= 2,			help="padding size for 3D reconstruction (default=2)")
	parser.add_option("--nsoft",    	type="int",    default= 1,			help="Use SHC in first phase of refinement iteration (default=1, to turn it off set to 0)")
	parser.add_option("--startangles",  action="store_true", default=False,	help="Use orientation parameters in the input file header to jumpstart the procedure")

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
	ali3d_options.pwreference = ""  #   It will have to be turned on after exhaustive done by setting to options.pwreference
	ali3d_options.fl     = 0.4
	ali3d_options.initfl = 0.4
	ali3d_options.aa     = 0.1

	if( ali3d_options.xr == "-1" ):  ali3d_options.xr = "2"


	mpi_init(0, [])



	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	#mpi_finalize()
	#exit()

	nxinit = -1  #int(280*0.3*2)

	mempernode = 4.0e9


	#  PARAMETERS OF THE PROCEDURE 
	#  threshold error
	thresherr = 0
	fq = 50 # low-freq resolution to which fuse ref volumes. [A] 

	# Get the pixel size, if none set to 1.0, and the original image size
	if(myid == main_node):
		a = get_im(orgstack)
		nnxo = a.get_xsize()
		if( nnxo%2 == 1 ):
			ERROR("Only even-dimensioned data allowed","sxmeridien",1)
		if ali3d_options.CTF:
			i = a.get_attr('ctf')
			pixel_size = i.apix
			fq = pixel_size/fq
		else:
			pixel_size = 1.0
			#  No pixel size, fusing computed as 5 Fourier pixels
			fq = 5.0/nnxo
		del a
	else:
		nnxo = 0
		pixel_size = 1.0
	pixel_size = bcast_number_to_all(pixel_size, source_node = main_node)
	nnxo = bcast_number_to_all(nnxo, source_node = main_node)
	fq   = bcast_number_to_all(fq, source_node = main_node)


	if(radi < 1):  radi = nnxo//2-2
	elif((2*radi+2)>nnxo):  ERROR("Particle radius set too large!","sxmeridien",1,myid)
	ali3d_options.ou = radi
	if(nxinit < 0):  nxinit = min(32, nnxo)
	else:
		if(nxinit%2 == 1): ERROR("Only even dimensions allowed","sxmeridien",1,myid)

	nxshrink = nxinit
	minshrink = 32.0/float(nnxo)
	shrink = max(float(nxshrink)/float(nnxo),minshrink)
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
	#  Do prealignment of 2D data using reference-free alignment
	if(not options.startangles):
		#  The maximum step size for 2D alignment
		nxrsteps = 4
		init2dir = os.path.join(masterdir,"2dalignment")
		doit, keepchecking = checkstep(init2dir, keepchecking, myid, main_node)
		if  doit:
			if( myid == main_node ):
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				print(line,"INITIALIZATION 2D")

			from applications import ali2d_base

			if(myid == main_node):
				import subprocess
				from logger import Logger, BaseLogger_Files
				#  Create output directory
				log2d = Logger(BaseLogger_Files())
				log2d.prefix = os.path.join(init2dir)
				cmd = "mkdir "+log2d.prefix
				outcome = subprocess.call(cmd, shell=True)
				log2d.prefix += "/"
			else:
				outcome = 0
				log2d = None
			txrm = (nnxo - (2*radi-1))//2
			if(txrm < 0):  			ERROR( "ERROR!!   Radius of the structure larger than the window data size permits   %d"%(radi), "sxmeridien",1, myid)
			if(txrm/nxrsteps>0):
				tss = ""
				txr = ""
				while(txrm/nxrsteps>0):
					tts=txrm/nxrsteps
					tss += "  %d"%tts
					txr += "  %d"%(tts*nxrsteps)
					txrm =txrm//2
			else:
				tss = "1"
				txr = "%d"%txrm

			params2d = ali2d_base(stack, init2dir, None, 1, radi, 1, txr, txr, tss, \
				False, 90.0, -1, 20, options.CTF, 1.0, False, \
				"ref_ali2d", "", log2d, \
				nproc, myid, main_node, MPI_COMM_WORLD,\
				write_headers = False)
			#  Convert 2d to 3D parameters
			if( myid == main_node ):
				for i in xrange(len(params2d)):
					params2d = params_2D_3D(params2d[0], params2d[1], params2d[2], int(params2d[3]))
				write_text_row(params2d,os.path.join(init2dir, "initial3Dshifts.txt"))

	#  Run exhaustive projection matching to get initial orientation parameters
	#  Estimate initial resolution
	initdir = os.path.join(masterdir,"main000")
	#  make sure the initial volume is not set to zero outside of a mask, as if it is it will crash the program
	if( myid == main_node and (not options.startangles)):
		viv = get_im(volinit)
		if(options.mask3D == None):  mask33d = model_circle(radi,nnxo,nnxo,nnxo)
		else:  mask33d = get_im(options.mask3D)
		st = Util.infomask(viv, mask33d, False)
		if( st[0] == 0.0 ):
			viv += (model_blank(nnxo,nnxo,nnxo,1.0) - mask33d)*model_gauss_noise(st[1]/1000.0,nnxo,nnxo,nnxo)
			viv.write_image(volinit)
		del mask33d, viv

	#  This is initial setting, has to be initialized here, we do not want it to run too long.
	#    INITIALIZATION THAT FOLLOWS WILL HAVE TO BE CHANGED SO THE USER CAN PROVIDE INITIAL GUESS OF RESOLUTION
	#  If we new the initial resolution, it could be done more densely
	xr = min(txrm,(nnxo - (2*radi-1))//2)
	ts = "1.0"

	delta = int(options.delta)
	if(delta <= 0.0):
		delta = "%f"%round(degrees(atan(1.0/float(radi))), 2)

	paramsdict = {	"stack":stack,"delta":"2.0", "ts":ts, "xr":"%f"%xr, "an":angular_neighborhood, \
					"center":options.center, "maxit":1, "local":False,\
					"lowpass":0.4, "initialfl":0.4, "falloff":0.1, "radius":radi, \
					"nsoft":0, "delpreviousmax":True, "shrink":1.0, "saturatecrit":1.0, "pixercutoff":2.0,\
					"refvol":volinit, "mask3D":options.mask3D}

	doit, keepchecking = checkstep(initdir, keepchecking, myid, main_node)
	if  doit:
		partids = os.path.join(masterdir, "ids.txt")
		partstack = os.path.join(init2dir, "initial3Dshifts.txt")

		if(options.startangles):

			if( myid == main_node ):
				cmd = "mkdir "+initdir
				cmdexecute(cmd)
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				print(line,"INITIALIZATION")
				cmd = "{} {}".format("sxheader.py --params=xform.projection  --export="+os.path.join(initdir,"params-chunk0.txt"), stack)
				cmdexecute(cmd)
				print(line,"Executed successfully: ","Imported initial parameters from the input stack")

		else:
			if( myid == main_node ):
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				print(line,"INITIALIZATION 3D")
				write_text_file(range(total_stack), partids)

			metamove(paramsdict, partids, partstack, initdir, 0, myid, main_node, nproc)
			if(myid == main_node):
				print(line,"Executed successfully: ","initialization ali3d_base_MPI")

			

		#  store params
		partids = [None]*2
		for procid in xrange(2):  partids[procid] = os.path.join(initdir,"chunk%01d.txt"%procid)
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(initdir,"params-chunk%01d.txt"%procid)
		from random import shuffle
		if(myid == main_node):
			#  split randomly
			params = read_text_row(os.path.join(initdir,"params-chunk0.txt"))
			assert(len(params) == total_stack)
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
			lowpass, falloff, currentres = get_resolution(vol, radi, nnxo, initdir)
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(  line,"Initial resolution %6.2f, filter cut-off %6.4f and fall-off %6.4f"%(currentres,lowpass,falloff))
			write_text_row([lowpass, falloff, currentres],os.path.join(initdir,"current_resolution.txt"))
		else:
			lowpass    = 0.0
			falloff    = 0.0
			currentres = 0.0
		lowpass    = bcast_number_to_all(lowpass, source_node = main_node)
		falloff    = bcast_number_to_all(falloff, source_node = main_node)
		currentres = bcast_number_to_all(currentres, source_node = main_node)
	else:
		if(myid == main_node): [lowpass, falloff, currentres] = read_text_row(os.path.join(initdir,"current_resolution.txt"))[0]
		else:
			lowpass    = 0.0
			falloff    = 0.0
			currentres = 0.0
		lowpass    = bcast_number_to_all(lowpass, source_node = main_node)
		lowpass    = round(lowpass,4)
		falloff    = bcast_number_to_all(falloff, source_node = main_node)
		falloff    = round(falloff,4)
		currentres = bcast_number_to_all(currentres, source_node = main_node)
		currentres = round(currentres,2)

	# set for the first iteration
	nxshrink = min(max(32, int( (lowpass+paramsdict["falloff"]/2.)*2*nnxo + 0.5)), nnxo)
	nxshrink += nxshrink%2
	shrink = float(nxshrink)/nnxo
	nsoft = options.nsoft
	paramsdict["initialfl"] = lowpass
	falloff = 0.2
	tracker = {"resolution":currentres,"lowpass":lowpass, "falloff":falloff, "initialfl":lowpass,  \
				"movedup":False,"eliminated-outliers":False,"PWadjustment":"","local":False,"nsoft":nsoft, \
				"nx":nxshrink, "shrink":shrink, "extension":0.0,"directory":"none"}
	history = [tracker.copy()]
	previousoutputdir = initdir
	#  MAIN ITERATION
	test_outliers = True
	mainiteration = 0
	keepgoing = 1
	paramsdict["xr"] = "3"
	paramsdict["ts"] = "1.0"
	while(keepgoing):
		mainiteration += 1

		#  prepare output directory
		history[-1]["directory"] = "main%03d"%mainiteration
		mainoutputdir = os.path.join(masterdir,history[-1]["directory"])

		if(myid == main_node):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"MAIN ITERATION  #%2d     shrink, nxshrink, currentres, lowpass, falloff "%mainiteration, shrink, nxshrink, currentres, lowpass, falloff)
			print_dict(history[-1],"TRACKER")

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

		# prepare names of input file names, they are in main directory, 
		#   log subdirectories contain outputs from specific refinements
		partids = [None]*2
		for procid in xrange(2):  partids[procid] = os.path.join(previousoutputdir,"chunk%01d.txt"%procid)
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(previousoutputdir,"params-chunk%01d.txt"%procid)

		mpi_barrier(MPI_COMM_WORLD)
		doit = bcast_number_to_all(doit, source_node = main_node)


		#mpi_finalize()
		#exit()

		#print("RACING  A ",myid)
		outvol = [os.path.join(previousoutputdir,"vol%01d.hdf"%procid) for procid in xrange(2)]

		if(myid == main_node):
			if  doit:
				vol = [ get_im(outvol[procid]) for procid in xrange(2) ]
				fuselowf(vol, fq)
				for procid in xrange(2):  vol[procid].write_image(os.path.join(mainoutputdir,"fusevol%01d.hdf"%procid) )
				del vol

		mpi_barrier(MPI_COMM_WORLD)

		#  Refine two groups at a current resolution
		lastring = int(shrink*radi + 0.5)
		if(lastring < 2):
			ERROR( "ERROR!!   lastring too small  %f    %f   %d"%(radi, shrink, lastring), "sxmeridien",1, myid)

		#  REFINEMENT
		#  Part "a"  SHC         <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		for procid in xrange(2):
			coutdir = os.path.join(mainoutputdir,"loga%01d"%procid)
			doit, keepchecking = checkstep(coutdir, keepchecking, myid, main_node)
			#  here ts has different meaning for standard and continuous
			delta = round(degrees(atan(1.0/lastring)), 2)
			subdict( paramsdict, { "delta":"%f"%delta , "an":angular_neighborhood, \
							"lowpass":lowpass, "falloff":falloff, "nsoft":nsoft, \
							"pixercutoff":get_pixercutoff(radi*shrink, delta, 0.5), \
							"delpreviousmax":True, "shrink":shrink, \
							"refvol":os.path.join(mainoutputdir,"fusevol%01d.hdf"%procid) } )
			if(len(history)>1):  old_nx = history[-2]["nx"]
			else:    old_nx = tracker["nx"]
			paramsdict["xr"] = "%s"%max(3,int(1.5*tracker["nx"]/float(old_nx) +0.5))
			if( paramsdict["nsoft"] > 0 ):
				if( float(paramsdict["an"]) == -1.0 ):
					paramsdict["saturatecrit"] = 0.75					
				else:
					paramsdict["saturatecrit"] = 0.90  # Shake and bake for local
				paramsdict["maxit"] = 1500
			else:
				if(paramsdict["local"]):
					paramsdict["saturatecrit"] = 0.95
					paramsdict["xr"] = "2.0"
					paramsdict["maxit"] = 5 #  ?? Lucky guess
				else:
					paramsdict["saturatecrit"] = 0.95
					paramsdict["maxit"] = 50 #  ?? Lucky guess

			if  doit:
				metamove(paramsdict, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)

		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(mainoutputdir, "loga%01d"%procid, "params-chunk%01d.txt"%procid)

		#  check it for the first, if it does not exist, run the program
		outvol = os.path.join(mainoutputdir,"loga0","vol0.hdf")
		doit, keepchecking = checkstep(outvol, keepchecking, myid, main_node)

		if  doit:
			xlowpass, xfalloff, xcurrentres = compute_resolution(stack, mainoutputdir, partids, partstack, radi, nnxo, ali3d_options.CTF, myid, main_node, nproc)
			del xlowpass, xfalloff, xcurrentres
			if( myid == main_node):
				# Move output to proper directories
				for procid in xrange(2):
					cmd = "{} {} {}".format("mv", os.path.join(mainoutputdir,"vol%01d.hdf"%procid), os.path.join(mainoutputdir,"loga%01d"%procid) )
					cmdexecute(cmd)
					cmd = "{} {} {}".format("mv", os.path.join(mainoutputdir,"within-fsc%01d.txt"%procid), os.path.join(mainoutputdir,"loga%01d"%procid) )
					cmdexecute(cmd)
				cmd = "{} {} {}".format("mv", os.path.join(mainoutputdir,"fsc.txt") , os.path.join(mainoutputdir,"afsc.txt"))
				cmdexecute(cmd)
				cmd = "{} {} {}".format("mv", os.path.join(mainoutputdir,"current_resolution.txt") , os.path.join(mainoutputdir,"acurrent_resolution.txt"))
				cmdexecute(cmd)

		#  fuse shc volumes to serve as starting point for the next, deterministic part.
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
				for procid in xrange(2):  vol.append(get_im(os.path.join(mainoutputdir,"loga%01d"%procid,"vol%01d.hdf"%procid) ))
				fuselowf(vol, fq)
				for procid in xrange(2):  vol[procid].write_image( os.path.join(mainoutputdir,"loga%01d"%procid,"fusevol%01d.hdf"%procid) )
				del vol
		else:  doit = 0
		mpi_barrier(MPI_COMM_WORLD)
		doit = bcast_number_to_all(doit, source_node = main_node)

		#  Part "b"  deterministic   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(mainoutputdir,"loga%01d"%procid,"params-chunk%01d.txt"%procid)

		for procid in xrange(2):
			coutdir = os.path.join(mainoutputdir,"logb%01d"%procid)
			doit, keepchecking = checkstep(coutdir, keepchecking, myid, main_node)
			if( nsoft > 0 and doit):  #  Only do finishing up when the previous step was SHC
				#  Run hard to finish up matching
				subdict(paramsdict, \
				{ "maxit":10, "nsoft":0, "saturatecrit":0.95, "delpreviousmax":True, \
				"refvol":os.path.join(mainoutputdir,"loga%01d"%procid,"fusevol%01d.hdf"%procid)} )

				if  doit:
					metamove(paramsdict, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)
			else:
				if( myid == main_node and doit):
					#  Simply copy from loga to logb the necessary stuff
					cmd = "{} {}".format("mkdir",coutdir)
					cmdexecute(cmd)
					cmd = "{} {} {}".format("cp -p", os.path.join(mainoutputdir,"loga%01d"%procid,"params-chunk%01d.txt"%procid) , os.path.join(coutdir,"params-chunk%01d.txt"%procid))
					cmdexecute(cmd)
					
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(mainoutputdir,"logb%01d"%procid,"params-chunk%01d.txt"%procid)

		#  Compute current resolution, store result in the main directory
		doit, keepchecking = checkstep(os.path.join(mainoutputdir,"current_resolution.txt"), keepchecking, myid, main_node)
		if doit:
			if( nsoft > 0 ):
				#  There was first soft phase, so the volumes have to be computed
				#  low-pass filter, current resolution
				lowpass, falloff, currentres = compute_resolution(stack, mainoutputdir, partids, partstack, radi, nnxo, ali3d_options.CTF, myid, main_node, nproc)
			else:
				#  Previous phase was hard, so the volumes exist
				vol = []
				for procid in xrange(2):  vol.append(get_im(os.path.join(mainoutputdir,"loga%01d"%procid,"vol%01d.hdf"%procid) ))
				newlowpass, newfalloff, currentres = compute_resolution(vol, mainoutputdir, partids, partstack, radi, nnxo, ali3d_options.CTF, myid, main_node, nproc)
				del vol
		else:
			if(myid == main_node):
				[newlowpass, newfalloff, currentres] = read_text_row( os.path.join(mainoutputdir,"current_resolution.txt") )[0]
			else:
				newlowpass = 0.0
				newfalloff = 0.0
				currentres = 0.0
			newlowpass = bcast_number_to_all(newlowpass, source_node = main_node)
			newlowpass = round(newlowpass,4)
			newfalloff = bcast_number_to_all(newfalloff, source_node = main_node)
			newfalloff = round(newfalloff,4)
			currentres = bcast_number_to_all(currentres, source_node = main_node)
			currentres = round(currentres,2)

		#  Here I have code to generate presentable results.  IDs and params have to be merged and stored and an overall volume computed.
		doit, keepchecking = checkstep(os.path.join(mainoutputdir,"volf.hdf"), keepchecking, myid, main_node)
		if  doit:
			if( myid == main_node ):
				pinids = map(int, read_text_file(partids[0]) ) + map(int, read_text_file(partids[1]) )
				params = read_text_row(partstack[0]) + read_text_row(partstack[1])

				assert(len(pinids) == len(params))

				for i in xrange(len(pinids)):
					pinids[i] = [ pinids[i], params[i] ]
				del params
				pinids.sort()

				write_text_file([pinids[i][0] for i in xrange(len(pinids))], os.path.join(mainoutputdir,"indexes.txt"))
				write_text_row( [pinids[i][1] for i in xrange(len(pinids))], os.path.join(mainoutputdir,"params.txt"))
			mpi_barrier(MPI_COMM_WORLD)
			ali3d_options.fl = currentres
			ali3d_options.ou = radi
			projdata = getindexdata(stack, os.path.join(mainoutputdir,"indexes.txt"), os.path.join(mainoutputdir,"params.txt"), myid, nproc)
			volf = do_volume(projdata, ali3d_options, mainiteration, mpi_comm = MPI_COMM_WORLD)
			if(myid == main_node): volf.write_image(os.path.join(mainoutputdir,"volf.hdf"))

		mpi_barrier(MPI_COMM_WORLD)

		#print("RACING  X ",myid)
		#  If the resolution stalled, try to eliminate outliers
		if(False):
		#if(currentres == tracker["resolution"]  and test_outliers):
			test_outliers = False

			for procid in xrange(2):
				coutdir = os.path.join(mainoutputdir,"logc%01d"%procid)
				doit, keepchecking = checkstep(coutdir, keepchecking, myid, main_node)

				if  doit:
					#  Do cross-check of the results
					subdict(paramsdict, \
						{ "maxit":1, "lowpass":lowpass, "xfalloff":falloff, "nsoft":0, "saturatecrit":0.95, \
									"delpreviousmax":True, "shrink":shrink, \
									"refvol":os.path.join(mainoutputdir,"vol%01d.hdf"%(1-procid)) } )
					#  The cross-check uses parameters from step "b" to make sure shifts are correct.  
					#  As the check is exhaustive, angles are ignored
					metamove(paramsdict, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)

			# identify bad apples
			doit, keepchecking = checkstep(os.path.join(mainoutputdir,"badapples.txt"), keepchecking, myid, main_node)
			if  doit:
				if(myid == main_node):
					from utilities import get_symt
					from pixel_error import max_3D_pixel_error
					ts = get_symt(ali3d_options.sym)
					badapples = []
					pixercutoff = get_pixercutoff(radi*shrink, float(paramsdict["delta"]), 0.5)
					total_images_now = 0
					for procid in xrange(2):
						bad = []
						ids  = map(int,read_text_file( partids[procid] ))
						total_images_now += len(ids)
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

							if(pixel_error > pixercutoff):
								bad.append(i)
						if(len(bad)>0):
							badapples += [ids[bad[i]] for i in xrange(len(bad))]
							for i in xrange(len(bad)-1,-1,-1):
								del oldp[bad[i]],ids[bad[i]]
						if(len(ids) == 0):
							ERROR("program diverged, all images have large angular errors, most likely the initial model is badly off","sxmeridien",1,myid)
						else:
							#  This generate new parameters, hopefully to be used as starting ones in the new iteration
							write_text_file(ids,os.path.join(mainoutputdir,"chunk%01d.txt"%procid))
							write_text_row(oldp,os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))
					if(len(badapples)>0):
						badapples.sort()
						write_text_file(badapples,os.path.join(mainoutputdir,"badapples.txt"))
						eli = 100*float(len(badapples))/float(total_images_now)
						line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
						print(line,"Elimination of outliers: %5.1f percent"%eli )
					else:  eli = 0.0
					del badapples, oldp,ids,bad,newp,ts
				else:  eli =0.0
				eli = bcast_number_to_all(eli, source_node = main_node)
				
				#  This part under MPI
				if(eli > 0.0):
					#  Compute current resolution
					depfilter, depfalloff, depres = compute_resolution(stack, mainoutputdir, \
							[os.path.join(mainoutputdir,"chunk%01d.txt"%procid) for procid in xrange(2)], \
							[os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid) for procid in xrange(2)], \
							radi, nnxo, ali3d_options.CTF, myid, main_node, nproc)
					#  AGAIN, BASED ON RESOLUTION!!
					if(depres < currentres):
						#  elimination of outliers decreased resolution, ignore the effort
						eliminated_outliers = False
					else:
						eliminated_outliers = True
						currentres = depres
						newlowpass = depfilter
						newfalloff = depfalloff
						"""
						#  It does not seem to be needed, as data is there, we just point to the directory
						for procid in xrange(2):
							#  set pointers to current parameters in main, which are for the reduced set stored above
							partids[procid]   = os.path.join(mainoutputdir,"chunk%01d.txt"%procid
							partstack[procid] = os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid)
						"""
				else:
					eliminated_outliers = False
		else:
			test_outliers = True
			eliminated_outliers = False
			newlowpass = lowpass
			newfalloff = falloff
		#  HERE the lowpass has the true meaning
		lowpass = newlowpass
		falloff = newfalloff

		if(myid == main_node and not eliminated_outliers and doit):  # I had to add here doit, otherwise during the restart it incorrectly copies the files.
			#  It still has a problem with the restart, has to be corrected.
			for procid in xrange(2):
				#  This is standard path, copy parameters to be used to the main
				cmd = "{} {} {}".format("cp -p", partids[procid] , os.path.join(mainoutputdir,"chunk%01d.txt"%procid))
				cmdexecute(cmd)
				cmd = "{} {} {}".format("cp -p", partstack[procid], os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))
				cmdexecute(cmd)

		keepgoing = 0
		if(angular_neighborhood == "-1" ):
			stepforward = 0.05
			increment   = 0.02
		else:
			stepforward = 0.02
			increment   = 0.01

		if(myid == main_node):
			print(" New resolution %6.3f   Previous resolution %6.3f"%(currentres , tracker["resolution"]))

		if( ( currentres > tracker["resolution"] ) or (eliminated_outliers and not tracker["eliminated-outliers"]) or mainiteration == 1):
			if(myid == main_node):
				if( currentres > tracker["resolution"]):  print("  Resolution improved, full steam ahead!")
				else:  print("  While the resolution did not improve, we eliminated outliers so we follow the _resolution_improved_ path.")
			if(currentres >= 0.45 ):
				print(" Resolution exceeded 0.45, i.e., approached Nyquist limit, program will terminate")
			else:
				# We need separate rules for each case
				if( currentres > tracker["resolution"] ):  tracker["movedup"] = True
				else:   tracker["movedup"] = False
				#  Exhaustive searches
				if(angular_neighborhood == "-1" and not tracker["local"]):
					tracker["extension"] = min(stepforward, 0.45 - currentres)  # lowpass cannot exceed 0.45
					paramsdict["initialfl"] = lowpass
					#  For exhaustive searches do sharp cut-off to prevent volumes size to grow too large
					paramsdict["falloff"]   = 0.2
					lowpass = currentres + tracker["extension"]
					shrink = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
					nxshrink = min(int(nnxo*shrink + 0.5),nnxo)
					nxshrink += nxshrink%2
					shrink = float(nxshrink)/nnxo
				#  Local discrete/gridding searches
				elif(angular_neighborhood != "-1" and not tracker["local"]):
					tracker["extension"] = min(stepforward, 0.45 - currentres)  # lowpass cannot exceed 0.45
					paramsdict["initialfl"] = lowpass
					paramsdict["falloff"]   = falloff
					lowpass = currentres# + tracker["extension"]
					shrink = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
					nxshrink = min(int(nnxo*shrink + 0.5),nnxo)
					nxshrink += nxshrink%2
					shrink = float(nxshrink)/nnxo
				#  Local/gridding  searches, move only as much as the resolution increase allows
				elif(tracker["local"]):
					tracker["extension"] =    0.0  # lowpass cannot exceed 0.45
					paramsdict["initialfl"] = lowpass
					paramsdict["falloff"]   = falloff
					lowpass = currentres + tracker["extension"]
					shrink = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
					nxshrink = min(int(nnxo*shrink + 0.5),nnxo)
					nxshrink += nxshrink%2
					shrink = float(nxshrink)/nnxo
				else:
					print(" Unknown combination of settings in improved resolution path",angular_neighborhood,tracker["local"])
					exit()  #  This will crash the program, but the situation is unlikely to occure
				tracker["nx"]                  = nxshrink
				tracker["resolution"]          = currentres
				tracker["lowpass"]             = lowpass
				tracker["falloff"]             = falloff
				tracker["initialfl"]           = paramsdict["initialfl"]
				tracker["eliminated-outliers"] = eliminated_outliers
				bestoutputdir = mainoutputdir
				keepgoing = 1

		elif( currentres < tracker["resolution"] ):
			#  The resolution decreased.  For exhaustive or local, backoff and switch to the next refinement.  For gridding, terminate
			if(not tracker["movedup"] and tracker["extension"] < increment and mainiteration > 1):
				if( angular_neighborhood == "-1" and not tracker["local"]):
					angular_neighborhood = options.an
					ali3d_options.pwreference = options.pwreference
					tracker["PWadjustment"] = ali3d_options.pwreference
					tracker["extension"] = stepforward + increment # so below it will be set to stepforward
					mainoutputdir = bestoutputdir
					keepgoing = 1
					if(myid == main_node):  print("  Switching to local searches with an %s"%angular_neighborhood)
				elif(angular_neighborhood != "-1" and not tracker["local"]):
					tracker["extension"] = stepforward + increment # so below it will be set to stepforward
					mainoutputdir = bestoutputdir
					keepgoing = 1
					if(myid == main_node):  print("  Switching to local/gridding searches")
				else:
					keepgoing = 0
					if(myid == main_node):  print("  Cannot improve resolution, the best result is in the directory %s"%bestoutputdir)
			else:
				if(not tracker["movedup"] and tracker["extension"] > 0.01 and mainiteration > 1):
					if(myid == main_node):  print("  Resolution decreased.  Will decrease target resolution and will fall back on the best so far:  main%03d"%tracker["bestsolution"])
					bestoutputdir = os.path.join(masterdir,"main%03d"%tracker["bestsolution"])
				elif( tracker["movedup"] and tracker["extension"] > 0.01 and mainiteration > 1):
					if(myid == main_node):  print("  Resolution decreased.  Will decrease target resolution and will try starting from previous stage:  main%03d"%(mainiteration - 1))
					bestoutputdir = os.path.join(masterdir,"main%03d"%(mainiteration-1))
				elif( mainiteration == 1):
					if(myid == main_node):  print("  Resolution decreased in the first iteration.  It is expected, not to worry")
					bestoutputdir = mainoutputdir
					tracker["extension"] = increment  # so it will become zero
				else:  # missing something here?
					ERROR(" Should not be here, ERROR 175!", "sxmeridien", 1, myid)
				if( bestoutputdir != mainoutputdir ):
					#  This is the key, we just reset the main to previous, so it will be eventually used as a starting in the next iteration
					mainoutputdir = bestoutputdir
					"""
					#  Set data from the main previous best to the current.
					for procid in xrange(2):
						partids[procid]   = os.path.join(bestoutputdir,"chunk%01d.txt"%procid)
						partstack[procid] = os.path.join(bestoutputdir,"params-chunk%01d.txt"%procid)
				"""
				if(myid == main_node):
					lowpass, falloff, currentres = read_text_row( os.path.join(bestoutputdir,"current_resolution.txt") )[0]
				currentres = bcast_number_to_all(currentres, source_node = main_node)
				currentres = round(currentres,2)
				#lowpass = bcast_number_to_all(lowpass, source_node = main_node)
				#lowpass = round(lowpass,4)
				#falloff = bcast_number_to_all(falloff, source_node = main_node)
				#falloff = round(falloff,4)
				tracker["extension"] -= increment
				lowpass = currentres + tracker["extension"]
				if(mainiteration > 1):
					#  Here to be consistent I would have to know what was shrink for this run
					k = -1
					for i in xrange(len(history)):
						if(history[i]["directory"] == bestoutputdir[-6:]):
							k = i
							break
					if(k == -1):
						print("  something wrong with bestoutputdir")
						exit()
					shrink                   = history[i]["shrink"]
					nxshrink                 = history[i]["nxshrink"]
					paramsdict["initialfl"]  = history[i]["initialfl"]
					paramsdict["falloff"]    = history[i]["falloff"]
					tracker["initialfl"]     = history[i]["initialfl"]
				tracker["resolution"]    = currentres
				tracker["lowpass"]       = lowpass
				tracker["falloff"]       = paramsdict["falloff"]
				tracker["eliminated-outliers"] = eliminated_outliers
				tracker["movedup"] = False
				keepgoing = 1

		elif( currentres == tracker["resolution"] ):
			# We need separate rules for each case
			if(myid == main_node):
				print("  Resolution did not improve, swith to the next move", angular_neighborhood, tracker["local"],currentres,lowpass)
			#  Exhaustive searches
			if(angular_neighborhood == "-1" and not tracker["local"]):
				if(myid == main_node):
					print("  Switching to local searches with an  %s"%options.an)
				falloff = 0.2
				lowpass = currentres
				tracker["extension"] = min(stepforward, 0.45 - currentres)  # lowpass cannot exceed 0.45
				paramsdict["initialfl"] = lowpass
				#  For exhaustive searches do sharp cut-off to prevent volumes size to grow too large
				paramsdict["falloff"]   = falloff
				#lowpass = currentres + tracker["extension"]
				shrink = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
				nxshrink = max( min(int(nnxo*shrink + 0.5),nnxo), nxshrink )  # shrink cannot decrease
				nxshrink += nxshrink%2
				shrink = float(nxshrink)/nnxo
				tracker["local"] = True

				angular_neighborhood = options.an
				ali3d_options.pwreference = options.pwreference
				tracker["PWadjustment"] = ali3d_options.pwreference
				keepgoing = 1
				if(myid == main_node):
					print("  Switching to local searches with an %s"%angular_neighborhood)
					if(tracker["PWadjustment"] != ""):
						print("  Turning on power spectrum adjustment %s"%tracker["PWadjustment"])
			#  Local discrete/gridding searches
			elif(angular_neighborhood != "-1" and not tracker["local"]):
				if(myid == main_node):
					print("  Switching to local searches without an")
				falloff = 0.2
				tracker["extension"] = min(stepforward, 0.45 - currentres)  # lowpass cannot exceed 0.45
				paramsdict["initialfl"] = lowpass
				paramsdict["falloff"]   = falloff
				lowpass = currentres# + tracker["extension"]
				shrink = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
				nxshrink = min(int(nnxo*shrink + 0.5),nnxo)
				nxshrink += nxshrink%2
				shrink = float(nxshrink)/nnxo
				tracker["local"] = True

				angular_neighborhood = options.an
				ali3d_options.pwreference = options.pwreference
				tracker["PWadjustment"] = ali3d_options.pwreference
				keepgoing = 1
				if(myid == main_node):
					print("  Switching to local searches with an %s"%angular_neighborhood)
					if(tracker["PWadjustment"] != ""):
						print("  Turning on power spectrum adjustment %s"%tracker["PWadjustment"])
			#  Local/gridding  searches, move only as much as the resolution increase allows
			elif(tracker["local"]):
				tracker["extension"] =    0.0  # lowpass cannot exceed 0.45
				paramsdict["initialfl"] = lowpass
				paramsdict["falloff"]   = falloff
				lowpass = currentres + tracker["extension"]
				shrink = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
				nxshrink = min(int(nnxo*shrink + 0.5),nnxo)
				nxshrink += nxshrink%2
				shrink = float(nxshrink)/nnxo
			else:
				print(" Unknown combination of settings in improved resolution path",angular_neighborhood,tracker["local"])
				exit()  #  This will crash the program, but the situation is unlikely to occur
			tracker["nx"]                  = nxshrink
			tracker["resolution"]          = currentres
			tracker["lowpass"]             = lowpass
			tracker["falloff"]             = falloff
			tracker["initialfl"]           = paramsdict["initialfl"]
			tracker["eliminated-outliers"] = eliminated_outliers
			bestoutputdir = mainoutputdir
			keepgoing = 1
			"""
			if( tracker["movedup"] ):
				if(myid == main_node):  print("The resolution did not improve. This is look ahead move.  Let's try to relax slightly and hope for the best")
				tracker["extension"]    = min(stepforward,0.45-currentres)
				tracker["movedup"]      = False
				tracker["initialfl"]    = lowpass
				paramsdict["initialfl"] = lowpass
				paramsdict["falloff"]   = falloff
				lowpass = currentres + tracker["extension"]
				shrink  = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
				nxshrink = min(int(nnxo*shrink + 0.5) + tracker["extension"],nnxo)
				nxshrink += nxshrink%2
				shrink = float(nxshrink)/nnxo
				if( tracker["nx"] == nnxo):
					keepgoing = 0
				else:
					tracker["resolution"] = currentres
					tracker["lowpass"]    = lowpass
					tracker["falloff"]    = falloff
					tracker["eliminated-outliers"] = eliminated_outliers
					keepgoing = 1
			else:
				#  The resolution is not moving up.  Check whether this is exhaustive search, 
				#    if yes switch to local searches by activating an, tun on PW adjustment, if given.
				#  
				if( angular_neighborhood == "-1" ):
					angular_neighborhood = options.an
					ali3d_options.pwreference = options.pwreference
					tracker["PWadjustment"] = ali3d_options.pwreference
					tracker["movedup"] = True
					nsoft = 0
					tracker["nsoft"] = 0
					paramsdoct["nsoft"] = 0
					if(myid == main_node):  print("  Switching to local searches with an %s and turning off SHC"%angular_neighborhood)
					keepgoing = 1
				elif(angular_neighborhood > 0.0 ):
					if( not paramsdict["local"] ):
						paramsdict["local"]     = True
						paramsdict["ts"]        = "2.0"
						paramsdict["falloff"]   = falloff
						tracker["local"]        = True
						tracker["movedup"]      = False
						tracker["initialfl"]    = lowpass
						paramsdict["initialfl"] = lowpass
						lowpass = currentres + tracker["extension"]
						shrink  = max(min(2*lowpass + paramsdict["falloff"], 1.0), minshrink)
						nxshrink = min(int(nnxo*shrink + 0.5) + tracker["extension"],nnxo)
						nxshrink += nxshrink%2
						shrink = float(nxshrink)/nnxo
						if( tracker["nx"] == nnxo):
							if(myid == main_node):  print("The resolution did not improve and image we reached the full image size.")
							keepgoing = 0
						else:
							tracker["resolution"] = currentres
							tracker["lowpass"]    = lowpass
							tracker["falloff"]    = falloff
							tracker["eliminated-outliers"] = eliminated_outliers
							tracker["movedup"] = False
							if(myid == main_node):  print("  Switching to local searches")
							#  We have to decrease angular error as these are "continuous" searches
							paramsdict["pixercutoff"] = get_pixercutoff(radi*shrink, degrees(atan(1.0/float(radi*shrink)))/4.0, 0.1)
							keepgoing = 1
					else:
						#  If the resolution did not improve for local, keep current parameters, but increase the image size to full.
						paramsdict["ts"]    = "2.0"
						paramsdict["falloff"]   = falloff
						tracker["local"]  = True
						tracker["movedup"]  = False
						tracker["initialfl"]    = lowpass
						paramsdict["initialfl"] = lowpass
						nxshrink = nnxo
						shrink = float(nxshrink)/nnxo
						if( tracker["nx"] == nnxo):
							if(myid == main_node):  print("The resolution did not improve and image we reached the full image size.")
							keepgoing = 0
						else:
							tracker["resolution"] = currentres
							tracker["lowpass"]    = lowpass
							tracker["falloff"]    = falloff
							tracker["eliminated-outliers"] = eliminated_outliers
							tracker["movedup"] = False
							if(myid == main_node):  print("  Resolution id not improve, do local searches at full size")
							#  We have to decrease angular error as these are "continuous" searches
							paramsdict["pixercutoff"] = get_pixercutoff(radi*shrink, degrees(atan(1.0/float(radi*shrink)))/4.0, 0.1)
							keepgoing = 1
				"""	
#			else:
#				if(myid == main_node):  print("The resolution did not improve.")
#				keepgoing = 0


		if( keepgoing == 1 ):
			if(myid == main_node):
				print("  New shrink and image dimension :",shrink,nxshrink)
				"""
				#  It does not look like it is necessary, we just have to point to the directory as the files should be there.
				#  Will continue, so update the params files
				for procid in xrange(2):
					#  partids ads partstack contain parameters to be used as starting in the next iteration
					if(not os.path.exists(os.path.join(mainoutputdir,"chunk%01d.txt"%procid))):
						cmd = "{} {} {}".format("cp -p", partids[procid] , os.path.join(mainoutputdir,"chunk%01d.txt"%procid))
						cmdexecute(cmd)
					if(not os.path.exists(os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))):
						cmd = "{} {} {}".format("cp -p", partstack[procid], os.path.join(mainoutputdir,"params-chunk%01d.txt"%procid))
						cmdexecute(cmd)
				"""
			previousoutputdir = mainoutputdir
			tracker["shrink"]     = shrink
			tracker["nx"]         = nxshrink
			history.append(tracker.copy())

		else:
			if(myid == main_node):
				print("  Terminating, the best solution is in the directory main%03d"%tracker["bestsolution"])
		mpi_barrier(MPI_COMM_WORLD)

	mpi_finalize()


if __name__=="__main__":
	main()

