#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
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

from global_def import *

def rec2D(  lines, idrange=None, snr=None ):
	"""Perform a 2D reconstruction on a set of 1D lines using nearest neighbouring reverse FFT algorithm.
	   Input: a set of 1D lines
	   Output: a 2D image
	"""
	from EMAN2 import Reconstructors

	assert len(lines) > 0


	size = lines[0].get_xsize();

        if snr is None:
	    params = {"size":size, "npad":4, "ndim":2}
        else: 
	    params = {"size":size, "npad":4, "ndim":2, "snr":snr}
        
	r = Reconstructors.get("nn4", params)
	r.setup()

        if idrange is None:
           idrange = xrange( len(lines) )

	t = Transform()
	for i in idrange:
		r.insert_slice( lines[i], t )

	return r.finish(True)


def insert_slices(reconstructor, proj):
	xforms = [ proj.get_attr("xform.projection") ]
	weights = [ proj.get_attr_default("weight", 1.0) ]
	ixform = 0
	while True:
		ixform += 1
		xform_proj = proj.get_attr_default("xform.projection" + str(ixform), None)
		if xform_proj == None:
			break
		# putting params in a list does not seem to be necessary, one could call reconstructor as one goes.
		xforms.append(xform_proj)
		#weights.append(proj.get_attr_default("weight" + str(ixform), 1.0))
		weights.append(1.0)
	for i in xrange(len(xforms)):
		reconstructor.insert_slice( proj, xforms[i], weights[i] )


def recons3d_4nn(stack_name, list_proj=[], symmetry="c1", npad=4, snr=None, weighting=1, varsnr=True, xysize=-1, zsize = -1):
	"""
	Perform a 3-D reconstruction using Pawel's FFT Back Projection algorithm.
	   
	Input:
	   stack_name - name of the file with projection data.
	   
	   list_proj -  list of projections to be used in the reconstruction

	   symmetry - Point group of the target molecule (defaults to "C1")
	   
	   npad - 

	   Angles and shifts are passed in the file header as set_attr. Keywords are phi, theta, psi, sx, sy

	   Return:  3D reconstructed volume image

	   Usage:
	     vol = recons3d_4nn(filepattern, list_proj, symmetry)
	"""
	import types
	from EMAN2 import Reconstructors

	if list_proj == []:
		if type(stack_name) == types.StringType: nima = EMUtil.get_image_count(stack_name)
		else : nima = len(stack_name)
		list_proj = xrange(nima) 
	# read first image to determine the size to use
	if type(stack_name) == types.StringType:
		proj = EMData()
		proj.read_image(stack_name, list_proj[0])
	else:    proj = stack_name[list_proj[0]].copy()

	size = proj.get_xsize()
	# sanity check -- image must be square
	if size != proj.get_ysize():
		ERROR("input data has to be square","recons3d_4nn",1)

	# reconstructor
	fftvol = EMData()
	weight = EMData()
	params = {"npad":npad, "symmetry":symmetry, "weighting":weighting, "fftvol":fftvol, "weight":weight}
	if ( xysize == -1 and zsize == -1 ):
		if snr is None:
			params["size"] = size
		else:
			params["size"] = size
			params["snr"] = snr
			params["varsnr"] = int(varsnr)
		r = Reconstructors.get("nn4", params)
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/size
			ry = float(xysize)/size
			rz = float(zsize)/size
		elif( xysize != -1):
			rx = float(xysize)/size
			ry = float(xysize)/size
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/size
							
		if snr is None:
			params["sizeprojection"] = size
			params["xratio"] = rx
			params["yratio"] = ry
			params["zratio"] = rz
		else:
			params["sizeprojection"] = size
			params["snr"] = snr
			params["varsnr"] = int(varsnr)
			params["xratio"] = rx
			params["yratio"] = ry
			params["zratio"] = rz	
		r = Reconstructors.get("nn4_rect", params)
		
	r.setup()

	if type(stack_name) == types.StringType:
		for i in xrange(len(list_proj)):
			proj.read_image(stack_name, list_proj[i])
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = proj.get_attr_default('active', 1)
			# if active == 1:
			# 	insert_slices(r, proj)
			insert_slices(r, proj)
	else:
		for i in list_proj:
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = stack_name[i].get_attr_default('active', 1)
			# if active == 1:
			# 	insert_slices(r, stack_name[i])
			insert_slices(r, stack_name[i])

	dummy = r.finish(True)
	return fftvol


def recons3d_4nn_MPI(myid, prjlist, symmetry="c1", info=None, npad=2, xysize=-1, zsize=-1, mpi_comm=None):
	from utilities  import reduce_EMData_to_root, pad
	from EMAN2      import Reconstructors
	from utilities  import iterImagesList
	from mpi        import MPI_COMM_WORLD
	import types

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if type(prjlist) == types.ListType:
		prjlist = iterImagesList(prjlist)

	if not prjlist.goToNext():
		ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist.image().get_xsize()
	if prjlist.image().get_ysize() != imgsize:
		imgsize = max(imgsize, prjlist.image().get_ysize())
		dopad = True
	else:
		dopad = False
	prjlist.goToPrev()

	fftvol = EMData()		
	weight = EMData()
	if (xysize == -1 and zsize == -1 ):
		params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
		r = Reconstructors.get( "nn4", params )
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = float(zsize)/imgsize
		elif( xysize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/imgsize
		params = {"sizeprojection":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol,"weight":weight,"xratio":rx,"yratio":ry,"zratio":rz}
		r = Reconstructors.get( "nn4_rect", params )
	r.setup()

	if not (info is None): nimg = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if(active == 1):
		# 	if dopad:
		# 		prj = pad(prj, imgsize,imgsize, 1, "circumference")
		# 	insert_slices(r, prj)
		# 	if( not (info is None) ):
		# 		nimg += 1
		# 		info.write("Image %4d inserted.\n" %(nimg) )
		# 		info.flush()

		if dopad:
			prj = pad(prj, imgsize,imgsize, 1, "circumference")
		insert_slices(r, prj)
		if( not (info is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

	if not (info is None): 
		info.write( "Begin reducing ...\n" )
		info.flush()

	reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if myid == 0:  dummy = r.finish(True)
	else:
		from utilities import model_blank
		if ( xysize == -1 and zsize == -1 ):
			fftvol = model_blank(imgsize, imgsize, imgsize)
		else:
			if zsize == -1:
				fftvol = model_blank(xysize, xysize, imgsize)
			elif xysize == -1:
				fftvol = model_blank(imgsize, imgsize, zsize)
			else:
				fftvol = model_blank(xysize, xysize, zsize)
	return fftvol

# secondrun
def secondrunrecons3d_4nnw_MPI(myid, prjlist, prevol, symmetry="c1", finfo=None, npad=2, mpi_comm=None):
	from utilities     import reduce_EMData_to_root, pad, get_params_proj
	from EMAN2         import Reconstructors
	from utilities     import iterImagesList, model_blank, model_circle, reshape_1d, read_text_file
	from fundamentals  import fft, rops
	from mpi           import MPI_COMM_WORLD
	import types

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if type(prjlist) == types.ListType:
		prjlist = iterImagesList(prjlist)

	if not prjlist.goToNext():
		ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist.image().get_xsize()
	bigsize = imgsize*npad
	bnx     = bigsize//2+1

	prjlist.goToPrev()

	fftvol = EMData()
	weight = EMData()
	print "   NEW"
	#t = read_text_file('fromrun8model.txt',4)
	#  GET FSC
	t = read_text_file('data_model.txt',4)
	from math import sqrt
	for i in xrange(len(t)):
		t[i] = max(t[i],0.0)
		#  This is what is used to get the SSNR
		t[i] = sqrt(2*t[i]/(1.0+t[i]))
	t = reshape_1d(t,len(t),npad*len(t))
	refvol = model_blank(2*bnx,1,1,0.5)
	for i in xrange(len(t)):  refvol.set_value_at(i,t[i])
	"""
	from math import tanh,pi
	fl = 0.15
	aa = 0.15
	for i in xrange(bnx):
		r = float(i)/bigsize
		refvol.set_value_at(i, 0.5*( tanh(pi*(r+fl)/2.0/fl/aa) - tanh(pi*(r-fl)/2.0/2.0/fl/aa) ) )
		print "  FILTER  ",i,refvol.get_value_at(i)
	"""
	#print " DONE refvol"
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "refvol":refvol, "weight":weight, "weighting":0, "snr":1.0}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()

	from projection import prep_vol, prgs
	from filter import filt_ctf
	#volft,kb = prep_vol(prevol)

	#mask2d = model_circle(imgsize//2-2, imgsize,imgsize)
	#maskbi = model_circle(imgsize//2-2, bigsize,bigsize)
	#  noise model of 2D data.
	models = [None]
	for ml in xrange(len(models)):
		temp = read_text_file('sigma2.txt',2)
		temp = reshape_1d(temp, len(temp), 2*len(temp))
		models[ml] = model_blank(len(temp)+10)
		for lm in xrange(len(temp)):  models[ml].set_value_at(lm,1.0/(temp[lm]*4*imgsize**2/npad))
		#from sys import exit
		#print "metadata/model-%04d.txt"%groupkeys[1][ml]
		#for lm in xrange(len(temp)):  print  lm,models[ml].get_value_at(lm)
		#exit()


	if not (finfo is None): nimg = 0
	ll = 0
	while prjlist.goToNext():
		prj = prjlist.image()

		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if(active == 1):
		if ll%100 == 0:  print "  moved  ",ll
		ll +=1
		prj.set_attr("sigmasq2", models[0])
		#if ll == 0:
		#	write_text_file([range(bigsize),[pqdif[i] for i in xrange(bigsize)] ],"pqdif.txt")
		#	ll+=1
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

		if ll%100 == 0:  print "  moved  ",ll
		ll +=1
		prj.set_attr("sigmasq2", models[0])
		#if ll == 0:
		#	write_text_file([range(bigsize),[pqdif[i] for i in xrange(bigsize)] ],"pqdif.txt")
		#	ll+=1
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

	if not (finfo is None):
		info.write( "Begin reducing ...\n" )
		info.flush()

	#del qdif, pqdif, mask2d, maskbi

	reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if myid == 0:
		print  "  STARTING FINISH"
		dummy = r.finish(True)
	else:
		from utilities import model_blank
		fftvol = model_blank(imgsize, imgsize, imgsize)
	return fftvol

#chc5
def recons3d_4nnw_MPI(myid, prjlist, prevol, symmetry="c1", finfo=None, npad=2, mpi_comm=None):
	from utilities     import reduce_EMData_to_root, pad, get_params_proj
	from EMAN2         import Reconstructors
	from utilities     import iterImagesList, model_blank, model_circle, reshape_1d, read_text_file
	from fundamentals  import fft, rops
	from mpi           import MPI_COMM_WORLD
	import types

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if type(prjlist) == types.ListType:
		prjlist = iterImagesList(prjlist)

	if not prjlist.goToNext():
		ERROR("empty input list","recons3d_4nn_MPI",1)

	imgsize = prjlist.image().get_xsize()
	bigsize = imgsize*npad
	bnx     = bigsize//2+1

	prjlist.goToPrev()

	fftvol = EMData()
	weight = EMData()
	"""
	if myid == 0:
		model_blank(bnx, bigsize, bigsize)
		temp = fft(pad(prevol,bigsize,bigsize,bigsize,0.0))
		temp.set_attr("is_complex",0)
		st = 0.5/(bigsize*bigsize)
		for kk in xrange(bigsize):
			for jj in xrange(bigsize):
				for ii in xrange(0,bnx,2):
					#print ii,jj,kk,temp.get_value_at(ii,jj,kk), temp.get_value_at(ii,jj,kk+1)
					refvol.set_value_at_fast(ii//2,jj,kk,st/((temp.get_value_at(ii,jj,kk))**2+(temp.get_value_at(ii+1,jj,kk))**2) )
					#refvol.set_value_at_fast(ii//2,jj,kk,1.0 )
		refvol.set_value_at_fast(0,0,0,0.0)
		del temp

		st = rops(pad(prevol,bigsize,bigsize,bigsize,0.0))*(bigsize**6)/4.
		from utilities import info
		from utilities import write_text_file
		#zizi = [st.get_value_at(i) for i in xrange(st.get_xsize())]
		#for i in xrange(st.get_xsize()):  st.set_value_at(i,1.0)#/st.get_value_at(i))
		#info(st,None,"refvol")
		refvol = model_blank(bigsize,1,1,1.0)
		for i in xrange(st.get_xsize()):  refvol.set_value_at(i,1.0/(211*st.get_value_at(i)))
	else:  refvol = EMData()
	"""
	print "   NEW"
	#t = read_text_file('fromrun8model.txt',4)
	t = read_text_file('../for-pawel/fsc-relion.txt',1)
	from math import sqrt
	for i in xrange(len(t)):
		t[i] = max(t[i],0.0)
		t[i] = sqrt(2*t[i]/(1.0+t[i]))
	t = reshape_1d(t,len(t),npad*len(t))
	refvol = model_blank(2*bnx,1,1,0.5)
	#for i in xrange(len(t)):  refvol.set_value_at(i,t[i])
	"""
	from math import tanh,pi
	fl = 0.15
	aa = 0.15
	for i in xrange(bnx):
		r = float(i)/bigsize
		refvol.set_value_at(i, 0.5*( tanh(pi*(r+fl)/2.0/fl/aa) - tanh(pi*(r-fl)/2.0/2.0/fl/aa) ) )
		print "  FILTER  ",i,refvol.get_value_at(i)
	"""
	#print " DONE refvol"
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "refvol":refvol, "weight":weight, "weighting":0, "snr":1.0}
	r = Reconstructors.get( "nn4_ctfw", params )
	r.setup()

	from projection import prep_vol, prgs
	from filter import filt_ctf
	#volft,kb = prep_vol(prevol)

	#mask2d = model_circle(imgsize//2-2, imgsize,imgsize)
	#maskbi = model_circle(imgsize//2-2, bigsize,bigsize)

	groupkeys = read_text_file("groupkeys.txt",-1)
	for ml in xrange(3):  groupkeys[ml] = map(int, groupkeys[ml])
	models = [None]*len(groupkeys[0])
	for ml in xrange(len(models)):
		temp = read_text_file("metadata/model-%04d.txt"%groupkeys[1][ml],-1)
		temp = reshape_1d(temp[2], len(temp[0]), 2*len(temp[0]))
		models[ml] = model_blank(len(temp)+10)
		for lm in xrange(len(temp)):  models[ml].set_value_at(lm,1.0/(temp[lm]*4*imgsize**2/npad))
		#from sys import exit
		#print "metadata/model-%04d.txt"%groupkeys[1][ml]
		#for lm in xrange(len(temp)):  print  lm,models[ml].get_value_at(lm)
		#exit()


	if not (finfo is None): nimg = 0
	ll = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		if ll%100 == 0:  print "  moved  ",ll
		ll +=1
		ml = prj.get_attr('groupindex')#int(prj.get_attr('data_path')[4:8])
		prj.set_attr("sigmasq2", models[groupkeys[1].index(ml)])
		insert_slices(r, prj)
		if( not (finfo is None) ):
			nimg += 1
			info.write("Image %4d inserted.\n" %(nimg) )
			info.flush()

	if not (finfo is None):
		info.write( "Begin reducing ...\n" )
		info.flush()

	#del qdif, pqdif, mask2d, maskbi

	reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if myid == 0:
		print  "  STARTING FINISH"
		dummy = r.finish(True)
	else:
		from utilities import model_blank
		fftvol = model_blank(imgsize, imgsize, imgsize)
	return fftvol


def recons3d_4nn_ctf(stack_name, list_proj = [], snr = 1.0, sign=1, symmetry="c1", verbose=0, npad=2, xysize = -1, zsize = -1 ):
	"""Perform a 3-D reconstruction using Pawel's FFT Back Projection algoritm.
	   
	   Input:
	    stack_name - name of the stack file on a disk,
	                 each image has to have the following attributes set:
			 psi, theta, phi, sx, sy, defocus, 
	    list_proj - list of images from stack_name to be included in the reconstruction
	    symmetry	 -- Point group of the target molecule (defaults to "C1")

	   Return:  3d reconstructed volume image

	   Usage:
	     
	     anglelist = getAngles("myangles.txt") # not yet written
	     vol = do_reconstruction(filepattern, start, end, anglelist, symmetry)
	"""
	import types
	from EMAN2     import Reconstructors
	from utilities import pad

	# read first image to determine the size to use
	if list_proj == []:	
		if type(stack_name) == types.StringType: nima = EMUtil.get_image_count(stack_name)
		else : nima = len(stack_name)
		list_proj = xrange(nima) 
	# read first image to determine the size to use
	if type(stack_name) == types.StringType:
		proj = EMData()
		proj.read_image(stack_name, list_proj[0])
	else:    proj = stack_name[list_proj[0]].copy()

	# convert angles to transform (rotation) objects
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# active = proj.get_attr_default('active', 1)
	size   = proj.get_xsize()
	if proj.get_ysize() != size:
		size = max(size, proj.get_ysize())
		dopad = True
	else:
		dopad = False

	# reconstructor
	fftvol = EMData()
	weight = EMData()
	params = {"npad":npad, "symmetry":symmetry, "snr":snr, "sign":sign, "fftvol":fftvol, "weight":weight}
	if ( xysize == -1 and zsize == -1 ):
		params["size"] = size
		r = Reconstructors.get("nn4_ctf", params)
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/size
			ry = float(xysize)/size
			rz = float(zsize)/size
		elif( xysize != -1):
			rx = float(xysize)/size
			ry = float(xysize)/size
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/size

		params["sizeprojection"] = size
		params["xratio"] = rx
		params["yratio"] = ry
		params["zratio"] = rz
		r = Reconstructors.get("nn4_ctf_rect", params)
	r.setup()

	if type(stack_name) == types.StringType:
		for i in xrange(len(list_proj)):
			proj.read_image(stack_name, list_proj[i])
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = proj.get_attr_default('active', 1)
			# if(active == 1):
			# 	if dopad: 
			# 		proj = pad(proj, size, size, 1, "circumference")
			# 	insert_slices(r, proj)
			if dopad: 
				proj = pad(proj, size, size, 1, "circumference")
			insert_slices(r, proj)
	else:
		for i in xrange(len(list_proj)):
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# active = stack_name[list_proj[i]].get_attr_default('active', 1)
			# if active == 1:
			# 	insert_slices(r, stack_name[list_proj[i]])
			insert_slices(r, stack_name[list_proj[i]])
	dummy = r.finish(True)
	return fftvol


def recons3d_4nn_ctf_MPI(myid, prjlist, snr = 1.0, sign=1, symmetry="c1", info=None, npad=2, xysize=-1, zsize=-1, mpi_comm=None, smearstep = 0.0):
	"""
		recons3d_4nn_ctf - calculate CTF-corrected 3-D reconstruction from a set of projections using three Eulerian angles, two shifts, and CTF settings for each projeciton image
		Input
			stack: name of the stack file containing projection data, projections have to be squares
			list_proj: list of projections to be included in the reconstruction or image iterator
			snr: Signal-to-Noise Ratio of the data 
			sign: sign of the CTF 
			symmetry: point-group symmetry to be enforced, each projection will enter the reconstruction in all symmetry-related directions.
	"""
	from utilities  import reduce_EMData_to_root, pad
	from EMAN2      import Reconstructors
	from utilities  import iterImagesList, set_params_proj
	from mpi        import MPI_COMM_WORLD
	import types

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if type(prjlist) == types.ListType:
		prjlist = iterImagesList(prjlist)
	if not prjlist.goToNext():
		ERROR("empty input list","recons3d_4nn_ctf_MPI",1)
	imgsize = prjlist.image().get_xsize()
	if prjlist.image().get_ysize() != imgsize:
		imgsize = max(imgsize, prjlist.image().get_ysize())
		dopad = True
	else:
		dopad = False
	prjlist.goToPrev()

	fftvol = EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in xrange(-ns,ns+1):
			if( j != 0):
				for i in xrange(-ns,ns+1):
					for k in xrange(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in xrange(-ns,ns+1):
			for k in xrange(-ns,ns+1):
				prj.append(i+k)
		for i in xrange(-2*ns,2*ns+1,1):
			 smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol.set_attr("smear", smear)

	weight = EMData()
	if (xysize == -1 and zsize == -1 ):
		params = {"size":imgsize, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
		r = Reconstructors.get( "nn4_ctf", params )
	else:
		if ( xysize != -1 and zsize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = float(zsize)/imgsize
		elif( xysize != -1):
			rx = float(xysize)/imgsize
			ry = float(xysize)/imgsize
			rz = 1.0
		else:
			rx = 1.0
			ry = 1.0
			rz = float(zsize)/imgsize
		#  There is an error here with sizeprojection  PAP 10/22/2014
		params = {"size":sizeprojection, "npad":npad, "snr":snr, "sign":sign, "symmetry":symmetry, "fftvol":fftvol, "weight":weight,"xratio":rx,"yratio":ry,"zratio":rz}
		r = Reconstructors.get( "nn4_ctf_rect", params )
	r.setup()

	#if not (info is None):
	nimg = 0
	while prjlist.goToNext():
		prj = prjlist.image()
		if dopad:
			prj = pad(prj, imgsize, imgsize, 1, "circumference")
		#if params:
		insert_slices(r, prj)
		if not (info is None):
			nimg += 1
			info.write(" %4d inserted\n" %(nimg) )
			info.flush()
	del pad
	if not (info is None): 
		info.write( "begin reduce\n" )
		info.flush()

	reduce_EMData_to_root(fftvol, myid, comm=mpi_comm)
	reduce_EMData_to_root(weight, myid, comm=mpi_comm)

	if not (info is None): 
		info.write( "after reduce\n" )
		info.flush()

	if myid == 0 :
		dummy = r.finish(True)
	else:
		from utilities import model_blank
		if ( xysize == -1 and zsize == -1 ):
			fftvol = model_blank(imgsize, imgsize, imgsize)
		else:
			if zsize == -1:
				fftvol = model_blank(xysize, xysize, imgsize)
			elif xysize == -1:
				fftvol = model_blank(imgsize, imgsize, zsize)
			else:
				fftvol = model_blank(xysize, xysize, zsize)
	return fftvol


def recons3d_nn_SSNR(stack_name,  mask2D = None, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0):

	"""
	Perform a 3-D reconstruction using nearest neighbor interpolation and 
	calculate 3D spectral signal-to-noise ratio (SSNR)	   
	Input : stack_name - Name of the file with projection data.
		CTF        - 
	        symmetry   - Point group of the target molecule (default "c1")
		npad       - Times of padding applied, default is 1
		sign       - Currently not used, may be used in the future
		w          - The thickness of the shell, default is 1
		filename   - The filename in which you can save the SSNR results 
	Return: reconstructed 3D SSNR volume
        Usage : vol = recons3d_nn_SSNR(stack_name, CTF, symmetry, npad, snr, sign, w, filename])
	CTF true:
	variance at one voxel  = Gamma^2d->3d [ |F_k^2D|^2   +  ctf^2*|P^2D->3D(F^3D)|^2 -
	          -2*Real(conj(F_k^2D)*ctf*P^2D->3D(F^3D))]	
	signal  at one voxel   = Gamma^2d->3d [ |F_k^2D|^2  ]
	SSNR =  sum_rot [ wght*signal/Kn ]/sum_rot[ wght*variance /(Kn(Kn-1))] -1
	Notice: wght is always turned on during SSNR calculation.
	"""
	import types
	from EMAN2 import Reconstructors

	# Yang add a safety on 05/22/07
	if type(stack_name) == types.StringType: nima = EMUtil.get_image_count(stack_name)
	else :                                   nima = len(stack_name)
	# read first image to determine the size to use
	if type(stack_name) == types.StringType:
		proj = EMData()
		proj.read_image(stack_name, 0)
	else:    
		proj = stack_name[0].copy()
	#active = proj.get_attr('active')
	size   = proj.get_xsize()
	# sanity check -- image must be square
	if size != proj.get_ysize(): ERROR("input data has to be square","recons3d_nn_SSNR",1)
	# reconstructor
	SSNR = EMData()
	fftvol = EMData()
	weight = EMData()
	weight2 = EMData()
	vol_ssnr = EMData()
	params = {"size":size, "npad":npad, "symmetry":symmetry, "SSNR":SSNR, "w":ring_width, "fftvol":fftvol, "weight":weight, "weight2":weight2, "vol_ssnr":vol_ssnr}
	if CTF:
		weight3 = EMData()
		params["sign"] = sign
		params["weight3"] = weight3
		r = Reconstructors.get("nnSSNR_ctf", params)
	else:
		r = Reconstructors.get("nnSSNR", params)
	r.setup()

	for i in xrange(nima):
		if type(stack_name) == types.StringType:
			proj.read_image(stack_name, i)
		else:
			proj = stack_name[i]
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = proj.get_attr_default('active', 1)
		# if(active == 1):
		if(random_angles  == 2):
			from  random import  random
			phi    = 360.0*random()
			theta  = 180.0*random()
			psi    = 360.0*random()
			xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
		elif(random_angles  == 3):
			from  random import  random
			phi    = 360.0*random()
			theta  = 180.0*random()
			psi    = 360.0*random()
			tx     = 6.0*(random() - 0.5)
			ty     = 6.0*(random() - 0.5)
			xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
		elif(random_angles  == 1):
			from  random import  random
			old_xform_proj = proj.get_attr( "xform.projection" )
			dict = old_xform_proj.get_rotation( "spider" )
			dict["psi"] = 360.0*random()
			xform_proj = Transform( dict )
		else:
			xform_proj = proj.get_attr( "xform.projection" )
	
		if mask2D:
			stats = Util.infomask(proj, mask2D, True)
			proj -= stats[0]
			proj *= mask2D
		r.insert_slice(proj, xform_proj)
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1  END

	dummy = r.finish(True)
	outlist = [[] for i in xrange(6)]
	nn = SSNR.get_xsize()
	for i in xrange(1,nn): outlist[0].append((float(i)-0.5)/(float(nn-1)*2))
	for i in xrange(1,nn):
		if(SSNR(i,1,0) > 0.0):
			outlist[1].append(max(0.0,(SSNR(i,0,0)/SSNR(i,1,0)-1.)))      # SSNR
		else:
			outlist[1].append(0.0)
	for i in xrange(1,nn): outlist[2].append(SSNR(i,1,0)/SSNR(i,2,0))	          # variance
	for i in xrange(1,nn): outlist[3].append(SSNR(i,2,0))				  # number of points in the shell
	for i in xrange(1,nn): outlist[4].append(SSNR(i,3,0))				  # number of added Fourier points
	for i in xrange(1,nn): outlist[5].append(SSNR(i,0,0))				  # square of signal
	return [outlist, vol_ssnr]

def recons3d_nn_SSNR_MPI(myid, prjlist, mask2D, ring_width=1, npad =1, sign=1, symmetry="c1", CTF = False, random_angles = 0, mpi_comm = None):
	from utilities import reduce_EMData_to_root
	from EMAN2 import Reconstructors
	from mpi import MPI_COMM_WORLD

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	if( len(prjlist) == 0 ):    ERROR("empty input list","recons3d_nn_SSNR_MPI",1)
	imgsize = prjlist[0].get_xsize()
	if prjlist[0].get_ysize() != imgsize:  ERROR("input data has to be square","recons3d_nn_SSNR_MPI",1)
	fftvol   = EMData()
	weight   = EMData()
	weight2  = EMData()
	SSNR     = EMData()
	vol_ssnr = EMData()
	params = {"size":imgsize, "npad":npad, "symmetry":symmetry, "SSNR":SSNR, "fftvol":fftvol, "weight":weight, "weight2":weight2, "vol_ssnr":vol_ssnr, "w":ring_width }
	if CTF:
		weight3  = EMData()
		params["sign"] = sign
		params["weight3"] = weight3
		r = Reconstructors.get("nnSSNR_ctf", params)
	else:
		r = Reconstructors.get("nnSSNR", params)
	r.setup()

	if prjlist[0].get_xsize() != imgsize or prjlist[0].get_ysize() != imgsize: ERROR("inconsistent image size","recons3d_nn_SSNR_MPI",1)
	for prj in prjlist:
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
		# active = prj.get_attr_default('active', 1)
		# if active == 1:
		if random_angles  == 2:
			from  random import  random
			phi	 = 360.0*random()
			theta    = 180.0*random()
			psi	 = 360.0*random()
			xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi} )
		elif random_angles  == 3:
			from  random import  random
			phi    = 360.0*random()
			theta  = 180.0*random()
			psi    = 360.0*random()
			tx     = 6.0*(random() - 0.5)
			ty     = 6.0*(random() - 0.5)
			xform_proj = Transform( {"type":"spider", "phi":phi, "theta":theta, "psi":psi, "tx":tx, "ty":ty} )
		elif random_angles  == 1:
			from  random import  random
			old_xform_proj = prj.get_attr( "xform.projection" )
			dict = old_xform_proj.get_rotation( "spider" )
			dict["psi"] = 360.0*random()
			xform_proj = Transform( dict )
		else:
			xform_proj = prj.get_attr( "xform.projection" )
		if mask2D:
			stats = Util.infomask(prj, mask2D, True)
			prj -= stats[0]
			prj *= mask2D
		r.insert_slice(prj, xform_proj )
		# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1 END
		
	#from utilities import info
	reduce_EMData_to_root(weight,  myid, 0, comm=mpi_comm)
	reduce_EMData_to_root(fftvol,  myid, 0, comm=mpi_comm)
	reduce_EMData_to_root(weight2, myid, 0, comm=mpi_comm)
	if CTF:
		reduce_EMData_to_root(weight3, myid, 0, comm=mpi_comm)
	if myid == 0 :
		dummy = r.finish(True)		
		outlist = [[] for i in xrange(6)]
		nn = SSNR.get_xsize()
		for i in xrange(1,nn): outlist[0].append((float(i)-0.5)/(float(nn-1)*2))
		for i in xrange(1,nn):
			if SSNR(i,1,0) > 0.0:
				outlist[1].append(max(0.0,(SSNR(i,0,0)/SSNR(i,1,0)-1.)))     # SSNR
			else:
				outlist[1].append(0.0)
		for i in xrange(1,nn): 
			if SSNR(i,2,0) > 0.0:
				outlist[2].append(SSNR(i,1,0)/SSNR(i,2,0))	          # variance
			else:
				outlist[2].append(0.0)
		for i in xrange(1,nn): outlist[3].append(SSNR(i,2,0))				  # number of points in the shell
		for i in xrange(1,nn): outlist[4].append(SSNR(i,3,0))				  # number of added Fourier points
		for i in xrange(1,nn): outlist[5].append(SSNR(i,0,0))				  # square of signal
		return [outlist, vol_ssnr]


class memory_store:
    def __init__(self, npad):
        self.m_npad = npad
        self.m_imgs = []

    def add_image(self, img):
        self.m_imgs.append(img)
    
    def get_image(self, id):
        return self.m_imgs[id]

def bootstrap_nn(proj_stack, volume_stack, list_proj, niter, media="memory", npad=4, symmetry="c1", output=-1, CTF=False, snr=1.0, sign=1, myseed=None ):
	from random import seed
	from random import randint
	from time   import time
	from sys    import stdout
	from utilities import set_ctf
	from EMAN2 import Reconstructors

	if(output == -1):
		import sys
		output=sys.stdout

	if not(myseed is None):
		seed(myseed) 

	nimages = len(list_proj)
	if nimages == 0 :
		print "empty list of projections input!"
		return None


	if media=="memory" :
		store = memory_store(npad)
	else :
		store = file_store(media,npad, 0, CTF)
		if not(output is None):
			output.flush()

	proj = EMData()
	proj.read_image(proj_stack,list_proj[0])

	size = proj.get_xsize()
	if size != proj.get_ysize():
		print "Image projections must be square!"
		return None

	#for j in xrange(nimages):
	#    proj = EMData()
	#    proj.read_image(proj_stack, j)
        #store.add_image( proj )

        #if (j+1) % 1000  == 0 :
        #    output.write( "%d image read\n" % (j+1) )
        #    output.flush()

	overall_start = time()
	for i in xrange(niter):
		iter_start = time()
		mults = nimages*[0]
		for j in xrange(nimages):
			imgid = randint(0,nimages-1)
			mults[imgid]=mults[imgid]+1

		if CTF:
			params = {"size":size, "npad":npad, "symmetry":symmetry, "snr":snr, "sign":sign}
			r = Reconstructors.get("nn4_ctf", params);
		else:
			params = {"size":size, "npad":npad, "symmetry":symmetry, "snr":snr}
			r = Reconstructors.get("nn4", params);

		r.setup()
		if not(output is None):
			output.write( "Bootstrap volume %8d " % i )
			output.flush()

		store.restart()

		if not(output is None):
			output.write( "Inserting images " )
			output.flush()

		for j in xrange(nimages):
			if mults[j] > 0 :
				img_j = EMData()
				store.get_image( j, img_j );
				phi_j = img_j.get_attr( "phi" )
				tht_j = img_j.get_attr( "theta" )
				psi_j = img_j.get_attr( "psi" )
				tra_j = Transform( {"type":"spider", "phi":phi_j, "theta":tht_j, "psi":psi_j} )

				if CTF:
					cs = img_j.get_attr( "Cs" )
					pixel   = img_j.get_attr( "Pixel_size" )
					defocus = img_j.get_attr( "defocus" )
					voltage = img_j.get_attr( "voltage" )
					ampcont = img_j.get_attr( "amp_contrast" )
					bfactor = 0.0
					set_ctf( img_j, [defocus, cs, voltage, pixel, bfactor, ampcont] )

				r.insert_slice(img_j, tra_j, mults[j])

        			#[mean,sigma,min,max]= Util.infomask(img_j, None, False)
        			#output.write( "img %4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % (j, mean, sigma, min, max, phi, theta, psi) )
        			#output.flush()

		if not(output is None):
			output.write( "Finishing... " )
			output.flush( )

		vol = r.finish(True)

		if not(output is None):
			output.write( "Writing... " )
			output.flush()

		vol.write_image(volume_stack,i)

		if not(output is None):
			output.write( " done!" )
			output.write( " time %15.3f %15.3f \n" % (time()-iter_start,time()-overall_start) )
			output.flush()


def recons3d_em(projections_stack, max_iterations_count = 100, radius = -1, min_avg_abs_voxel_change = 0.01, use_weights = False, symmetry = "c1"):
	"""
	Reconstruction algorithm basing on the Expectation Maximization method
		projections_stack            -- file or list with projections
		max_iterations_count         -- stop criterion 
		min_avg_abs_voxel_change     -- stop criterion 
		use_weights                  -- true == multiply projections by extra weights
		symmetry                     -- type of symmetry
	#
	"""
	from time import clock
	from utilities import model_blank, model_circle, model_square
	from morphology import threshold_to_minval
	import types
	min_allowed_divisor = 0.0001
	
	if type(projections_stack) is types.StringType:
		projections = EMData.read_images(projections_stack)
	else:
		projections = projections_stack
	
	if len(projections) == 0:
		ERROR("Stack of projections cannot be empty", "recons3d_em")
	
	nx = projections[0].get_xsize()
	if (projections[0].get_ysize() != nx) or (projections[0].get_zsize() != 1):
		ERROR("This procedure works only for square images", "recons3d_em")

	if radius < 0:  radius = nx // 2 - 1
	sphere2D = model_circle(radius, nx, nx)   
	sphere3D = model_circle(radius, nx, nx, nx)
	solution = model_blank(nx, nx, nx)
	a = model_blank(nx, nx, nx) # normalization volume
	e2D = model_square(nx, nx, nx)
	sphere3D_volume = model_blank(nx,nx,nx).cmp("lod",sphere3D,{"negative":0,"normalize":0})
	#print "Parameters:  size=%d  radius=%d  projections_count=%d  max_iterations_count=%d min_avg_abs_voxel_change=%f" % (
	#					nx, radius, len(projections), max_iterations_count, min_avg_abs_voxel_change )
	
	# ----- create initial solution, calculate weights and normalization image (a)
	projections_angles = []  # list of lists of angles
	projections_data   = []  # list of lists of projections' images with weights
	for proj in projections:
		angles = [] # list of angles
		data = []   # list of projections' images with weights
		RA = proj.get_attr( "xform.projection" )
		Util.mul_img( proj, sphere2D )
		for j in range(RA.get_nsym(symmetry)):
			angdict = RA.get_sym(symmetry,j).get_rotation("spider") 
			angles.append( [angdict["phi"], angdict["theta"], angdict["psi"]] )
			chao_params = {"anglelist":angles[j],"radius":radius}
			Util.add_img( solution, proj.backproject("chao", chao_params) )
			Util.add_img( a, e2D.backproject("chao", chao_params) )
			if use_weights:
				proj3Dsphere = sphere3D.project("chao", chao_params)
				Util.mul_scalar( proj3Dsphere, 1.0 / Util.infomask(proj3Dsphere, None, True)[3] )
				Util.mul_img( proj, proj3Dsphere )
			data.append(proj)
		projections_angles.append(angles)
		projections_data.append(data)
	a = threshold_to_minval(a, min_allowed_divisor)  # make sure that voxels' values are not too small (image a is divisior)
	Util.mul_img( solution, sphere3D )
	Util.div_img( solution, a )
	#print "Projections loading COMPLETED"
	# ----- iterations
	prev_avg_absolute_voxel_change = 999999999.0
	time_projection = 0.0
	time_backprojection = 0.0
	time_iterations = clock()
	for iter_no in xrange(max_iterations_count):
		q = model_blank(nx, nx, nx)
		for i in range(len(projections_angles)):
			for j in range(len(projections_angles[i])):
				chao_params = {"anglelist":projections_angles[i][j],"radius":radius}
				time_start = clock()
				w = solution.project("chao", chao_params)
				time_projection += clock() - time_start
				p = projections_data[i][j] / threshold_to_minval(w, min_allowed_divisor)
				time_start = clock()
				q += p.backproject("chao", chao_params)
				time_backprojection += clock() - time_start
		Util.div_img( q, a )
		Util.mul_img( q, solution ) # q <- new solution  
		avg_absolute_voxel_change = q.cmp("lod",solution,{"mask":sphere3D,"negative":0,"normalize":0}) / sphere3D_volume
		if avg_absolute_voxel_change > prev_avg_absolute_voxel_change:
			#print "Finish and return last good solution"
			break
		prev_avg_absolute_voxel_change = avg_absolute_voxel_change
		solution = q
		#print "Iteration ", iter_no, ",  avg_abs_voxel_change=", avg_absolute_voxel_change 
		if min_avg_abs_voxel_change > avg_absolute_voxel_change:
			break
	time_iterations = clock() - time_iterations
	# ----- return solution and exit
	#print "Times: iterations=", time_iterations, "  project=", time_projection, "  backproject=", time_backprojection
	return solution


def recons3d_em_MPI(projections_stack, output_file, max_iterations_count = 100, radius = -1, min_norm_absolute_voxel_change = 0.01, use_weights = False, symmetry = "c1", min_norm_squared_voxel_change = 0.0001):
	"""
	Reconstruction algorithm basing on the Expectation Maximization method.
		projections_stack              -- file or list with projections
		max_iterations_count           -- stop criterion 
		min_norm_absolute_voxel_change -- stop criterion (set -1 to switch off) 
		min_norm_squared_voxel_change  -- stop criterion (set -1 to switch off)
		use_weights                    -- true == multiply projections by extra weights
		symmetry                       -- type of symmetry
	#
	"""
	from time import clock
	from utilities import model_blank, model_circle, model_square, circumference
	from morphology import threshold_to_minval
	import types
	from string import replace
	from mpi import mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
	from utilities import reduce_EMData_to_root, bcast_EMData_to_all, bcast_number_to_all, send_EMData, recv_EMData
	min_allowed_divisor = 0.0001
	
	mpi_n = mpi_comm_size(MPI_COMM_WORLD)
	mpi_r = mpi_comm_rank(MPI_COMM_WORLD)
	
	# ----- read projections 
	if type(projections_stack) is types.StringType:
		all_projs_count = EMUtil.get_image_count(projections_stack)
	else:
		all_projs_count = len(projections_stack)
	
	if all_projs_count < mpi_n:
		ERROR("Number of projections cannot be less than number of MPI processes", "recons3d_em")
	
	projs_begin = (mpi_r * all_projs_count) // mpi_n
	projs_end = ((mpi_r+1) * all_projs_count) // mpi_n
	
	if type(projections_stack) is types.StringType:
		projections = EMData.read_images(projections_stack, range(projs_begin,projs_end))
	else:
		#projections = projections_stack[projs_begin:projs_end]
		projections = projections_stack
	# ----------------------------------------------
	
	nx = projections[0].get_xsize()
	if (projections[0].get_ysize() != nx) or (projections[0].get_zsize() != 1):
		ERROR("This procedure works only for square images", "recons3d_em")

	if radius < 0: radius = nx // 2 - 1
	sphere2D = model_circle(radius, nx, nx)   
	sphere3D = model_circle(radius, nx, nx, nx)
	solution = model_blank(nx, nx, nx)
	a = model_blank(nx, nx, nx) # normalization volume
	e2D = model_square(nx, nx, nx)
	if mpi_r == 0:
		print "MPI processes: ", mpi_n
		print "Parameters:  size=%d  radius=%d  projections_count=%d  max_iterations_count=%d min_norm_absolute_voxel_change=%f" % (
						nx, radius, all_projs_count, max_iterations_count, min_norm_absolute_voxel_change )	
	
	# ----- create initial solution, calculate weights and normalization image (a)
	projections_angles = []  # list of lists of angles
	projections_data   = []  # list of lists of projections' images with weights
	for proj in projections:
		angles = [] # list of angles
		data = []   # list of projections' images with weights
		RA = proj.get_attr( "xform.projection" )
		Util.mul_img( proj, sphere2D )
		for j in range(RA.get_nsym(symmetry)):
			angdict = RA.get_sym(symmetry,j).get_rotation("spider") 
			angles.append( [angdict["phi"], angdict["theta"], angdict["psi"]] )
			chao_params = {"anglelist":angles[j],"radius":radius}
			Util.add_img( solution, proj.backproject("chao", chao_params) )
			Util.add_img( a, e2D.backproject("chao", chao_params) )
			if use_weights:
				proj3Dsphere = sphere3D.project("chao", chao_params)
				Util.mul_scalar( proj3Dsphere, 1.0 / Util.infomask(proj3Dsphere, None, True)[3] )
				Util.mul_img( proj, proj3Dsphere )
			data.append(proj)
		projections_angles.append(angles)
		projections_data.append(data)
	# reduce_scatter(solution)
	reduce_EMData_to_root(solution, mpi_r)
	bcast_EMData_to_all  (solution, mpi_r)
	# reduce_scatter(a)
	reduce_EMData_to_root(a, mpi_r)
	bcast_EMData_to_all  (a, mpi_r)
	# ------------------------
	a = threshold_to_minval(a, min_allowed_divisor)  # make sure that voxels' values are not too small (image a is divisior)
	Util.mul_img( solution, sphere3D )
	Util.div_img( solution, a )
	if mpi_r == 0: print "Projections loading COMPLETED"
	# ----- iterations
	prev_avg_absolute_voxel_change = 999999999.0
	time_projection = 0.0
	time_backprojection = 0.0
	time_iterations = clock()
	for iter_no in xrange(max_iterations_count):
		q = model_blank(nx, nx, nx)
		for i in range(len(projections_angles)):
			for j in range(len(projections_angles[i])):
				chao_params = {"anglelist":projections_angles[i][j],"radius":radius}
				time_start = clock()
				w = solution.project("chao", chao_params)
				time_projection += clock() - time_start
				p = projections_data[i][j] / threshold_to_minval(w, min_allowed_divisor)
				time_start = clock()
				q += p.backproject("chao", chao_params)
				time_backprojection += clock() - time_start
		# reduce_scatter(q)
		reduce_EMData_to_root(q, mpi_r)
		bcast_EMData_to_all  (q, mpi_r)
		# ----------------------
		Util.div_img( q, a )
		Util.mul_img( q, solution ) # q <- new solution  
		norm_absolute_voxel_change = q.cmp("lod",solution,{"mask":sphere3D,"negative":0,"normalize":0}) / q.cmp("lod",model_blank(nx,nx,nx),{"mask":sphere3D,"negative":0,"normalize":0})
		norm_squared_voxel_change  = q.cmp("sqEuclidean",solution,{"mask":sphere3D}) / q.cmp("sqEuclidean",model_blank(nx,nx,nx),{"mask":sphere3D})
		if norm_absolute_voxel_change > prev_avg_absolute_voxel_change:
			if mpi_r == 0: print "Finish and return last good solution"
			break
		prev_avg_absolute_voxel_change = norm_absolute_voxel_change
		solution = q
		solution = circumference(solution, radius-2, radius)
		if (iter_no+1)%5 == 0 and mpi_r == 0:
			solution.write_image(replace(output_file, ".hdf", "_%03d.hdf"%(iter_no+1)))
		if mpi_r == 0: print "Iteration ", iter_no+1, ",  norm_abs_voxel_change=", norm_absolute_voxel_change, ",  norm_squared_voxel_change=", norm_squared_voxel_change 
		if min_norm_absolute_voxel_change > norm_absolute_voxel_change or min_norm_squared_voxel_change > norm_squared_voxel_change:
			break
	time_iterations = clock() - time_iterations
	# ----- return solution and exit
	if mpi_r == 0: print "Times: iterations=", time_iterations, "  project=", time_projection, "  backproject=", time_backprojection
	return solution


def recons3d_sirt(stack_name, list_proj, radius, lam=1.0e-4, maxit=100, symmetry="c1", tol=0.001):
	"""
	SIRT
		
		tol   -- convergence tolerance
		lam   -- damping parameter
		maxit -- maximum number of iterations
	#
	"""
	from math import sqrt
	from utilities import model_circle
	#  analyze the symmetries Phil's code has all symmetries ready...
	nsym=1

	#  get image size from the first image
	data = EMData()
	data.read_image(stack_name,list_proj[0])
	nx = data.get_xsize()
	mask2d=model_circle(radius,nx,nx)  # SIRT works for squares only!
	mask2d = 1.0 - mask2d  # invert the mask to get average in corners
	nangles = len(list_proj)
	#
	mask3d=model_circle(radius,nx,nx,nx) # a 3D mask for error calculation
	#
	# create a volume to hold the reconstruction 
	#
	xvol = EMData()
	xvol.set_size(nx,nx,nx)
	xvol.to_zero()
	#
	# create a volume to hold trans(P)*P*xvol
	#
	pxvol = xvol.copy()
	#  array of symmetrized angles
	symangles=3*[0.0]
	angles = []

	# start iterating
	iter  = 1
	while iter <= maxit:
		if (iter == 1):
			#
			# backproject 2D images first to create the right-hand side of the
			# the normal equation
			#
			bvol = EMData()
			bvol.set_size(nx,nx,nx)
			bvol.to_zero()
			for i in xrange(nangles):
				# read projections and do initial  backprojection
				data.read_image(stack_name,list_proj[i])
				stat = Util.infomask(data, mask2d, False)
				data = data-stat[0]   # subtract the background average in the corners
				
				RA = data.get_attr( "xform.projection" )

				angles.append(RA)
				#ATTENTION
				#for transform in Symmetry3D.get_symmetries(symmetry):
					#Tf = transform*RA
					# though why do you go from 1 to nysm? why not 0 to nysm-1 ? It should be
					# equivalent unless I am missing something
					#angdict = Tf.get_params("spider")
					# then continue as before
				for ns in xrange(1,nsym+1):
					# multiply myangles by symmetry using Phil's Transform class
					Tf=RA.get_sym(symmetry,ns) #
					angdict = Tf.get_rotation("spider")
					#   Chao - please check the order of phi, theta, psi
					symangles[0] = angdict["phi"]
					symangles[1] = angdict["theta"]
					symangles[2] = angdict["psi"]
					myparams = {"anglelist":symangles, "radius":radius}
					bvol += data.backproject("chao", myparams)
			old_rnorm = bnorm = sqrt(bvol.cmp("dot",bvol,{"mask":mask3d,"negative":0}))
			grad  = bvol
		else:
			#  Insert your favorite MPI here
			pxvol.to_zero() 
			for i in xrange(nangles):
				# just a single slice of phi, theta, psi
				RA = angles[i]
				for ns in xrange(1,nsym+1):
					# multiply myangles by symmetry using Phil's Transform class
					Tf = RA.get_sym(symmetry,ns)#Tf.get_rotation()
					angdict = Tf.get_rotation("spider")
					#				    Chao - please check the order of phi, theta, psi
					symangles[0] = angdict["phi"]
					symangles[1] = angdict["theta"]
					symangles[2] = angdict["psi"]
					myparams = {"anglelist":symangles, "radius":radius}
					data  = xvol.project("chao", myparams) 
					pxvol += data.backproject("chao",myparams)
			grad  = bvol - pxvol

		rnorm = sqrt(grad.cmp("dot",grad,{"mask":mask3d,"negative":0}))
		print 'iter = %3d,  rnorm = %6.3f,  rnorm/bnorm = %6.3f' % (iter,rnorm,rnorm/bnorm)
		if (rnorm < tol or rnorm > old_rnorm): break
		old_rnorm = rnorm
		xvol = xvol + lam*grad
		iter = iter + 1

	return  xvol

def recons3d_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1", radius=None): 
	"""
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		list_proj  - list of projections to be used in the reconstruction
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	""" 
	import types
	from utilities import get_im

	if type(stack_name) == types.StringType:
		B = EMData()
		B.read_image(stack_name,list_proj[0])
	else : B = stack_name[list_proj[0]].copy()

	ny = B.get_ysize()  # have to take ysize, because xsize is different for real and fft images

	if radius == None: radius = (ny - 1) // 2

	CUBE = EMData()
	CUBE.set_size(ny, ny, ny)
	CUBE.to_zero()

	nsym = Transform.get_nsym(symmetry)
	nimages = len(list_proj)

	ss = [0.0]*(6*nsym*nimages)
	symmetry_transforms = [ [None for i in xrange(nsym)] for j in xrange(nimages) ] # list of nimages lists of nsym elements
	for iProj in xrange(nimages):
		if type(stack_name) == types.StringType:
			B.read_image(stack_name,list_proj[iProj], True)
		else:
			B = stack_name[list_proj[iProj]]
		transform = B.get_attr("xform.projection")
		for iSym in xrange(nsym):
			symmetry_transforms[iProj][iSym] = transform.get_sym(symmetry, iSym)
			d = symmetry_transforms[iProj][iSym].get_params("spider")
			DMnSS = Util.CANG(d["phi"], d["theta"], d["psi"])
			ss[ (iProj*nsym+iSym)*6 : (iProj*nsym+iSym+1)*6 ] = DMnSS["SS"]

	if method=="exact":    
		const = int(const)

	for iProj in xrange(nimages):
		proj = get_im(stack_name, list_proj[iProj])
		for iSym in xrange(nsym):
			B = proj.copy()
			B.set_attr("xform.projection", symmetry_transforms[iProj][iSym])
			if   method=="general":  Util.WTF(B, ss, const, iProj*nsym+iSym+1)  # counting in WTF start from 1!
			elif method=="exact"  :  Util.WTM(B, ss, const, iProj*nsym+iSym+1)  # counting in WTM start from 1!
			Util.BPCQ(B, CUBE, radius)

	return CUBE

      
def recons3d_vwbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1",outstack="bdb:temp"): 
	"""
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		list_proj  - list of projections to be used in the reconstruction
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object

		WARNING - symmetries not implemented!!!!!!!!!
	""" 
	import types
	from utilities import get_im

	if type(stack_name) == types.StringType:
		B = EMData()
		B.read_image(stack_name,list_proj[0])
	else : B = stack_name[list_proj[0]].copy()

	ny = B.get_ysize()  # have to take ysize, because xsize is different for real and fft images

	nsym = Transform.get_nsym(symmetry)
	nimages = len(list_proj)

	ss = [0.0]*(6*nsym*nimages)
	symmetry_transforms = [ [None for i in xrange(nsym)] for j in xrange(nimages) ] # list of nimages lists of nsym elements
	for iProj in xrange(nimages):
		if type(stack_name) == types.StringType:
			B.read_image(stack_name,list_proj[iProj], True)
		else:
			B = stack_name[list_proj[iProj]]
		transform = B.get_attr("xform.projection")
		for iSym in xrange(nsym):
			symmetry_transforms[iProj][iSym] = transform.get_sym(symmetry, iSym)
			d = symmetry_transforms[iProj][iSym].get_params("spider")
			DMnSS = Util.CANG(d["phi"], d["theta"], d["psi"])
			ss[ (iProj*nsym+iSym)*6 : (iProj*nsym+iSym+1)*6 ] = DMnSS["SS"]

	if method=="exact":    
		const = int(const)

	for iProj in xrange(nimages):
		if(iProj%100 == 0):  print "BPCQ  ",iProj
		CUBE = EMData()
		CUBE.set_size(ny, ny, ny)
		CUBE.to_zero()
		proj = get_im(stack_name, list_proj[iProj])
		for iSym in xrange(nsym):
			B = proj.copy()
			B.set_attr("xform.projection", symmetry_transforms[iProj][iSym])
			if   method=="general":  Util.WTF(B, ss, const, iProj*nsym+iSym+1)  # counting in WTF start from 1!
			elif method=="exact"  :  Util.WTM(B, ss, const, iProj*nsym+iSym+1)  # counting in WTM start from 1!
			from filter import filt_tanl
			B = filt_tanl(B, 0.3, 0.1)
			Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)
		CUBE.write_image(outstack, iProj)
	return CUBE


def prepare_wbp(stack_name, list_proj, method = "general", const=1.0E4, symmetry="c1"):
	"""
		Prepare auxiliary arrays dm and ss.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		list_proj  - list of projections to be used in the reconstruction
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	"""
	import types

	if type(stack_name) == types.StringType:
		B = EMData()
		B.read_image(stack_name,list_proj[0])
	else : B = stack_name[list_proj[0]].copy()

	nx = B.get_xsize()

	RA = Transform()
	nsym = RA.get_nsym(symmetry)

	nimages = len(list_proj)
	ntripletsWnsym = nsym*nimages
	dm=[0.0]*(9*ntripletsWnsym)
	ss=[0.0]*(6*ntripletsWnsym)
	count = 0
	from utilities import get_params_proj
	for i in xrange(nimages):
		if type(stack_name) == types.StringType:
			B.read_image(stack_name,list_proj[i], True)
			PHI, THETA, PSI, s2x, s2y = get_params_proj( B )
		else:  
			PHI, THETA, PSI, s2x, s2y = get_params_proj( stack_name[list_proj[i]] )
		DMnSS = Util.CANG(PHI,THETA,PSI)
		dm[(count*9) :(count+1)*9] = DMnSS["DM"]
		ss[(count*6) :(count+1)*6] = DMnSS["SS"]
		count += 1
	return dm,ss


def recons3d_swbp(A, transform, L, ss, method = "general", const=1.0E4, symmetry="c1"):
	"""
	    Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, and backproject.
		The projection number is L counting from zero.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	"""
	B = A.copy()
	nx = B.get_xsize()
	if(method=="exact"  ):    const = int(const)
	nsym = 1
	CUBE = EMData()
	CUBE.set_size(nx, nx, nx)
	CUBE.to_zero()

	org_transform = B.get_attr("xform.projection")

	count = 0
	for j in xrange(1):
	  	count += 1   # It has to be there as counting in WTF and WTM start from 1!
	  	if   (method=="general"):    Util.WTF(B, ss, const, L+1)
	  	elif (method=="exact"  ):    Util.WTM(B, ss, const, L+1)

		B.set_attr("xform.projection", transform)
		Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)
		
	B.set_attr("xform.projection", org_transform)
	return CUBE, B

def weight_swbp(A, L, ss, method = "general", const=1.0E4, symmetry="c1"):
	"""
	    Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, return weighted projection.
		The projection number is L counting from zero.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	"""

	if(method=="exact"  ):    const = int(const)
	nsym = 1
	B = A.copy()
	count = 0
	for j in xrange(1):
	  	count += 1   # It has to be there as counting in WTF and WTM start from 1!
	  	if   (method=="general"):    Util.WTF(B, ss, const, L+1)
	  	elif (method=="exact"  ):    Util.WTM(B, ss, const, L+1)

	return B

def backproject_swbp(B, transform = None, symmetry="c1"): 
	"""
	    Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, and backproject.
		The projection number is L counting from zero.
		Weigthed back-projection algorithm.
		stack_name - disk stack with projections or in-core list of images
		method - "general" Radermacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	""" 

	ny = B.get_ysize()
	CUBE = EMData()
	CUBE.set_size(ny, ny, ny)
	CUBE.to_zero()

	org_transform = B.get_attr("xform.projection")
	if transform != None:
		B.set_attr("xform.projection", transform)
	Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)
	B.set_attr("xform.projection", org_transform)

	return CUBE

def one_swbp(CUBE, B, transform = None, symmetry="c1"): 
	"""
	        Take one projection, but angles form the entire set.  Build the weighting function for the given projection taking into account all,
		apply it, and backproject.
		The projection number is L counting from zero.
		method - "general" Rademacher's Gaussian, "exact" MvHs triangle
		const  - for "general" 1.0e4 works well, for "exact" it should be the diameter of the object
		symmetry - point group symmetry of the object
	""" 
	org_transform = B.get_attr("xform.projection")
	if transform != None:
		B.set_attr("xform.projection", transform)
	Util.BPCQ(B, CUBE, (B.get_ysize()-1)//2)  
	B.set_attr("xform.projection", org_transform)

def prepare_recons(data, symmetry, myid, main_node_half, half_start, step, index, finfo=None, npad = 2, mpi_comm=None):
	from random     import randint
	from utilities  import reduce_EMData_to_root
	from mpi        import mpi_barrier, MPI_COMM_WORLD
	from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	nx = data[0].get_xsize()

	fftvol_half = EMData()
	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4", half_params )
	half.setup()

	group = -1
	for i in xrange(half_start, len(data), step):
		if(index >-1 ):  group = data[i].get_attr('group')
		if(group == index):
			# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
			# if( data[i].get_attr_default('active',1) == 1):
			# 	xform_proj = data[i].get_attr( "xform.projection" )
			# 	half.insert_slice(data[i], xform_proj )
			xform_proj = data[i].get_attr( "xform.projection" )
			half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = randint(0, 1000000)
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:  return fftvol_half_file, weight_half_file

	return None, None

'''  Not used anywhere?  07/29/2015  PAP
def prepare_recons_ctf_fftvol(data, snr, symmetry, myid, main_node_half, pidlist, finfo=None, npad = 2, mpi_comm=None):
	from utilities import reduce_EMData_to_root
	from EMAN2 import Reconstructors
	from mpi import MPI_COMM_WORLD

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	nx = data[0].get_xsize()

	fftvol_half = EMData()
	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in pidlist:
		xform_proj = data[i].get_attr( "xform.projection" )
		half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	return fftvol_half, weight_half
'''

def prepare_recons_ctf(nx, data, snr, symmetry, myid, main_node_half, half_start, step, finfo=None, npad = 2, mpi_comm=None, smearstep = 0.0):
	from random     import randint
	from utilities  import reduce_EMData_to_root
	from mpi        import mpi_barrier, MPI_COMM_WORLD
	from EMAN2 import Reconstructors

	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD

	fftvol_half = EMData()

	if( smearstep > 0.0 ):
		#if myid == 0:  print "  Setting smear in prepare_recons_ctf"
		ns = 1
		smear = []
		for j in xrange(-ns,ns+1):
			if( j != 0):
				for i in xrange(-ns,ns+1):
					for k in xrange(-ns,ns+1):
						smear += [i*smearstep,j*smearstep,k*smearstep,1.0]
		# Deal with theta = 0.0 cases
		prj = []
		for i in xrange(-ns,ns+1):
			for k in xrange(-ns,ns+1):
				prj.append(i+k)
		for i in xrange(-2*ns,2*ns+1,1):
			 smear += [i*smearstep,0.0,0.0,float(prj.count(i))]
		#if myid == 0:  print "  Smear  ",smear
		fftvol_half.set_attr("smear", smear)

	weight_half = EMData()
	half_params = {"size":nx, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol_half, "weight":weight_half}
	half = Reconstructors.get( "nn4_ctf", half_params )
	half.setup()

	for i in xrange(half_start, len(data), step):
		xform_proj = data[i].get_attr( "xform.projection" )
		half.insert_slice(data[i], xform_proj )

	if not(finfo is None):
		finfo.write( "begin reduce half\n" )
		finfo.flush()

	reduce_EMData_to_root(fftvol_half, myid, main_node_half, mpi_comm)
	reduce_EMData_to_root(weight_half, myid, main_node_half, mpi_comm)

	if not(finfo is None):
		finfo.write( "after reduce half\n" )
		finfo.flush()

	if myid == main_node_half:
		tmpid = randint(0, 1000000) 
		fftvol_half_file = ("fftvol_half%d.hdf" % tmpid)
		weight_half_file = ("weight_half%d.hdf" % tmpid)
		fftvol_half.write_image(fftvol_half_file)
		weight_half.write_image(weight_half_file)
	mpi_barrier(mpi_comm)

	fftvol_half = None
	weight_half = None

	if myid == main_node_half:
		return fftvol_half_file, weight_half_file

	return None,None


def recons_from_fftvol(size, fftvol, weight, symmetry, npad = 2):
	from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "symmetry":symmetry, "fftvol":fftvol, "weight":weight}
	r = Reconstructors.get("nn4", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol


def recons_ctf_from_fftvol(size, fftvol, weight, snr, symmetry, weighting=1, npad = 2):
	from EMAN2 import Reconstructors

	params = {"size":size, "npad":npad, "snr":snr, "sign":1, "symmetry":symmetry, "fftvol":fftvol, "weight":weight, "weighting":weighting}
	r = Reconstructors.get("nn4_ctf", params)
	r.setup()
	dummy = r.finish(True)
	return fftvol


def get_image_size( imgdata, myid ):
	from mpi import mpi_gather, mpi_bcast, MPI_COMM_WORLD, MPI_INT
	nimg = len(imgdata)

        nimgs = mpi_gather( nimg, 1, MPI_INT, 1, MPI_INT, 0, MPI_COMM_WORLD )

        if myid==0:
		src = -1
		for i in xrange( len(nimgs) ):
			if int(nimgs[i]) > 0 :
				src = i
				break
		if src==-1:
			return 0
	else:
		src = -1

	size_src = mpi_bcast( src, 1, MPI_INT, 0, MPI_COMM_WORLD )

	if myid==int(size_src[0]):
		assert nimg > 0
		size = imgdata[0].get_xsize()
	else:
		size = -1

	nx = mpi_bcast( size, 1, MPI_INT, size_src[0], MPI_COMM_WORLD )
	return int(nx[0])


def rec3D_MPI(data, snr = 1.0, symmetry = "c1", mask3D = None, fsc_curve = None, \
		myid = 0, main_node = 0, rstep = 1.0, odd_start=0, eve_start=1, finfo=None, \
		index=-1, npad = 2, mpi_comm=None, smearstep = 0.0):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept 
	  in the memory, computes reconstruction and through odd-even, in order to get the resolution
	'''
	import os
	from statistics import fsc_mask
	from utilities  import model_blank, model_circle, get_image, send_EMData, recv_EMData
	from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD
	
	nproc = mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005


	if index != -1 :
		grpdata = []
		for i in xrange(len(data)):
			if data[i].get_attr('group') == index:
				grpdata.append(data[i])
		imgdata = grpdata
	else:
		imgdata = data

	nx = get_image_size(imgdata, myid)
	if nx == 0:
		ERROR("Warning: no images were given for reconstruction, this usually means there is an empty group, returning empty volume", "rec3D", 0)
		return model_blank( 2, 2, 2 ), None

	fftvol_odd_file,weight_odd_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_odd, odd_start, 2, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	fftvol_eve_file,weight_eve_file = prepare_recons_ctf(nx, imgdata, snr, symmetry, myid, main_node_eve, eve_start, 2, finfo, npad, mpi_comm=mpi_comm, smearstep = smearstep)
	del imgdata

	if nproc == 1:
		fftvol = get_image(fftvol_odd_file)
		weight = get_image(weight_odd_file)
		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		fftvol = get_image(fftvol_eve_file)
		weight = get_image(weight_eve_file)
		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
		fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3d

		fftvol = get_image( fftvol_odd_file )
		fftvol_tmp = get_image(fftvol_eve_file)
		fftvol += fftvol_tmp
		fftvol_tmp = None

		weight = get_image( weight_odd_file )
		weight_tmp = get_image(weight_eve_file)
		weight += weight_tmp
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )

		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			weight = get_image( weight_odd_file )
			volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			
			if( not mask3D ):
				nx = volodd.get_xsize()
				ny = volodd.get_ysize()
				nz = volodd.get_zsize()
				mask3D = model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)
			fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve, mask3D
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			weight = get_image( weight_eve_file )
			voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			fftvol_tmp = recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			fftvol += fftvol_tmp
			fftvol_tmp = None

			weight = get_image( weight_odd_file )
			weight_tmp = recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			weight += weight_tmp
			weight_tmp = None

			volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )

			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			weight = get_image( weight_eve_file )
			send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)
			send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file )
			return model_blank(nx,nx,nx), None

	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = get_image( fftvol_odd_file )
		send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = get_image( weight_odd_file )
		send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		del fftvol, weight
		voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)

		if( not mask3D ):
			nx = volodd.get_xsize()
			ny = volodd.get_ysize()
			nz = volodd.get_zsize()
			mask3D = model_circle(min(nx,ny,nz)//2 - 2, nx,ny,nz)

		fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve, mask3D
		volall = recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file )
		return volall, fscdat

	if myid == main_node_eve:
		ftmp = recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = get_image( fftvol_eve_file )
		Util.add_img( ftmp, fftvol )
		send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = get_image( weight_eve_file )
		send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_ctf_from_fftvol(nx, fftvol, weight, snr, symmetry, npad = npad)
		send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return model_blank(nx,nx,nx),None

	return model_blank(nx,nx,nx),None


def rec3D_MPI_noCTF(data, symmetry = "c1", mask3D = None, fsc_curve = None, myid = 2, main_node = 0, \
		rstep = 1.0, odd_start=0, eve_start=1, finfo=None, index = -1, npad = 2, mpi_comm=None):
	'''
	  This function is to be called within an MPI program to do a reconstruction on a dataset kept in the memory 
	  Computes reconstruction and through odd-even, in order to get the resolution
	  if index > -1, projections should have attribute group set and only those whose group matches index will be used in the reconstruction
	    this is for multireference alignment
	'''
	import os
	from statistics import fsc_mask
	from utilities  import model_blank, get_image,send_EMData, recv_EMData
	from mpi        import mpi_comm_size, MPI_COMM_WORLD
	
	if mpi_comm == None:
		mpi_comm = MPI_COMM_WORLD
	
	nproc = mpi_comm_size(mpi_comm)

	if nproc==1:
		assert main_node==0
		main_node_odd = main_node
		main_node_eve = main_node
		main_node_all = main_node
	elif nproc==2:
		main_node_odd = main_node
		main_node_eve = (main_node+1)%2
		main_node_all = main_node

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002
	else:
		#spread CPUs between different nodes to save memory
		main_node_odd = main_node
		main_node_eve = (int(main_node)+nproc-1)%int(nproc)
		main_node_all = (int(main_node)+nproc//2)%int(nproc)

		tag_voleve     = 1000
		tag_fftvol_eve = 1001
		tag_weight_eve = 1002

		tag_fftvol_odd = 1003
		tag_weight_odd = 1004
		tag_volall     = 1005

	nx = data[0].get_xsize()

	fftvol_odd_file,weight_odd_file = prepare_recons(data, symmetry, myid, main_node_odd, odd_start, 2, index, finfo, npad, mpi_comm=mpi_comm)
	fftvol_eve_file,weight_eve_file = prepare_recons(data, symmetry, myid, main_node_eve, eve_start, 2, index, finfo, npad, mpi_comm=mpi_comm) 

	if nproc == 1:
		fftvol = get_image( fftvol_odd_file )
		weight = get_image( weight_odd_file )
		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fftvol = get_image( fftvol_eve_file )
		weight = get_image( weight_eve_file )
		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)

		fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve

		fftvol = get_image( fftvol_odd_file )
		Util.add_img( fftvol, get_image(fftvol_eve_file) )

		weight = get_image( weight_odd_file )
		Util.add_img( weight, get_image(weight_eve_file) )

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
		return volall,fscdat

	if nproc == 2:
		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			weight = get_image( weight_odd_file )
			volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)
			fscdat = fsc_mask( volodd, voleve, mask3D, rstep, fsc_curve)
			del  volodd, voleve
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			weight = get_image( weight_eve_file )
			voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)

		if myid == main_node_odd:
			fftvol = get_image( fftvol_odd_file )
			fftvol_tmp = recv_EMData( main_node_eve, tag_fftvol_eve, mpi_comm)
			Util.add_img( fftvol, fftvol_tmp )
			fftvol_tmp = None

			weight = get_image( weight_odd_file )
			weight_tmp = recv_EMData( main_node_eve, tag_weight_eve, mpi_comm)
			Util.add_img( weight, weight_tmp )
			weight_tmp = None
			volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
			os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
			return volall,fscdat
		else:
			assert myid == main_node_eve
			fftvol = get_image( fftvol_eve_file )
			send_EMData(fftvol, main_node_odd, tag_fftvol_eve, mpi_comm)

			weight = get_image( weight_eve_file )
			send_EMData(weight, main_node_odd, tag_weight_eve, mpi_comm)
			os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );
			return model_blank(nx,nx,nx), None
	# cases from all other number of processors situations
	if myid == main_node_odd:
		fftvol = get_image( fftvol_odd_file )
		send_EMData(fftvol, main_node_eve, tag_fftvol_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("fftvol odd sent\n")
			finfo.flush()

		weight = get_image( weight_odd_file )
		send_EMData(weight, main_node_all, tag_weight_odd, mpi_comm)

		if not(finfo is None):
			finfo.write("weight odd sent\n")
			finfo.flush()

		volodd = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		del fftvol, weight
		voleve = recv_EMData(main_node_eve, tag_voleve, mpi_comm)
		fscdat = fsc_mask(volodd, voleve, mask3D, rstep, fsc_curve)
		del  volodd, voleve
		volall = recv_EMData(main_node_all, tag_volall, mpi_comm)
		os.system( "rm -f " + fftvol_odd_file + " " + weight_odd_file );
		return volall,fscdat

	if myid == main_node_eve:
		ftmp = recv_EMData(main_node_odd, tag_fftvol_odd, mpi_comm)
		fftvol = get_image( fftvol_eve_file )
		Util.add_img( ftmp, fftvol )
		send_EMData(ftmp, main_node_all, tag_fftvol_eve, mpi_comm)
		del ftmp

		weight = get_image( weight_eve_file )
		send_EMData(weight, main_node_all, tag_weight_eve, mpi_comm)

		voleve = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		send_EMData(voleve, main_node_odd, tag_voleve, mpi_comm)
		os.system( "rm -f " + fftvol_eve_file + " " + weight_eve_file );

		return model_blank(nx,nx,nx), None


	if myid == main_node_all:
		fftvol = recv_EMData(main_node_eve, tag_fftvol_eve, mpi_comm)
		if not(finfo is None):
			finfo.write( "fftvol odd received\n" )
			finfo.flush()

		weight = recv_EMData(main_node_odd, tag_weight_odd, mpi_comm)
		weight_tmp = recv_EMData(main_node_eve, tag_weight_eve, mpi_comm)
		Util.add_img( weight, weight_tmp )
		weight_tmp = None

		volall = recons_from_fftvol(nx, fftvol, weight, symmetry, npad)
		send_EMData(volall, main_node_odd, tag_volall, mpi_comm)

		return model_blank(nx,nx,nx),None


	return model_blank(nx,nx,nx),None
