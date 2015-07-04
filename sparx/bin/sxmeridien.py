#!/usr/bin/env python
#
#  06/27/2015
#  New version.  
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

global cushion
cushion = 8


def AI( icurrentres, Tracker, HISTORY ):
	#  
	#  Possibilities we will consider:
	#    1.  resolution improved: keep going with current settings.
	#    2.  resolution stalled and no pwadjust: turn on pwadjust
	#    3.  resolution stalled and pwadjust: move to the next phase
	#    4.  resolution decreased: back off and move to the next phase
	#    5.  All phases tried and nxinit < nnxo: set nxinit == nnxo and run local searches.
	from sys import exit
	reset_data = False
	Tracker["delpreviousmax"] = False
	if(Tracker["mainiteration"] == 1):
		nxinit                 = Tracker["nxinit"]
		while( icurrentres + cushion > nxinit//2 ): nxinit += Tracker["nxstep"]
		Tracker["nxinit"]      = min(nxinit,Tracker["constants"]["nnxo"])
		Tracker["icurrentres"] = icurrentres
		Tracker["nsoft"]       = 0
		Tracker["local"]       = False
		Tracker["zoom"]        = True
		if not Tracker["applyctf"] :  reset_data  = True
		Tracker["applyctf"]    = True
		#  Develop something intelligent
		Tracker["xr"] = "6 3"
		Tracker["ts"] = "2 1"
		Tracker["an"]          = "-1"
		Tracker["anger"]       = -1.0
		Tracker["shifter"]     = -1.0
		Tracker["constants"]["best"] = Tracker["mainiteration"]
		keepgoing = 1
	elif(Tracker["mainiteration"] == 2):
		#  Go back to initial window size and exhaustive search to improve centering of the data.
		nxinit = HISTORY[0]["nxinit"]
		Tracker["nxinit"]       = nxinit
		Tracker["icurrentres"] = min(icurrentres, nxinit//2-3)
		Tracker["nsoft"]       = 0
		Tracker["local"]       = False
		Tracker["zoom"]        = True
		if not Tracker["applyctf"] :  reset_data  = True
		Tracker["applyctf"]    = True
		#  Go back to initial settings
		Tracker["xr"] , Tracker["ts"] = stepali(Tracker["nxinit"] , Tracker["constants"]["nnxo"], Tracker["constants"]["radius"])
		Tracker["an"]          = "-1"
		Tracker["anger"]       = -1.0
		Tracker["shifter"]     = -1.0
		Tracker["constants"]["best"] = Tracker["mainiteration"]
		keepgoing = 1
	elif(Tracker["mainiteration"] == 3):
		#  Here we do soft with the initial window size
		nxinit                 = Tracker["nxinit"]
		Tracker["icurrentres"] = min(icurrentres, nxinit//2-3)
		Tracker["nsoft"]       = 1
		Tracker["delpreviousmax"] = True
		Tracker["local"]       = False
		Tracker["zoom"]        = False
		if not Tracker["applyctf"] :  reset_data  = True
		Tracker["applyctf"]    = True
		#  This is for soft
		Tracker["xr"] = "3"
		Tracker["ts"] = "1"
		Tracker["an"]          = "-1"
		Tracker["anger"]       = -1.0
		Tracker["shifter"]     = -1.0
		Tracker["constants"]["best"] = Tracker["mainiteration"]
		keepgoing = 1
	#  First four steps were preset mainly to establish shifts and decent initial structure.
	#    Now the program has to start making decisions.
	elif(Tracker["mainiteration"] > 3):
		#  Each case needs its own settings.  We arrived here after soft at initial window size.
		#  Possibilities we will consider:
		#    1.  resolution improved: keep going with current settings.
		#    2.  resolution stalled and no pwadjust: turn on pwadjust
		#    3.  resolution stalled and pwadjust: move to the next phase
		#    4.  resolution decreased: back off and move to the next phase
		#    5.  All phases tried and nxinit < nnxo: set nxinit == nnxo and run local searches.
		if( Tracker["state"] == "INITIAL" ):
			move_up_phase = True
		else:
			#  For all other states make a decision based on resolution change
			direction = icurrentres - Tracker["icurrentres"]
			if Tracker["movedback"] :
				# previous move was back, but the resolution did not improve
				if direction <= 0 :
					if Tracker["constants"]["pwreference"] :
						if Tracker["PWadjustment"] :
							keepgoing = 0
						else:
							Tracker["PWadjustment"] = Tracker["constants"]["pwreference"]
							Tracker["movedback"] = False
							keepgoing = 1
					else:
						keepgoing = 0
				else:
					#  Resolution improved
					Tracker["movedback"] = False
					#  And then continue with normal moves

			move_up_phase = False
			#  Resolution improved, adjust window size and keep going
			if(  direction > 0 ):
				keepgoing = 1
				if(Tracker["state"] == "FINAL2"):
					if( icurrentres + 4 > Tracker["constants"]["nnxo"]): keepgoing = 0
				else:
					nxinit = Tracker["nxinit"]
					while( icurrentres + cushion > nxinit//2 ): nxinit += Tracker["nxstep"]
					Tracker["nxinit"] = min(nxinit,Tracker["constants"]["nnxo"])
				Tracker["icurrentres"] = icurrentres
				Tracker["constants"]["best"] = Tracker["mainiteration"]

				if(Tracker["state"] == "EXHAUSTIVE"):
					#  Move up if changes in angles are less than 30 degrees (why 30??)
					if(Tracker["anger"]   < 30.0 ):  move_up_phase = True
					else:
						Tracker["xr"] = "%d"%(int(Tracker["shifter"]*float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))+1)
				elif(Tracker["state"] == "RESTRICTED"):
					#switch to the next phase if restricted searches make no sense anymore
					#  The next sadly repeats what is in the function that launches refinement, but I do nto know how to do it better.
					sh = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
					rd = int(Tracker["constants"]["radius"] * sh +0.5)
					rl = float(Tracker["icurrentres"])/float(Tracker["nxinit"])
					dd = degrees(atan(0.5/rl/rd))
					if( Tracker["anger"]  < dd ):  move_up_phase = True
					else:
						Tracker["xr"] = "%d"%(int(Tracker["shifter"]*float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))+1)
				else:
					Tracker["anger"]   = -1.0
					Tracker["shifter"] = -1.0

			#  Resolution stalled
			elif( direction == 0 ):
				if Tracker["constants"]["pwreference"] :
					if Tracker["PWadjustment"] :
						# move up with the state
						move_up_phase = True
					else:
						# turn on pwadjustment
						Tracker["PWadjustment"] = Tracker["constants"]["pwreference"]
						if(Tracker["state"] == "EXHAUSTIVE"):
								xr = int(Tracker["shifter"]*float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))+1
								Tracker["zoom"] = True
								Tracker["xr"] = "%d   %d"%(2*xr, xr)
								Tracker["ts"] =  "2 1"
						elif(Tracker["state"] == "RESTRICTED"):
								xr = int(Tracker["shifter"]*float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))+1
								Tracker["zoom"] = True
								Tracker["xr"] = "%d   %d"%(2*xr, xr)
								Tracker["ts"] =  "2 1"
								Tracker["an"] =  "%6.2f  %6.2f"(2*Tracker["anger"],2*Tracker["anger"])					
						keepgoing = 1
				elif( Tracker["state"] == "FINAL2"):  keepgoing = 0
				else:
					# move up with the state
					move_up_phase = True
				Tracker["constants"]["best"] = Tracker["mainiteration"]
	
			# Resolution decreased
			elif( direction < 1 ):
				# Come up with rules
				lb = -1
				for i in xrange(len(HISTORY)):
					if( HISTORY[i]["mainiteration"] == Tracker["constants"]["best"] ):
						lb = i
						break
				if( lb == -1 ):
					ERROR("No best solution in HISTORY, cannot be","sxmeridien",1)
					exit()
				#  Here we have to jump over the current state
				#  However, how to avoid cycling between two states?
				Tracker["movedback"] = True
				stt = [Tracker["state"], Tracker["PWadjustment"], Tracker["mainiteration"] ]
				Tracker = HISTORY[lb].copy()
				Tracker["state"]          = stt[0]
				Tracker["PWadjustment"]   = stt[1]
				Tracker["mainiteration"]  = stt[2]
				icurrentres               = Tracker["icurrentres"]
				#  This will set previousoutputdir to the best parames back then.
				move_up_phase = True
				
	
	
		if move_up_phase:

			#  INITIAL
			if(Tracker["state"] == "INITIAL"):
				#  Switch to EXHAUSTIVE
				nxinit = Tracker["nxinit"]
				while( icurrentres + cushion > nxinit//2 ): nxinit += Tracker["nxstep"]
				Tracker["nxinit"] = min(nxinit,Tracker["constants"]["nnxo"])
				Tracker["icurrentres"] = icurrentres
				Tracker["nsoft"]       = 0
				Tracker["local"]       = False
				Tracker["zoom"]        = False
				if not Tracker["applyctf"] :  reset_data  = True
				Tracker["applyctf"]    = True
				#  Switch to exhaustive
				Tracker["upscale"]     = 0.5
				Tracker["state"]       = "EXHAUSTIVE"
				#  Develop something intelligent
				Tracker["xr"] = "%d"%(int(Tracker["shifter"]*float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))+1)
				Tracker["ts"] = "1"
				keepgoing = 1
			#  Exhaustive searches
			elif(Tracker["state"] == "EXHAUSTIVE"):
				#  Switch to RESTRICTED
				nxinit = Tracker["nxinit"]
				while( icurrentres + cushion > nxinit//2 ): nxinit += Tracker["nxstep"]
				Tracker["nxinit"] = min(nxinit,Tracker["constants"]["nnxo"])
				Tracker["icurrentres"] = icurrentres
				Tracker["nsoft"]       = 0
				Tracker["local"]       = False
				Tracker["zoom"]        = False
				if Tracker["applyctf"] :  reset_data  = True
				Tracker["upscale"]     = 0.5
				Tracker["applyctf"]    = True
				Tracker["an"]          = "30"
				Tracker["state"]       = "RESTRICTED"
				Tracker["maxit"]       = 50
				Tracker["xr"] = "%d"%(int(Tracker["shifter"]*float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"]))+1)
				Tracker["ts"] = "1"
				keepgoing = 1
			#  Restricted searches
			elif(Tracker["state"] == "RESTRICTED"):
				#  Switch to LOCAL
				nxinit = Tracker["nxinit"]
				while( icurrentres + cushion > nxinit//2 ): nxinit += Tracker["nxstep"]
				Tracker["nxinit"] = min(nxinit,Tracker["constants"]["nnxo"])
				Tracker["icurrentres"] = icurrentres
				Tracker["nsoft"]       = 0
				Tracker["local"]       = True
				Tracker["zoom"]        = False
				if Tracker["applyctf"] :  reset_data  = True
				Tracker["upscale"]     = 0.5
				Tracker["applyctf"]    = False
				Tracker["an"]          = "-1"
				Tracker["state"]       = "LOCAL"
				Tracker["maxit"]       = 1
				Tracker["xr"] = "2"
				Tracker["ts"] = "2"
				keepgoing = 1
			#  Local searches
			elif(Tracker["state"] == "LOCAL"):
				#  Switch to FINAL1
				Tracker["nxinit"] = Tracker["constants"]["nnxo"]
				Tracker["icurrentres"] = icurrentres
				Tracker["nsoft"]       = 0
				Tracker["local"]       = True
				Tracker["zoom"]        = False
				if Tracker["applyctf"] :  reset_data  = True
				Tracker["applyctf"]    = False
				Tracker["upscale"]     = 0.5
				Tracker["an"]          = "-1"
				Tracker["state"]       = "FINAL1"
				Tracker["maxit"]       = 10
				Tracker["xr"] = "2"
				Tracker["ts"] = "2"
				keepgoing = 1
			elif(Tracker["state"] == "FINAL1"):
				#  Switch to FINAL2
				Tracker["nxinit"] = Tracker["constants"]["nnxo"]
				Tracker["icurrentres"] = icurrentres
				Tracker["nsoft"]       = 0
				Tracker["local"]       = True
				Tracker["zoom"]        = False
				if Tracker["applyctf"] :  reset_data  = True
				Tracker["upscale"]     = 1.0
				Tracker["applyctf"]    = False
				Tracker["an"]          = "-1"
				Tracker["state"]       = "FINAL2"
				Tracker["maxit"]       = 10
				Tracker["xr"] = "2"
				Tracker["ts"] = "2"
				keepgoing = 1
			elif(Tracker["state"] == "FINAL2"):
				keepgoing = 0
			else:
				ERROR(" Unknown phase","sxmeridien",1)
				exit()  #  This will crash the program, but the situation cannot occur

	Tracker["previousoutputdir"] = Tracker["directory"]
	return keepgoing, reset_data, Tracker


def threshold_params_changes(currentdir, previousdir, th = 0.95, sym = "c1"):
	#  Indexes contain list of images processed - sorted integers, subset of the full range.
	#  params - contain parameters associated with these images
	#  Both lists can be of different sizes, so we have to find a common subset
	from utilities    import getang3
	from pixel_error  import max_3D_pixel_error
	from EMAN2        import Vec2f
	import sets

	cids    = read_text_file(os.path.join(currentdir,"indexes.txt"))
	pids    = read_text_file(os.path.join(previousdir,"indexes.txt"))
	cparams = read_text_row(os.path.join(currentdir,"params.txt"))
	pparams = read_text_row(os.path.join(previousdir,"params.txt"))
	u = list(set(cids) & set(pids))
	u.sort()
	#  Extract common subsets of parameters
	cp = []
	pp = []
	i = 0
	for q in u:
		l = cids.index(q,i)
		cp.append(cparams[l])
		i = l+1
	i = 0
	for q in u:
		l = pids.index(q,i)
		pp.append(pparams[l])
		i = l+1
	del pparams,cparams
	n = len(u)
	anger       = [360.0]*n
	shifter     = [1.0e23]*n
	if(sym == "c1"):
		for i in xrange(n):
			anger[i]       = getang3([cp[i][0], cp[i][1]],[pp[i][0], pp[i][1]])
			shifter[i]     = max(abs(cp[i][3] - pp[i][3]),abs(cp[i][4] - pp[i][4]))
			#pixel_error[i] = max_3D_pixel_error(t1, t2, radius)
	else:
		#from utilities import get_symt
		#ts = get_symt(sym)
		for i in xrange(n):
			#t1 = Transform({"type":"spider","phi":pp[i][0],"theta":pp[i][1],"psi":pp[i][2]})
			#t1.set_trans(Vec2f(-pp[i][3], -pp[i][4]))
			t2 = Transform({"type":"spider","phi":cp[i][0],"theta":cp[i][1],"psi":cp[i][2]})
			t2.set_trans(Vec2f(-cp[i][3], -cp[i][4]))
			ts = t2.get_sym_proj(sym)
			shifter[i] = max(abs(cp[i][3] - pp[i][3]),abs(cp[i][4] - pp[i][4]))
			for kts in ts:
				# we do not care which position minimizes the error
				du = kts.get_params("spider")
				qt = getang3([pp[i][0], pp[i][1]], [du["phi"], du["theta"]])
				if(qt < anger[i]):   anger[i] = qt

	anger.sort()
	shifter.sort()
	la = min(int(th*n + 0.5), n-1)
	#  Returns error thresholds under which one has th fraction of images
	#  The shifter is given in the full scale displacement
	return anger[la], shifter[la]
	"""

	h1,h2 = hist_list(anger,32)
	#print h1
	#print h2
	h3 = h2[:]
	for i in xrange(1,len(h3)):
		h3[i] = h3[i-1]+h2[i]
	u1,u2 = hist_list(shifter,32)
	#print u1
	#print u2
	u3 = u2[:]
	for i in xrange(1,len(u3)):
		u3[i] = u3[i-1]+u2[i]
	return h1,h2,h3,u1,u2,u3
	"""



def subdict(d,u):
	# substitute values in dictionary d by those given by dictionary u
	for q in u:  d[q] = u[q]


def stepali(nxinit, nnxo, irad, nxrsteps = 3):
	txrm = (nxinit - 2*(int(irad*float(nxinit)/float(nnxo) + 0.5)+1))//2
	if (txrm < 0): ERROR("ERROR!! Shift value ($d) is too large for the mask size"%txrm)
	
	if (txrm/nxrsteps>0):
		tss = ""
		txr = ""
		while(txrm/nxrsteps>0):
			tts=txrm/nxrsteps
			tss += "%d  "%tts
			txr += "%d  "%(tts*nxrsteps)
			txrm -= nxrsteps
	else:
		txr = "%d"%txrm
		tss = "1"
	return txr, tss


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


# NOTE: 2015/06/12 Toshio Moriya
# This function seems to be not used
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

def read_fsc(fsclocation, lc, myid, main_node, comm = -1):
	# read fsc and fill it with zeroes pass lc location
	from utilities import bcast_list_to_all, read_text_file
	if comm == -1 or comm == None: comm = MPI_COMM_WORLD
	if(myid == main_node):
		f = read_text_file(fsclocation,1)
		n = len(f)
		if(n > lc+1 ):  f = f[:lc+1] +[0.0 for i in xrange(lc+1,n)]
	else: f = 0.0
	mpi_barrier(comm)
	f = bcast_list_to_all(f, myid, main_node)
	return f


def get_resolution_mrk01(vol, radi, nnxo, fscoutputdir, mask_option):
	# this function is single processor
	#  Get updated FSC curves, user can also provide a mask using radi variable
	import types
	if(type(radi) == int):
		if(mask_option is None):  mask = model_circle(radi,nnxo,nnxo,nnxo)
		else:                           mask = get_im(mask_option)
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


def get_pixel_resolution(vol, mask, nnxo, fscoutputdir):
	# this function is single processor
	nx = vol[0].get_xsize()
	nfsc = fsc(vol[0]*mask,vol[1]*mask, 1.0 )
	if(nx<nnxo):
		for i in xrange(3):
			for k in xrange(nx//2+1,nnxo/2+1):
				nfsc[i].append(0.0)
		for i in xrange(nnxo/2+1):
			nfsc[0][i] = float(i)/nnxo
	write_text_file( nfsc, os.path.join(fscoutputdir,"fsc.txt") )
	ns = len(nfsc[1])
	currentres = -1
	'''
	#  This is actual resolution, as computed by 2*f/(1+f)
	for i in xrange(1,ns-1):
		if ( nfsc[1][i] < 0.333333333333333333333333):
			currentres = nfsc[0][i-1]
			break
	'''
	#  0.5 cut-off
	for i in xrange(1,ns-1):
		if ( nfsc[1][i] < 0.5):
			currentres = i
			break
	if(currentres < 0):
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
	#lowpass, falloff = fit_tanh1(nfsc, 0.01)
	lowpass = nfsc[0][currentres]
	falloff = 0.2

	return  round(lowpass,4), round(falloff,4), currentres


def compute_resolution(stack, outputdir, partids, partstack, org_radi, nnxo, CTF, mask_option, sym, myid, main_node, nproc, pixel=1.0):
	#  while the code pretends to accept volumes as input, it actually does not.
	import types
	vol = [None]*2
	fsc = [None]*2
	if( type(stack[0]) == list ):
		nx = stack[0][0].get_xsize()
		nz = stack[0][0].get_zsize()
	else:
		nz = 1
	if(mask_option is None):  mask = model_circle(org_radi,nnxo,nnxo,nnxo)
	else:                     mask = get_im(mask_option)

	projdata = []
	for procid in xrange(2):
		if(type(stack[0]) == str or ( nz == 1 )):
			if(type(stack[0]) == str):
				projdata.append(getindexdata(stack, partids[procid], partstack[procid], myid, nproc))
			else:
				projdata.append(None)
				projdata[procid] = stack[procid]
			if( procid == 0 ):
				nx = projdata[procid][0].get_xsize()
				if( nx != nnxo):
					mask = Util.window(rot_shift3D(mask,scale=float(nx)/float(nnxo)),nx,nx,nx)

			if CTF:
				from reconstruction import rec3D_MPI
				vol[procid],fsc[procid] = rec3D_MPI(projdata[procid], symmetry = sym, \
					mask3D = mask, fsc_curve = None, \
					myid = myid, main_node = main_node, odd_start = 1, eve_start = 0, finfo = None, npad = 2)
			else:
				from reconstruction import rec3D_MPI_noCTF
				vol[procid],fsc[procid] = rec3D_MPI_noCTF(projdata[procid], symmetry = sym, \
					mask3D = mask, fsc_curve = None, \
					myid = myid, main_node = main_node, odd_start = 1, eve_start = 0, finfo = None, npad = 2)

			if(type(stack) == str):  del projdata
		else:
			#  Volumes
			vol[procid] = stack[procid]
			nx = vol[0].get_xsize()
			if( nx != nnxo ):
				mask = Util.window(rot_shift3D(mask,scale=float(nx)/float(nnxo)),nx,nx,nx)

		if( myid == main_node):
			from fundamentals import fpol
			fpol(vol[procid], nnxo, nnxo, nnxo).write_image(os.path.join(outputdir,"vol%01d.hdf"%procid))
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(  line,"Generated vol #%01d  using  image size %d "%(procid, nx))


	lowpass    = 0.0
	falloff    = 0.0
	icurrentres = 0

	if(myid == main_node):
		if(type(stack) == str or ( nz == 1 )):
			if(nx<nnxo):
				for procid in xrange(2):
					for i in xrange(3):
						for k in xrange(nx/2+1, nnxo/2+1):
							fsc[procid][i].append(0.0)
					for k in xrange(nnxo/2+1):
						fsc[procid][0][k] = float(k)/nnxo
			for procid in xrange(2):
				#  Compute adjusted within-fsc as 2*f/(1+f)
				fsc[procid].append(fsc[procid][1][:])
				for k in xrange(len(fsc[procid][1])):  fsc[procid][-1][k] = 2*fsc[procid][-1][k]/(1.0+fsc[procid][-1][k])
				write_text_file( fsc[procid], os.path.join(outputdir,"within-fsc%01d.txt"%procid) )

		lowpass, falloff, icurrentres = get_pixel_resolution(vol, mask, nnxo, outputdir)
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(  line,"Current resolution  %6.2f  %6.2f A  (%d), low-pass filter cut-off %6.2f and fall-off %6.2f"%\
						(icurrentres/float(nnxo),pixel*float(nnxo)/float(icurrentres),icurrentres,lowpass,falloff))
		write_text_row([[lowpass, falloff, icurrentres]],os.path.join(outputdir,"current_resolution.txt"))
	#  Returns: low-pass filter cutoff;  low-pass filter falloff;  current resolution
	icurrentres = bcast_number_to_all(icurrentres, source_node = main_node)
	lowpass     = bcast_number_to_all(lowpass, source_node = main_node)
	falloff     = bcast_number_to_all(falloff, source_node = main_node)
	return round(lowpass,4), round(falloff,4), icurrentres



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



def get_shrink_data(Tracker, nxinit, partids, partstack, myid, main_node, nproc, preshift = False):
	# The function will read from stack a subset of images specified in partids
	#   and assign to them parameters from partstack with optional CTF application and shifting of the data.
	# So, the lengths of partids and partstack are the same.
	#  The read data is properly distributed among MPI threads.
	if( myid == main_node ):
		print("    ")
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(  line, "Reading data  onx: %3d, nx: %3d, CTF: %s, applyctf: %s, preshift: %s."%(Tracker["constants"]["nnxo"], Tracker["nxinit"], Tracker["constants"]["CTF"], Tracker["applyctf"], preshift) )
		print("                       stack: %s\n                       partids: %s\n                       partstack:%s\n"%(Tracker["constants"]["stack"], partids, partstack) )
	if( myid == main_node ): lpartids = read_text_file(partids)
	else:  lpartids = 0
	lpartids = wrap_mpi_bcast(lpartids, main_node)
	ndata = len(lpartids)
	if( myid == main_node ):  partstack = read_text_row(partstack)
	else:  partstack = 0
	partstack = wrap_mpi_bcast(partstack, main_node)
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
	#  Preprocess the data
	mask2D  = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
	nima = image_end - image_start
	oldshifts = [[0.0,0.0]]#*nima
	data = [None]*nima
	shrinkage = nxinit/float(Tracker["constants"]["nnxo"])
	radius = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	#  Note these are in Fortran notation for polar searches
	#txm = float(nxinit-(nxinit//2+1) - radius -1)
	#txl = float(2 + radius - nxinit//2+1)
	txm = float(nxinit-(nxinit//2+1) - radius)
	txl = float(radius - nxinit//2+1)
	for im in xrange(nima):
		data[im] = get_im(Tracker["constants"]["stack"], lpartids[im])
		phi,theta,psi,sx,sy = partstack[im][0], partstack[im][1], partstack[im][2], partstack[im][3], partstack[im][4]
		if( Tracker["constants"]["CTF"] and Tracker["applyctf"] ):
			ctf_params = data[im].get_attr("ctf")
			st = Util.infomask(data[im], mask2D, False)
			data[im] -= st[0]
			data[im] = filt_ctf(data[im], ctf_params)
			data[im].set_attr('ctf_applied', 1)
		if preshift:
			data[im] = fshift(data[im], sx, sy)
			set_params_proj(data[im],[phi,theta,psi,0.0,0.0])
		#oldshifts[im] = [sx,sy]
		#  resample will properly adjusts shifts and pixel size in ctf
		data[im] = resample(data[im], shrinkage)
		#  We have to make sure the shifts are within correct range, shrinkage or not
		set_params_proj(data[im],[phi,theta,psi,max(min(sx*shrinkage,txm),txl),max(min(sy*shrinkage,txm),txl)])
		#  For local SHC set anchor
		#if(nsoft == 1 and an[0] > -1):
		#  We will always set it to simplify the code
		set_params_proj(data[im],[phi,theta,psi,0.0,0.0], "xform.anchor")
	assert( nxinit == data[0].get_xsize() )  #  Just to make sure.
	#oldshifts = wrap_mpi_gatherv(oldshifts, main_node, MPI_COMM_WORLD)
	return data, oldshifts


def metamove(projdata, oldshifts, Tracker, partids, partstack, outputdir, procid, myid, main_node, nproc):
	from applications import slocal_ali3d_base, sali3d_base
	#  Takes preshrunk data and does the refinement as specified in Tracker
	#
	#  Will create outputdir
	#  Will write to outputdir output parameters: params-chunk0.txt and params-chunk1.txt
	from utilities  import get_input_from_string
	shrinkage = float(Tracker["nxinit"])/float(Tracker["constants"]["nnxo"])
	if(myid == main_node):
		#  Create output directory
		log = Logger(BaseLogger_Files())
		log.prefix = os.path.join(outputdir)
		cmd = "mkdir "+log.prefix
		cmdexecute(cmd)
		log.prefix += "/"
		ref_vol = get_im(Tracker["constants"]["refvol"])
		nnn = ref_vol.get_xsize()
		if(Tracker["nxinit"] != nnn ):
			# Good enough?
			ref_vol = Util.window(rot_shift3D(ref_vol,scale=shrinkage),Tracker["nxinit"],Tracker["nxinit"],Tracker["nxinit"])
			#from fundamentals import resample
			#ref_vol = resample(ref_vol, shrinkage)
	else:
		log = None
		ref_vol = model_blank(Tracker["nxinit"], Tracker["nxinit"], Tracker["nxinit"])
	mpi_barrier(MPI_COMM_WORLD)
	bcast_EMData_to_all(ref_vol, myid, main_node)
	#  
	#  Compute current values of some parameters.
	Tracker["radius"] = int(Tracker["constants"]["radius"] * shrinkage +0.5)
	if(Tracker["radius"] < 2):
		ERROR( "ERROR!!   lastring too small  %f    %f   %d"%(Tracker["radius"], Tracker["constants"]["radius"]), "sxmeridien",1, myid)
	Tracker["lowpass"] = float(Tracker["icurrentres"])/float(Tracker["nxinit"])
	if( Tracker["state"] == "LOCAL" or Tracker["state"][:-1] == "FINAL"):
		Tracker["delta"] = "2.0"
		Tracker["ts"]    = "2.0"
	else:
		delta = "%f  "%min(round(degrees(atan(0.5/Tracker["lowpass"]/Tracker["radius"])), 2), 3.0)
		Tracker["delta"] = ""
		for i in xrange(len(get_input_from_string(Tracker["xr"]))):  Tracker["delta"] += delta
	Tracker["pixercutoff"] = get_pixercutoff(Tracker["radius"], float(get_input_from_string(Tracker["delta"])[0]), 0.5)


	if(Tracker["delpreviousmax"]):
		for i in xrange(len(projdata)):
			try:  projdata[i].del_attr("previousmax")
			except:  pass

	if(myid == main_node):
		print_dict(Tracker,"METAMOVE parameters")
		if Tracker["PWadjustment"] :
			print("                    =>  PW adjustment       :  ",Tracker["PWadjustment"])
		print("                    =>  partids             :  ",partids)
		print("                    =>  partstack           :  ",partstack)

	#  Run alignment command
	if Tracker["local"] : params = slocal_ali3d_base(projdata, get_im(Tracker["refvol"]), \
									Tracker, mpi_comm = MPI_COMM_WORLD, log = log, chunk = 1.0)
	else: params = sali3d_base(projdata, ref_vol, Tracker, mpi_comm = MPI_COMM_WORLD, log = log )

	del log
	#  store params
	if(myid == main_node):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"Executed successfully: ","sali3d_base_MPI, nsoft = %d"%Tracker["nsoft"],"  number of images:%7d"%len(params))
		for i in xrange(len(params)):
			params[i][3] = params[i][3]/shrinkage# + oldshifts[i][0]
			params[i][4] = params[i][4]/shrinkage# + oldshifts[i][1]
		write_text_row(params, os.path.join(Tracker["directory"],"params-chunk%01d.txt"%procid) )
	return  Tracker

def print_dict(dict,theme):
	line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	print(line,theme)
	spaces = "                    "
	for key, value in sorted( dict.items() ):
		if(key != "constants"):  print("                    => ", key+spaces[len(key):],":  ",value)


# 
# - "Tracker" (dictionary) object
#   Keeps the current state of option settings and dataset 
#   (i.e. particle stack, reference volume, reconstructed volume, and etc)
#   Each iteration is allowed to add new fields/keys
#   if necessary. This happes especially when type of 3D Refinement or metamove changes.
#   Conceptually, each iteration will be associated to a specific Tracker state.
#   Therefore, the list of Tracker state represents the history of process.
#   (HISTORY is doing this now)
#   This can be used to restart process from an arbitrary iteration.
#   The program will store the HISTORY in the form of file.
#   
#
def main():

	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	import user_functions
	from applications import MPI_start_end
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc
	from development import do_volume_mrk01
	from logger import Logger, BaseLogger_Files
	import sys
	import os
	import time
	import socket
	
	# ------------------------------------------------------------------------------------
	# PARSE COMMAND OPTIONS
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  [output_directory]  initial_volume  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=an  --center=center_type --fl --aa --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	#parser.add_option("--ir",      		type= "int",   default= 1,			help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--radius",      		type= "int",   default= -1,			help="Outer radius [in pixels] for rotational correlation < int(nx/2)-1 (Please set to the radius of the particle)")
	##parser.add_option("--rs",      		type= "int",   default= 1,			help="step between rings in rotational correlation >0  (set to 1)" ) 
	#parser.add_option("--xr",      		type="string", default= "-1",		help="range for translation search in x direction, search is +/xr (default 0)")
	#parser.add_option("--yr",      		type="string", default= "-1",		help="range for translation search in y direction, search is +/yr (default = same as xr)")
	#parser.add_option("--ts",      		type="string", default= "1",		help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	#parser.add_option("--delta",   		type="string", default= "-1",		help="angular step of reference projections during initialization step (default automatically selected based on radius of the structure.)")
	parser.add_option("--an",      		type="string", default= "-1",		help="angular neighborhood for local searches (phi and theta) (Default exhaustive searches)")
	parser.add_option("--center",  		type="float",  default= -1,			help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	#parser.add_option("--maxit",   		type="int",  	default= 400,		help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--outlier_percentile",type="float",    default= 95,	help="percentile above which outliers are removed every iteration")
	parser.add_option("--iteration_start",type="int",    default= 0,		help="starting iteration for rviper, 0 means go to the most recent one (default).")
	parser.add_option("--CTF",     		action="store_true", default=False,	help="Use CTF (Default no CTF correction)")
	#parser.add_option("--snr",     		type="float",  default= 1.0,		help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",   		type="string", default= "S",		help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",     		type="string", default= "c1",		help="Point-group symmetry of the refined structure")
	#parser.add_option("--npad",    		type="int",    default= 2,			help="padding size for 3D reconstruction (default=2)")
	parser.add_option("--nsoft",    	type="int",    default= 0,			help="Use SHC in first phase of refinement iteration (default=0, to turn it on set to 1)")
	parser.add_option("--startangles",  action="store_true", default=False,	help="Use orientation parameters in the input file header to jumpstart the procedure")

	#options introduced for the do_volume function
	#parser.add_option("--fl",			type="float",	default=0.12,		help="cut-off frequency of hyperbolic tangent low-pass Fourier filter (default 0.12)")
	#parser.add_option("--aa",			type="float",	default=0.1,		help="fall-off of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--inires",		type="float",	default=25.,		help="Resolution of the initial_volume volume (default 25A)")
	parser.add_option("--pwreference",	type="string",	default="",			help="text file with a reference power spectrum (default no power spectrum adjustment)")
	parser.add_option("--mask3D",		type="string",	default=None,		help="3D mask file (default a sphere with radius (nx/2)-1)")


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
	# ------------------------------------------------------------------------------------
	# Initialize MPI related variables
	mpi_init(0, [])

	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	# ------------------------------------------------------------------------------------
	#  INPUT PARAMETERS
	global_def.BATCH = True


	#  Constant settings of the project
	Constants				= {}
	Constants["stack"]        = args[0]
	Constants["rs"]           = 1
	Constants["radius"]       = options.radius
	Constants["an"]           = options.an
	Constants["maxit"]        = 50
	sym = options.sym
	Constants["sym"]          = sym[0].lower() + sym[1:]
	Constants["npad"]         = 2
	Constants["center"]       = options.center
	Constants["pwreference"]  = options.pwreference
	Constants["CTF"]          = options.CTF
	Constants["ref_a"]        = options.ref_a
	Constants["snr"]          = 1.0
	Constants["mask3D"]       = options.mask3D
	Constants["nnxo"]         = -1
	Constants["pixel_size"]   = 1.0
	Constants["refvol"]       = volinit
	Constants["masterdir"]    = masterdir
	Constants["best"]         = 0
	Constants["states"]       = ["INITIAL", "EXHAUSTIVE", "RESTRICTED", "LOCAL", "FINAL1", "FINAL2"]
	#Constants["mempernode"]   = 4.0e9
	#  The program will use three different meanings of x-size
	#  nnxo         - original nx of the data, will not be changed
	#  nxinit       - window size used by the program during given iteration, 
	#                 will be increased in steps of 32 with the resolution
	#
	#  nxstep       - step by wich window size increases
	#
	# Create and initialize Tracker Dictionary with input options
	Tracker					  = {}
	Tracker["constants"]      = Constants
	Tracker["maxit"]          = Tracker["constants"]["maxit"]
	Tracker["radius"]         = Tracker["constants"]["radius"]
	Tracker["xr"]             = ""
	Tracker["yr"]             = "-1"  # Do not change!
	Tracker["ts"]             = 1
	Tracker["an"]             = "-1"
	Tracker["delta"]          = "2.0"
	Tracker["zoom"]           = True
	Tracker["nsoft"]          = options.nsoft
	Tracker["local"]          = False
	Tracker["PWadjustment"]   = ""
	Tracker["upscale"]        = 0.5
	Tracker["applyctf"]       = True  #  Should the data be premultiplied by the CTF.  Set to False for local continues.
	Tracker["refvol"]         = None

	Tracker["nxinit"]         = 64
	Tracker["nxstep"]         = 32
	Tracker["icurrentres"]    = -1
	Tracker["lowpass"]        = 0.4
	Tracker["falloff"]        = 0.2
	Tracker["inires"]         = options.inires  # Now in A, convert to absolute before using
	Tracker["fuse_freq"]      = 50  # Now in A, convert to absolute before using
	Tracker["delpreviousmax"] = False
	Tracker["anger"]          = -1.0
	Tracker["shifter"]        = -1.0
	Tracker["saturatecrit"]   = 0.95
	Tracker["pixercutoff"]    = 2.0
	Tracker["directory"]      = ""
	Tracker["previousoutputdir"] = ""
	#Tracker["eliminated-outliers"] = False
	Tracker["mainiteration"]  = 0
	Tracker["movedback"]      = False
	Tracker["state"]          = Tracker["constants"]["states"][0]

	# ------------------------------------------------------------------------------------

	# Get the pixel size; if none, set to 1.0, and the original image size
	if(myid == main_node):
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		print(line,"INITIALIZATION OF MERIDIEN")



		a = get_im(orgstack)
		nnxo = a.get_xsize()
		if( Tracker["nxinit"] > nnxo ):
			ERROR("Image size less than minimum permitted $d"%Tracker["nxinit"],"sxmeridien",1)
			nnxo = -1
		else:
			if Tracker["constants"]["CTF"]:
				i = a.get_attr('ctf')
				pixel_size = i.apix
				fq = pixel_size/Tracker["fuse_freq"]
			else:
				pixel_size = 1.0
				#  No pixel size, fusing computed as 5 Fourier pixels
				fq = 5.0/nnxo
			del a
	else:
		nnxo = 0
		fq = 0.0
		pixel_size = 1.0
	nnxo = bcast_number_to_all(nnxo, source_node = main_node)
	if( nnxo < 0 ):
		mpi_finalize()
		exit()
	pixel_size = bcast_number_to_all(pixel_size, source_node = main_node)
	fq         = bcast_number_to_all(fq, source_node = main_node)
	Tracker["constants"]["nnxo"]         = nnxo
	Tracker["constants"]["pixel_size"]   = pixel_size
	Tracker["fuse_freq"]    = fq
	del fq, nnxo, pixel_size

	if(Tracker["constants"]["radius"]  < 1):
		Tracker["constants"]["radius"]  = Tracker["constants"]["nnxo"]//2-2
	elif((2*Tracker["constants"]["radius"] +2) > Tracker["constants"]["nnxo"]):
		ERROR("Particle radius set too large!","sxmeridien",1,myid)


	# ------------------------------------------------------------------------------------
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



	if(myid == main_node):
		print_dict(Tracker["constants"], "Permanent settings of meridien")


	#  create a vstack from input stack to the local stack in masterdir
	#  Stack name set to default
	Tracker["constants"]["stack"] = "bdb:"+masterdir+"/rdata"
	# Initialization of stacks
	if(myid == main_node):
		if keepchecking:
			if(os.path.exists(os.path.join(masterdir,"EMAN2DB/rdata.bdb"))):  doit = False
			else:  doit = True
		else:  doit = True
		if  doit:
			if(orgstack[:4] == "bdb:"):	cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+Tracker["constants"]["stack"])
			else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, Tracker["constants"]["stack"])
			cmdexecute(cmd)
			cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", Tracker["constants"]["stack"])
			cmdexecute(cmd)
			keepchecking = False
		total_stack = EMUtil.get_image_count(Tracker["constants"]["stack"])
	else:
		total_stack = 0

	total_stack = bcast_number_to_all(total_stack, source_node = main_node)

	# ------------------------------------------------------------------------------------
	#  INITIALIZATION
	initdir = os.path.join(masterdir,"main000")

	# Create first fake directory main000 with parameters filled with zeroes or copied from headers.  Copy initial volume in.
	doit, keepchecking = checkstep(initdir, keepchecking, myid, main_node)
	if  doit:
		partids = os.path.join(initdir, "indexes.txt")
		if( myid == main_node ):
			cmd = "mkdir "+initdir
			cmdexecute(cmd)
			write_text_file(range(total_stack), partids)
		mpi_barrier(MPI_COMM_WORLD)


		#  store params
		partids = [None]*2
		for procid in xrange(2):  partids[procid] = os.path.join(initdir,"chunk%01d.txt"%procid)
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(initdir,"params-chunk%01d.txt"%procid)
		from random import shuffle
		if(myid == main_node):
			#  split randomly
			ll = range(total_stack)
			shuffle(ll)
			l1 = ll[:total_stack//2]
			l2 = ll[total_stack//2:]
			del ll
			l1.sort()
			l2.sort()
			write_text_file(l1,partids[0])
			write_text_file(l2,partids[1])
			
			if(options.startangles):
				tp_list = EMUtil.get_all_attributes(Tracker["constants"]["stack"], "xform.projection")
				for i in xrange(len(tp_list)):
					dp = tp_list[i].get_params("spider")
					tp_list[i] = [dp["phi"], dp["theta"], dp["psi"], -dp["tx"], -dp["ty"]]
				write_text_row(tp_list, os.path.join(initdir,"params.txt"))
				write_text_row([tp_list[i] for i in l1], partstack[0])
				write_text_row([tp_list[i] for i in l2], partstack[1])
				line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
				print(line,"Executed successfully: Imported initial parameters from the input stack")
				
			else:
				write_text_row([[0,0,0,0,0] for i in xrange(len(l1))], partstack[0])
				write_text_row([[0,0,0,0,0] for i in xrange(len(l2))], partstack[1])
				write_text_row([[0,0,0,0,0] for i in xrange(len(l1)+len(l2))], os.path.join(initdir,"params.txt"))
			
			del l1, l2
				
			# Create reference models for each particle group
			# make sure the initial volume is not set to zero outside of a mask, as if it is it will crash the program
			for procid in xrange(2):
				# make a copy of original reference model for this particle group (procid)
				file_path_viv = os.path.join(initdir,"vol%01d.hdf"%procid)
				cmd = "{} {} {}".format("cp -p", volinit, file_path_viv)
				cmdexecute(cmd)
			    # add small noise to the reference model
				viv = get_im(file_path_viv)
				if(options.mask3D == None):  mask33d = model_circle(Tracker["constants"]["radius"],Tracker["constants"]["nnxo"],\
													Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
				else:                        mask33d = get_im(options.mask3D)
				st = Util.infomask(viv, mask33d, False)
				if( st[0] == 0.0 ):
					viv += (model_blank(Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],1.0) - mask33d)*\
							model_gauss_noise(st[1]/1000.0,Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
					viv.write_image(file_path_viv)
				del mask33d, viv

		mpi_barrier(MPI_COMM_WORLD)


	#  This is initial setting for the first iteration.
	#   It will be exhaustive search with zoom option to establish good shifts.
	Tracker["inires"] = Tracker["constants"]["pixel_size"]/Tracker["inires"]  # This is in full size image units.
	Tracker["icurrentres"] = int(Tracker["constants"]["nnxo"]*Tracker["inires"]+0.5)
	#  Make sure nxinit is at least what it was set to as an initial size.
	Tracker["nxinit"] =  max(Tracker["icurrentres"]*2 + cushion , Tracker["nxinit"])
	if( Tracker["nxinit"] > Tracker["constants"]["nnxo"] ):
			ERROR("Resolution of initial volume at the range of Nyquist frequency for given window and pixel sizes","sxmeridien",1, myid)

	#  Here we need an algorithm to set things correctly
	Tracker["xr"] , Tracker["ts"] = stepali(Tracker["nxinit"] , Tracker["constants"]["nnxo"], Tracker["constants"]["radius"])
	Tracker["previousoutputdir"] = initdir
	subdict( Tracker, {"zoom":True} )

	#  remove projdata, if it existed, initialize to nonsense
	projdata = [[model_blank(1,1)], [model_blank(1,1)]]
	HISTORY = []
	oldshifts = [[],[]]
	
	# ------------------------------------------------------------------------------------
	#  MAIN ITERATION

	mainiteration = 0
	keepgoing = 1

	while(keepgoing):
		mainiteration += 1
		Tracker["mainiteration"] = mainiteration
		#  prepare output directory,  the settings are awkward
		Tracker["directory"]     = os.path.join(masterdir,"main%03d"%Tracker["mainiteration"])

		if(myid == main_node):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"MAIN ITERATION  #%2d   %s  nxinit, icurrentres, resolution  %d   %d  %f"%\
				(Tracker["mainiteration"], Tracker["state"],Tracker["nxinit"],  Tracker["icurrentres"], \
				Tracker["constants"]["pixel_size"]*Tracker["constants"]["nnxo"]/float(Tracker["icurrentres"])))
			print(line,"  mainoutputdir  previousoutputdir  ",Tracker["directory"], Tracker["previousoutputdir"])

			if keepchecking:
				if(os.path.exists(Tracker["directory"] )):
					doit = 0
					print("Directory  ",Tracker["directory"] ,"  exists!")
				else:
					doit = 1
					keepchecking = False
			else:
				doit = 1

			if doit:
				cmd = "{} {}".format("mkdir", Tracker["directory"])
				cmdexecute(cmd)

		mpi_barrier(MPI_COMM_WORLD)
		
		# prepare names of input file names, they are in main directory,
		#   log subdirectories contain outputs from specific refinements
		partids = [None]*2
		for procid in xrange(2):  partids[procid] = os.path.join(Tracker["previousoutputdir"],"chunk%01d.txt"%procid)
		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(Tracker["previousoutputdir"],"params-chunk%01d.txt"%procid)

		mpi_barrier(MPI_COMM_WORLD)
		doit = bcast_number_to_all(doit, source_node = main_node)


		#mpi_finalize()
		#exit()

		#print("RACING  A ",myid)
		outvol = [os.path.join(Tracker["previousoutputdir"],"vol%01d.hdf"%procid) for procid in xrange(2)]

		doit, keepchecking = checkstep(outvol[1], keepchecking, myid, main_node)
		if(myid == main_node):
			if  doit:
				vol = [ get_im(outvol[procid]) for procid in xrange(2) ]
				fuselowf(vol, Tracker["fuse_freq"])
				for procid in xrange(2):  vol[procid].write_image(os.path.join(Tracker["directory"], "fusevol%01d.hdf"%procid) )
				del vol

		mpi_barrier(MPI_COMM_WORLD)

		#  REFINEMENT   ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
		for procid in xrange(2):
			coutdir = os.path.join(Tracker["directory"], "loga%01d"%procid)
			doit, keepchecking = checkstep(coutdir, keepchecking, myid, main_node)
			#  here ts has different meaning for standard and continuous
			Tracker["refvol"] = os.path.join(Tracker["directory"], "fusevol%01d.hdf"%procid)\

			if  doit:
				mpi_barrier(MPI_COMM_WORLD)
				if( Tracker["nxinit"] != projdata[procid][0].get_xsize() ):
					projdata[procid] = []
					projdata[procid], oldshifts[procid] = get_shrink_data(Tracker, Tracker["nxinit"],\
						partids[procid], partstack[procid], myid, main_node, nproc, preshift = False)

				# METAMOVE
				Tracker = metamove(projdata[procid], oldshifts[procid], Tracker, partids[procid], partstack[procid], coutdir, procid, myid, main_node, nproc)

		# Update HISTORY
		HISTORY.append(Tracker.copy())


		#  REFINEMENT  ASSESSMENT  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

		partstack = [None]*2
		for procid in xrange(2):  partstack[procid] = os.path.join(Tracker["directory"], "params-chunk%01d.txt"%procid)

		#  Due to efficiency reasons, we have to check it for the volf, if it does not exist, estimate the resolution and compute volf
		doit, keepchecking = checkstep(os.path.join(Tracker["directory"] ,"volf.hdf"), keepchecking, myid, main_node)

		if  doit:
			#  Projections have to be read
			repeat = True
			nxinit = Tracker["nxinit"]
			while  repeat:
				if( nxinit != projdata[procid][0].get_xsize() ):
					projdata = [[],[]]
					for procid in xrange(2):
						projdata[procid], oldshifts[procid] = get_shrink_data(Tracker, nxinit,\
									partids[procid], partstack[procid], myid, main_node, nproc, preshift = False)

				xlowpass, xfalloff, xcurrentres = compute_resolution(projdata, \
													Tracker["directory"], partids, partstack, \
													Tracker["constants"]["radius"], Tracker["constants"]["nnxo"], \
													Tracker["constants"]["CTF"], Tracker["constants"]["mask3D"], \
													Tracker["constants"]["sym"], myid, main_node, nproc, Tracker["constants"]["pixel_size"])
				if( xcurrentres > (nxinit-cushion)/2 and nxinit < Tracker["constants"]["nnxo"] ):
					nxinit = Tracker["constants"]["nnxo"]
					projdata = [[model_blank(1,1)], [model_blank(1,1)]]
				else:  repeat = False

			#  Make sure these variables are not carried
			del xlowpass, xfalloff, xcurrentres, nxinit
			if( myid == main_node):
				# Carry over chunk information
				for procid in xrange(2):
					cmd = "{} {} {}".format("cp -p", os.path.join(Tracker["previousoutputdir"],"chunk%01d.txt"%procid), \
											os.path.join(Tracker["directory"],"chunk%01d.txt"%procid) )
					cmdexecute(cmd)


			#  PRESENTABLE RESULT                <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
			#  Here I have code to generate presentable results.  IDs and params have to be merged and stored and the overall volume computed.
			if( myid == main_node ):
				pinids = read_text_file(partids[0]) + read_text_file(partids[1])
				params = read_text_row(partstack[0]) + read_text_row(partstack[1])

				assert(len(pinids) == len(params))

				for i in xrange(len(pinids)):
					pinids[i] = [ pinids[i], params[i] ]
				del params
				pinids.sort()

				write_text_file([pinids[i][0] for i in xrange(len(pinids))], os.path.join(Tracker["directory"] ,"indexes.txt"))
				write_text_row( [pinids[i][1] for i in xrange(len(pinids))], os.path.join(Tracker["directory"] ,"params.txt"))
			mpi_barrier(MPI_COMM_WORLD)
			nfsc = read_fsc(os.path.join(Tracker["directory"] ,"fsc.txt"),Tracker["constants"]["nnxo"], myid, main_node)
			Tracker["lowpass"], Tracker["falloff"] = fit_tanh([[float(i)/Tracker["constants"]["nnxo"] for i in xrange(len(nfsc))],nfsc], 0.01)
			del nfsc

			# The next will have to be decided later, i.e., under which circumstances we should recalculate full size volume.
			#  Most likely it should be done only when the program terminates.
			"""
			if( nxinit == Tracker["constants"]["nnxo"] ):
				projdata = getindexdata(Tracker["constants"]["stack"], os.path.join(Tracker["directory"] ,"indexes.txt"), \
						os.path.join(Tracker["directory"] ,"params.txt"), myid, nproc)
			"""
			if( myid == main_node ):\
				volf = 0.5*(get_im(os.path.join(Tracker["directory"] ,"vol0.hdf"))+get_im(os.path.join(Tracker["directory"] ,"vol0.hdf")))
			else:  volf = model_blank(Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"],Tracker["constants"]["nnxo"])
			volf = do_volume_mrk01(volf, Tracker, mainiteration, mpi_comm = MPI_COMM_WORLD)
			if(myid == main_node):
				fpol(volf, Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"], Tracker["constants"]["nnxo"]).write_image(os.path.join(Tracker["directory"] ,"volf.hdf"))



		doit, keepchecking = checkstep(os.path.join(Tracker["directory"] ,"error_thresholds.txt"), keepchecking, myid, main_node)
		if  doit:
			#  ANALYZE CHANGES IN OUTPUT PARAMETERS WITH RESPECT TO PREVIOUS INTERATION  <><><><><><><><><><><><><><><><><><><><><><><><><><><>
			if(myid == main_node):
				anger, shifter = threshold_params_changes( Tracker["directory"], Tracker["previousoutputdir"], 0.95, Tracker["constants"]["sym"])
				write_text_row( [[anger, shifter]], os.path.join(Tracker["directory"] ,"error_thresholds.txt") )
			else:
				anger   = 0.0
				shifter = 0.0
		else:
			if(myid == main_node):
				[anger, shifter] = read_text_row( os.path.join(Tracker["directory"] ,"error_thresholds.txt") )[0]
			else:
				anger   = 0.0
				shifter = 0.0

		if( myid == main_node):
			line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
			print(line,"Maximum displacements at 0.95 level for angular directions  %6.2f  and shifts %6.1f"%(anger, shifter) )

		anger   = bcast_number_to_all(anger,   source_node = main_node)
		shifter = bcast_number_to_all(shifter, source_node = main_node)
		Tracker["anger"]   = anger
		Tracker["shifter"] = shifter


		if( myid == main_node):
			[newlowpass, newfalloff, icurrentres] = read_text_row( os.path.join(Tracker["directory"],"current_resolution.txt") )[0]
		else:
			newlowpass = 0.0
			newfalloff = 0.0
			icurrentres = 0
		newlowpass = bcast_number_to_all(newlowpass, source_node = main_node)
		newlowpass = round(newlowpass,4)
		newfalloff = bcast_number_to_all(newfalloff, source_node = main_node)
		Tracker["falloff"] = round(newfalloff,4)  # For the time being
		icurrentres = bcast_number_to_all(icurrentres, source_node = main_node)


		#mpi_barrier(MPI_COMM_WORLD)
		#mpi_finalize()
		#exit()
		if myid == main_node:  print("   >>> AI  <<<  ",Tracker["mainiteration"] ,icurrentres,Tracker["icurrentres"])

		keepgoing, reset_data, Tracker = AI( icurrentres, Tracker, HISTORY)

		if( keepgoing == 1 ):
			if reset_data :  projdata = [[model_blank(1,1)],[model_blank(1,1)]]
			if(myid == main_node):
				print("  New image dimension :", Tracker["nxinit"])
		else:
			if(myid == main_node):
				print("  Terminating  , the best solution is in the directory  %s"%Tracker["constants"]["best"])
			mpi_barrier(MPI_COMM_WORLD)
			mpi_finalize()


if __name__=="__main__":
	main()

