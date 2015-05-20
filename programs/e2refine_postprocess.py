#!/usr/bin/env python

#
# Author: Steve Ludtke 06/20/2013 (sludtke@bcm.edu)
# Copyright (c) 2013- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
import time


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	This program performs the post-processing steps for e2refine_easy. This gives an easy way to re-run the post-processing if necessary.

"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--even", dest="even", type=str,default=None, help="The filename of the map from the even 1/2 of the data")
	parser.add_argument("--odd", dest="odd", type=str,default=None, help="The filename of the map from the odd 1/2 of the data")
	parser.add_argument("--output", dest="output", type=str,default=None, help="Filename for the final averaged/filtered result.")
	parser.add_argument("--mass", default=0, type=float,help="The rough mass of the particle in kilodaltons, used to run normalize.bymass. Due to resolution effects, not always the true mass.")
	parser.add_argument("--restarget", default=5, type=float,help="The specified target resolution to avoid underfiltering artifacts.")
	parser.add_argument("--setsf",type=str,help="Force the structure factor to match a 'known' curve prior to postprocessing (<filename>, auto or none). default=none",default="none")
	parser.add_argument("--iter", dest = "iter", type = int, default=6, help = "Iteration number to generate FSC filenames")
	parser.add_argument("--align",action="store_true",default=False,help="Will do o to e alignment and test for handedness flips. Should not be repeated as it overwrites the odd file with the aligned result.")
	parser.add_argument("--m3dpostprocess", type=str, default=None, help="Default=none. An arbitrary post-processor to run after all other automatic processing.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--automaskexpand", default=-1, type=int,help="Default=boxsize/20. Specify number of voxels to expand mask before soft edge. Only used if automask3d not specified." )
	parser.add_argument("--automask3d", default=None, type=str,help="Default=auto. Specify as a processor, eg - mask.auto3d:threshold=1.1:radius=30:nshells=5:nshellsgauss=5." )
	parser.add_argument("--automask3d2", default=None, type=str,help="Default=None. Specify as a processor. This will be applied to the mask produced by the first automask." )
	parser.add_argument("--underfilter",action="store_true",default=False,help="This will shift the computed Wiener filter to be about 10%% more resolution than has been achieved.")
	parser.add_argument("--sym", dest="sym", type=str,default="c1", help="Symmetry so we can decide how to align the particle.")

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv,options.ppid)

	if options.setsf.lower()=="none" : options.setsf=None

	if options.underfilter: underfilter=":sscale=1.1"
	else: underfilter=""

	if options.m3dpostprocess==None or len(options.m3dpostprocess.strip())==0 : m3dpostproc=""
	else : m3dpostproc="--process "+options.m3dpostprocess

	lpfilt=1.15/max(15.0,options.restarget)	# low-pass for alignment

	### Post-processing ###
	### Even/Odd Alignment
	evenfile=options.even
	oddfile=options.odd
	combfile=options.output
	path=os.path.dirname(combfile)+"/"
	if (path=="/") : path="./"
	if options.align :
		if options.sym.lower() in ("icos","tet","oct") or options.sym[0].lower()=="d" : align="" 	# no alignment with higher symmetries
		elif options.sym[0].lower()=="c" and options.sym[1]!="1" : align=align=" --ralignzphi={path}tmp0.hdf".format(path=path)		# z/phi alignment only
		else: align="--alignref={path}tmp0.hdf --align=refine_3d".format(path=path)	# full 3-D alignment for C1

		cmd="e2proc3d.py {evenfile} {path}tmp0.hdf --process=filter.lowpass.gauss:cutoff_freq={lpfilt}".format(path=path,evenfile=evenfile,lpfilt=lpfilt)
		run(cmd)
		cmd="e2proc3d.py {oddfile} {path}tmp1.hdf --process=filter.lowpass.gauss:cutoff_freq={lpfilt} {align}".format(path=path,oddfile=oddfile,align=align,lpfilt=lpfilt)
		run(cmd)
		# in case the handedness got swapped due to too much randomization, we need to double-check the inverted handedness in the alignment
		cmd="e2proc3d.py {oddfile} {path}tmp2.hdf --process=filter.lowpass.gauss:cutoff_freq={lpfilt} --process=xform.flip:axis=z {align}".format(path=path,oddfile=oddfile,align=align,lpfilt=lpfilt)
		run(cmd)
		# now we have to check which of the two handednesses produced the better alignment
		# Pick the best handedness
		a=EMData("{path}tmp1.hdf".format(path=path),0)
		b=EMData("{path}tmp2.hdf".format(path=path),0)
		c=EMData("{path}tmp0.hdf".format(path=path),0)
		ca=c.cmp("ccc",a)
		cb=c.cmp("ccc",b)
		if ca<cb :
			try: ali=a["xform.align3d"]
			except: ali=Transform()
			o=EMData(oddfile,0)
			print "correct hand detected ",ali
		else :
			try: ali=b["xform.align3d"]
			except: ali=Transform()
			o=EMData(oddfile,0)
			o.process_inplace("xform.flip",{"axis":"z"})
			print "handedness flip required",ali
		o.transform(ali)

		os.unlink(oddfile)
		o.write_image(oddfile)
		os.unlink("{path}tmp1.hdf".format(path=path))
		os.unlink("{path}tmp2.hdf".format(path=path))
		os.unlink("{path}tmp0.hdf".format(path=path))

	### Unmasked FSC
	cmd="e2proc3d.py {evenfile} {path}fsc_unmasked_{itr:02d}.txt --calcfsc={oddfile}".format(path=path,itr=options.iter,evenfile=evenfile,oddfile=oddfile)
	run(cmd)

	### Filtration & Normalize
	even=EMData(evenfile,0)
	odd=EMData(oddfile,0)
	combined=even+odd
	try: combined["ptcl_repr"]=even["ptcl_repr"]+odd["ptcl_repr"]
	except: pass
	combined.write_image(combfile,0)

	if options.setsf and options.setsf!="none" : 
		setsf="--setsf "+options.setsf
	else: 
		setsf=""

	nx,ny,nz=combined["nx"],combined["ny"],combined["nz"]
	run("e2proc3d.py {combfile} {combfile} {setsf} --process=filter.wiener.byfsc:fscfile={path}fsc_unmasked_{itr:02d}.txt:snrmult=2:maxfreq={maxfreq} --process=normalize.bymass:thr=1:mass={mass}".format(
		path=path,itr=options.iter,mass=options.mass,setsf=setsf,combfile=combfile,maxfreq=1.0/options.restarget))

	# combined2 is the merged file after setsf & initial filtration
	combined2=EMData(combfile,0)
	sigmanz=combined2["sigma_nonzero"]
	apix=combined["apix_x"]

	### Masking
	if options.automask3d==None or len(options.automask3d.strip())==0 :
		# This loop runs automatic masking with real parameters to test if the mask is extending to the edge of the box
		# if it is, it adjusts the threshold and seed parameters to make the mask smaller
		# This is a bit of a hack, and unfortunately a slow way of handling the problem
		radav=1.0
		seeds=24
		itr=0
		while radav>0.001:
			sigmanz*=1.1
			seeds=max(0,seeds-8)
			vol=EMData(combfile,0)
			vol.process_inplace("mask.auto3d",{"threshold":sigmanz*.85,"radius":nx/10,"nshells":int(nx*0.05+.5),"nshellsgauss":int(options.restarget*1.5/apix),"nmaxseed":seeds,"return_mask":1})
			dis=vol.calc_radial_dist(vol["nx"]/2,0,1,0)
			radav=sum(dis[-4:])
			itr+=1
			if itr>5 :
				print "WARNING: failed to achieve a properly isolated volume, FSC artifacts may occur. Often this is caused by an incorrect A/pix value or specifying too large a mask. It could also indicate that the box size is too small."
				sigmanz=combined2["sigma_nonzero"]*1.1
				seeds=0
				break
		
		if options.automaskexpand<0 : options.automaskexpand=int(nx*0.05+0.5)
		# Final automasking parameters
		amask3d="--process mask.auto3d:threshold={thresh}:radius={radius}:nshells={shells}:nshellsgauss={gshells}:nmaxseed={seed}".format(
			thresh=sigmanz*.75,radius=nx/10,shells=options.automaskexpand,gshells=int(options.restarget*1.5/apix),seed=seeds)
		amask3dtight="--process mask.auto3d:threshold={thresh}:radius={radius}:nshells={shells}:nshellsgauss={gshells}:nmaxseed={seed}".format(
			thresh=sigmanz,radius=nx/10,shells=int(options.restarget*1.2/apix),gshells=int(options.restarget*1.2/apix),seed=seeds)
	else:
		amask3d="--process "+options.automask3d
		amask3dtight=amask3d

	combined2=0
	if options.automask3d2==None or len(options.automask3d2.strip())==0 : amask3d2=""
	else : amask3d2="--process "+options.automask3d2

	# this is a terrible hack. mask.auto3d needs the actual data to generate a mask, but the other mask. processors don't in general, and don't have the return_mask option
	if amask3d.split(":")[0]=="--process mask.auto3d" : maskopt=":return_mask=1"
	else: maskopt=" --inputto1"
	
	run("e2proc3d.py {cfile} {path}mask.hdf {mask}{maskopt} {amask3d2}".format(path=path,cfile=combfile,mask=amask3d,amask3d2=amask3d2,maskopt=maskopt))
	run("e2proc3d.py {cfile} {path}mask_tight.hdf {mask}{maskopt} {amask3d2}".format(path=path,cfile=combfile,mask=amask3dtight,amask3d2=amask3d2,maskopt=maskopt))
	
	### Masked tight FSC
	run("e2proc3d.py {evenfile} {path}tmp_even.hdf --multfile {path}mask_tight.hdf".format(evenfile=evenfile,path=path))
	run("e2proc3d.py {oddfile} {path}tmp_odd.hdf --multfile {path}mask_tight.hdf".format(oddfile=oddfile,path=path))

	# New FSC between the two masked volumes, which we will use for the final filter
	cmd="e2proc3d.py {path}tmp_even.hdf {path}fsc_maskedtight_{itr:02d}.txt --calcfsc {path}tmp_odd.hdf".format(path=path,itr=options.iter)
	run(cmd)

	### Masked FSC, and tmp_ files for later filtration
	run("e2proc3d.py {evenfile} {path}tmp_even.hdf --multfile {path}mask.hdf".format(evenfile=evenfile,path=path))
	run("e2proc3d.py {oddfile} {path}tmp_odd.hdf --multfile {path}mask.hdf".format(oddfile=oddfile,path=path))

	# New FSC between the two masked volumes, which we will use for the final filter
	cmd="e2proc3d.py {path}tmp_even.hdf {path}fsc_masked_{itr:02d}.txt --calcfsc {path}tmp_odd.hdf".format(path=path,itr=options.iter)
	run(cmd)

	# _unmasked volumes are filtered
	run("e2proc3d.py {evenfile} {path}threed_even_unmasked.hdf {setsf} --process filter.wiener.byfsc:fscfile={path}fsc_masked_{itr:02d}.txt:snrmult=2{underfilter}:maxfreq={maxfreq}".format(evenfile=evenfile,path=path,itr=options.iter,mass=options.mass,setsf=setsf,underfilter=underfilter,maxfreq=1.0/options.restarget))
	run("e2proc3d.py {oddfile} {path}threed_odd_unmasked.hdf {setsf} --process filter.wiener.byfsc:fscfile={path}fsc_masked_{itr:02d}.txt:snrmult=2{underfilter}:maxfreq={maxfreq}".format(oddfile=oddfile,path=path,itr=options.iter,mass=options.mass,setsf=setsf,underfilter=underfilter,maxfreq=1.0/options.restarget))

	# Technically snrmult should be 1 here, but we use 2 to help speed convergence
	cmd="e2proc3d.py {path}tmp_even.hdf {evenfile} {setsf} --process filter.wiener.byfsc:fscfile={path}fsc_masked_{itr:02d}.txt:snrmult=2{underfilter}:maxfreq={maxfreq} --multfile {path}mask.hdf --process normalize.bymass:thr=1:mass={mass} {postproc}".format(
	evenfile=evenfile,path=path,itr=options.iter,mass=options.mass,setsf=setsf,underfilter=underfilter,maxfreq=1.0/options.restarget,postproc=m3dpostproc)
	run(cmd)

	cmd="e2proc3d.py {path}tmp_odd.hdf {oddfile} {setsf} --process filter.wiener.byfsc:fscfile={path}fsc_masked_{itr:02d}.txt:snrmult=2{underfilter}:maxfreq={maxfreq} --multfile {path}mask.hdf --process normalize.bymass:thr=1:mass={mass} {postproc}".format(
	oddfile=oddfile,path=path,itr=options.iter,mass=options.mass,setsf=setsf,underfilter=underfilter,maxfreq=1.0/options.restarget,postproc=m3dpostproc)
	run(cmd)

	### Refilter/mask
	combined.write_image(combfile,0)	# write the original average back to disk

	# Note that the snrmult=4 below should really be 2 (due to the averaging of the 2 maps), the 4 is a somewhat arbitrary compensation for the .143 cutoff being a bit low
	nx,ny,nz=combined["nx"],combined["ny"],combined["nz"]
	
	# we impose the symmetry in real-space, since this is what people expect
	if options.sym=="c1" : symopt=""
	else: symopt="--sym {}".format(options.sym)
	
	run("e2proc3d.py {combfile} {combfile} {setsf} --process filter.wiener.byfsc:fscfile={path}fsc_masked_{itr:02d}.txt:snrmult=2{underfilter}:maxfreq={maxfreq} --multfile {path}mask.hdf --process normalize.bymass:thr=1:mass={mass} {symopt} {postproc}".format(
		combfile=combfile,path=path,itr=options.iter,mass=options.mass,setsf=setsf,postproc=m3dpostproc,symopt=symopt,underfilter=underfilter,maxfreq=1.0/options.restarget))

	try:
		os.unlink("{path}tmp_even.hdf".format(path=path))
		os.unlink("{path}tmp_odd.hdf".format(path=path))
		#os.system("gzip {path}mask.hdf".format(path=path)
		#os.system("gzip {path}mask_tight.hdf".format(path=path)
	except:
		pass

	E2end(logid)

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print "{}: {}".format(time.ctime(time.time()),command)
#	append_html("<p>{}: {}</p>".format(time.ctime(time.time()),command))

	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print "Error running: ",command
		sys.exit(1)

	return

if __name__ == "__main__":
    main()
