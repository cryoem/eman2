#!/usr/bin/env python
#
# Author: Steve Ludtke (sludtke@bcm.edu)
# Copyright (c) 2000-2012 Baylor College of Medicine


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
import os, re
from EMAN2 import *
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <stackname1> <stackname2> ...
	This will take selected images from stacks of per-micrograph particles and build "sets", which are pseudo-stack files referencing the images in the original files.
	This supports the old convention where BDB files are used for both particles and stacks as well as the new convention where particles are in HDF files and stacks are in LST format.
	Inputs are a list of micrograph base names, and must be in a "particles" directory. For example, particles/jj1234_ctf_flip.hdf would be specified simply as jj1234. Some
	attempts will be made to correct improper specifications.

	e2buildsets.py dh1234 dh2318 dh7965 --excludebad --setname=myset

	This will look in the particles directory for files such as:
	dh1234_ptcls
	dh1234_ptcls_ctf_flip
	...

	and build a separate set for each type (_ctf_flip, _wienesr, ...)
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="stack_files",help="List of micrograph names", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, nosharedb=True)
#	parser.add_header(name="buildheader", help='Options below this label are specific to e2buildsets', title="### e2buildsets options ###", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--setname",type=str,help="Name of the stack to build", default='my_stack', guitype='strbox',row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--excludebad",action="store_true",help="Exclude bad particles.",default=False, guitype='boolbox',row=4,col=0,rowspan=1,colspan=1)
	parser.add_argument("--allparticles",action="store_true",help="Will process all particle stacks stored in the particles subdirectory (if specified, list of files will be ignored)",default=False, guitype='boolbox',row=1, col=0)
	parser.add_argument("--withflipped",action="store_true",help="Only include images with phase-flipped counterparts!",default=False,guitype='boolbox', row=4, col=1, rowspan=1, colspan=1)
	parser.add_argument("--minptcl",type=int,help="Files with fewer than the specified number of particles will be skipped",default=0,guitype='intbox', row=5, col=0)
	parser.add_argument("--minqual",type=int,help="Files with a quality value lower than specified will be skipped",default=0,guitype='intbox', row=5, col=1)
	parser.add_argument("--mindf",type=float,help="Files with a defocus lower than specified will be skipped",default=0,guitype='floatbox', row=6, col=0)
	parser.add_argument("--maxdf",type=float,help="Files with a defocus higher than specified will be skipped",default=20.0,guitype='floatbox', row=6, col=1)
	parser.add_argument("--minlosnr",type=float,help="Integrated SNR from 1/200-1/20 1/A must be larger than this",default=0,guitype='floatbox', row=8, col=0)
	parser.add_argument("--minhisnr",type=float,help="Integrated SNR from 1/10-1/4 1/A must be larger than this",default=0,guitype='floatbox', row=8, col=1)
	parser.add_argument("--minbfactor",type=float,help="Files with a B-factor lower than specified will be skipped",default=0.0,guitype='floatbox', row=7, col=0)
	parser.add_argument("--maxbfactor",type=float,help="Files with a B-factor higher than specified will be skipped",default=5000.0,guitype='floatbox', row=7, col=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	# If allparticles, list of files is ignored !
	if options.allparticles:
		args=[base_name("particles/"+i) for i in os.listdir("particles") if "__" not in i and i[0]!="." and ".hed" not in i ]
		args.sort()
		print "%d particle stacks identified"%len(args)
	else :
		# refactor the arguments in case someone gave us a full specification
		for i in range(len(args)):
			args[i]=base_name(args[i]).replace("bdb:","").split("_ctf")[0]	# This should give us the true base name

	try:
		nimg = EMUtil.get_image_count("particles/{}__ptcls.hdf".format(args[0]))
		basetype="__ptcls"
	except :
		try:
			nimg = EMUtil.get_image_count("particles/{}_ptcls.hdf".format(args[0]))
			basetype="_ptcls"
		except:
			nimg = EMUtil.get_image_count("particles/{}.hdf".format(args[0]))
			basetype=""


	# remove any files that don't have enough particles from the list
	if options.minptcl>0 :
		print "Filtering by particle count"
		args=[i for i in args if imcount("particles/{}{}.hdf".format(i,basetype))>=options.minptcl]
		if options.verbose: print "{} stacks after minptcl filter".format(len(args))

	# remove files with quality too low
	if options.minqual>0 :
		print "Filtering by quality"
		outargs=[]
		for i in args:
			try:
				if js_one_key(info_name(i+".hdf"),"quality")>=options.minqual : outargs.append(i)
			except:
				traceback.print_exc()
				print "Unknown quality for {},{}, including it".format(i,info_name(i+".hdf"))
				outargs.append(i)

		args=outargs

	# remove files without phase flipped particles
	if options.withflipped :
		print "Insuring that files have phase flipped particles"
		ptcls=[i for i in os.listdir("particles") if i[0]!="."]
		args=[i for i in args if i+"__ctf_flip_hp.hdf"in ptcls or i+"__ctf_flip.hdf" in ptcls]	# Not super-efficient, but functional

	print "Filtering by Defocus and B-factor"
	ctfmsg=0
	outargs=[]
	errsets=[]
	for i in args:
		try:
			ctf=js_one_key(info_name(i+".hdf"),"ctf")[0]
			r1=int(floor(1.0/(200.0*ctf.dsbg)))
			r2=int(ceil(1.0/(20.0*ctf.dsbg)))
			r3=int(floor(1.0/(10.0*ctf.dsbg)))
			r4=int(ceil(1.0/(4.0*ctf.dsbg)))
			losnr=sum(ctf.snr[r1:r2])/(r2-r1)
			hisnr=sum(ctf.snr[r3:r4])/(r4-r3)
			if ctf.defocus>=options.mindf and ctf.defocus<=options.maxdf and ctf.bfactor>=options.minbfactor and ctf.bfactor<=options.maxbfactor and losnr>options.minlosnr and hisnr>options.minhisnr : outargs.append(i)
			if options.verbose > 1: print "{}<{}<{}   {}<{}<{}   {}>{}   {}>{}".format(options.mindf,ctf.defocus,options.maxdf,options.minbfactor,ctf.bfactor,options.maxbfactor, losnr,options.minlosnr,hisnr,options.minhisnr)
		except:
			if options.verbose>2 : traceback.print_exc()
			errsets.append(i)
			ctfmsg+=1
			outargs.append(i)

	if len(errsets)>0 : 
		print "Warning, ",len(errsets)," images were missing CTF information, and were included irrespective of specified CTF limits."
		if options.verbose>1 : print "Missing CTF images were: ",errsets

	if ctfmsg: print "Included {} images with undefined CTF".format(ctfmsg)

	args=outargs
	if len(args)==0 :
		print "ERROR: No images left to include after applying filters!"
		sys.exit(1)

	print "%d files to include in processing after filters"%len(args)

	logid=E2init(sys.argv)

	# identify particle groups present
	ptcls=[i for i in os.listdir("particles") if i[0]!="."]

	groups=None
	for f in args:
		group=set([i.replace(f,"").rsplit(".",1)[0] for i in ptcls if f in i])		# files matching current name with i removed
		# groups includes only those filetypes common to ALL specified files
		if groups==None: groups=group
		else: groups.intersection_update(group)

	print "Making sets for the following types: ",
	for i in groups: print "'{}' ".format(i),
	print ""

	totptcl=0
	totbad=0
	lsx={}
	for t in groups:
		try: os.unlink("sets/{}{}.lst".format(options.setname,t))		# we remake each set from scratch
		except: pass
		lsx[t]=LSXFile("sets/{}{}.lst".format(options.setname,t))

	for f in args:
		nimg=EMUtil.get_image_count("particles/{}{}.hdf".format(f,basetype))

		if options.excludebad :
			try : bad=set(js_one_key(info_name(f+".hdf"),"sets")["bad_particles"])
			except :
				if options.verbose : print "No badlist for ",f
				bad=set()
		else : bad=set()
		if options.verbose>1 : print "File: {} -> {} particles - {} bad".format(f,nimg,len(bad))
		totptcl+=nimg-len(bad)
		totbad+=len(bad)

		for t in groups:
			for i in xrange(nimg):
				if i not in bad : lsx[t].write(-1,i,"particles/{}{}.hdf".format(f,t))

	print "Done - {} particles included".format(totptcl),
	if totbad>0 : print ". {} excluded as bad.".format(totbad)
	else: print ""

	E2end(logid)

def imcount(fsp):
	try: return EMUtil.get_image_count(fsp)
	except: return 0

if __name__ == "__main__":
	main()
