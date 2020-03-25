#!/usr/bin/env python
# This program performs simple processing of .LST files

# Author: Steven Ludtke, 10/06/14 (sludtke@bcm.edu), modified: May 15, 2017 (Jesus GalazMontoya)
# Copyright (c) 2014- Baylor College of Medicine
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

from past.utils import old_div
from builtins import range
from EMAN2 import *
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\nproclst.py [options] <lst 1> <lst 2> ... \nSimple manipulations of LST files. If your goal is to produce an actual image file rather than the
sort of virtual stack represented by .lst files, use e2proc2d.py or e2proc3d.py instead. Those other programs will treat LST files as normal image files for input.\n."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
#	parser.add_argument("--average", action="store_true", help="Averages all input images (without alignment) and writes a single output image")
	
	parser.add_argument("--create", type=str, default=None, help="to use this option, the input files should be image files. Specify an .lst or .lsx file to create here (e.g., --create mylst.lst) with references to all of the images in the inputs.")
	parser.add_argument("--eosplit", action="store_true", help="Will generate _even and _odd .lst files for each specified input .lst file")

	parser.add_argument("--dereforig", type=str, default=None, help="Extract the data_source and data_n parameters from each image in the file and create a new .lst file referencing the original image(s)")

	parser.add_argument("--exclude", type=str, default=None, help="only works if --create is supplied. comma-separated list of indexes from the input file(s) to EXCLUDE from the created .lst file.")

	#parser.add_argument("--first", type=int, default=0, help="Default=0 (first image index in input(s)). This will be the first particle index in the input images to put in the output lsx/lst file.")

	parser.add_argument("--include", type=str, default=None, help="only works if --create is supplied. comma-separated list of indexes to take from the input file(s) to INCLUDE in the created .lst file. if you have the list of indexes to include in a .txt file, you can provide it through --list.")
	parser.add_argument("--inplace", action="store_true", default=False, help="only works with --create. if the stack specified in --create already exists, this will prevent appending to it. rather, the file will be modified in place.")

	#parser.add_argument("--last", type=str, default=-1, help="Default=-1 (last image index in input (s)). This will be the first particle index in the input images to put in the output lsx/lst file.")
	parser.add_argument("--list", type=str, default=None, help="only works if --create is supplied. .txt file with a list of indexes (one per line/row) to take from the input file(s) to INCLUDE in the created .lst file.")

	parser.add_argument("--merge", type=str, default=None, help="Specify the output name here. This will concatenate all of the input .lst files into a single output")
	parser.add_argument("--mergesort", type=str, default=None, help="Specify the output name here. This will merge all of the input .lst files into a single (resorted) output")
	parser.add_argument("--mergeeo", action="store_true", default=False, help="Merge even odd lst.")
	parser.add_argument("--minhisnr", type=float, help="Integrated SNR from 1/10-1/4 1/A must be larger than this",default=-1,guitype='floatbox', row=8, col=1)
	parser.add_argument("--minlosnr", type=float, help="Integrated SNR from 1/200-1/20 1/A must be larger than this",default=-1,guitype='floatbox', row=8, col=0)
	parser.add_argument("--mindf", type=float, help="Minimum defocus",default=-1,guitype='floatbox', row=8, col=1)
	parser.add_argument("--maxdf", type=float, help="Maximum defocus",default=-1,guitype='floatbox', row=8, col=0)
	
	parser.add_argument("--numaslist", type=str, default=None, help="extract the particle indexes (numbers) only from an lst file into a text file (one number per line).")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--range", type=str, default=None, help="Range of particles to use. Works only with --create option. Input of 0,10,2 means range(0,10, step=2).")
	parser.add_argument("--retype", type=str, default=None, help="If a lst file is referencing a set of particles from particles/imgname__oldtype.hdf, this will change oldtype to the specified string in-place (modifies input files)")
	parser.add_argument("--refile", type=str, default=None, help="similar to retype, but replaces the full filename of the source image file with the provided string")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higher number means higher level of verboseness",default=1)


	(options, args) = parser.parse_args()

	if len(args)<1 :
		parser.error("At least one lst file required")
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	#if options.numaslist != None:
	if options.eosplit:
		for inp in args:
			if inp[-4:].lower()!=".lst" : continue
			lin=LSXFile(inp,True)
			ename="{}_even.lst".format(inp[:-4])
			oname="{}_odd.lst".format(inp[:-4])
			try: os.unlink(ename)
			except: pass
			try: os.unlink(oname)
			except: pass
			loute=LSXFile(ename,False)
			louto=LSXFile(oname,False)
			for i in range(len(lin)):
				imt=lin.read(i)
				if i%2: louto.write(-1,imt[0],imt[1],imt[2])
				else: loute.write(-1,imt[0],imt[1],imt[2])
		print("Generated: ",ename,oname)
		sys.exit(0)

	if options.numaslist:
		out=open(options.numaslist,"w")

		for f in args:
			lst=LSXFile(f,True)
			for i in range(len(lst)):
				out.write("{}\n".format(lst[i][0]))

	if options.dereforig:
		newlst=LSXFile(options.dereforig)

		for f in args:
			n=EMUtil.get_image_count(f)
			for i in range(n):
				im=EMData(f,i,True)
				# It shouldn't be possible for this to go infinitely, or there would have been a problem on the previous line
				while im["data_source"][-4:]==".lst" : im=EMData(im["data_source"],im["data_n"],True)
				newlst.write(-1,im["data_n"],im["data_source"])
				if options.verbose>1 : print("{},{} -> {},{}".format(f,i,im["data_source"],im["data_n"]))
		
		print("exiting after --dereforig")
		sys.exit(0)

	if options.create:

		if '.lst' not in options.create and '.lsx' not in options.create:
			print("\nERROR: the extension of the output file in --create must be .lst or .lsx")
			sys.exit(1)

		lst=LSXFile(options.create,False)
		
		if options.mergeeo:
			print("Merging two image stacks...")
			if len(args)!=2:
				print("Error: Need two inputs...")
				exit()
			n0=EMUtil.get_image_count(args[0])
			n1=EMUtil.get_image_count(args[1])
			n=max(n0,n1)

			if args[0].endswith(".lst"):
				lste=LSXFile(args[0],True)
				lsto=LSXFile(args[1],True)
				fromlst=True
			else:
				fromlst=False

			for i in range(n):
				if fromlst:
					if i<n0:
						ln=lste.read(i)
						lst.write(-1,ln[0],ln[1],ln[2])
					if i<n1:
						ln=lsto.read(i)
						lst.write(-1,ln[0],ln[1],ln[2])
				else:
					if i<n0:
						lst.write(-1,i,args[0])
					if i<n1:
						lst.write(-1,i,args[1])
			lst=None
			sys.exit(1)

		else:
			for f in args:
				n=EMUtil.get_image_count(f)
				if f.endswith(".lst"):
					lstin=LSXFile(f,True)
					fromlst=True
				else:
					fromlst=False
				
				

				indxsinclude = list(range(n)) #by default, assume all particles in input file will be part of output lsx; otherwise, modify indexes to include according to options

				if options.range:
					indxsinclude = eval("range({})".format(options.range))
		
				elif options.exclude:
					indxs=set(range(n))
					indxsexclude = set(options.exclude.split(','))
					indxsinclude = [int(j) for j in indxs-indxsexclude]

				elif options.include:
					indxsinclude = [int(j) for j in options.include.split(',')]
					
				elif options.list:
					ff=open(options.list)
					lines=ff.readlines()
					ff.close()
					
					indxsinclude=[]
					k=0
					for line in lines:
						if line:	#check that the line is not empty
							indxsinclude.append( int(line.replace('\n','')))
						else:
							print("\nWARNING, line {} in {} seems to be empty!".format(k,options.list)) 
						k+=1
				
				if options.verbose :
					print("Processing {} images in {}".format(len(indxsinclude),f))
					

				kk=0
				for i in indxsinclude:
				
					if options.range and i>=n: 
						break
					
					if fromlst:
						ln = lstin.read(i)
						if options.inplace:
							lst.write(kk,ln[0],ln[1],ln[2])
						else:
							lst.write(-1,ln[0],ln[1],ln[2])
					else:
						if options.inplace:
							lst.write(kk,i,f)
						else:
							lst.write(-1,i,f)
					kk+=1

		sys.exit(0)


	if options.retype != None:
		if options.minlosnr>0 or options.minhisnr>0 or options.mindf>0 or options.maxdf>0 :
			print("ERROR: --minlosnr and --minhisnr not compatible with --retype")
			sys.exit(1)

		# if the user provided the leading __ for us, we strip it off and add it back later
		if options.retype[:2]=="__" :
			options.retype=options.retype[2:]

		for f in args:
			if options.verbose : print("Processing ",f)
			lst=LSXFile(f,True)

			a=lst.read(0)
			if a[1][:10]!="particles/" and a[1][:12]!="particles3d/" :
				print("To use the --retype option, the .lst file must reference image files in particles/*")

			if options.verbose>1 :
				b=base_name(a[1])
				print("{} -> {}".format(a[1],b+"__"+options.retype+".hdf"))

			# loop over the images in the lst file
			for i in range(len(lst)):
				im=lst.read(i)
				if "3d" in a[1] : outname="particles3d/{}__{}.hdf".format(base_name(im[1]),options.retype)
				else: outname="particles/{}__{}.hdf".format(base_name(im[1]),options.retype)
				lst.write(i,im[0],outname,im[2])

			lst.normalize()			# clean up at the end

			if options.verbose>1 : print(len(lst)," particles adjusted")

		if options.verbose : print("Done processing {} files".format(len(args)))

	if options.refile != None:
		if options.minlosnr>0 or options.minhisnr>0 or options.mindf>0 or options.maxdf>0 :
			print("ERROR: --minlosnr and --minhisnr not compatible with --refile")
			sys.exit(1)

		for f in args:
			if options.verbose : print("Processing ",f)
			lst=LSXFile(f,True)

			# loop over the images in the lst file
			for i in range(len(lst)):
				im=lst.read(i)
				lst.write(i,im[0],options.refile,im[2])

			lst.normalize()			# clean up at the end

			if options.verbose>1 : print(len(lst)," particles adjusted")

		if options.verbose : print("Done processing {} files".format(len(args)))


	if options.merge!=None:

		if options.minlosnr>0 or options.minhisnr>0 or options.mindf>0 or options.maxdf>0 :
			print("ERROR: --minlosnr and --minhisnr not compatible with --merge. Please use --mergesort instead.")
			sys.exit(1)

		# create/update output lst
		lsto=LSXFile(options.merge)
		ntot=0

		# loop over input files
		for f in args:
			lst=LSXFile(f,True)
			ntot+=len(lst)

			for i in range(len(lst)):
				im=lst.read(i)
				lsto.write(-1,im[0],im[1],im[2])

		if options.verbose : print("{} particles added to {}".format(ntot,options.merge))

	if options.mergesort!=None:
		# create/update output lst
		lsto=LSXFile(options.mergesort)
		ntot=0

		# loop over input files
		ptcls=[]
		pfiles=set()
		for f in args:
			lst=LSXFile(f,True)
			ntot+=len(lst)

			for i in range(len(lst)):
				im=lst.read(i)
				ptcls.append((im[1],im[0],im[2]))
				pfiles.add(im[1])

		ptcls.sort()

		# remove particles in files not meeting our criteria
		if options.minlosnr>0 or options.minhisnr>0 or options.mindf>0 or options.maxdf>0 :
			# the list conversion here is so we are iterating over a copy and not modifying the set while we iterate over it
			for pfile in list(pfiles):
				js=js_open_dict(info_name(pfile))
				ctf=js["ctf"][0]
				js.close()
				r1=int(floor(1.0/(200.0*ctf.dsbg)))		# lowsnr is 200-20 A
				r2=int(ceil(1.0/(20.0*ctf.dsbg)))
				r3=int(floor(1.0/(10.0*ctf.dsbg)))		# hisnr is 10 to 4 A
				r4=int(ceil(1.0/(4.0*ctf.dsbg)))
				losnr=old_div(sum(ctf.snr[r1:r2]),(r2-r1))
				hisnr=old_div(sum(ctf.snr[r3:r4]),(r4-r3))
				if losnr<options.minlosnr or hisnr<options.minhisnr or (options.mindf>0 and ctf.defocus<options.mindf) or (options.maxdf>0 and ctf.defocus>options.maxdf) :
					pfiles.remove(pfile)
					if options.verbose: print(pfile," removed due to SNR or defocus limits")

		nwrt=0
		for i in ptcls:
			if i[0] in pfiles :
				lsto.write(-1,i[1],i[0],i[2])
				nwrt+=1

		if options.verbose :
			if nwrt==ntot : print("{} particles in {}".format(ntot,options.mergesort))
			else : print("{} of {} particles written to {}".format(nwrt,ntot,options.mergesort))

	E2end(logid)


if __name__ == "__main__":
	main()
