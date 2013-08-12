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
	usage = """prog [options] <classes file>
	This program will go through a standard classes file (containing EMAN2-generated class-averages), and extract the
particles associated with one or more class-averages. A similar task can be performed graphically with e2evalparticles.py.
This program is used by e2refinemulti to extract particles associated with each of the output volumes for subsequent
single-model refinement.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--refinemulti",action="store_true",help="Extracts particles based on the model_id header value in each class-average, normally produced by e2refinemulti",default=False)
	parser.add_argument("--setname",type=str,help="Name of the stack to build", default=None)
	parser.add_argument("--classlist",type=str,help="Filename of a text file containing a (comma or whitespace separated) list of class average numbers to operate on. ", default=None)
	parser.add_argument("--excludebad",action="store_true",help="Excludes the particles from the generated set(s). They are included by default.",default=False)
	parser.add_argument("--noderef",action="store_true",help="If particle file was .lst, normally the output .lst will reference the original image file. With this option, the output will reference the .lst file instead, creating a lst pointing to another lst.",default=False)
	parser.add_argument("--sort",action="store_true",help="If set, output .lst file will be sorted. The default is to leave the output grouped by class-average.",default=False)
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if not options.refinemulti and options.classlist==None :
		print "Please specify one of --refinemulti or --classlist=<listfile>"
		sys.exit(1)
		
	try:
		ncls=EMUtil.get_image_count(args[0])
		if ncls<1 : raise Exception
	except:
		print "Error, no class-averages found"
		traceback.print_exc()
		sys.exit(1)


	logid=E2init(sys.argv)

	if options.refinemulti:
		# possible that the first few images may not exist. We look for the first good one:
		for i in xrange(ncls):
			try :
				hdr=EMData(args[0],i,True)
				break
			except : continue
		
		# find the existing set/stack containing the particle data used to make the averages
		inset=hdr["class_ptcl_src"]		# theoretically this could be different for different class-averages, but in practice no program does that
		
		if inset.lower()[:4]=="bdb:" :
			print "Sorry, this program only works with EMAN2.1+, and cannot deal with BDB style sets"
			sys.exit(1)
		
		# This seems a bit odd, as after this point, inset could be either a string or an LSXFile object, but it is useful later
		if not options.noderef and inset[-4:]==".lst" : inset=LSXFile(inset)
		
		outlst={}
		for c in xrange(ncls):
			try : h=EMData(args[0],c,True)
			except:
				if options.verbose>0 : print "Bad class-average: ",c
				continue
			
			# this is a list of all of the particle indices from the input set (inset). 
			# We may need to dereference these to the original file if inset is already a .lst file (as it normally will be)
			ptcl=[]
			try : ptcl.extend(h["class_ptcl_idxs"])
			except: pass
			if not options.excludebad :
				try: ptcl.extend(h["exc_class_ptcl_idxs"])
				except: pass
				
			if len(ptcl)==0 :
				if options.verbose>0 : print "No particles in class-average: ",c
				continue

			# this one is fatal, since we should only have gotten here with good averages
			try : mdl=h["model_id"]
			except:
				print "No model_id in class average {}. Was this classes file created with e2refinemulti.py ?".format(c)
				sys.exit(1)
				
			if not outlst.has_key(mdl) :
				if options.setname!=None :
					if options.setname[-4:]==".lst" : fsp="{}_m{}.lst".format(options.setname[:-4],mdl)
					else: fsp="{}_m{}.lst".format(options.setname,mdl)
				else :
					try : fsp="{}_m{}.lst".format(inset.rsplit(".",1)[0],mdl)
					except : fsp="{}_m{}.lst".format(inset.path.rsplit(".",1)[0],mdl)
				try: os.unlink(fsp)		# make sure we start with a new file, since we're appending
				except: pass
				outlst[mdl]=LSXFile(fsp)

			#### This is where we actually generate the new sets
			for p in sorted(ptcl):
				if isinstance(inset,str) :
					outlst[mdl].write(-1,p,inset)
				else :
					nextf,extf,com=inset.read(p)	# read the original information for this particle from the old set
					outlst[mdl].write(-1,nextf,extf)
					
		if options.sort :
			for k in outlst:
				ptcls=[]
				for i in xrange(outlst[k].n):
					nextf,extf,com=outlst[k].read(i)
					ptcls.append((extf,nextf))
				ptcls.sort()
				
				for i,v in enumerate(ptcls):
					outlst[k].write(i,v[1],v[0])
					
		if options.verbose>0 :
			print "Output files:"
			for k in sorted(outlst.keys()) : print "model_id = {} : {} ({})".format(k,outlst[k].path,outlst[k].n)
					
					
	else :
		print "Sorry, this mode is not yet complete. Please email sludtke@bcm.edu"
		sys.exit(1)

	E2end(logid)


if __name__ == "__main__":
	main()
