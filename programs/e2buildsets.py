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

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <stackname1> <stackname2> ...
	This will take selected images from stacks of per-micrograph particles and build "sets", which are pseudo-stack files referencing the images in the original files.
	This supports the old convention where BDB files are used for both particles and stacks as well as the new convention where particles are in HDF files and stacks are in LST format.
	Inputs are a list of micrograph base names, and must be in a "particles" directory. For example, particles/jj1234_ctf_flip.hdf would be specified simply as jj1234. Some
	attempts will be made to correct improper specifications.

	e2buildsets.py dh1234 dh2318 dh7965 --excludebad

	This will look in the particles directory for files such as:
	dh1234_ptcls
	dh1234_ptcls_ctf_flip
	...

	and build a separate set for each type (_ctf_flip, _wienesr, ...)
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="stack_files",help="List of micrograph names", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, nosharedb=True)
	parser.add_header(name="buildheader", help='Options below this label are specific to e2buildsets', title="### e2buildsets options ###", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--setname",type=str,help="Name of the stack to build", default='my_stack', guitype='strbox',row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--filetype",help="Type of file",default='lst',guitype='combobox',choicelist='["lst","bdb"]',row=3,col=0,rowspan=1,colspan=1)
	parser.add_argument("--excludebad",action="store_true",help="Exclude bad particles.",default=False, guitype='boolbox',row=4,col=0,rowspan=1,colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()


	if options.filetype.lower() not in ("bdb","lst") :
		print "Only BDB and LST filetypes are accepted"
		sys.exit(1)

	# refactor the arguments in case someone gave us a full specification
	for i in range(len(args)):
		args[i]=base_name(args[i]).replace("bdb:","").split("_ctf")[0]	# This should give us the true base name

	print args

	logid=E2init(sys.argv)

	# identify particle groups present
	if options.filetype=="bdb" : ptcls=db_list_dicts("bdb:particles")
	else : ptcls=[i for i in os.listdir("particles") if i[0]!="."]

	groups=None
	for f in args:
		group=set([i.replace(f,"").rsplit(".",1)[0] for i in ptcls if f in i])		# files matching current name with i removed
		# groups includes only those filetypes common to ALL specified files
		if groups==None: groups=group
		else: groups.intersection_update(group)

	print "Making sets for the following types: ",
	for i in groups: print "'{}' ".format(i),
	print ""

	if options.filetype.lower()=="bdb":

		for n in range(0,len(args),10):			# process in groups of 10 for efficiency
			for t in groups:
				cmd="e2bdb.py "
				for f in args[n:n+10]:
					if "_ctf" in f : f=f.split("_ctf")[0]						# if the user provided ctf modified names, we clean up the mess for them
					if options.excludebad: cmd+="bdb:particles#%s%s?exclude.%s "%(f,t,f)
					else : cmd+="bdb:particles#%s%s "%(f,t)

				cmd+="--appendvstack=bdb:sets#%s%s"%(options.setname,t)

				launch_childprocess(cmd)
				#print cmd

	elif options.filetype.lower()=="lst" :
		totptcl=0
		lsx={}
		for t in groups:
			try: os.unlink("sets/{}{}.lst".format(options.setname,t))		# we remake each set from scratch
			except: pass
			lsx[t]=LSXFile("sets/{}{}.lst".format(options.setname,t))

		for f in args:
			try:
				nimg = EMUtil.get_image_count("particles/{}__ptcls.hdf".format(f))
				basetype="__ptcls"
			except :
				try:
					nimg = EMUtil.get_image_count("particles/{}_ptcls.hdf".format(f))
					basetype="_ptcls"
				except:
					nimg = EMUtil.get_image_count("particles/{}.hdf".format(f))
					basetype=""

			if options.excludebad :
				try : bad=set(js_open_dict(info_name(f))["sets"]["bad_particles"])
				except :
					if options.verbose : print "No badlist for ",f
					bad=set()
			else : bad=set()
			if options.verbose>1 : print "File: {} -> {} particles - {} bad".format(f,nimg,len(bad))
			totptcl+=nimg-len(bad)

			for t in groups:
				for i in xrange(nimg):
					if i not in bad : lsx[t].write(-1,i,"particles/{}{}.hdf".format(f,t))

		print "Done - {} particles total".format(totptcl)
	E2end(logid)


if __name__ == "__main__":
	main()
