#!/usr/bin/env python

#
# Author: Steven Ludtke, 05/23/2012 (sludtke@bcm.edu)
# Copyright (c) 2000-2013 Baylor College of Medicine
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

from math import *
import time
import os
import sys
import re
import traceback

from EMAN2 import *

def main():
	global debug
	progname = os.path.basename(sys.argv[0])
	usage = """prog <dest folder> 
	
This program will update an EMAN2.0x project to EMAN2.1. Since this is a one-way process, rather than modifying your project in-place, \
this script will create a new complete copy of your existing project in a new folder. This insures you always have the original in case \
the script fails, or there are other bugs in EMAN 2.1. It does mean you will want to insure you have sufficient disk space. 

This program must be run from within an existing EMAN2.0x project directory. It takes a single argument, the name of a new folder in which \
to place the EMAN2.1 copy.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--yes",action="store_true",help="This will skip the 'are you sure' question, and proceed with the conversion",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	(options, args) = parser.parse_args()

	if not os.path.exists("EMAN2DB/project.bdb") :
		print "ERROR: This does not appear to be a valid project directory !"
		sys.exit(1)
	
	if len(args)!=1 :
		print """ERROR: please specify a path to the folder where you want the EMAN2.1
duplicate project to be stored"""
		sys.exit(1)
	
	try:
		if os.path.samefile(".",args[0]) :
			print """ERROR: The specified destination is the same as the current folder. You must specify
the name of a new folder to be created with a copy of this project."""
			sys.exit(1)
	except: pass
	
	print """Warning:
This tool will convert an EMAN2.0x project to be as compatible as possible with the new standards in 
EMAN2.1. Rather than modify the project in-place, it will duplicate the entire project 
including data, with the new project follwing the 2.1 standards. This will leave 
you with the original 2.0x project intact, in case something goes wrong (as it may well 
with this early version of this tool). This means you need to have enough disk space
on the target drive for a copy of this project.

This program will:
- Convert all of the metadata in BDB databases in e2boxercache and the project EMAN2DB directory to corresponding info/* files.
- Convert all images in BDB databases in subdirectories into HDF format
- convert BDB sets into LST files
- NOTE - this program will not currently copy any other files you may have created outside the normal EMAN2 structure

It is strongly suggested that you run 'e2bdb.py -c' prior to running this program.
"""

	if not options.yes :
		a=raw_input("Are you sure you want to proceed (enter YES) ? ")
		if a!="YES" : sys.exit(2)


	dest=os.path.abspath(args[0])
	try: os.makedirs(dest)
	except: pass
	
	try: 
		os.mkdir(dest+"/info")
		os.mkdir(dest+"/sets")
	except: pass

	# The basic 'project database'
	try:
		if options.verbose>0 : print "Converting project database"
		prj =db_open_dict("bdb:.#project",ro=True)
		prjo=js_open_dict("{}/info/project.json".format(dest))
		prjo.update(prj)
		prj.close()
		prjo.close()
	except:
		traceback.print_exc()
		print "Could not convert project database"
	
	# CTF
	try:
		if options.verbose : print "Converting CTF"
		ctf=db_open_dict("bdb:.#e2ctf.parms",ro=True)
		ctfbg=db_open_dict("bdb:.#e2ctf.bg2d",ro=True)
		ctffg=db_open_dict("bdb:.#e2ctf.im2d",ro=True)
		for k in ctf.keys():
			try:
				if options.verbose>1 : print "\t",k
				js=js_open_dict("{}/{}".format(dest,info_name("bdb:particles#"+k)))
				if options.verbose>2 : print "Generating {}/{}".format(dest,info_name("bdb:particles#"+k))
				js.setval("quality",int(ctf[k][-1]),True)
				c=EMAN2Ctf()
				c.from_string(ctf[k][0])
				cfg=ctffg[k]
				cfg.del_attr("micrograph_id")	# unicode problems  :^(
				cbg=ctfbg[k]
				cbg.del_attr("micrograph_id")
				js["ctf"]=[c,ctf[k][1],ctf[k][2],cfg,cbg]
			except:
				print "Warning: CTF conversion failed for ",k
	except:
		traceback.print_exc()
		print "\n\nUnable to convert projects without CTF information(%s). If this is a major problem, please contact sludtke@bcm.edu"%k
		sys.exit(1)
	
	# Boxes
	if options.verbose : print "Converting Boxes"
	boxes=db_open_dict("bdb:e2boxercache#boxes",ro=True)
	try:
		for k in boxes.keys():
			try:
				if options.verbose>1 : print "\t",k
				js=js_open_dict("{}/{}".format(dest,info_name(k)))
				js["boxes"]=boxes[k]
			except:
				print "Warning: Error converting boxes for ",k
	except:
		print "Note: No box locations found to convert"
	
	dl=[i for i in os.listdir(".") if os.path.isdir(i) and i not in ("sets","e2boxercache","EMAN2DB")]

	# handle sets specially
	if options.verbose : print "Converting sets:"
	dcts=db_list_dicts("bdb:sets")
	for dct in dcts:
		try:
			lst=LSXFile("{}/sets/{}.lst".format(dest,dct))
			src=db_open_dict("bdb:sets#{}".format(dct))
			if options.verbose>1 : print "\t",dct
			for i in xrange(len(src)):
				attr=src.get_header(i)
				lst.write(i,attr["data_n"],"particles/{}.hdf".format(attr["data_source"].split("#")[-1].replace("_ctf","__ctf")))
		except:
			print "Unable to convert set: ",dct

	# This handles generic conversion of subdirectories
	for d in dl:
		# First we handle BDBs in subdirectories
		try: os.mkdir("{}/{}".format(dest,d))
		except: pass

		dcts=db_list_dicts("bdb:"+d)
		for dct in dcts:
			try:
				if options.verbose>0 : print "Processing {}/{}".format(d,dct)
				n=EMUtil.get_image_count("bdb:{}#{}".format(d,dct))
				tmp2=db_open_dict("bdb:{}#{}".format(d,dct),ro=True)
				if n==0:
					tmp=js_open_dict("{}/{}/{}.json".format(dest,d,dct))
					tmp.update(tmp2)
				else:
					if d=="particles" : dct2=dct.replace("_ctf","__ctf")
					else : dct2=dct
					for k in xrange(n):
						if options.verbose>1 : 
							print "  {}/{}  \r".format(k,n),
							sys.stdout.flush()
						im=tmp2[k]
						im.write_image("{}/{}/{}.hdf".format(dest,d,dct2),k)
			except:
				print "Error in converting ",dct," (skipping)"

		# Now we handle any regular files that may exist
		fls=[i for i in os.listdir(d) if i!="EMAN2DB" and os.path.isfile(i)]
		for i in fls: 
			try: copyfile("{}/{}".format(d,i),"{}/{}/{}".format(dest,d,i))
			except: print "Unable to copy {}/{}".format(d,i)

	# files in the project directory
	fls=[i for i in os.listdir(".") if os.path.isfile(i)]
	for i in fls: copyfile(i,"{}/{}".format(dest,i))

def copyfile(a,b):
	"""There really should be an os.copy()"""
	fin=file(a,"rb")
	fout=file(b,"wb")
	while True:
		data=fin.read(16777216)
		if len(data)==0 : break
		fout.write(data)

if __name__ == "__main__":
	main()
