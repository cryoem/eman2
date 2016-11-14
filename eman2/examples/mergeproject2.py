#!/usr/bin/env python

#
# Author: Steven Ludtke, 01/09/2013 (sludtke@bcm.edu)
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

from EMAN2 import *
import sys
import time

print "WARNING: this program was designed for EMAN2.07. It will not work properly with EMAN2.1+ projects. With EMAN2.1 projects, you can simply combine the files from micrographs, particles and info into a single project directory. No specific programs need to be used"

if len(sys.argv)<2 : 
	print "Usage: mergeproject.py <project dir> ... [--minquality=<0-9>]\n\nThis is a new version of mergeproject desinged to work with e2projectmanager projects rather than e2workflow project. \n\nRun this program from the target project directory, and specify the path to another project whose particles you wish to incorporate.\nDeals only with boxed particles & CTF. New sets will need to be made in the current project.\nIf minquality is provided, only particles with a quality at least as large as the specified value will be copied."

minq=0
for i in sys.argv[1:]:
	if i[:13]=="--minquality=" :
		minq=int(i[13:])
	elif i[0]=="-" :
		print "Usage: mergeproject2.py <project dir> ... [--minquality=<0-9>]"
		print "Unknown option: ",i
		sys.exit(1)

for source in sys.argv[1:]:
	if source[0]=="-" : continue
	print "================= Processing project : ",source

	dest=db_open_dict("bdb:.#project")
	src=db_open_dict("bdb:%s#project"%source,ro=True)

	# get a list of particle dictionaries
	srcptcl=db_list_dicts("bdb:%s/particles"%source)
	destptcl=db_list_dicts("bdb:particles")

	# This means we're in a new project directory, so we copy a bit more
	if len(destptcl)==0 : 
		for k in ("global.apix","global.microscope_cs","global.microscope_voltage","global.num_cpus"):
			dest[k]=src[k]

	# filter out images with low quality if requested
	if minq>0 :
		src2=db_open_dict("bdb:%s#e2ctf.parms"%source,ro=True)
		goodkeys=[]
		for k in src2.keys():
			if src2[k][3]<minq :
				print "%s excluded with quality %d"%(k,src2[k][3])
				srcptcl=[i for i in srcptcl if not k in i]
			else: goodkeys.append(k)

	# now copy the particle databases from the source project
	for f in srcptcl:
		print "process ",f

		# This copies the image files
		src=db_open_dict("bdb:%s/particles#%s"%(source,f),ro=True)
		dest=db_open_dict("bdb:particles#%s"%f)
		for i in xrange(len(src)): dest[i]=src[i]

	# now copy CTF parameters
	print "Copy e2ctf.bg2d"
	src=db_open_dict("bdb:%s#e2ctf.bg2d"%source,ro=True)
	dest=db_open_dict("bdb:.#e2ctf.bg2d")
	for k in goodkeys: 
		print k
		dest[k]=src[k]

	print "Copy e2ctf.im2d"
	src=db_open_dict("bdb:%s#e2ctf.im2d"%source,ro=True)
	dest=db_open_dict("bdb:.#e2ctf.im2d")
	for k in goodkeys: 
		print k
		dest[k]=src[k]

	print "Copy e2ctf.parms"
	src=db_open_dict("bdb:%s#e2ctf.parms"%source,ro=True)
	dest=db_open_dict("bdb:.#e2ctf.parms")
	for k in goodkeys: 
		print k
		dest[k]=src[k]
	
	# Copy bad particles list
	try:
		print "Copy bad particles lists"
		src=db_open_dict("bdb:%s#select"%source,ro=True)
		dest=db_open_dict("bdb:.#select")
		for k in goodkeys: 
			print k
			dest[k]=src[k]
	except:
		print "Error in copying bad particle list. No bad particles marked ?"
