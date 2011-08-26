#!/usr/bin/env python

#
# Author: Steven Ludtke, 08/26/2011 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

if len(sys.argv)<2 : 
	print "Usage: mergeproject.py <project dir>\n\nRun this program from the target project directory, and specify the path to another project whose particles you wish to incorporate.\nDeals only with boxed particles & CTF. New sets will need to be made in the current project."

dest=db_open_dict("bdb:.#project")
src=db_open_dict("bdb:%s#project"%sys.argv[1],ro=True)

# get a dictionary containing dictionaries of particle paths
srcptcl=src["global.spr_ptcls_dict"]
destptcl=dest["global.spr_ptcls_dict"]

# This means we're in a new project directory, so we copy a bit more
if destptcl==None : 
	destptcl={}
	for k in ("global.apix","global.microscope_cs","global.microscope_voltage","global.num_cpus"):
		dest[k]=src[k]

# now copy the particle databases from the source project
for f in srcptcl:
	print "process ",f
	destptcl[f]=srcptcl[f]

	# f is a dictionary containing type:path pairs for each original source frame
	for tp in srcptcl[f]:
		# tp iterates over the keys of f
		src2path="bdb:%s/%s"%(sys.argv[1],srcptcl[f][tp][4:])
		src2=db_open_dict(src2path,ro=True)

		dst2=db_open_dict(srcptcl[f][tp])

		print "Copy %s to %s"%(src2path,srcptcl[f][tp])
		# copy the images
		for i in range(len(src2)):
			dst2[i]=src2[i]


dest["global.spr_ptcls_dict"]=destptcl

# now copy CTF parameters
print "Copy e2ctf.bg2d"
src=db_open_dict("bdb:%s#e2ctf.bg2d"%sys.argv[1],ro=True)
dest=db_open_dict("bdb:.#e2ctf.bg2d")
for k in src.keys(): 
	print k
	dest[k]=src[k]

print "Copy e2ctf.im2d"
src=db_open_dict("bdb:%s#e2ctf.im2d"%sys.argv[1],ro=True)
dest=db_open_dict("bdb:.#e2ctf.im2d")
for k in src.keys(): 
	print k
	dest[k]=src[k]

print "Copy e2ctf.parms"
src=db_open_dict("bdb:%s#e2ctf.parms"%sys.argv[1],ro=True)
dest=db_open_dict("bdb:.#e2ctf.parms")
for k in src.keys(): 
	print k
	dest[k]=src[k]

