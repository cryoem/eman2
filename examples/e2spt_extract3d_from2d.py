#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/17/2026 (sludtke@bcm.edu)
# Copyright (c) 2000-2026 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds

from EMAN3 import *
import sys,os

# This program will extract the 3-D particles corresponding to a subset of 2-D particles in SPT
# it will also generate a new 2-D file with renumbered 3-D references
# usage program.py <2d substack lst> <3d stack lst> <new 2d out> <new 3d out>

try: os.unlink(sys.argv[3])
except: pass
try: os.unlink(sys.argv[4])
except: pass

in2d=LSXFile(sys.argv[1])
in3d=LSXFile(sys.argv[2])
out2d=LSXFile(sys.argv[3])
out3d=LSXFile(sys.argv[4])


pt3did={int(x[2]["ptcl3d_id"]) for x in in2d}
pt3dlst=sorted(list(pt3did))
map=dict()

for i,j in enumerate(pt3dlst):
	out3d[i]=in3d[j]
	map[j]=i

for num,fsp,meta in in2d:
	meta["ptcl3d_id"]=map[meta["ptcl3d_id"]]
	out2d[-1]=num,fsp,meta

