#!/usr/bin/env python

# extractfsc.py
# Author: Steven Ludtke, 04/06/2011 (sludtke@bcm.edu)
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

import sys
from EMAN2 import *
from EMAN2db import db_open_dict

if len(sys.argv)==1 :
	print "Please provide the name of a refine_xx directory to extract FSC curves from"
	exit(1)

db=db_open_dict("bdb:%s#convergence.results"%sys.argv[1],ro=True)

for k in db.keys():
	curve=db[k]
	out=file("fsc_%s_%s.txt"%(sys.argv[1].rsplit("_",1)[-1],k),"w")
	for i in xrange(len(curve[0])): out.write("%1.5f\t%1.4f\n"%(curve[0][i],curve[1][i]))
