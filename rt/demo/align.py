#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

# usage: align.py referencefile [inputfile outputfile]

from EMAN2 import *
import sys

if (len(sys.argv) != 1 and len(sys.argv) != 4):
	print "usage: align.py [referencefile inputfile outputfile]"
	sys.exit(1)

if len(sys.argv) == 1:
	reffile = TestUtil.get_debug_image("samesize1.mrc")
	inputfile = TestUtil.get_debug_image("samesize2.mrc")
	outputfile = "align.mrc"
	
else:
	reffile = sys.argv[1]
	inputfile = sys.argv[2]
	outputfile = sys.argv[3]
	
ref = EMData()
ref.read_image(reffile)

e = EMData()
e.read_image(inputfile)

e.align("rtf_best", {"to": ref, "maxshift": 8})

e.write_image(outputfile)
