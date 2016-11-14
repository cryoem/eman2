#!/usr/bin/env python

#
# Author: Steven Ludtke, 12/20/2011 (sludtke@bcm.edu)
# Copyright (c) 2011- Baylor College of Medicine
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

# This will read a series of images, translationally align them, average them
# together, and optionally iterate. Translational alignment only.
# transalignavg2.py <infile> ...

from EMAN2 import *
import sys
from math import *


ref0=EMData(sys.argv[1])
ref0.process_inplace("normalize.edgemean")
ref0.translate(5,0,0)							# small shift to avoid nozero issue
nx,ny=ref0["nx"],ref0["ny"]

for i in xrange(2):
	print "Iter ",i
	avgr=Averagers.get("mean", {"ignore0":True})
	ref0.write_image("seq.hdf",-1)
	for fsp in sys.argv[1:]:
		im=EMData(fsp,0)
		im.process_inplace("normalize.edgemean")
		if im["nx"]!=nx or im["ny"]!=ny :
			im=im.get_clip(Region(-(nx-im["nx"])/2,-(ny-im["ny"])/2,nx,ny))
	
		im.write_image("seq.hdf",-1)
		ima=im.align("translational",ref0,{"nozero":1,"maxshift":ref0["nx"]/4.0},"ccc",{})
		ima.write_image("seq.hdf",-1)
		print fsp,ima["xform.align2d"],ima.cmp("ccc",ref0)
		ima.process_inplace("normalize.toimage",{"to":ref0,"ignore_zero":1})
		avgr.add_image(ima)
	

	ref0=avgr.finish()

ref0.write_image("result.hdf",0)
display(ref0)
