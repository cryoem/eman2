#!/usr/bin/env python

#
# Author: John Flanagan May 13th 2011  
# ported and refactored using EMAN2 by David Woolford October 6th 2008
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

def main():
	usage = """%prog inputmap to unwrap
	Unwrap a 3D map into cyliderical coordinates. Useful for inspecting the interior of radially symmetric molecules
	"""
	
	map3d = EMData()
	map3d.read_image(sys.argv[1])
	
	unwraped3d = None
	for z in xrange(map3d.get_zsize()):
		t = Transform({"type":"eman","tz":(map3d.get_zsize()/2 - z)})
		slice2d = EMData(map3d.get_xsize(), map3d.get_ysize())
		slice2d.cut_slice(map3d, t, 0)
		polarslice = slice2d.unwrap(0,-1,180,0,0,1,0)
		
		if unwraped3d == None:
			unwraped3d = EMData(polarslice.get_xsize(), polarslice.get_ysize(), map3d.get_zsize())
		polarslice.uncut_slice(unwraped3d,t)
	
	unwraped3d.process_inplace("xform.flip",{"axis":"z"}) # A cheat
	unwraped3d.write_image(sys.argv[2])
	
if __name__ == "__main__":
	main()

