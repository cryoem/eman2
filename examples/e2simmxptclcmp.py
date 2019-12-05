#!/usr/bin/env python
from __future__ import division

#
# Author: Steven Ludtke, 03/16/2012 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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
from builtins import range
import os
from EMAN2 import *
from numpy import *
from eman2_gui.emapplication import EMApp

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2simmxptclcmp.py [options] <simmx file> ...
Plots the set of similarity quality values for a single particle from a set of simmx files.        
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--ptcl",type=int,help="particle number",default=0)
#        parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
        #parser.add_argument("--output",type=str,help="Output text file",default="zvssim.txt")
        #parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
        #parser.add_argument("--inimgs",type=str,help="Input image file",default=None)
        #parser.add_argument("--outimgs",type=str,help="Output image file",default="imgs.hdf")
        #parser.add_argument("--filtimgs",type=str,help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
        #parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
        
	(options, args) = parser.parse_args()

	app = EMApp()

	from eman2_gui.emplot2d import EMPlot2DWidget
	plotw=EMPlot2DWidget(application=app)
	
	plts=[]
	for simmx in args:
		mx=EMData(simmx,0)
		plts.append(mx.get_clip(Region(0,options.ptcl,mx["nx"],1)))
		plts[-1].process_inplace("normalize")
		plotw.set_data(plts[-1],simmx)

	plotw.setWindowTitle("ptcl vs similarity")
	plotw.show()
	try: plotw.raise_()
	except: pass
	app.exec_()

	#out=file("simmxplot.txt","w")
	#for i in range(plts[0]["nx"]):
		#out.write(str(i))
		#for p in plts: out.write("\t{}".format(p[i]))
		#out.write("\n")
	#out.close()
		
	#os.system("e2display.py --plot simmxplot.txt")

if __name__ == "__main__":  main()

