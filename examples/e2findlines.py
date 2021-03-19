#!/usr/bin/env python

#
# Author: Steve Ludtke, 2/09/20 
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

from future import standard_library
standard_library.install_aliases()
from math import *
import os
import sys
from EMAN2 import *
import queue
import numpy as np

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """e2findlines sets/img.lst
	
	** EXPERIMENTAL **
	this program looks for ~ straight line segments in images, such as wrinkles in graphene oxide films or possible C-film edges

	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--threshold", type=float, help="Threshold for separating particles, default=3", default=3.0)
	parser.add_argument("--newsets",default=False,action="store_true",help="Split lines/nolines into 2 new sets")
	#parser.add_argument("--output",type=str,help="Output filename (text file)", default="ptclplot.txt")
	parser.add_argument("--gui",default=False, action="store_true",help="show histogram of values")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on the local computer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--invar",default=False, action="store_true",help="create the invar file for the newsets. The newsets option must be used.")

	(options, args) = parser.parse_args()
	
	if (len(args)<1 ): parser.error("Please specify an input stack/set to operate on")
	
	E2n=E2init(sys.argv, options.ppid)
	
	options.threads+=1		# one extra thread for storing results

	im0=EMData(args[0],0)	# first image
	r2=im0["ny"]/4	# outer radius

	# we build up a list of 'Z scores' which should be larger for images containing one or more parallel lines.
	# if 2 lines aren't parallel the number may be lower, even if the lines are strong, but should still be higher
	# than images without lines in most cases
	n=EMUtil.get_image_count(args[0])
	step=max(n//500,1)
	Z=[]
	im2d=[]
	for i in range(n):
		im=EMData(args[0],i)
		a=im.do_fft().calc_az_dist(60,-88.5,3,4,r2)
		d=np.array(a)
		Z.append((d.max()-d.mean())/d.std())
		if i%step==0: 
			im["zscore"]=(d.max()-d.mean())/d.std()
			im2d.append(im)

	if options.gui:
		# GUI display of a histogram of the Z scores
		from eman2_gui.emhist import EMHistogramWidget
		from eman2_gui.emimagemx import EMImageMXWidget
		from eman2_gui.emapplication import EMApp
		app = EMApp()
		histw=EMHistogramWidget(application=app)
		histw.set_data(Z)
		app.show_specific(histw)
		imd=EMImageMXWidget(application=app)
		im2d.sort(key=lambda x:x["zscore"])
		imd.set_data(im2d)
		app.show_specific(imd)
		app.exec_()

"""
		if options.newsets:
			lstin=LSXFile(args[0])

			# output containing images with lines
			linesfsp=args[0].rsplit(".",1)[0]+"_lines.lst"
			try: os.unlink(linesfsp)
			except: pass
			lstlines=LSXFile(linesfsp)	

			# output containin images without lines
			nolinesfsp=args[0].rsplit(".",1)[0]+"_nolines.lst"
			try: os.unlink(nolinesfsp)
			except: pass
			lstnolines=LSXFile(nolinesfsp)	

			for i,z in enumerate(Z):
				if z>options.threshold: lstlines[-1]=lstin[i]
				else: lstnolines[-1]=lstin[i]
"""

	if options.newsets and not options.invar:
		lstin=LSXFile(args[0])

		# output containing images with lines
		linesfsp=args[0].split("__",1)[0]+"_lines__"+args[0].split("__",1)[1]
		try: os.unlink(linesfsp)
		except: pass
		lstlines=LSXFile(linesfsp)	

		# output containin images without lines
		nolinesfsp=args[0].split("__",1)[0]+"_nolines__"+args[0].split("__",1)[1]
		try: os.unlink(nolinesfsp)
		except: pass
		lstnolines=LSXFile(nolinesfsp)	

		for i,z in enumerate(Z):
			if z>options.threshold: lstlines[-1]=lstin[i]
			else: lstnolines[-1]=lstin[i]

	if options.newsets and options.invar:
		lstin=LSXFile(args[0])

		# output containing images with lines
		fnamemod=input("Type the filename modifier:   ")
		linesfsp=args[0].split("__",1)[0]+f"_lines_{fnamemod}__"+args[0].split("__",1)[1]
		try: os.unlink(linesfsp)
		except: pass
		lstlines=LSXFile(linesfsp)	

		# output containing images without lines
		nolinesfsp=args[0].split("__",1)[0]+f"_nolines_{fnamemod}__"+args[0].split("__",1)[1]
		try: os.unlink(nolinesfsp)
		except: pass
		lstnolines=LSXFile(nolinesfsp)	
		
		# output copy of nolines folder
		invarfsp = nolinesfsp.rsplit("_",1)[0] + "_invar.lst"
		try: os.unlink(invarfsp)
		except: pass
		lstinvar=LSXFile(invarfsp)	
		

		for i,z in enumerate(Z):
			if z>options.threshold: lstlines[-1]=lstin[i]
			else: 
				lstnolines[-1]=lstin[i]
				lstinvar[-1]=lstin[i]
        
		"""invarfsp = nolinesfsp.rsplit("_",1)[0] + "_invar.lst"
		try: os.unlink(nolinesfsp)
		except: pass
		lstinvar=LSXFile(invarfsp)	
		print(f"len(lstinvar)={lstinvar.__len__()}")
		
		for i in lstnolines: lstinvar[-1]=lstnolines[i]
		print(f"len(lstinvar)={len(lstinvar)}. len(lstnolines)={len(lstnolines)}")
		print(lstinvar[-1], lstnolines[-1])

		#print(f"copying {nolinesfsp} to {invarfsp}.")
		#os.system(f"cp {nolinesfsp} {lstinvar}")"""
		print(f"running: e2proclst.py {invarfsp} --retype ctf_flip_invar" )
		os.system(f"e2proclst.py {invarfsp} --retype ctf_flip_invar" )
		print("invar file created.")
	
		
	E2end(E2n)
	



if __name__ == "__main__":
    main()
