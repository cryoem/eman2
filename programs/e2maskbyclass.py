#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# Author: Steven Ludtke, 09/04/2017 (sludtke@bcm.edu)
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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import time
import os
import threading
import queue
from sys import argv,exit

def maskfile(jsd,n,fsp,classes,masks,clsmap,options):

	fspout=fsp.rsplit(".",1)[0].split("__")[0]+"__ctf_flip_masked.hdf"
	fspbout=fsp.rsplit(".",1)[0].split("__")[0]+"__ctf_flip_bispec.hdf"

	for i in range(len(clsmap)):
		ptcl=EMData(fsp,i)
		# if the particle isn't in any classes we put the unmasked image in the output file
		if clsmap[i]!=-1 :
			# align the projection to the particle
			alip=classes[clsmap[i]].align("rotate_translate_tree",ptcl)

			# mask in the same orientation
			alim=masks[clsmap[i]].process("xform",{"transform":alip["xform.align2d"]})
			ptcl.mult(alim)

		ptcl.write_image(fspout,i)

		if options.redobispec:
			bspec=ptcl.process("math.bispectrum.slice",{"fp":bispec_invar_parm[1],"size":bispec_invar_parm[0]})
			bspec.write_image(fspbout,i)

	# signal that we're done
	jsd.put(n)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2maskbyclass.py [options]
This program will generate a masked stack of particles for each micrograph based on the class-average each particle belongs to. While this will run with any EMAN2.2 generated class-average stack, it is designed to work on averages produced from __ctf_flip_fullres particles.

In short:
- read a class-average stack
- convert each class-average into a 2-D mask
- find the original phase-flipped particle for each particle in the class-average
- mask each particle on a per-micrograph basis

once complete, bispectra can be recomputed based on the masked particles, or the masked particles can be used directly, potentially producing better alignments and averages.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--classes",type=str,default=None,help="Path to a class-average file (must be EMAN2 HDF averages)")
	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--nofullresok",action="store_true",help="Overrides the requirement that the class-averages be made from _fullres particle images.",default=False)
	parser.add_argument("--redobispec",action="store_true",help="Recomputes bispectra from masked particles",default=False)
	parser.add_argument("--gui",action="store_true",help="Permits interactive adjustment of mask parameters",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="tuning[True]")
# 	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
# 	parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
# 	parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
# 	parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls.hdf) containing the aligned subtomograms.",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	NTHREADS=max(options.threads+1,2)		# we have one controlling thread

	logid=E2init(sys.argv, options.ppid)

	classes=EMData.read_images(options.classes)		# read all class-averages
	nx=classes[0]["nx"]

	if options.gui:
		parm=maskparmgui(classes)
	else: parm=(0.03,8,nx//8,0.3)

	# Make a mask for each class-average
	for c in classes: c.process_inplace("filter.lowpass.gauss",{"cutoff_freq":parm[0]})
	masks=[c.process("mask.auto2d",{"nmaxseed":parm[1],"nshells":parm[2],"radius":old_div(nx,10),"return_mask":1,"sigma":parm[3]}) for c in classes]
	for i in masks: i.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.03})

	# Find all of the particles
	ptcls={}
	lsx=LSXFile(classes[0]["class_ptcl_src"])
	for cn,cl in enumerate(classes):
		ptclis=cl["class_ptcl_idxs"]
		try: ptclis.extend(cl["exc_class_ptcl_idxs"])
		except: pass
		for p in ptclis:
			orig=lsx[p]							# nextfile,extfile,comment
			try: ptcls[orig[1]][orig[0]]=cn		# dictionary of filenames, value is a list of class numbers for each of the n particles in the file
			except:
				n=EMUtil.get_image_count(orig[1])
				ptcls[orig[1]]=[-1]*n			# initialize with -1 for all class numbers
				ptcls[orig[1]][orig[0]]=cn

	if not options.nofullresok :
		for fsp in ptcls:
			if "ctf_flip_fullres" not in fsp:
				print("ERROR: This program is meant to be used with full resolution phase flipped particles (__ctf_flip_fullres). ",fsp," does not appear to be this type of file. You can override this behavior with --nofullresok")

	if options.verbose:
		print(sum([i.count(-1) for i in list(ptcls.values())])," missing particles in classes. They will be unmasked")

# 	import pprint
# 	pprint.pprint(ptcls)

	jsd=queue.Queue(0)

	n=-1
	thrds=[(jsd,i,k,classes,masks,ptcls[k],options) for i,k in enumerate(ptcls)]

	# here we run the threads and save the results, no actual alignment done here
	print(len(thrds)," threads")
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch]=threading.Thread(target=maskfile,args=thrds[thrtolaunch])
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)

		# no return other than the thread that finished
		while not jsd.empty():
			n=jsd.get()
			thrds[n].join()
			thrds[n]=None

	if options.verbose: print("Finished processing ",len(ptcls), "particle files")
	E2end(logid)

def maskparmgui(classes):
	try:
		from eman2_gui.emapplication import EMApp
		from PyQt4 import QtCore, QtGui, QtOpenGL
		from OpenGL import GL,GLUT
		from eman2_gui.valslider import ValSlider,CheckBox
		from eman2_gui.emimagemx import EMImageMXWidget
		
	except:
		print("Error: PyQt4 must be usable to use the --gui option")
		sys.exit(1)


	class GUImask(QtGui.QWidget):
		def __init__(self,app,classes):
			"""Effectively a modal dialog for selecting masking parameters interactively
			"""
			self.app=app
			QtGui.QWidget.__init__(self,None)
			nx=classes[0]["nx"]
			
			self.classes=classes
			self.classview=EMImageMXWidget(self,classes)
			
			self.vbl = QtGui.QVBoxLayout(self)
			self.vbl.addWidget(self.classview)
			
			self.hbl = QtGui.QHBoxLayout()
			
			self.cmode=CheckBox(self,"orig",value=1)
			self.hbl.addWidget(self.cmode)

			self.slpres=ValSlider(self,(0.001,0.2),"Low-pass Filter:",0.03,90)
			self.hbl.addWidget(self.slpres)
			
			self.snmax=ValSlider(self,(0,20),"NMax:",5,90)
			self.snmax.intonly=1
			self.hbl.addWidget(self.snmax)
			
			self.sshells=ValSlider(self,(0,40),"NShells:",nx//8,90)
			self.sshells.intonly=1
			self.hbl.addWidget(self.sshells)
			
			self.ssigma=ValSlider(self,(0,2),"Sigma:",0.333,90)
			self.hbl.addWidget(self.ssigma)
			
			self.bok=QtGui.QPushButton("OK")
			self.hbl.addWidget(self.bok)

			self.vbl.addLayout(self.hbl)

			self.cmode.valueChanged.connect(self.newParm)
			self.slpres.valueChanged.connect(self.newParm)
			self.snmax.valueChanged.connect(self.newParm)
			self.sshells.valueChanged.connect(self.newParm)
			self.ssigma.valueChanged.connect(self.newParm)
			self.bok.clicked[bool].connect(self.close)
	
			self.newParm()

		def quit(self):
			self.app.close_specific(self)
			
		def newParm(self):
			if self.cmode.getValue(): 
				self.classview.set_data(self.classes)
				return
			self.masked=[i.process("filter.lowpass.gauss",{"cutoff_freq":self.slpres.value}) for i in self.classes]
			nx=self.masked[0]["nx"]
			for i,im in enumerate(self.masked):
				im.process_inplace("mask.auto2d",{"nmaxseed":int(self.snmax.value),"nshells":int(self.sshells.value),"radius":old_div(nx,10),"return_mask":1,"sigma":self.ssigma.value})
				im.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.03})
				im.mult(self.classes[i])
				
			self.classview.set_data(self.masked)
	app=EMApp()
	gui=GUImask(app,classes)
	gui.show()
	gui.raise_()
	app.exec_()
	
	return((gui.slpres.value,gui.snmax.value,gui.sshells.value,gui.ssigma.value))

if __name__ == "__main__":
	main()

