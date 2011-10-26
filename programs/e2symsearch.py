#!/usr/bin/env python

#
# Author: John Flanagan Oct 20th 2011 (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine
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
import math
import os
from EMAN2db import EMTask

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program aligns a paricle to its symmetry axis. There are tow algorithmic modes. A course search followed by simplex
	minizization OR monte carlo course search followed by simplex minizization. The Goal is to align the paricle to its 
	symmetry axis so symmetry can be applied for avergaing and for alignment spped up(it is only necessary to search over the
	assymetric unit!
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_header(name="refineheader", help='Options below this label are specific to e2refine', title="### e2refine options ###", row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of input volume", guitype='filebox', row=0, col=0, rowspan=1, colspan=2)
	parser.add_argument("--output", dest="output", default=None,type=str, help="The name of the output volume", guitype='filebox', filecheck=False, row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--sym", dest = "sym", default="c1", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction omit this option or specify c1.", guitype='symbox', row=3, col=0, rowspan=1, colspan=2)
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. For speed purposes. Default=0, no shrinking", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--steps", dest="steps", type = int, default=10, help="Number of steps (for the MC)", guitype='intbox', row=4, col=1, rowspan=1, colspan=1)
	parser.add_argument("--symmetrize", default=True, action="store_true", help="Symmetrize volume after alignment.", guitype='boolbox', row=5, col=0, rowspan=1, colspan=1)
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used in comparing the symmtrized object to unsymmetrized", default="ccc", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\', True)', row=6, col=0, rowspan=1, colspan=2)
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value>",default=None, guitype='strbox', row=7, col=0, rowspan=1, colspan=2)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	
	# Initialize parallelism if being used
	from EMAN2PAR import EMTaskCustomer
	if options.parallel :
		etc=EMTaskCustomer(options.parallel)
	else:
		etc=EMTaskCustomer("thread:1")
		
	logid=E2init(sys.argv,options.ppid)
	
	volume = EMData(options.input)
	if options.shrink:
		volume.process_inplace("math.meanshrink",{"n":options.shrink})
		
	symalgorithm = SymALignStrategy(volume, options.sym, options.steps, options.cmp, etc)
	symxform = symalgorithm.execute()
	
	print "Found Best alignment, now writting output..."
	if options.shrink:
		trans = symxform.get_trans()
		symxform.set_trans(trans[0]*options.shrink, trans[1]*options.shrink, trans[2]*options.shrink)
		output = EMData(options.input)
		output.process_inplace('xform',{'transform':symxform})
		if options.symmetrize:
			output = output.process('xform.applysym',{'sym':options.sym})
	else:
		output = volume.process('xform',{'transform':symxform})
		if options.symmetrize:
			output = output.process('xform.applysym',{'sym':options.sym})
	
	output.write_image(options.output)
	
	E2end(logid)

# Use strategy pattern here. Any new stategy needs to inherit this
class Strategy:
	def __init__(self, volume, sym, steps, comp, etc):
		self.volume = volume
		self.sym = sym
		self.steps = steps
		self.cmp=parsemodopt(comp)
		self.etc = etc
		
	def execute(self):
		raise NotImplementedError("Subclass must implement abstract method")

class SymALignStrategy(Strategy):
	def __init__(self, volume, sym, steps, comp, etc):
		Strategy.__init__(self, volume, sym, steps, comp, etc)
		
	def execute(self):
		Util.set_randnum_seed(Util.get_randnum_seed())
		tasks=[]
		for i in xrange(self.steps):
			az = Util.get_frand(0,360) 					# theta
			alt  = math.degrees(math.acos(2*Util.get_frand(0,1) - 1))	# phi
			phi = Util.get_frand(0,360)					# kappa
			t = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})
			# Now do the simplex alignment
			tasks.append(SymAlignTask(self.volume, self.sym, self.cmp, t))
		tids=self.etc.send_tasks(tasks)
		
		solns = []
		# Hang till tasks are done
		while 1:
			time.sleep(5)
			proglist=self.etc.check_task(tids)
			
			for i,prog in enumerate(proglist):
				if prog==100:
					print "Finished a MC trial"
					r=self.etc.get_results(tids[i])
					solns.append(r[1]["symalign"])
					
			tids=[j for i,j in enumerate(tids) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
			if len(tids)==0: break
		
		print "\nFinished alignments...\n"
		# Search for the best scoring
		bestscore = 0
		bestxform = Transform()
		for i,soln in enumerate(solns):
			print "score for MC trial %d is %f"%(i,soln.get_attr('score'))
			if soln.get_attr('score') < bestscore:
				bestscore = soln.get_attr('score')
				bestxform = soln.get_attr('xform.align3d')
				
		return bestxform
		
		
class SymAlignTask(EMTask):
	def __init__(self, volume, sym, comp, xform):
		data = {"volume":volume}
		EMTask.__init__(self,"CmpTilt",data,{},"")
		
		self.sym = sym
		self.cmp=comp
		self.xform=xform
		
	def execute(self,callback=None):
		symalign = self.data['volume'].align('symalignquat',self.data['volume'],{"sym":self.sym,"xform.align3d":self.xform},self.cmp[0],self.cmp[1])
		return {"symalign":symalign}
if __name__ == "__main__":
    main()