#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/11/15 
# Copyright (c) 2015- Baylor College of Medicine
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
from copy import deepcopy
import os
import sys
import random
from random import choice
import traceback

from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]

	The goal of this program is to reduce the heterogeneity of a reconstruction by splitting a single map
	into two maps, each more homogeneous. You must run e2refine_easy to completion before using this program.
	It will take the class-averaging results from the final iteration, and split the particles from each 
	class-average into 2 groups, producing 2 class-averages for each. The program then attempts to construct
	a maximally self-consistent grouping of these pairs of class averages into 2 3-D maps. 
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path", default=None, type=str,help="The name of an existing refine_xx folder, where e2refine_easy ran to completion")
	parser.add_argument("--parallel", default="thread:2", help="Standard parallelism option. Default=thread:2")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path==None:
		paths=[i for i in os.listdir(".") if "refine_" in i and len(i)==9]
		paths.sort()
		options.path=paths[-1]

	# check the specified path for the files we need
	try:
		olddb = js_open_dict(options.path+"/0_refine_parms.json")
		last_map=olddb["last_map"]
		last_iter=int(last_map.split("_")[-1][:2])
		try: ptcls=olddb["inputavg"]
		except: ptcls=olddb["input"]
		
		if options.verbose : print "Found iteration {} in {}, using {}".format(last_iter,options.path," & ".join(ptcls))
	except:
		traceback.print_exc()
		print "Error: Cannot find necessary files in ",options.path
		sys.exit(1)
		
	logger=E2init(sys.argv,options.ppid)

	# classmx is a list with 2 elements. Each element is a list of EMData from the corresponding cls_result file
	classmx=[]
	classmx.append(EMData.read_images("{}/cls_result_{:02d}_even.hdf".format(options.path,last_iter)))
	classmx.append(EMData.read_images("{}/cls_result_{:02d}_odd.hdf".format(options.path,last_iter)))
	ncls=max(int(classmx[0][0]["maximum"])+1,int(classmx[0][0]["maximum"])+1)

	# Rearrange the info in classmx
	classlists=[[] for i in xrange(ncls)]	# empty list for each class
	
	# This will produce a list of particles with Transforms for each class
	for eo in (0,1):
		for y in xrange(classmx[eo][0]["ny"]):
			ptcl=[eo,y,Transform({"type":"2d","tx":classmx[eo][2][0,y],"ty":classmx[eo][3][0,y],"alpha":classmx[eo][4][0,y],"mirror":int(classmx[eo][5][0,y])})]
			classlists[int(classmx[eo][0][0,y])].append(ptcl)
	

	# Initialize parallelism
	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel)

	# prepare tasks
	tasks=[]
	gc=0
	ns=[classmx[eo][0]["ny"] for eo in (0,1)]
	for c,cl in enumerate(classlists):
		if len(cl)<6 : continue
		tasks.append(ClassSplitTask(ptcls,ns,cl,gc,options.verbose-1))
		gc+=1

	# execute task list
	taskids=etc.send_tasks(tasks)
	alltaskids=taskids[:]

	while len(taskids)>0 :
		curstat=etc.check_task(taskids)
		for i,j in enumerate(curstat):
			if j==100 :
				rslt=etc.get_results(taskids[i])
				rsltd=rslt[1]
				cls=rslt[0].options["classnum"]
				
#					print rsltd["avg"],cls
				rsltd["avg"].write_image("test.hdf",cls)
				for ii,i in enumerate(rsltd["basis"]): i.write_image("basis.hdf",cls*5+ii)
				
		taskids=[j for i,j in enumerate(taskids) if curstat[i]!=100]

		if options.verbose and 100 in curstat :
			print "%d/%d tasks remain"%(len(taskids),len(alltaskids))
		if 100 in curstat :
			E2progress(logger,1.0-(float(len(taskids))/len(alltaskids)))

		time.sleep(3)


		if options.verbose : print "Completed all tasks"


	print "Class averaging complete"
	E2end(logger)

class ClassSplitTask(JSTask):
	"""This task will create a single task-average"""

	def __init__(self,ptclfiles,ns,ptcls,nc,verbose):
		"""ptclfiles is a list of 2 (even/odd) particle stacks. ns is the number of particles in each of ptcfiles. ptcls is a list of lists containing [eo,ptcl#,Transform]"""
		data={"particles1":["cache",ptclfiles[0],(0,ns[0])],"particles2":["cache",ptclfiles[1],(0,ns[1])]}
		JSTask.__init__(self,"ClassSplit",data,{},"")

		self.options={"particles":ptcls,"classnum":nc,"verbose":verbose}

	def execute(self,callback=None):
		"""This does the actual class-averaging, and returns the result"""
		options=self.options

		if options["verbose"]>0 : print "Start averaging class {} with {} particles ".format(options["classnum"],len(options["particles"]))
		
		files=self.data["particles1"][1],self.data["particles2"][1]
		
		ptcls=options["particles"]			# just a shortcut

		#if len(options["particles"])<5 :
			#z=EMData(str(files[ptcls[0][0]]),ptcls[0][1]).to_zero()
			#return {"avg":z,"basis":[z,z,z,z,z]}
		 		
#		print files,ptcls[0]
		# read in all particles and append each to element in ptcls
		avgr=Averagers.get("mean")
		for p in ptcls: 
			p.append(EMData(str(files[p[0]]),p[1]).process("xform",{"transform":p[2]}).process("normalize.circlemean",{"radius":-6}).process("mask.soft",{"outer_radius":-8,"width":4}))
			avgr.add_image(p[3])
		
		avg=avgr.finish()
		if options["verbose"]>0: print "averaging class {}".format(options["classnum"])

		pca=Analyzers.get("pca_large",{"nvec":5})
		for p in ptcls: 
			pca.insert_image(p[3]-avg)
		
		basis=pca.analyze()

		
		#for p in ptcls: 
			#avgr.add_image(p[3].process("xform",{"transform":p[2]}))
		
		if options["verbose"]>0: print "Finish averaging class {}".format(options["classnum"])
		if callback!=None : callback(100)
		return {"avg":avg,"basis":basis}

jsonclasses["ClassSplitTask"]=ClassSplitTask.from_jsondict


if __name__ == "__main__":
    main()
