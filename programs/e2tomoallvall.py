#!/usr/bin/env python

#
# Author: David Woolford 04/16/2009 (woolford@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

from optparse import OptionParser
from EMAN2 import file_exists,EMData,E2init,E2progress,E2end,EMANVERSION,check_eman2_type_string,EMUtil,remove_file
import EMAN2
from EMAN2db import EMTask

def check_options(options,args):
	error = []
	if len(args) < 2:
		error.append("Error - to average you must supply atleast two images")
	else:
		for arg in args:
			if not file_exists(arg):
				error.append("%s does not exist" %arg)

	if options.align == None:
		error.append("Error - you have to supply the align option")
	else:
		e = check_eman2_type_string(options.align,EMAN2.Aligners,"Aligners")
		if e != None:
			error.append(e)
			
	if options.ralign != None: # not strictly necessary
		e = check_eman2_type_string(options.ralign,EMAN2.Aligners,"Aligners")
		if e != None:
			error.append(e)
	
	if options.cmp == None:
		error.append("Error - you have to supply the cmp option")
	else:
		e = check_eman2_type_string(options.cmp,EMAN2.Cmps,"Cmps")
		if e != None:
			error.append(e)
	
	if file_exists(options.output):
		if not options.force:
			error.append("Error - output file exists. Remove it or supply the force argument")
		else:
			remove_file(options.output)
	
	return error
	
class EMAllVAll:
	'''
	This class breaks the jobs of the boot-strapped average generation procedure
	so that they can be run in parallel. 
	'''
	def __init__(self,files,options,logger):
		'''
		@param options - the options returned by the call to (options, args) = parser.parse_args() 
		@param args - a list of image names - args is that which is returned by the call to (options, args) = parser.parse_args()
		'''
		self.options = options
		self.files = files
		self.logger = logger
		self.using_cuda = EMUtil.cuda_available() and options.cuda
		
		self.output_images = None # these will be a list of images to write the output to
		self.init_output()
		
	def init_output(self):
		
		self.output_images = []
		e = EMData(len(self.files),len(self.files))
		e.to_zero()
		
		# need to store correlation score (or cmp score)
		self.output_images.append(e)
		
		# need to store translation x
		self.output_images.append(e.copy())
		# need to store translation y
		self.output_images.append(e.copy())
		# need to store translation z
		self.output_images.append(e.copy())
		
		# need to store translation az
		self.output_images.append(e.copy())
		# need to store translation alt
		self.output_images.append(e.copy())
		# need to store translation phi
		self.output_images.append(e.copy())
		
		
	def execute(self):
		alignment_jobs = []
		for i in range(len(self.files)):
			for j in range(i+1,len(self.files)):
				alignment_jobs.append([i,j])
		from e2tomoaverage import EMTomoAlignments
		alignments_manager = EMTomoAlignments(self.options)
		alignments_manager.execute(alignment_jobs, self.files,self)
		
#		options = self.options
#		align_data = EMAN2.parsemodopt(options.align)
#		cmp_data = EMAN2.parsemodopt(options.cmp)
#		ralign_data = None
#		if options.ralign != None: ralign_data = EMAN2.parsemodopt(options.ralign)
#		
#		if self.options.parallel and len(self.options.parallel) > 2 and self.options.parallel[:2] == "dc":
#			task_customers = []
#			tids = []
#			
#			for i in range(len(self.files)):
#				for j in range(i+1,len(self.files)):
#					data = {}
#					data["probe"] = ("cache",self.files[i],0)
#					data["target"] = ("cache",self.files[j],0)
#					data["target_idx"] = j
#					data["probe_idx"] = i
#
#			
#					task = EMTomoAlignTaskDC(data=data,align_data=align_data,cmp_data=cmp_data,ralign_data=ralign_data,using_cuda=self.using_cuda)
#				
#					from EMAN2PAR import EMTaskCustomer
#					etc=EMTaskCustomer(self.options.parallel)
#					#print "Est %d CPUs"%etc.cpu_est()
#					tid=etc.send_task(task)
#					#print "Task submitted tid=",tid
#					
#					task_customers.append(etc)
#					tids.append(tid)
#					
#			self.dc_monitor(task_customers,tids)
#		else:
#			n = len(self.files)
#			n = n*(n-1)/2 # geometris sequence
#			p = 0.0
#			for i in range(len(self.files)):
#				for j in range(i+1,len(self.files)):
#					probe = EMData(self.files[i])
#					target = EMData(self.files[j])
#					data = {}
#					data["target"] = target
#					data["probe"] = probe
#					data["target_idx"] = j
#					data["probe_idx"] = i
#
#			
#					task = EMTomoAlignTask(data=data,align_data=align_data,cmp_data=cmp_data,ralign_data=ralign_data,using_cuda=self.using_cuda)
#				
#					rslts = task.execute(self.progress_callback)
#					self.write_output(rslts)
#					
#					p += 1.0
#					E2progress(self.logger,p/n)
					
		self.finalize_writing()

			
#	def progress_callback(self,percent):
#		'''
#		Need this function in order to use a uniform interface for parallel and non parallel
#		'''
#		pass
	
#	def dc_monitor(self,task_customers,tids):
#		'''
#		This program runs a loop that ends only when all of the argument tasks are complete
#		'''
#		import time
#		n = len(task_customers)
#		while 1:
#			if len(task_customers) == 0: break
#			print len(task_customers),"tomo averaging tasks left in main loop"
#			st_vals = task_customers[0].check_task(tids)
#			for i in xrange(len(task_customers)-1,-1,-1):
#				st = st_vals[i]
#				if st==100:
#					task_customer = task_customers[i]
#					tid = tids[i] 
#					
#					task_customers.pop(i)
#					tids.pop(i)
#
#					rslts = task_customer.get_results(tid)
#					
#					self.write_output(rslts[1])
#
#					if self.logger != None:
#						E2progress(self.logger,1.0-len(task_customers)/(n))
#			
#			time.sleep(5)
			
	def write_output(self,results):
		
		cmp = results["cmp"]
		ali = results["ali"]
		target_idx = results["target_idx"]
		probe_idx = results["probe_idx"]
		
		a = ali.get_params("eman")
		self.output_images[0].set(target_idx,probe_idx,cmp)
		self.output_images[1].set(target_idx,probe_idx,a["tx"])
		self.output_images[2].set(target_idx,probe_idx,a["ty"])
		self.output_images[3].set(target_idx,probe_idx,a["tz"])
		self.output_images[4].set(target_idx,probe_idx,a["az"])
		self.output_images[5].set(target_idx,probe_idx,a["alt"])
		self.output_images[6].set(target_idx,probe_idx,a["phi"])
		
		
		# here we put the inverted alignment as well, yes this is redundant but it 
		# seams like the most logical way to fit with the current eman2 philosophies
		ali_inv = ali.inverse()
		ai = ali_inv.get_params("eman")
		
		self.output_images[0].set(probe_idx,target_idx,cmp)
		self.output_images[1].set(probe_idx,target_idx,ai["tx"])
		self.output_images[2].set(probe_idx,target_idx,ai["ty"])
		self.output_images[3].set(probe_idx,target_idx,ai["tz"])
		self.output_images[4].set(probe_idx,target_idx,ai["az"])
		self.output_images[5].set(probe_idx,target_idx,ai["alt"])
		self.output_images[6].set(probe_idx,target_idx,ai["phi"])
		
	def finalize_writing(self):
		'''
		Writes the output images to disk
		'''
		
		for i,image in enumerate(self.output_images):
			image.set_attr("file_names",self.files)
			image.write_image(self.options.output,i)
		
		


class EMTomoAlignTask:
	'''
	A class the knows how to align two 3D volumes using the current program context
	'''
	def __init__(self,command="e2tomoallvall",data=None,options=None,align_data=["rt.3d.grid",{}],cmp_data=["dot",{}],ralign_data=None,using_cuda=False):
		self.data = data
		self.align_data = align_data
		self.ralign_data = ralign_data
		self.cmp_data = cmp_data
		self.using_cuda = using_cuda
	
	def execute(self,progress_callback):
		'''
		Called to perform class averaging 
		May boot strap the original average, iteratively refines averages, aligns final average to ref 
		'''
		from EMAN2db import db_open_dict
		progress_callback(0)
		
		probe = self.data["probe"]
		target = self.data["target"]
		
		return self.align(probe,target,progress_callback)
		
	def align(self,probe,target,progress_callback):
		
		if self.using_cuda:
			target.set_gpu_rw_current()
			target.cuda_lock() # locking it prevents if from being overwritten
			probe.set_gpu_rw_current()
			probe.cuda_lock()
		
		progress = 0.0
		max_progress = 3
		progress += 1.0
		progress_callback(int(100*(progress/float(max_progress))))

		ali = probe.align(self.align_data[0],target,self.align_data[1])
		
		progress += 1.0
		progress_callback(int(100*(progress/float(max_progress))))
		if self.ralign_data != None:
			xform_align3d = ali.get_attr("xform.align3d")
			self.ralign_data[1]["xform.align3d"] = xform_align3d
			ali = probe.align(self.ralign_data[0],target,self.ralign_data[1])
		
		if self.using_cuda:
			# this is just for consistency. The EMData destructor would have done this automatically anyway
			target.cuda_unlock()
			probe.cuda_unlock()
		
		progress += 1.0
		progress_callback(int(100*(progress/float(max_progress))))
	
		results = {}
		results["cmp"] = ali.cmp(self.cmp_data[0],target,self.cmp_data[1])
		results["ali"] = ali.get_attr("xform.align3d")
		results["target_idx"] = self.data["target_idx"]
		results["probe_idx"] = self.data["probe_idx"]
		
		return results
	  	
			
class EMTomoAlignTaskDC(EMTask):
	'''
	A class the knows how to align two 3D volumes using parallel DC framework
	'''
	def __init__(self,command="e2tomoallvall",data=None,options=None,align_data=["rt.3d.grid",{}],cmp_data=["dot",{}],ralign_data=None,using_cuda=False):
		EMTask.__init__(self,command,data,options)
		# need this import for parallelism to work - somebody fixme?
		from e2tomoallvall import EMTomoAlignTask
		self.align_task = EMTomoAlignTask(data=data,options=options,align_data=align_data,cmp_data=cmp_data,ralign_data=ralign_data,using_cuda=using_cuda)
#		self.align_data = align_data
#		self.ralign_data = ralign_data
#		self.cmp_data = cmp_data
	
	def execute(self,progress_callback):
		'''
		Called to perform class averaging 
		May boot strap the original average, iteratively refines averages, aligns final average to ref 
		'''
		from EMAN2db import db_open_dict
		progress_callback(0)
		
		cache_name=self.data["probe"][1]
		cache=db_open_dict(cache_name)
		probe = cache[self.data["probe"][2]]
		
		cache_name=self.data["target"][1]
		cache=db_open_dict(cache_name)
		target = cache[self.data["target"][2]]
		
		return self.align_task.align(probe, target, progress_callback)
		
	  	
import os,sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image1> <image2> <image3> <image4> ....
	
Boot straps an initial probe doing all versus all alignment of the input images

"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--align",type="string",help="The aligner and its parameters. e.g. --align=rt.3d.grid:ralt=180:dalt=10:dphi=10:rphi=180:search=5", default="rt.3d.grid")
	parser.add_option("--cmp",type="string",help="The comparator used to obtain the final similarity", default="dot")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_option("--output",type="string",default="e2tomoallvall.hdf",help="The output image which will store the results matrix")
	parser.add_option("--force",action="store_true",default=False,help="Force overwrite the output file if it already exists")
	parser.add_option("--parallel",type="string",default=None,help="Use parallelism")
	if EMUtil.cuda_available():
		parser.add_option("--cuda",action="store_true",help="GPU acceleration using CUDA. Experimental", default=False)
 
	(options, args) = parser.parse_args()
	
	error_messages = check_options(options,args)
	if len(error_messages) != 0:
		msg = "\n"
		for error in error_messages:
			msg += error +"\n"
		parser.error(msg)
		exit(1)
	
	logger=E2init(sys.argv)

	module = EMAllVAll(args,options,logger)
	module.execute()
	
	E2end(logger)
	
	
# If executed as a program
if __name__ == '__main__':
	main() 
	