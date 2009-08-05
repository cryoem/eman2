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
from EMAN2 import file_exists,EMData,E2init,E2progress,E2end,EMANVERSION,check_eman2_type_string,numbered_bdb,Transform,EMUtil
import EMAN2
from EMAN2db import EMTask,db_open_dict
from e2tomoallvall import EMTomoAlignTask,EMTomoAlignTaskDC,check_tomo_options

tomo_ave_path_root = "tomo_ave" # this  string is used for making output directories automatically

def check_options(options,args):
	'''
	A way to check the options, arg as returned by parser.parse_args()  in e2tomoaverage
	'''
	error = []
	
	error.extend(check_tomo_options(options)) # there is a big bunch of generic options
	
	if len(args) < 2:
		error.append("Error - to average you must supply atleast two images")
	else:
		for arg in args:
			if not file_exists(arg):
				error.append("%s does not exist" %arg)
	
	if options.path != None:
		if not os.path.exists(options.path):
			error.append( "Error: the path %s does not exist" %options.path)
	else:
		options.path = EMAN2.numbered_path(tomo_ave_path_root,True)
			
	return error
	
class EMBootStrappedAverages:
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
		self.images = None
		self.using_cuda = EMUtil.cuda_available() and options.cuda
		self.current_files = None
		self.jobs_completed = 0
		self.total_jobs = 0 # you can never know how many jobs will need to be completed
		self.current_progress = 0.0
		
	def get_all_v_all_cmd(self):
		'''
		A function for get the allvall command
		'''
		options = self.options
		cmd = "e2tomoallvall.py"
		for file in self.files: cmd += " " + file
		cmd += " --align=" + options.align
		cmd += " --aligncmp=" + options.aligncmp
		cmd += " --cmp=" + options.cmp
		if options.ralign != None:
			cmd += " --ralign=" + options.ralign
			cmd += " --raligncmp=" + options.raligncmp

		cmd += " --nsoln=" + str(options.nsoln)
		
		if options.parallel:
			cmd += " --parallel="+options.parallel
		return cmd
	
	def get_all_v_all_output(self):
		'''
		A function for getting the name of an output file
		'''
		return numbered_bdb("bdb:"+self.options.path+"#all_v_all")
		#return numbered_bdb("bdb:tomo_ave#all_v_all")
	
	def get_couples(self,cmp_data):
		'''
		A function for finding the couples in a matrix of similarity data
		Only considers the top half of the triangle - i.e. assumes similarity matrix is from an all
		versus all program and therefore the diagonal is redundant etc.
		The approach taken here is from Mike Schmid and is best explained inductively - first find the best couple and record
		it (there is always at least one couple). Put the indices of this couple in list of
		'taken' indices, the importance of which will become imminently obvious - Now find the next best couple - one of two scenarios occurs:
		1. Both of the indices are not present in the 'taken' list - Make them a couple and add both indices to the 'taken' list
		2. One of the indices is in the 'taken indices' list, in which case add the index which was not already in the 'taken' list 
		into it.
		This process finds all nicely isolated pairs and is guaranteed to find at least one pair.
		@returns a list of index pairs which are the located couples
		'''
		couples = []
		taken = []
		
		cmp_data_copy = cmp_data.copy()
		cmp_max = cmp_data_copy["maximum"] + 0.1 # just worst then the min
		
		# set the redundant part of the matrix
		# to a small value so it can't be identified
		# as a maximum / can be ignored automatically
		for i in range(cmp_data.get_ysize()):
			for j in range(i+1):
				cmp_data_copy.set(j,i,cmp_max)
		
		# the best match - we always get atleast one
		best = cmp_data_copy.calc_min_location() # remember it's a min because EMAN2 is designed for the result of all comparison to be better if they're smaller
		couples.append([best[0],best[1]])
		taken.extend(best[:2])
		cmp_data_copy.set(best[0],best[1],cmp_max)
		
		# this loop sorts the couples out from the singles
		n = cmp_data.get_xsize()
		while True:
			best = cmp_data_copy.calc_min_location() # remember it's a min because EMAN2 is designed for the result of all comparison to be better if they're smaller
			
			i = best[0]
			j = best[1]
			yes = True
			# if any of the indices have already been encountered then they can not form a couple
			try:
				(val for val in taken if val == i ).next()
				yes = False
			except: pass
			try:
				(val for val in taken if val == j ).next()
				yes = False
			except: pass
			
			if yes:
				couples.append([best[0],best[1]])
				taken.extend(best[:2])
			else:
				# it can't be made into a couple so add any indices into the taken list that are not already there
				for idx in [i,j]:
					try:  (val for val in taken if idx == val ).next()
					except: taken.append(idx)
			
			cmp_data_copy.set(best[0],best[1],cmp_max)
			
			if len(taken) == n: break
			
		return couples
	
	def save_to_workflow_db(self,output_name):
		'''
		Used by the workflow to automatically add the names of averages to a list
		that is displayed in a form
		'''
		if self.options.dbls:
			pdb = db_open_dict("bdb:project")
			tmp_data = pdb.get(self.options.dbls, dfl={})
			s = {}
			s["Original Data"] = output_name
			tmp_data[output_name]= s
			# global.spr_ref_free_class_aves
			pdb[self.options.dbls] = tmp_data
	
	def execute(self):
		'''
		The main function - executes the job of performing all v all boot strapped probe generation
		'''
		if self.logger: E2progress(self.logger,0.0)
		all_v_all_cmd = self.get_all_v_all_cmd()
		all_v_all_output = self.get_all_v_all_output()
		
		# NOTE: calling the allvall program is probably not strictly necessary, seeing
		# as there is a generic framework for generating and executing alignment jobs
		# implemented below that would be easily adaptable to this - however I left it
		# because doing it this way is absolutely equivalent and has the same cost. 
		all_v_all_cmd += " --output="+all_v_all_output
		print "executing",all_v_all_cmd
		if self.logger:	E2progress(self.logger,0.01)
		if ( os.system(all_v_all_cmd) != 0 ):
			print "Failed to execute %s" %all_v_all_cmd
			sys.exit(1)
		if self.logger:	E2progress(self.logger,0.02)
		
		images = []
		images.append(EMData(all_v_all_output,0))
		images.append(EMData(all_v_all_output,1))
		images.append(EMData(all_v_all_output,2))
		images.append(EMData(all_v_all_output,3))
		images.append(EMData(all_v_all_output,4))
		images.append(EMData(all_v_all_output,5))
		images.append(EMData(all_v_all_output,6))
		
		# keep tracks of the names of the new files
		big_n = images[0].get_xsize()*(images[0].get_xsize()-1)/2.0
		
		iter = 1
		current_files = self.files
		# this loop 
		while True:
			couples = self.get_couples(images[0])
			taken = range(images[0].get_xsize())
			#print len(couples),len(taken)
			new_files = []

			# write the averages of the couples to disk, store the new names
			for i,j in couples:
				image_1 = EMData(current_files[j],0)
				image_2 = EMData(current_files[i],0)
				
				d = {}
				d["type"] = "eman"
				d["tx"] = images[1].get(i,j)
				d["ty"] = images[2].get(i,j)
				d["tz"] = images[3].get(i,j)
				d["az"] = images[4].get(i,j)
				d["alt"] = images[5].get(i,j)
				d["phi"] = images[6].get(i,j)
				t = Transform(d)
				image_1.process_inplace("math.transform",{"transform":t})
				image_2 += image_1
				image_2.mult(.5)
				image_2.set_attr("src_image",current_files[j]) # so we can recollect how it was created
				image_2.set_attr("added_src_image",current_files[i]) # so we can recollect how it was created
				image_2.set_attr("added_src_transform",t) # so we can recollect how it was created
				image_2.set_attr("added_src_cmp",images[0](i,j)) # so we can recollect how it was created
				output_name = numbered_bdb("bdb:"+self.options.path+"#tomo_ave_0"+str(iter))
				image_2.write_image(output_name,0)
				if self.options.dbls: self.save_to_workflow_db(output_name)
				new_files.append(output_name)
				taken.remove(i)
				taken.remove(j)
			
			num_new = len(new_files) # the number of averages produced
			new_n = len(new_files) + len(taken)
			new_images = []
			e = EMData(new_n,new_n)
			e.to_zero()
			new_images.append(e)
			for j in range(6): new_images.append(e.copy())
			
			for i,idxi in enumerate(taken):
				new_files.append(current_files[idxi])
				for j,idxj in enumerate(taken):
					if i == j: continue
					else:
						new_images[0].set(num_new+i,num_new+j,images[0].get(idxi,idxj))
						new_images[1].set(num_new+i,num_new+j,images[1].get(idxi,idxj))
						new_images[2].set(num_new+i,num_new+j,images[2].get(idxi,idxj))
						new_images[3].set(num_new+i,num_new+j,images[3].get(idxi,idxj))
						new_images[4].set(num_new+i,num_new+j,images[4].get(idxi,idxj))
						new_images[5].set(num_new+i,num_new+j,images[5].get(idxi,idxj))
						new_images[6].set(num_new+i,num_new+j,images[6].get(idxi,idxj))
			
			alignment_jobs = []# a list of comparisons to be performed
			for i in range(num_new):
				for j in range(i+1,len(new_files)):
					alignment_jobs.append([i,j])
					
			if self.logger: 
				E2progress(self.logger,1.0-len(alignment_jobs)/big_n)
					
			self.register_current_images(new_images)
			self.register_current_files(new_files)
			alignments_manager = EMTomoAlignments(self.options)
			alignments_manager.execute(alignment_jobs, new_files,self)
			
			if len(couples) == 1 and len(taken) == 2: break
			
			self.write_current_images(new_files)
			current_files = new_files
			images = new_images
			iter += 1
			
		if self.logger: E2progress(self.logger,1.0)
					
			
	
	def register_current_images(self,images): 
		'''
		This is just a formalization of doing self.images = images
		'''
		self.images = images
	def register_current_files(self,files):
		'''
		This is just a formalization of doing self.current_files = files
		'''
		self.current_files = files
	
	def write_current_images(self,files):
		'''
		Writes the image results to disk
		'''
		output = self.get_all_v_all_output()
		for i,image in enumerate(self.images):
			image.set_attr("file_names",files)
			image.write_image(output,i)
	
	def process_output(self,results):
		'''
		Stores results (alignment, scores etc) in  memory
		'''
		from e2tomoallvall import EMTomoOutputWriter
		output_writer = EMTomoOutputWriter()
		output_writer.write_ouptut(results, self.images, self.current_files)
		
		
class EMTomoAlignments:
	'''
	A class for performing many alignments, takes care of parallel considerations automatically
	This class is used extensively, in e2tomoallvall, e2tomoaverage, and e2tomohunter
	'''
	def __init__(self,options):
		self.options = options
		self.using_cuda = EMUtil.cuda_available() and options.cuda
		self.nsoln = options.nsoln
		
	def execute(self,alignment_jobs,files,caller):
		'''
		The main function
		@param alignment_jobs a list of alignment pair indices like this [[0,1],[2,1],[2,3],[0,5],...] etc the indices pair represent images to be aligned and correspond to the order of the files argument
		@param files a list of filenames - used to read image based on the indices present in alignment_jobs
		@param caller - the calling object - it needs to have a function called process_output that takes a dictionary as the argument 
		'''
		options = self.options
		align_data = EMAN2.parsemodopt(options.align)
		align_cmp_data = EMAN2.parsemodopt(options.aligncmp)
		cmp_data = EMAN2.parsemodopt(options.cmp)
		ralign_data = None
		if options.ralign != None: 
			ralign_data = EMAN2.parsemodopt(options.ralign)
			ralign_cmp_data = EMAN2.parsemodopt(options.raligncmp)
			
		
		data = {}
		data["align"] = align_data
		data["aligncmp"] = align_cmp_data
		data["cmp"] = cmp_data
		if ralign_data:
			data["ralign"] = ralign_data
			data["raligncmp"] = ralign_cmp_data
			
		data["using_cuda"] = self.using_cuda
		data["nsoln"] = self.nsoln
			
		if self.options.parallel and len(self.options.parallel) > 2 and self.options.parallel[:2] == "dc":
			task_customers = []
			tids = []

			for i,j in alignment_jobs:
				data["probe"] = ("cache",files[i],0)
				data["target"] = ("cache",files[j],0)
				data["target_idx"] = j
				data["probe_idx"] = i

				task = EMTomoAlignTaskDC(data=data)
				
				from EMAN2PAR import EMTaskCustomer
				etc=EMTaskCustomer(self.options.parallel)
				#print "Est %d CPUs"%etc.cpu_est()
				tid=etc.send_task(task)
				#print "Task submitted tid=",tid
				
				task_customers.append(etc)
				tids.append(tid)
			
			self.dc_monitor(task_customers,tids,caller)
		else:
			n = len(alignment_jobs)
			p = 0.0
			for i,j in alignment_jobs:
				probe = EMData(files[i],0)
				target = EMData(files[j],0)
				data["target"] = target
				data["probe"] = probe
				data["target_idx"] = j
				data["probe_idx"] = i

		
				task = EMTomoAlignTask(data=data)
				rslts = task.execute(self.progress_callback)
				caller.process_output(rslts)
				
				p += 1.0
				
	def progress_callback(self,val):
		'''
		The function needs to be supplied in order to make the task execution interface consistent
		DC tasks are given a progress_callback function that does something - tasks that 
		are executed in the current thread call this function and currently it does nothing
		'''
		pass
			
	def dc_monitor(self,task_customers,tids,caller):
		'''
		This program runs a loop that ends only when all of the argument tasks are complete
		'''
		import time
		n = len(task_customers)
		while 1:
			if len(task_customers) == 0: break
			print len(task_customers),"tomo averaging tasks left in main loop"
			st_vals = task_customers[0].check_task(tids)
			for i in xrange(len(task_customers)-1,-1,-1):
				st = st_vals[i]
				if st==100:
					task_customer = task_customers[i]
					tid = tids[i] 
					
					task_customers.pop(i)
					tids.pop(i)

					rslts = task_customer.get_results(tid)
					
					caller.process_output(rslts[1])
			
			time.sleep(5)			


import os,sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image1> <image2> <image3> <image4> ....
	
Currently only supports bootstrapping an initial probe doing all versus all alignment of the input images

"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--align",type="string",help="The aligner and its parameters. e.g. --align=rt.3d.grid:ralt=180:dalt=10:dphi=10:rphi=180:search=5", default="rt.3d.grid")
	parser.add_option("--aligncmp",type="string",help="The comparator used for determing the best initial alignment", default="dot.tomo:threshold=0")
	parser.add_option("--cmp",type="string",help="The comparator used to obtain the final similarity", default="dot.tomo:threshold=0")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparator used for determing the refined alignment", default="dot.tomo:threshold=0")
	parser.add_option("--bootstrap",action="store_true",default=True,help="Boot strap alignment")
	parser.add_option("--output",type="string",default="e2tomoave.hdf",help="The output image which will store the results matrix")
	parser.add_option("--parallel",type="string",default=None,help="Use parallelism")
	parser.add_option("--path", default=None, type="string",help="The name of a directory where results are placed. If unspecified will generate one automatically of type tomo_ave_??.")
	parser.add_option("--nsoln", default=1, type="int",help="If supplied and greater than 1, the nsoln-best alignments will be written to a text file. This is useful for debug but may be left unspecified")
	parser.add_option("--dbls",type="string",help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
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
	
	if options.bootstrap:
		module = EMBootStrappedAverages(args,options,logger)
		module.execute()
	else:
		print "boot strap only supported technique"
	E2end(logger)
	
	
# If executed as a program
if __name__ == '__main__':
	main() 
	