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
from EMAN2 import file_exists,EMData,E2init,E2progress,E2end,EMANVERSION,check_eman2_type_string,numbered_bdb
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
		
	def get_all_v_all_cmd(self):
		options = self.options
		cmd = "e2tomoallvall.py"
		for file in self.files: cmd += " " + file
		cmd += " --align=" + options.align
		cmd += " --cmp=" + options.cmp
		if options.ralign != None:
			cmd += " --ralign=" + options.ralign

		
		
		if options.parallel:
			cmd += " --parallel="+options.parallel
		return cmd
	
	def get_all_v_all_output(self):
		return numbered_bdb("bdb:tomo_bootstrap#all_v_all")
	def execute(self):
		
		all_v_all_cmd = self.get_all_v_all_cmd()
		all_v_all_output = self.get_all_v_all_output()
		
		all_v_all_cmd += " --output="+all_v_all_output
		print "executing",all_v_all_cmd
		if ( os.system(all_v_all_cmd) != 0 ):
			print "Failed to execute %s" %all_v_all_cmd
			sys.exit(1)
			
		print "file exists?",file_exists(all_v_all_output)
		
		cmp_data = EMData(all_v_all_output,0)
		tx_data = EMData(all_v_all_output,1)
		ty_data = EMData(all_v_all_output,2)
		tz_data = EMData(all_v_all_output,3)
		az_data = EMData(all_v_all_output,4)
		alt_data = EMData(all_v_all_output,5)
		phi_data = EMData(all_v_all_output,6)
		
		couples = []
		taken = []
		
		
		cmp_data_copy = cmp_data.copy()
		cmp_min = cmp_data_copy["minimum"] - 0.1 # just worst then the min
		
		# make the redundant part of the matrix all bad
		for i in range(cmp_data.get_ysize()):
			for j in range(i+1):
				cmp_data_copy.set(j,i,cmp_min)
		
		best = cmp_data_copy.calc_max_location()
		couples.append([best[0],best[1]])
		taken.extend(best[:2])
		cmp_data_copy.set(best[0],best[1],cmp_min)
		
		n = cmp_data.get_xsize()
		
		while True:
			best = cmp_data_copy.calc_max_location()
			
			
			i = best[0]
			j = best[1]
			yes = True
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
				for val in [i,j]:
					try: (val for val in taken if val == val ).next()
					except: taken.extend(val)
			
			cmp_data_copy.set(best[0],best[1],cmp_min)
			print couples, taken
			
			if len(taken) == n: break
		
		print couples
		print taken
		
		
		
		

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
	parser.add_option("--boostrap",action="store_true",default=True,help="Boot strap alignment")
	parser.add_option("--output",type="string",default="e2tomoave.hdf",help="The output image which will store the results matrix")
	parser.add_option("--parallel",type="string",default=None,help="Use parallelism")
	
	(options, args) = parser.parse_args()
	
	error_messages = check_options(options,args)
	if len(error_messages) != 0:
		msg = "\n"
		for error in error_messages:
			msg += error +"\n"
		parser.error(msg)
		exit(1)
	
	logger=E2init(sys.argv)
	
	if options.boostrap:
		module = EMBootStrappedAverages(args,options,logger)
		module.execute()
	else:
		print "boot strap only supported technique"
	E2end(logger)
	
	
# If executed as a program
if __name__ == '__main__':
	main() 
	