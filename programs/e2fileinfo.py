#!/usr/bin/env python

#
# Author: Steven Ludtke 07/10/2007 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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
#

from EMAN2 import *
import sys
import re
import os

bfactor_expressions = ["bf", "bfactor", "bfactors", "bfac"]
defocus_expressions = ["df", "def", "defocus"]
ac_expressions = ["ac", "ampc", "ampcon", "ampcont", "ampcontrast", "acon", "acont", "acontrast" ]
amp_expressions = ["amp", "A", "amplitude" ]
ee_expressions = ["expenv","ee"]

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <input file> 
	Experimental program, not for general use. See e2iminfo.py and e2bdb.py for general image information queries"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_argument("--auto","-A",type=str,action="append",help="Autobox using specified method: circle, ref, grid",default=[])
#	parser.add_argument("--threshold","-T",type=float,help="Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
#	parser.add_argument("--maxbad","-M",type=int,help="Maximumum number of unassigned helices",default=2)
#	parser.add_argument("--minhelix","-H",type=int,help="Minimum residues in a helix",default=6)
#	parser.add_argument("--apix","-P",type=float,help="A/Pixel",default=1.0)
	
	parser.add_argument("--getinfo",type=str,help="getinfo from file (either defocus, ac (amplitude contrast), or bfactor)",default="")
	parser.add_argument("--remove",type=str,help="getinfo from file (either defocus, ac (amplitude contrast), or bfactor)",default="")
	parser.add_argument("--op",type=str,help="getinfo from file (either defocus, ac (amplitude contrast), or bfactor)",default="")
	parser.add_argument("--outfile",type=str,help="The output file name, may not be required",default="")
	parser.add_argument("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")
	print """WARNING: Experimental program, not for general use. See e2iminfo.py and e2bdb.py for general image information queries"""

	logid=E2init(sys.argv, options.ppid)

	if ( options.remove != "" ):
		checkoutput(options)
		a = parsemodopt_logical( options.remove )
		fileinfo_remove(args[0], a, options.outfile)
		E2end(logid)
		return
	elif ( options.op != "" ):
		checkoutput(options)
		a = parsemodopt_operation( options.op )
		fileinfo_op(args[0], a, options.outfile)
		E2end(logid)
		return
	elif ( options.getinfo != "" ):
		fileinfo_output(args[0],options.getinfo)
	else:
		fileinfo(args[0])
	
	E2end(logid)

def checkoutput(options):
	outfile = options.outfile
	if outfile == '':
		print "Error, must specify the --outfile argument"
		exit(0)
	if os.path.exists(outfile):
		if not options.force:
			print "Error, outfile",outfile, "already exists"
			exit(0)
		else:
			remove_file(options.outfile)

def getidx(paramstring):
	if ( paramstring in defocus_expressions ):
		idx = 0
	elif ( paramstring in bfactor_expressions ):
		idx = 1
	elif ( paramstring in amp_expressions ):
		idx = 2
	elif ( paramstring in ac_expressions ):
		idx = 3
	elif (paramstring in ee_expressions ):
		idx = 4
	else:
		print "error, cannot handle", paramstring
		exit(0)

	return idx
def fileinfo_op(filename,info,outfile):
	
	if ( len(info) != 3 ):
		print "ERROR - logical expression must be a single expression"
		print "Could not process the following: "
		print info
		exit(1)
		
	if ( info[1] not in ["+=", "-=", "*=", "/=", "%="]):
		print "ERROR: could not extract logical expression"
		print "Must be one of", "+=", "-=", "*=", "/=", "%="
		print info
		exit(1)
	
	checkInfoType(info[0])
	
	n=EMUtil.get_image_count(filename)

	idx = getidx(info[0])

	for i in xrange(0,n):
		d=EMData()
		d.read_image(filename,i)
	
		try:
			expr = d.get_attr("IMAGIC.label")
		except RuntimeError:
			print "ERROR: the image has no \"IMAGIC.label\" attribute"
			exit(1)
			
		#print expr
		vals = re.findall("\S*[\w*]", expr)
		
		if ( len(vals) < 4 ):
			print "ERROR: the CTF params were inconsistent with what was expected"
			print "I am examining image number %d, and its ctf params are as follows:" %(i+1)
			print vals
			exit(1)
		
		
		if ( idx == 0 ):
			f = re.findall("\d.*\d*", vals[0])
			score = f[0]
		else:
			score = vals[idx]

		score = float(score)
		
		op_value = float(info[2])
		
		if ( info[1] == "+=" ):
			score += op_value
		elif ( info[1] == "-=" ):
			score -= op_value
		elif ( info[1] == "*=" ):
			score *= op_value
		elif ( info[1] == "/=" ):
			score /= op_value
		elif ( info[1] == "%=" ):
			score %= op_value

		if idx != 0:
			vals[idx] = str(score)
		else:
			vals[idx] = '!--'+str(score)
			
		n = len(vals)
		ilabel = ''
		for j in range(0,n):
			ilabel += vals[j]
			if j != n-1:
				ilabel += ' '
		
		d.set_attr("IMAGIC.label",ilabel)
		d.write_image(outfile,-1)

def fileinfo_remove(filename, info,outfile):
	
	if ( len(info) != 3 ):
		print "ERROR - logical expression must be a single expression"
		print "Could not process the following: "
		print info
		exit(1)
		
	if ( info[1] not in ["==", "<=", ">=", "!=", "~=", "<", ">"] ):
		print "ERROR: could not extract logical expression"
		print "Must be one of \"==\", \"<=\", \">=\", \"<\", \">\" "
		print info
		exit(1)
		
	checkInfoType(infotype)
		
	n=EMUtil.get_image_count(filename)
	t=EMUtil.get_imagetype_name(EMUtil.get_image_type(filename))
	idx = getidx(info[0])
	#os.unlink("cleaned.hed")
	#os.unlink("cleaned.img")

	total_removed = 0

	for i in xrange(0,n):
		d=EMData()
		d.read_image(filename,i,True)
	
		try:
			expr = d.get_attr("IMAGIC.label")
		except RuntimeError:
			print "ERROR: the image has no \"IMAGIC.label\" attribute"
			exit(1)
			
		#print expr
		vals = re.findall("\S*[\w*]", expr)
		
		if ( len(vals) < 4 ):
			print "ERROR: the CTF params were inconsistent with what was expected"
			print "I am examining image number %d, and its ctf params are as follows:" %(i+1)
			print vals
			exit(1)
			
		if ( idx == 0 ):
			f = re.findall("\d.*\d*", vals[0])
			score = f[0]
		else:
			score = vals[idx]

		
		score = float(score)
		comparison_value = float(info[2])
		
		
		write_image = True
		if ( info[1] == "==" ):
			if ( score == comparison_value ):
				write_image = False
		if ( info[1] == "!=" or info[1] == "~="):
			if ( score != comparison_value ):
				write_image = False
		if ( info[1] == ">=" ):
			if ( score >= comparison_value ):
				write_image = False
		if ( info[1] == "<=" ):
			if ( score <= comparison_value ):
				write_image = False
		if ( info[1] == ">" ):
			if ( score > comparison_value ):
				write_image = False
		if ( info[1] == "<" ):
			if ( score < comparison_value ):
				write_image = False
				
		if write_image:
			dd=EMData()
			# now read the image data as well as the header
			dd.read_image(filename,i)
			dd.write_image(outfile, -1 )
		else:
			total_removed += 1
	
	print "Of a total of %d images %d were removed" %(n,total_removed)

def checkInfoType(infotype):
	if ( infotype not in defocus_expressions and infotype not in bfactor_expressions and infotype not in ac_expressions and infotype not in amp_expressions and infotype not in ee_expressions):
		print "Error, infotype %s must be in the following sets:" %infotype
		print bfactor_expressions
		print defocus_expressions
		print ac_expressions
		print amp_expressions
		print ee_expressions
		exit(1)

def fileinfo_output(filename, infotype):
	
	checkInfoType(infotype)
	
	#l=[len(i) for i in filenames]
	#l=max(l)
	
	n=EMUtil.get_image_count(filename)
	t=EMUtil.get_imagetype_name(EMUtil.get_image_type(filename))
	
	idx = getidx(infotype)
	for i in xrange(0,n):
		d=EMData()
		d.read_image(filename,i,True)
	
		try:
			expr = d.get_attr("IMAGIC.label")
		except RuntimeError:
			print "ERROR: the image has no \"IMAGIC.label\" attribute"
			exit(1)
			
					#print expr
		vals = re.findall("\S*[\w*]", expr)
		
		if idx == 0:
			f = re.findall("\d.*\d*", vals[0])
			defocus = f[0]
			print "%f" %float(defocus)
		else:
			print "%f" %float(vals[idx])


def fileinfo(filenames):
	if isinstance(filenames,str) : filenames=[filenames]
	
	l=[len(i) for i in filenames]
	l=max(l)
	
	for i in filenames:
		n=EMUtil.get_image_count(i)
		t=EMUtil.get_imagetype_name(EMUtil.get_image_type(i))
		d=EMData()
		d.read_image(i,0,True)
		if d.get_zsize()==1:
			s="%%-%ds%%s\t%%d\t%%d x %%d"%(l+2)
			print s%(i,t,n,d.get_xsize(),d.get_ysize())
		else:
			s="%%-%ds%%s\t%%d\t%%d x %%d x %%d"%(l+2)
			print s%(i,t,n,d.get_xsize(),d.get_ysize(),d.get_zsize())
		
# If executed as a program
if __name__ == '__main__':
	main()	
