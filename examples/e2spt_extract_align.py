#!/usr/bin/env python
#
# Author: Jesus Galaz-Montoya, March 2014; last update by Jesus Galaz-Montoya on March/20/2014
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

from EMAN2 import *
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This program takes a .json file produced by e2spt_classaverage.py and a stack of raw,
	unaligned 3-D images (or "subtomograms") and applies the transforms indicated in the .json file
	to average the particles in the stack.
	Alternatively, if the aligned stack is supplied and there's alignment information in
	the headers of the particles, a .json file with the alignment parameters can be produced.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", default='',type=str, help="The name of the hdf stack of volumes to process.")
	parser.add_argument("--output", default="avg.hdf",type=str, help="The name of the output average volume.")
	parser.add_argument("--rotationtype", default="eman",type=str, help="Valid options are: eman,imagic,mrc,spider,quaternion,sgirot,spi,xyz")
	
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean")
	
	parser.add_argument("--sym", type=str, default='', help = "Symmetry to impose - choices are: c<n>, d<n>, h<n>, tet, oct, icos")

	parser.add_argument("--path",default='',type=str,help="Name of directory where to save the output file.")
	parser.add_argument("--alifile",default='',type=str,help=".json file with alingment parameters, if raw stack supplied via --input.")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--saveali",action="store_true", default=False,help="""If set, will save the 
		aligned particle volumes in class_ptcl.hdf. Overwrites existing file.""")

	(options, args) = parser.parse_args()


	if options.averager: 
		options.averager=parsemodopt(options.averager)

	if not options.input:
		parser.print_help()
		sys.exit(0)

	if ".hdf" not in options.output and ".mrc" not in options.output:
		print "ERROR. The output must contain a valid format ending, for example '.hdf.' TERMINATING!"
		sys.exit()
	
	if ".hdf" not in options.input:
		print "ERROR. HDF is the only supported input format."
		sys.exit()
	
	logid=E2init(sys.argv,options.ppid)
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptextractali')
	
	n = EMUtil.get_image_count( os.getcwd() + '/' + options.input )
	
	if not options.alifile:
		a=Transform({"type":"eman","alt":1.0})
		#k=list(a.get_rotation(sys.argv[2]).keys())
	
		k=list(a.get_rotation( options.rotationtype ).keys())

		k.remove("type")
		if len(k)==3: 
			print "#{},{},{}".format(*k)
		else: 
			print "#{},{},{},{}".format(k)

		for i in range(n):
			im=EMData( options.input ,i)
			#xf=im["spt_ali_param"]
			
			xf=im['xform.align3d']
			r=xf.get_rotation( options.rotationtype )
			print "{}".format(i),
			for j in k: 
				print ", {}".format(r[j]),
			print ""
	
	elif options.alifile:
		preOrientationsDict = js_open_dict(options.alifile)
			
		avgr=Averagers.get(options.averager[0], options.averager[1])
	
		for i in range(n):			
			a=EMData( options.input, i)
			
			ptcl = a.copy()
			ID=unicode(i)
			t = preOrientationsDict[0][ID]
			
			#t=''
			#for key in ts.keys():
			#	print "Key is", key,type(key)
			#	if ID in key:
			#		t = ts[unicode(0)][unicode(key]
					
			
			print "Transform is",t
			
			
			
			
			ptcl['origin_x'] = 0
			ptcl['origin_y'] = 0
			ptcl['origin_z'] = 0
			ptcl['xform.align3d'] = Transform()
			
			if t:
				ptcl.transform(t)
			
			avgr.add_image(ptcl)
			
			if options.saveali:
				print "\nAligned particle should be saved",i
				pass
				
		avg=avgr.finish()
		if options.sym and options.sym is not 'c1' and options.sym is not 'C1':
			avg=avg.process('xform.applysym',{'sym':options.sym})
		
		avg.write_image( options.path + '/' + options.output, 0 )
		
		preOrientationsDict.close()			
	
	E2end(logid)
	
	return


if __name__ == '__main__':
	main()