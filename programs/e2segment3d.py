#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/02/2010 (ludtke@bcm.edu)
# Copyright (c) 2000-2010 Baylor College of Medicine
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
from optparse import OptionParser
import random
from math import *
import os
import sys
from e2simmx import cmponetomany
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <input volume> [options] 
	This program provides access to various algorithms for segmenting a 3-D volume into multiple pieces automatically.
	Note that you MUST have sufficient RAM to hold at least two copies of the volume in memory. Some segmentation algorithms
	may require more. The actual segmentation is performed using one of the segment.* processors. 'e2help.py processors |grep segment'
	for more information (-v 1 will give even more)."""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_option("--process", metavar="processor_name:param1=value1:param2=value2", type="string",default="segment.kmeans:ampweight=1:nseg=50:thr=0.8",
						help="The name and parameters of a processor to perform the segmentation. 'e2help.py processors -v' for a full list. Default=segment.kmeans:ampweight=1:nseg=50:thr=0.8 ")
	parser.add_option("--output", default=None, type="string",help="Name of output file for segmented volume")
	parser.add_option("--chimeraout", default=None, type="string",help="Name of file to write center of segments in UCSF Chimera marker format.")
	parser.add_option("--pdbout", default=None, type="string",help="Name of file to write center of segments in PDB format.")
	parser.add_option("--txtout", default=None, type="string",help="Name of file to write center of segments in text format (n\tx\ty\tz, with coordinates in pixels, 0,0,0 in the corner)")
	parser.add_option("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()

	
	if options.process[:8]!="segment." :
		print "You must specify a segment.* processor with any necessary parameters in the form segment.xxxx:parm=value:parm=value"
		sys.exit(1)

	E2n=E2init(sys.argv,options.ppid)

	if options.verbose>0: print "Reading volume"
	volume=EMData(args[0],0)
	
	if options.verbose>0: print "Executing segmentation"
	(processorname, param_dict) = parsemodopt(options.process)
	seg=volume.process(processorname,param_dict)
	seg["apix_x"]=volume["apix_x"]
	seg["apix_y"]=volume["apix_y"]
	seg["apix_z"]=volume["apix_z"]
	
	if options.verbose>0: print "Writing output"
	if options.output!=None : seg.write_image(options.output,0)
	
	# make a list of 3-tuples for the center locations
	centers=seg["segment_centers"]
	centers=[(centers[i],centers[i+1],centers[i+2]) for i in xrange(0,len(centers),3)]
	
	# write output
	if options.chimeraout : write_chimera_markers(options.chimeraout,centers,seg["apix_x"],seg["apix_y"],seg["apix_z"],seg["apix_x"]*3)
	if options.pdbout : write_pdb_markers(options.pdbout,centers,seg["apix_x"],seg["apix_y"],seg["apix_z"])
	if options.txtout :
		out=file(options.txtout,"w")
		for n,i in enumerate(centers): out.write("%d\t%1.3f\t%1.3f\t%1.3f\n"%(n,i[0],i[1],i[2]))
		out.close()
	E2end(E2n)
		
def write_chimera_markers(filename,centers,apix_x,apix_y,apix_z,marker_size=3.0) :
	"""Writes a set of coordinates (without connectors) in Chimera marker XML format"""
	try: 
		out=file(filename,"w")
		out.write('<marker_set name="%s">\n'%filename.split(".")[0])
		for j,c in enumerate(centers) :
			out.write('<marker id="%d" x="%1.3f" y="%1.3f" z="%1.3f" radius="%1.3f"/>\n'%(j,c[0]*apix_x,c[1]*apix_y,c[2]*apix_z,marker_size))
		out.write('</marker_set>\n')
		out.close()
	except:
		traceback.print_exc()
		print "\n---------------\nFailed to write Chimera output file, check permissions, etc"
		
	
def write_pdb_markers(filename,centers,apix_x,apix_y,apix_z):
	"""Writes a set of coordinates into a pseudo-PDB file"""
	try: 
		out=file(filename,"w")
		for j,c in enumerate(centers) :
			out.write("ATOM  %5d  CA  ALA  %4d    %8.3f%8.3f%8.3f  1.00%6.2f      S_00  0 \n"%(j,j,c[0]*apix_x,c[1]*apix_y,c[2]*apix_z,1.0))
		out.close()
	except:
		traceback.print_exc()
		print "\n---------------\nFailed to write PDB output file, check permissions, etc"

if __name__ == "__main__":
	main()
