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
import random
from math import *
import os
import sys
from e2simmx import cmponetomany
import traceback

def read_helix(filename,sx,sy,sz,ax,ay,az):
	print "Reading helix atoms from pdb file..."
	points = []
	pdbfile = open(filename, "r")
	lines = pdbfile.readlines()
	pdbfile.close()
	nhlx=0
	# Read atoms
	for line in (i for i in lines if i.startswith("HELIX  ")):
		atomid=[int(line[21:27].strip()), int(line[33:38].strip())]
		for atline in (i for i in lines if i.startswith("ATOM  ")):
			cn=int(atline[22:30].strip())
			if cn>=atomid[0] and cn<=atomid[1]:
				if atline[13:15]=="CA":
					pos = (float(atline[30:38].strip())/ax+sx/2, float(atline[38:46].strip())/ay+sy/2, float(atline[46:54].strip())/az+sz/2)#,nhlx)
					points.append( pos)
		nhlx+=1
	
	return points

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <input volume> [options] 
	This program provides access to various algorithms for segmenting a 3-D volume into multiple pieces automatically.
	Note that you MUST have sufficient RAM to hold at least two copies of the volume in memory. Some segmentation algorithms
	may require more. The actual segmentation is performed using one of the segment.* processors. 'e2help.py processors |grep segment'
	for more information (-v 1 will give even more)."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--process", metavar="processor_name:param1=value1:param2=value2", type=str,default="segment.kmeans:ampweight=1:nseg=50:thr=0.8",
						help="The name and parameters of a processor to perform the segmentation. 'e2help.py processor segment -v 1' for a full list. Default=segment.kmeans:ampweight=1:nseg=50:thr=0.8 ")
	parser.add_argument("--output", default=None, type=str,help="Name of output file for segmented volume")
	parser.add_argument("--chimeraout", default=None, type=str,help="Name of file to write center of segments in UCSF Chimera marker format.")
	parser.add_argument("--pdbout", default=None, type=str,help="Name of file to write center of segments in PDB format.")
	parser.add_argument("--txtout", default=None, type=str,help="Name of file to write center of segments in text format (n\tx\ty\tz, with coordinates in pixels, 0,0,0 in the corner)")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--shifttocenter", action="store_true", help="Shift the output pdb to center of the density map")
	parser.add_argument("--helixfile", default=None, type=str, help="Start with existing secondary structure.")
	parser.add_argument("--edgefile", default=None, type=str, help="Write an edge file for pathwalker.py. Only avaliable while using existing secondary structures.")
	
	
	(options, args) = parser.parse_args()

	
	if options.process[:8]!="segment." :
		print "You must specify a segment.* processor with any necessary parameters in the form segment.xxxx:parm=value:parm=value"
		sys.exit(1)

	E2n=E2init(sys.argv,options.ppid)

	if options.verbose>0: print "Reading volume"
	volume=EMData(args[0],0)
	
	if options.shifttocenter:
		sx=volume.get_xsize()
		sy=volume.get_ysize()
		sz=volume.get_zsize()
	else:
		sx=0
		sy=0
		sz=0
	print sx,sy,sz
	
	if options.helixfile!=None:
		helix=read_helix(options.helixfile,sx,sy,sz,volume["apix_x"],volume["apix_y"],volume["apix_z"])
		ps=options.process
		p1=ps.find("nseg")
		p2=ps.find(":",p1)
		if p2==-1: 
			p2=len(ps)
		num=int(ps[p1+5:p2])
		nps=ps[:p1+5]+str(num-len(helix))+ps[p2:]
		options.process=nps
		
		print "Existing helix length: "+ str(len(helix))
		#print helix[0]
		
	
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
	
	if options.helixfile!=None:
		for h in helix:
			centers.append(h)
		out=file(options.edgefile,"w")
		for i in range(1,len(helix)):
			d=sqrt((helix[i-1][0]-helix[i][0])**2+(helix[i-1][1]-helix[i][1])**2+(helix[i-1][2]-helix[i][2])**2)
			if d<3.8:
				out.write("%d\t%d\n"%(i-1+num-len(helix),i+num-len(helix)))
		
		out.close()
		#print centers
		
		
	# write output
	if options.chimeraout : write_chimera_markers(options.chimeraout,centers,seg["apix_x"],seg["apix_y"],seg["apix_z"],seg["apix_x"]*3)
	if options.pdbout : write_pdb_markers(options.pdbout,centers,seg["apix_x"],seg["apix_y"],seg["apix_z"],sx,sy,sz)
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
		
	
def write_pdb_markers(filename,centers,apix_x,apix_y,apix_z,sx,sy,sz):
	"""Writes a set of coordinates into a pseudo-PDB file"""
	print sx,sy,sz
	try: 
		out=file(filename,"w")
		for j,c in enumerate(centers) :
			out.write("ATOM  %5d  CA  ALA  %4d    %8.3f%8.3f%8.3f  1.00%6.2f      S_00  0 \n"%(j,j,(c[0]-sx/2)*apix_x,(c[1]-sy/2)*apix_y,(c[2]-sz/2)*apix_z,1.0))
		out.close()
	except:
		traceback.print_exc()
		print "\n---------------\nFailed to write PDB output file, check permissions, etc"

if __name__ == "__main__":
	main()
