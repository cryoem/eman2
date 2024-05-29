#!/usr/bin/env python
#
# Author: Steven Ludtke  4/6/22 
# Copyright (c) 2022- Baylor College of Medicine
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
# Imports tilt series movies collected in one canoncial form into tilt series stacks

from EMAN2 import *
import os

def main():
	
	usage="""Imports a set of MRC movie stacks comprising a tilt series following a canonical naming scheme where
the filename contains the tilt angle in square brackets, eg-
Position_1_001[30.00]-847511.mrc
Position_1_001[27.00]-847411.mrc
...

identified files in the specified target path will be averaged, stacked, and copied to tiltseries/
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path", type=str,help="specify path to seach for appropriate numbered files", default=".")
	parser.add_argument("--apix", type=float,help="override A/pix in header", default=-1.0)
	parser.add_argument("--skipframes", type=int,help="Number of frames to skip at the beginning of each movie",default=0)
	parser.add_argument("--compressbits", type=int,help="Number of bits to preserve in the output tilt series",default=6)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higher number means higher level of verboseness",default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	logid=E2init(sys.argv,options.ppid)

	try: os.mkdir("tiltseries")
	except: pass

	# files of interest
	files=[f for f in os.listdir(options.path) if "[" in f and "]" in f]

	# base names, each is a tilt series
	bases={f.split("[")[0].rsplit("_",1)[0] for f in files}
	if len(bases)==0 : 
		print("No appropriate files identified in target folder. Aborting.")
		sys.exit(1)
	
	for base in bases:
		# files associated with this base name
		bfiles=[(float(f.split("[")[1].split("]")[0]),int(f.split("[")[0].rsplit("_",1)[1]),f) for f in files if f.startswith(base)]
		bfiles.sort()
		try: os.unlink(f"tiltseries/{base}.hdf")
		except: pass
	
		print(f"{base} with {len(bfiles)} tilt movies")
		for ang,seq,f in bfiles:
			if options.verbose>1: print("\t",ang)
			try: h=EMData(f"{options.path}/{f}",0,True)
			except:
				print("Error reading ",f)
				sys.exit(1)
			if options.apix>0: apix=options.apix
			else: apix=h["apix_x"]

			if h["nz"]>1:
				# movie stored as a volume
				mov=[EMData(f"{options.path}/{f}",0,False,Region(0,0,i,h["nx"],h["ny"],1)) for i in range(options.skipframes,h["nz"])]
			else:
				# movies in a stack
				mov=EMData.read_images(f"{options.path}/{f}")[options.skipframes:]

			img=sum(mov)
			img["apix_x"]=apix
			img["apix_y"]=apix
			img["apix_z"]=apix
			img["tilt_angle"]=ang
			img["tilt_seq"]=seq
			
			img.write_compressed(f"tiltseries/{base}.hdf",-1,options.compressbits,nooutliers=True)
	
if __name__ == "__main__":
	main()
