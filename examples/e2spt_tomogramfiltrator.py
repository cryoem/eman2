#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya 05/2010. Last modification: 05/2013
# Copyright (c) 2011 Baylor College of Medicine
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
import os
import sys


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2spt_tomogramfiltrator.py
	This program pads a tomogram to filter it and then crops it back to the original size to avoid aliasing effects that occur when filtering straight with e2proc3d.py
	"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--output", type=str, help="""File name for the filtered tomogram.""", default=None)

	parser.add_argument("--tomogram", type=str, help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None, guitype='filebox', browser='EMSubTomosTable(withmodal=True,multiselect=False)', row=0, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")	

	(options, args) = parser.parse_args()

	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
	
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)

	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)

	a=EMData()
	a.read_image(options.tomogram) 
	x = a["nx"]
	y = a["ny"]
	z = a["nz"]
	print "The dimensions of the tomogram are", x,y,z

	pady = int(y*0.30)
	padx = int(x*0.30)
	padz = int(z*0.30)

	print "I will pad this much in x,y,z", padx,pady,padz
	final_y = y + 2*pady
	final_x = x + 2*padx
	final_z = z + 2*padz

	print "So the final size in x,y,z will b", final_x,final_y,final_z

	r1 = Region(-padx,-pady,-padz, x + padx, y + pady, z + padz)
	a.clip_inplace(r1)
	test1 = a.get_value_at(0,0,y+pady)
	test2 = a.get_value_at(0,0,pady -1)
	print "test values are", test1, test2
	x2 = a['nx']
	y2 = a['ny']
	z2 = a['nz']
	print "The new dimensions are ", x2,y2,z2
	#factor = 1.0/float(res)
	#b=a.process("filter.lowpass.gauss",{"cutoff_freq":factor})

	'''
	#Preprocess, lowpass and/or highpass
	'''
	if options.preprocess:
		a.process_inplace(options.preprocess[0],options.preprocess[1])
	
	if options.lowpass:
		a.process_inplace(options.lowpass[0],options.lowpass[1])
	
	if options.highpass:
		a.process_inplace(options.highpass[0],options.highpass[1])

	r2 = Region(padx,pady,padz,x,y,z)
	a.clip_inplace(r2)

	tomo = options.tomogram
	outname = tomo.replace('.','_LP.')
	if options.output:
		outname = options.output
	#else:
	#	if '.mrc' in outname:
	#		outname = outname.replace('.mrc','_LP.mrc')
	#	if '.hdf' in outname:
	#		outname = outname.replace('.hdf','_LP.hdf')
	#	if '.rec' in outname:
	#		outname = outname.replace(".rec","_LP.rec")

	#filtered_tomo = filtered_tomo.replace(".","_new.")
	x3 = a['nx']
	y3 = a['ny']
	z3 = a['nz']
	print "The filtered tomogram has been clipped presumably to the original size", x3,y3,z3
	print "The name of the filtered tomogram is", outname
	a.write_image(outname)

	return()

#tomogram = argv[1]
#resolution = float(argv[2])

#tomogram_filtrator(tomogram, resolution)


	
if __name__ == "__main__":
    main()

