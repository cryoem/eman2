#!/usr/bin/env python
from __future__ import division

#
# Author: Michael Bell, 6/19/2018 (jmbell@bcm.edu)
# Copyright (c) 2018-2022 Baylor College of Medicine
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

from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """
	prog <tomograms> [options]
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tomograms",help="List the tomograms you want to process.", default="", guitype='filebox', browser="EMTomoTable(withmodal=True,multiselect=True)",  filecheck=False, row=0, col=0,rowspan=1, colspan=3, mode='proc')
	parser.add_argument("--invert",action="store_true",help="Invert the contrast of the tomograms in output files (default false). If specified _proctag tomograms will be generated.",default=False, guitype='boolbox', row=1, col=1, rowspan=1, colspan=1, mode='proc[False]')
	parser.add_argument("--proctag",help="Tag added to the name of each tomogram when using the proc options. Default is 'preproc'.",default="preproc", guitype='strbox', row=1, col=2, rowspan=1, colspan=1, mode='proc["seg"]')
	parser.add_argument("--proc1",help="If specified _proctag files will be generated. Typical = filter.lowpass.gauss:cutoff_abs=.25",default="filter.lowpass.gauss:cutoff_abs=.25", guitype='strbox', row=2, col=0, rowspan=1, colspan=3, mode='proc["filter.lowpass.gauss:cutoff_abs=.25"]')
	parser.add_argument("--proc2",help="If specified _proctag tomograms will be generated. Typical = filter.highpass.gauss:cutoff_pixels=5",default="filter.highpass.gauss:cutoff_pixels=5", guitype='strbox', row=3, col=0, rowspan=1, colspan=3, mode='proc["filter.highpass.gauss:cutoff_pixels=5"]')
	parser.add_argument("--proc3",help="If specified _proctag tomograms will be generated.",default="normalize", guitype='strbox', row=4, col=0, rowspan=1, colspan=3, mode='proc[""]')
	parser.add_argument("--proc4",help="If specified _proctag tomograms will be generated. Typical = threshold.clampminmax.nsigma:nsigma=3",default="threshold.clampminmax.nsigma:nsigma=3", guitype='strbox', row=5, col=0, rowspan=1, colspan=3, mode='proc["threshold.clampminmax.nsigma:nsigma=3"]')
	parser.add_argument("--proc5",help="If specified _proctag tomograms will be generated.",default="", guitype='strbox', row=14, col=0, rowspan=1, colspan=3, mode='proc[""]')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	(options, args) = parser.parse_args()

	proc1 = parsemodopt(options.proc1)
	proc2 = parsemodopt(options.proc2)
	proc3 = parsemodopt(options.proc3)
	proc4 = parsemodopt(options.proc4)
	proc5 = parsemodopt(options.proc5)

	if not options.invert and proc1[0] == None and proc2[0] == None and proc3[0] == None and proc4[0] == None and proc5[0] == None:
		print("No processing specified. Exiting.")
		sys.exit(1)

	for i,filename in enumerate(args):
		name = os.path.basename(filename).split(".")[0]
		procout=["tomograms/{}_{}.hdf".format(name,options.proctag)]
		
		sys.stdout.write("\r{}/{}".format(i+1,len(args)))
		sys.stdout.flush()
		
		if proc1[0] != None: procout.append(proc1)
		if proc2[0] != None: procout.append(proc2)
		if proc3[0] != None: procout.append(proc3)
		if proc4[0] != None: procout.append(proc4)
		if proc5[0] != None: procout.append(proc5)

		im = EMData(filename,0)
		if options.invert: im.mult(-1.0)

		if len(procout) > 1:
			for op in procout[1:]: # we take a sequence of processor option 2-tuples
				if op[0] != None: 
					im.process_inplace(op[0],op[1])
		
		im.write_image(procout[0],0)

if __name__ == "__main__":
	main()
