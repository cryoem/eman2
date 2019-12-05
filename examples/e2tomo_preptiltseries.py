#!/usr/bin/env python
'''
====================
Author: Jesus Galaz - 10/August/2018, Last update: 10/August/2018
====================
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
'''
from past.utils import old_div
from EMAN2 import *
from EMAN2_utils import *
import os

def main():

	usage = """Program to automatically prefilter tiltseries to facilitate their alignment. 
			Should ideally be run after xray correction in IMOD."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--input",type=float,default=None,help="""Default=None. Sampling size (angstroms per pixel) to set on the header of the cylindrical mask.""")
	
	parser.add_argument("--highpass",type=str,default='filter.highpass.gauss:cutoff_pixels=4',help="Default=filter.highpass.gauss:cutoff_pixels=4 (zeros out the first four Fourier pixels)_. A highpass filtering processor (as in e2proc3d.py) to be applied to the tiltseries.")

	parser.add_argument("--lowpass",type=str,default='filter.lowpass.tanh:cutoff_freq=0.02',help="Default=filter.lowpass.tanh:cutoff_freq=0.02 (filters to 50 angstroms resolution). A lowpass filtering processor (as in e2proc3d.py) to be applied to the tiltseries.")

	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, used for cross platform PPID""",default=-1)

	parser.add_argument("--replace",action='store_true',default=False,help="""Default=False. This will save the filtered tiltseries to the filename of the original tiltseries. The original tiltseries will be backed up.""")
	
	parser.add_argument("--verbose", "-v", help="""verbose level [0-9], higner number means higher level of verboseness. Default=0.""",dest="verbose", action="store", metavar="n",type=int, default=0)

	(options, args) = parser.parse_args()

	tiltseries=args[0]

	checkinput(options)

	logger = E2init(sys.argv, options.ppid)
	
	extension = os.path.splitext(os.path.basename(tiltseries))[-1]
	if options.verbose:
		print("\nextension is {}".format(extension))
	
	print("\nextension is {}".format(extension))
	
	tiltseriesout = tiltseries.replace(extension,'_lp50_hpp4' + extension)
	
	lowpass = 'filter.lowpass.tanh:cutoff_freq=0.02'
	if options.lowpass:
		lowpass = options.lowpass
		frequency = float(options.lowpass.split('cutoff_freq=')[-1].split(':')[0])
		resolution = str(int(round(1.0/frequency)))
		tiltseriesout = tiltseries.replace(extension,'_lp' + resolution + extension)
	
	if options.verbose:
		print("\lowpass filter to apply is {}".format(lowpass))

	highpass = 'filter.highpass.gauss:cutoff_pixels=4'
	if options.highpass:
		highpass = options.highpass
		pixels = options.highpass.split('cutoff_pixels=')[-1].split(':')[0]
		if options.lowpass:
			tiltseriesout = tiltseriesout.replace(extension,'_hpp' + pixels + extension)
		else:
			tiltseriesout = tiltseries.replace(extension,'_hpp' + pixels + extension)

	if options.verbose:
		print("\highpass filter to apply is {}".format(highpass))

	cmd = 'e2proc2d.py ' + tiltseries + ' ' + tiltseriesout + ' --process ' + lowpass + ' --process ' + highpass
	if options.replace:
		cmd += ' && cp ' + tiltseries +' backup.' + tiltseries + ' && cp ' + tiltseriesout +  ' ' + tiltseries
	
	runcmd(options,cmd)
	
	E2end(logger)
	
	return


if __name__ == '__main__':
	main()	
	

