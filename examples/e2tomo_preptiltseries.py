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

from EMAN2 import *


def main():

	usage = """Program to generate a cylindrical mask. It can also create a cylindrical shell if
			you specify the --height_inner and --radius_inner parameters, in addition to the
			required outer --height and --radius.
			"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--input",type=float,default=None,help="""Default=None. Sampling size (angstroms per pixel) to set on the header of the cylindrical mask.""")
	
	parser.add_argument("--highpass",type=str,default='',help="Default=filter.highpass.gauss:cutoff_pixels=4 (zeros out the first four Fourier pixels)_. A highpass filtering processor (as in e2proc3d.py) to be applied to the tiltseries.")

	parser.add_argument("--lowpass",type=str,default='',help="Default=filter.lowpass.tanh:cutoff_freq=0.02 (filters to 50 angstroms resolution). A lowpass filtering processor (as in e2proc3d.py) to be applied to the tiltseries.")

	(options, args) = parser.parse_args()

	tiltseries=args[0]

	logger = E2init(sys.argv, options.ppid)

	E2end(logger)
	return


if __name__ == '__main__':
	main()	
	

