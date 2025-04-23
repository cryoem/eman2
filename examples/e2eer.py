#!/usr/bin/env python
#
# Author: Steven Ludtke, 11/9/21 (sludtke@bcm.edu)
# Copyright (c) 2021- Baylor College of Medicine
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
import os
import sys

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <single image file> ...

This program will process eer files in various convenient ways.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="input",help="Input EER stack file", default="")
	
	
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=None, guitype='floatbox', row=3, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getAPIX()']")
	parser.add_argument("--constbfactor",type=float,help="Set B-factor to fixed specified value, negative value autofits",default=-1.0, guitype='floatbox', row=8, col=0, rowspan=1, colspan=1, mode='eval[-1.0]')
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=None, guitype='floatbox', row=3, col=1, rowspan=1, colspan=1, mode="eval['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=None, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--phaseplate",action="store_true",help="Include phase/amplitude contrast in CTF estimation. For use with hole-less phase plates.",default=False, guitype='boolbox', row=3, col=2, rowspan=1, colspan=1, mode='filter[False]')
	parser.add_argument("--astigmatism",action="store_true",help="Includes astigmatism in automatic fitting",default=False, guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='eval')
	parser.add_argument("--box",type=int,help="Forced box size in grid mode. Overrides any previous setting. ",default=-1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--usefoldername",action="store_true",help="If you have the same image filename in multiple folders, and need to import into the same project, this will prepend the folder name on each image name",default=False,guitype='boolbox',row=5, col=1, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv,options.ppid)
