#!/usr/bin/env python

#
# Author: Steven Ludtke (sludtke@bcm.edu) 11/07/2012
# Copyright (c) 2000-2012 Baylor College of Medicine
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
from EMAN2fsc import *
from EMAN2db import db_open_dict, db_close_dict, db_check_dict
from math import *
import os
import sys
import traceback
import e2refine

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	
This program combines standard refinement as embodied by e2refine.py with robust resolution testing implemented in e2refine_evenodd.py. The 
even/odd test is run simultaneously with the refinement. In addition the complicated gamut of options provided by e2refine.py are largely
selected automatically. The user needs to provide basic experimental parameters or goals.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the set containing the particle data", browser='EMSetsTable(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--model", dest="model", type=str,default="threed.0a.mrc", help="The name of the 3D image that will seed the refinement", guitype='filebox', browser='EMModelsTable(withmodal=True,multiselect=False)', row=5, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--targetres",type=float,default=25.0,help="Targeted resolution in angstroms (eg - 8, not 0.125). This is used to determine many different parameters. Avoid over-optimism in early runs.")
	parser.add_argument("--path", default=None, type=str,help="The name of a directory where results are placed. If unspecified will generate one automatically of type refine_??.")
	parser.add_argument("--mass", default=0, type=float,help="The mass of the particle in kilodaltons, used to run normalize.bymass. If unspecified (set to 0) nothing happens. Requires the --apix argument.", guitype='floatbox', row=2, col=1, rowspan=1, colspan=1, mode="refinement['self.pm().getMass()']")
	parser.add_argument("--apix", default=0, type=float,help="The angstrom per pixel of the input particles. This argument is required if you specify the --mass argument. \n If unspecified (set to 0), the convergence plot is generated using either the project apix, or if not an apix of 1.", guitype='floatbox', row=2, col=0, rowspan=1, colspan=1, mode="refinement['self.pm().getAPIX()']")
	parser.add_argument("--sym", dest = "sym", default="c1", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos. For asymmetric reconstruction omit this option or specify c1.", guitype='symbox', row=10, col=0, rowspan=1, colspan=3, mode="refinement")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>:<value> EX thread:4",default=None, guitype='strbox', row=3, col=0, rowspan=1, colspan=2, mode="refinement")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	
	