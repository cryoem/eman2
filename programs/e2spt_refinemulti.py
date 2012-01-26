#!/usr/bin/env python

#
# Author: Jesus Galaz  21/1/2012 (rewritten)
# Copyright (c) 2011- Baylor College of Medicine
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

from sys import argv
import os
from EMAN2 import *

current = os.getcwd()
filesindir = os.listdir(current)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <stack>

	WARNING: This program still under development.
	
	Refinement of a 3D-volume stack against multiple models. The initial models can be provided in a stack, OR generated from the data itself.
	
	When no reference is provided, you define a number of models greater than 2 to use (for example, 4).
	[If you want to refine the data against ONE model, use e2classaverage3d.py]
	The data set is divded into that same number of groups (4). 
	An initial model will be generated with the data in each group.
	Then, the entire data set will be refined against all 4 initial models.
	
	You can increase the number of references used for each iteration by specifying the --addmodel parameter.
	This will take the "best initial model" (the one that most particles prefered) and include it as an initial model for the next round of refinement.
	For exampe, if you start with two references A and B, two averages will come out of aligning the data against them, avgA and avgB.
	So if --addmodel is on, instead of only using avgA and avgB as references for the next refinement round, the best of A and B will also be used,
	which means you will refine the data against 3 models in the next round, not just 2.
	
	If you supply a single reference/model then --addmodel MUST be supplied too; otherwise, to refine a data set against a single model use
	e2classaverage3d.py
	 """

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ncls", type=int, help="...", default=2)
	parser.add_argument("--nbasis", type=int, help="Basis vectors to use", default=3)

	parser.add_argument("--input", type=str, help="The name of the volumes stack that HAVE BEEN ALIGNED to a common reference", default=None)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")			
	parser.add_argument("--filter",type=int,help="Resolution (integer, in Angstroms) at which you want to apply a gaussian lowpass filter to the tomogram prior to loading it for boxing",default=0, guitype='intbox', row=3, col=2, rowspan=1, colspan=1)
	parser.add_argument("--preprocess",type=str,help="""A processor (as in e2proc3d.py) to be applied to the tomogram before opening it. \nFor example, a specific filter with specific parameters you might like. \nType 'e2proc3d.py --processors' at the commandline to see a list of the available processors and their usage""",default=None)
		
	(options, args) = parser.parse_args()
