#!/usr/bin/env python

#
# Author: David Woolford 04/16/2009 (woolford@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

from optparse import OptionParser

tomo_db_name = "bdb:e2tomoboxercache#"

def check(args,options):
	error = False


class EMBootStrappedAverages:
	'''
	This class breaks the jobs of the boot-strapped average generation procedure
	so that they can be run in parallel. 
	'''
	def __init__(self,options,args):
		'''
		@param options - the options returned by the call to (options, args) = parser.parse_args() 
		@param args - a list of image names - args is that which is returned by the call to (options, args) = parser.parse_args()
		'''
		self.options = options
		self.args = args


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image1> <image2> <image3> <image4> ....
	
Boot straps an initial probe doing all versus all alignment of the input images

"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=64)
	parser.add_option("--writeoutput",action="store_true",default=False,help="Uses coordinates stored in the tomoboxer database to write output")
	parser.add_option("--align",type="string",help="The aligner and its parameters. e.g. --align=rot.3d.grid:ralt=180:dalt=10:dphi=10:rphi=180:search=5", default=None)
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	
	(options, args) = parser.parse_args()
	
	error = check(args,options)
	
	
	