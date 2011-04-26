#!/usr/bin/env python
#
# Author: Pawel A.Penczek and Edward H. Egelman 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
# Copyright (c) 2008-Forever The University of Virginia
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

def main():
	import os
	import sys
	from optparse import OptionParser
	from global_def import SPARXVERSION
	import global_def
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + " stack --params='parm1 parm2 parm3 ...' --zero --one --randomize --rand_alpha --import=file --export=file --print --backup --suffix --restore --delete"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--params",	   type="string",       default=None,    help="parameter list")
	parser.add_option("--zero",	   action="store_true", default=False, help="set parameter to zero")
	parser.add_option("--one",	   action="store_true", default=False, help="set parameter to one")
	parser.add_option("--randomize",   action="store_true", default=False, help="set parameter to randomized value")
	parser.add_option("--rand_alpha",  action="store_true", default=False, help="set all angles to randomized value")
	parser.add_option("--import",	   type="string", dest="fimport", default=None, help="import parameters from file")
	parser.add_option("--export",	   type="string", dest="fexport", default=None, help="export parameters to file")
	parser.add_option("--print",	   action="store_true", dest="fprint", default=False, help="print parameters")
	parser.add_option("--backup",	   action="store_true", default=False, help="backup parameters")
	parser.add_option("--suffix",	   type="string",       default="_backup",    help="suffix for xform name in backup")
	parser.add_option("--restore",   action="store_true", default=False, help="restore parameters")
	parser.add_option("--delete",    action="store_true", default=False, help="delete parameters")

	(options,args) = parser.parse_args( arglist[1:] )

	if len(args) != 1 :
		print usage
		sys.exit(-1)

	if options.params == None:
		print "Error: no parameters given"
		exit(-1)

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
        from applications import header
	header(args[0], options.params, options.zero, options.one, options.randomize, options.rand_alpha, options.fimport, options.fexport, options.fprint, options.backup, options.suffix, options.restore, options.delete)

if __name__ == "__main__":
	main()
