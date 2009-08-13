#!/usr/bin/env python

#
# Author: David Woolford 08/13/08 (woolford@bcm.edu
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from emapplication import EMStandAloneApplication
from emimage3dsym import EM3DSymViewerModule,EMSymInspector
from e2eulerxplor import InputEventsManager
import os,sys
from EMAN2 import *
from optparse import OptionParser
	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog  <simmx file> <projection file>
	
Simmx xplor for comparator evaluation
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	(options, args) = parser.parse_args()
	
	
	logid=E2init(sys.argv)
	
	em_app = EMStandAloneApplication()
	window = EMSimmxExplorer(application=em_app)

	em_app.show()
	em_app.execute()
	
	E2end(logid)


class EMSimmxExplorer(EM3DSymViewerModule):
	def get_desktop_hint(self): return "image"
	
	def __init__(self,application=None,ensure_gl_context=True,application_control=True,projection_file="",simmx_file=""):
		self.init_lock = True # a lock indicated that we are still in the __init__ function
		self.au_data = None # This will be a dictionary, keys will be refinement directories, values will be something like available iterations for visual study	
		EM3DSymViewerModule.__init__(self,application,ensure_gl_context=ensure_gl_context,application_control=application_control)
		#InputEventsManager.__init__(self)
	
		self.project_file = projection_file
		self.simmx_file = simmx_file
		
	def set_projection_file(self,file): self.project_file = file
	def set_simmx_file(self,file): self.simmx_file = file
	
	
if __name__ == '__main__':
	main()