#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

import sys
import os
try:
	import wx
	from wx import py
except:
	print "WxPython 2.6 or above is required to run the EMAN2 GUI environment"
	sys.exit(1)

from optparse import OptionParser
#from panels import *
from EMAN2 import *
from wxgui import *

class EMAN2(wx.App):
	def OnInit(self):
		wx.InitAllImageHandlers()
		
#		self.paramdialog = ParamDialog(None, -1, "")
#		self.SetTopWindow(self.paramdialog)
#		self.paramdialog.Show()
#		self.frame = py.crust.CrustFrame()
		self.frame = py.shell.ShellFrame()
		self.frame.SetSize((800, 600)
		self.frame.Show()
		self.SetTopWindow(self.frame)
		self.timer=E2Timer()
		self.timer.Start(500)
#		wx.EVT_TIMER(self, self.OnTimer)
		return True
		
class E2Timer(wx.Timer):
	def Notify(self):
		for i in EMImage.allim.keys():
#			try: 
#				print i.data.get_attr("changecount"),i.changec
				if i.rdata.get_attr("changecount")!=i.changec :
					i.setdata(i.rdata)
#			except:pass
	
if __name__ == "__main__":
	global options
	
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] <image>
	
Local interface to EMEN2."""

	parser = OptionParser(usage=usage,version="0.1")

#	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
#	parser.add_option("--box","-B",type="int",help="Box size in pixels",default=-1)
#	parser.add_option("--web","-W",type="string",help="Access the database on a remote machine via XMLRPC. Provide the URL of the database.",default=None)
		
	(options, args) = parser.parse_args()

	app = EMAN2(0)
	sys.app=app
	app.MainLoop()
