#!/bin/env python

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

class EMAN2(wx.App):
	def OnInit(self):
		wx.InitAllImageHandlers()
		
#		self.paramdialog = ParamDialog(None, -1, "")
#		self.SetTopWindow(self.paramdialog)
#		self.paramdialog.Show()
		self.frame = py.crust.CrustFrame()
		self.frame.SetSize((800, 600))
		self.frame.Show()
		self.SetTopWindow(self.frame)
		return True
			
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
