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
from emimage import *

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
		self.timer=E2Timer()
		self.timer.Start(500)
#		wx.EVT_TIMER(self, self.OnTimer)
		return True
		
class E2Timer(wx.Timer):
	def Notify(self):
		for i in EMImage.allim:
			try: 
#				print i.data.get_attr("changecount"),i.changec
				if i.data.get_attr("changecount")!=i.changec : i.changed()
			except: pass
	
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
