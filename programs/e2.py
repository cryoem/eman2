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
from optparse import OptionParser
import EMAN2
from emimage import image_update
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
import threading
from IPython.Shell import *
from emapplication import EMStandAloneApplication

class IPShellQt4a(threading.Thread):
	"""Run a Qt event loop in a separate thread.

	Python commands can be passed to the thread where they will be executed.
	This is implemented by periodically checking for passed code using a
	Qt timer / slot."""

	TIMEOUT = 100 # Millisecond interval between timeouts.

	def __init__(self,argv=None,user_ns=None,user_global_ns=None,
				debug=0,shell_class=MTInteractiveShell):

		from PyQt4 import QtGui

		#class newQApplication2:
			#def __init__( self ):
				#self.QApplication = QtGui.QApplication

			#def __call__( *args, **kwargs ):
				#return QtGui.qApp

			#def exec_loop( *args, **kwargs ):
				#pass

			#def __getattr__( self, name ):
				#return getattr( self.QApplication, name )

		self.app = EMStandAloneApplication()
	
		QtGui.QApplication = self.app.get_app()

		# Allows us to use both Tk and QT.
		self.tk = get_tk()

		self.IP = make_IPython(argv,user_ns=user_ns,
							user_global_ns=user_global_ns,
							debug=debug,
							shell_class=shell_class,
							on_kill=[QtGui.qApp.exit])

		# HACK: slot for banner in self; it will be passed to the mainloop
		# method only and .run() needs it.  The actual value will be set by
		# .mainloop().
		self._banner = None

		threading.Thread.__init__(self)
    
	def get_app(self):
		return self.app
	
	def run(self):
		self.IP.runlines("from EMAN2 import *\nfrom emimage import *\n")
		self.IP.mainloop(self._banner)
		self.IP.kill()

	def mainloop(self,sys_exit=0,banner=None):

		from PyQt4 import QtCore, QtGui

		self._banner = banner

		if QtGui.QApplication.startingUp():
			a = QtGui.QApplication.QApplication(sys.argv)
		self.timer = QtCore.QTimer()
		QtCore.QObject.connect( self.timer, QtCore.SIGNAL( 'timeout()' ), self.on_timer )

		self.start()
		self.timer.start( self.TIMEOUT )
		while True:
			if self.IP._kill: break
			QtGui.qApp.exec_()
		self.join()

	def on_timer(self):
		update_tk(self.tk)
		result = self.IP.runcode()
		image_update()
		self.timer.start( self.TIMEOUT )
		return result


if __name__ == "__main__":
	
	EMAN2.GUIMode=True
	
	sh=IPShellQt4a()
	EMAN2.app = sh.get_app()
#	app = get_app()
	sh.mainloop(banner="Welcome to EMAN2\nPrompt provided by IPython\nEnter '?' for ipython help\n")
	
	sys.exit(0)
