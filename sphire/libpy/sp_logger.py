
from __future__ import print_function
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
class BaseLogger_Print(object):

	def logLine(self, prefix, line, *args, **kwargs):
		sp_global_def.sxprint(line)


class BaseLogger_Files(object):

	def logLine(self, prefix, line, file_name, do_print):
		if do_print:
			line = sp_global_def.sxprint(line)
		else:
			t = sp_global_def.get_timestamp()
			f = sys._getframe(1).f_code.co_name
			m = t + " " + f + " => " + line
		fstr = open(prefix + file_name, 'a+')
		fstr.write(line + "\n")
		fstr.close()


class Logger(object):
	
	base_logger = None
	prefix = ""

	def __init__(self, base_logger=BaseLogger_Print(), base_logger2=None, prefix="", file_name="log.txt", do_print=True):
		self.prefix = prefix
		self.base_logger = base_logger
		self.base_logger2 = base_logger2
		self.file_name = file_name
		self.do_print = do_print

	def add(self, param1=None, param2=None, param3=None, param4=None, param5=None, param6=None, param7=None, param8=None, param9=None, param10=None):
				
		if self.base_logger == None:
			return
		
		params = [param1, param2, param3, param4, param5, param6, param7, param8, param9, param10]
		line = ""
		for p in params:
			if p != None:
				line += " " + str(p)

		#print(self.prefix)
		self.base_logger.logLine(self.prefix, line, self.file_name, self.do_print)
		if self.base_logger2:
			self.base_logger2.logLine(self.prefix, line, self.file_name, self.do_print)

	def sublog(self, add_prefix):
		logger = Logger()
		logger.base_logger = self.base_logger
		logger.prefix = self.prefix + add_prefix
		return logger

#--------------------------------------------------------[ Simple Logger ][NEW]

"""
The simple logger does the same as the default logger above. It is simpler
though, easier to use, and more flexible. It is not clear whether replacing the
original logger is safe to do though so, for now, we keep both.
"""
class SimpleLogger( object ):
	"""
	Simple class to log lines during execution. Makes use of sxprint to add
	timestamps and caller information to each line.

	Constructor args:

		filename (string): If provided the log is also printed to a file of the
			given name. Otherwise the log is only printed to the screen.
			[Default: ""]
	
		print_log (bool): Set to False to not print the log lines to the screen.
			NOTE: This only makes sense when a file name is given, otherwise the
			log is not printed at all.
	"""
	def __init__( self, filename="", print_log=True ):
		self.filename  = filename
		self.print_log = print_log
		self.log( "Starting log" )

	def log( self, *args, **kwargs ):

		if self.filename != "":
			sxprint( *args, filename=self.filename, **kwargs )
		else:
			sxprint( *args, **kwargs )
from builtins import object

import time

import sys
import sp_global_def

#-----------------------------------------------------------[ original Logger ]

