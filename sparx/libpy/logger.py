
from __future__ import print_function
from builtins import object

import time

from global_def import sxprint

#-----------------------------------------------------------[ original Logger ]

class BaseLogger_Print(object):

	def logLine(self, prefix, line, *args, **kwargs):
		sxprint(line)


class BaseLogger_Files(object):

	def logLine(self, prefix, line, file_name):
		line = sxprint(line)
		fstr = open(prefix + file_name, 'a+')
		fstr.write(line + "\n")
		fstr.close()


class Logger(object):
	
	base_logger = None
	prefix = ""

	def __init__(self, base_logger=BaseLogger_Print(), base_logger2=None, prefix="", file_name="log.txt"):
		self.prefix = prefix
		self.base_logger = base_logger
		self.base_logger2 = base_logger2
		self.file_name = file_name

	def add(self, param1=None, param2=None, param3=None, param4=None, param5=None, param6=None, param7=None, param8=None, param9=None, param10=None):
				
		if self.base_logger == None:
			return
		
		params = [param1, param2, param3, param4, param5, param6, param7, param8, param9, param10]
		line = ""
		for p in params:
			if p != None:
				line += " " + str(p)

		#print(self.prefix)
		self.base_logger.logLine(self.prefix, line, self.file_name)
		if self.base_logger2:
			self.base_logger2.logLine(self.prefix, line, self.file_name)

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
