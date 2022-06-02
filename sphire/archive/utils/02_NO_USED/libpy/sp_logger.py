



































































































from __future__ import print_function
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
			sp_global_def.sxprint( *args, filename=self.filename, **kwargs )
		else:
			sp_global_def.sxprint( *args, **kwargs )









