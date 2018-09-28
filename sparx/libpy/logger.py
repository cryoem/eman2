
from __future__ import print_function
from builtins import object

class BaseLogger_Print(object):
	
	def logLine(self, prefix, line, *args, **kwargs):
		print(line)


class BaseLogger_Files(object):
	
# 	filename_prefix = "log_" 
# 	
# 	def __init__(self, prefix=None):
# 		if prefix == None:
# 			from time import localtime, strftime
# 			prefix = strftime("%Y%m%d%H%M%S", localtime())
# 		self.filename_prefix = prefix
	
	def logLine(self, prefix, line, file_name):
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
		from time import localtime, strftime
		
		if self.base_logger == None:
			return
		
		params = [param1, param2, param3, param4, param5, param6, param7, param8, param9, param10]
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
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

