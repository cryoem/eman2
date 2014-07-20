

class BaseLogger_Print:
	
	def logLine(self, prefix, line):
		print (prefix + " " + line)


class BaseLogger_Files:
	
# 	filename_prefix = "log_"
# 	
# 	def __init__(self, prefix=None):
# 		if prefix == None:
# 			from time import localtime, strftime
# 			prefix = strftime("%Y%m%d%H%M%S", localtime())
# 		self.filename_prefix = prefix
	
	def logLine(self, prefix, line):
		fstr = open(prefix + "log.txt", 'a+')
		fstr.write(line + "\n")
		fstr.close()


class Logger:
	
	base_logger = None
	prefix = ""

	def __init__(self, base_logger=BaseLogger_Print(), prefix=""):
		self.prefix = prefix
		self.base_logger = base_logger

	def add(self, param1=None, param2=None, param3=None, param4=None, param5=None, param6=None, param7=None, param8=None, param9=None, param10=None):
		from time import localtime, strftime
		
		if self.base_logger == None:
			return
		
		params = [param1, param2, param3, param4, param5, param6, param7, param8, param9, param10]
		line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
		for p in params:
			if p != None:
				line += " " + str(p)

		self.base_logger.logLine(self.prefix, line)

	def sublog(self, add_prefix):
		logger = Logger()
		logger.base_logger = self.base_logger
		logger.prefix = self.prefix + add_prefix
		return logger

