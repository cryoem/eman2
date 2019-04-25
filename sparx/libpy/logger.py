
from __future__ import print_function
from builtins import object
#
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Please do not copy or modify this file without written consent of the author.
# Copyright (c) 2000-2019 The University of Texas - Houston Medical School
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
#

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

