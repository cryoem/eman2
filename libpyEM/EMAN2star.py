#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/17/14 (sludtke@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

import os
import os.path
import re
import traceback

#from libpyEMData2 import EMData
#from libpyUtils2 import EMUtil

######
# This module implements access to STAR files (includes files used by Relion 
# and other software, as well as CIF/MMCIF. It does not implement the full STAR specification.
# We do try to cover the CIF specification reasonably well, though:
# http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
# 
# 
# The STAR file is represented as a dictionary-like abstraction of the physical file on disk. 
# Changes to the abstract object will be syncronized with the file when only explicitly requested.
# The entire file is kept in RAM, and the file is only re-read from disk if 'readfile()' is
# explicitly called (which will overwrite any changes in memory). 
#
# There is no support for schema, or constrained datatypes for values. Values may be int, float or string.
#
# The JSON file mechanism is  much more robust and greatly preferred for EMAN2 applications. 
# This module is mainly used for inter-software compatibility.
#
# keys have the leading "_" stripped off
#
# loop values are represented as a python list. keys from the same loop should have an identical number
# of elements. loops are identified internally as a list of lists (self.loop) independent of the
# actual data storage.
######

def goodval(vals): 
	val=max(vals)
	try: val=int(val)
	except:
		try: val=float(val)
		except: pass
	return val

class StarFile(dict):
	
	def __init__(self,filename):
		dict.__init__(self)
		self.filename=filename
		self.loops=[]
		
		if os.path.isfile(filename) :
			self.readfile()
	def _nextline(self):
		"""Used internally when parsing a star file to emulate readline"""
		self.lineptr+=1
		return self.lines[self.lineptr-1]
	
	def readfile(self):
		"""This parses the STAR file, replacing any previous contents in the dictionary"""
		
		self.loops=[]
		self.clear()
		
		matcher=re.compile("""("[^"]+")|('[^']+')|([^\s]+)""")
		
		# read the entire file into a buffer, this dramatically simplifies the logic, even if it eats a chunk of RAM
		self.lines=[i for i in file(self.filename,"r") if len(i.strip())!=0 and i[0]!="#"]
		self.lineptr=0
		
		while 1:
			try: line=self._nextline().strip()
			except: break
		
			if line[0]=="_" :				# A single key/value pair
				spl=line.split(None,1)		# split on whitespace
				key=spl[0][1:]
				
				if len(spl)==2:				# value on the same line
					if spl[1][0] in ("'",'"') : self[key]=spl[1:-1]		# we assume the last non-whitespace character is the ending delimeter
					else:
						try: val=int(spl[1])
						except: 
							try: val=float(spl[1])
							except: val=spl[1]			# if not an int or a float, must be a simple value string
					
						self[key]=val
				else:						# value starts on next line
					line2=self._nextline()
					if line2[0] in ("'",'"') :
						self[key]=line2.strip()[1:-1]
					elif line2[0]==";" :
						val=[line2[1:]]
						while 1:
							try: line2=self._nextline()
							except: raise Exception,"StarFile: Error found parsing multi-line string value for %s"%key
							if line2[0]==';' : break
							val.append(line2)
						val[-1]=val[-1].rstrip()		# remove trailing whitespace on the last line
						val="".join(val)
						self[key]=val
					else: raise Exception,"StarFile: Key-value pair error. Matching value for %s not found."%key
			elif line[:5].lower()=="data_":
				if len(self)>0 :
					print "WARNING: second data_ block encountered in ",self.filename,". Cannot deal with this at present. Second block ignored"
					return
				self.dataname=line[5:]
			elif line[:5].lower()=="loop_":
				loop=[]
				self.loops.append(loop)				# add it to the list of loops immediately then update it as we go
				# First we read the parameter names for the loop
				while 1:
					line2=self._nextline().strip()
					if line2[0]=="_":
						loop.append(line2.split()[0][1:])
						self[loop[-1]]=[]			# this will hold the data values when we read them
					else: break
				self.lineptr-=1
				
				# Now we read the actual loop data elements
				vals=[]
				while 1:
					try: line2=self._nextline().strip()
					except: break
				
					if line2[0]=="_" or line2.lower()[:5]=="loop_" : break
					elif line2[0]==";" :
						val=line2[0][1:]
						while 1:
							line2=self._nextline()
							if line2[0]==";":
								break
							val+=line2
						vals.append(val)
					else:
						vals.extend([goodval(i) for i in matcher.findall(line2)])
						if len(vals)<len(loop) :			# we may need to read multiple lines to get enough values
							continue
						if len(vals)>len(loop) : 
							print "mismatch"
							print line2
							print len(loop),loop
							print len(vals),vals
							break
						for i in range(len(vals)): self[loop[i]].append(vals[i])
						vals=[]
				self.lineptr-=1
			else:
				print "StarFile: Unknown content on line :",line
				break

				
	def writefile(self,filename=None):
		"""Writes the contents of the current dictionary back to disk using either the existing filename, or an alternative name passed in"""
		
		print "Sorry, writing not implemented yet"
					
			
			