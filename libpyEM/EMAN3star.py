#!/usr/bin/env python
#
# Author: Steven Ludtke, 06/25/2023 (sludtke@bcm.edu)
# Copyright (c) 2023- Baylor College of Medicine
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

from builtins import range
import os
import os.path
import re
import traceback
import numpy as np

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
# Changes to the abstract object will be synchronized with the file when only explicitly requested.
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
	"""This is a more complete Star file implementation (than EMAN2 StarFile) which also supports formats like mmCIF.

	This will appear as a dictionary of dictionaries. The outermost level of the dictionary will be the "data_" elements.
	if no "data_" element is present, the file will be keyed as "data_default". Within the "data_" element, all keys
	are at the same level. loop_ block columns are aggregated as lists/nparrays at the top level """

	def __init__(self,filename):
		self.readfile(filename)

	def readfile(self,filename):
		self.filename=filename
		self.loops={}			# A dictionary keyed in the same way as the main dict, each with a list of lists with loop_ parameter names
		self.clear()			# the dictionary superclass, which is keyed by data_ elements in the file

		# default in case the file lacks a data_ definition. We erase this element if still empty at the end
		self.curdata="default"
		self.curdict={}
		self[self.curdata]=self.curdict
		self.loops[self.curdata]=[]

		matcher=re.compile(r"""("[^"]+")|('[^']+')|([^\s]+)""")

		fin=open(self.filename,"r")

		line=self.__readline(fin)
		while len(line)>0:
			if line[0]=="#":
				line=self.__readline(fin)
				continue
			# single element
			if line[0]=="_":
				key,val=self.__read_value(fin,line)
				self.curdict[key]=val
				line=self.__readline(fin)
			elif line[:5]=="data_":
				self.curdata=line.strip()[5:]
				self.curdict={}
				self[self.curdata]=self.curdict
				self.loops[self.curdata]=[]
				line=self.__readline(fin)
			elif line[:5]=="loop_":
				# first read the list of column identifiers
				line=self.__readline(fin)
				lsts=[]						# self.curdict keyed by column name, this is a list of the same lists to make it easier to populate
				parms=[]					# list of column headers
				while len(line)>0 and line[0]=="_":
					parms.append(line.split()[0][1:])	# add new column header to list
					self.curdict[parms[-1]]=[]		# an empty list for the new column
					lsts.append(self.curdict[parms[-1]])
					line=self.__readline(fin)		# read lines at the end rather than beginning of the loop

				self.loops[self.curdata].append(parms)		# store the list of column headers for potential later use

				while len(line)>0 and line[0]!="_" and line[:5]!="data_" and line[:5]!="loop_":
					for i,v in enumerate(matcher.findall(line)): lsts[i].append(goodval(v))
					line=self.__readline(fin)

				# try and convert all of the non string lists into numpy arrays
				p=self.loops[self.curdata][-1]		# list of column names we just created
				for k in p:
					try:
						x=float(self.curdict[k][0])		# if we can convert the first item to a float then we do the array conversion, not perfect, but will block most strings from conversion
						self.curdict[k]=np.array(self.curdict[k])
					except: pass

		if len(self["default"])==0: del self["default"]

	def __readline(self,fileobj):
		r=fileobj.readline()
		if len(r)==0: return r
		while r[0]=="#" or len(r.strip())==0:
			r=fileobj.readline()
			if len(r)==0: return r
		return r.strip()

	def __read_value(self,fileobj,line):
		"""reads the value either from the current line or additional lines if necessary, returns key,value tuple"""
		spl=line.split(None,1)
		key=spl[0][1:]
		if len(spl)==2:
			if spl[1][0] in ("'",'"') : return key,spl[1:-1]		# we assume the last non-whitespace character is the ending delimiter
			else:
				try: return key,int(spl[1])
				except:
					try: return key,float(spl[1])
					except: return key,spl[1]			# if not an int or a float, must be a simple value string

		line2=self.__readline(fileobj)
		if line2[0] in ("'",'"') :
			return key,line2.strip()[1:-1]
		elif line2[0]==";" :
			val=[line2[1:]]
			while 1:
				try: line2=self.__readline(fileobj)
				except: raise Exception("StarFile: Error found parsing multi-line string value ",key)
				if line2[0]==';' : break
				val.append(line2)
			val[-1]=val[-1].rstrip()		# remove trailing whitespace on the last line
			val="".join(val)
			return key,val

		raise Exception("StarFile: Key-value pair error. Matching value not found. ",key)



