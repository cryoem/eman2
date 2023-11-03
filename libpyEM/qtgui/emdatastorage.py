#!/usr/bin/env python
#
# Author: this is a direct copy of emen2/Database/datastorage.py
# Date Copied: Nov 7 2008
# Copyright (c) 2000-2008 Baylor College of Medicine
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

from past.utils import old_div
from collections.abc import MutableMapping as DictMixin

from math import *
import time
import re
# validation/conversion for booleans
def boolconv(x):
	try:
		x=int(x)
		if x: return 1
		else: return 0
	except:
		if x[0] in ("T","t","Y","y") : return 1
		if x[0] in ("F","f","N","n") : return 0
		raise Exception("Invalid boolean %s"%str(x))

# validation/conversion for text, permitting unicode
def textconv(x):
	try: return str(x)
	except: return str(x)

def datetimeconv(x):
	return str(x)
	
def timeconv(x):
	return str(x)

def dateconv(x):
	return str(x)

valid_vartypes={
	"int":("d",lambda x:int(x)),			# 32-bit integer
	"longint":("d",lambda x:int(x)),		# not indexed properly this way
	"float":("f",lambda x:float(x)),		# double precision
	"longfloat":("f",lambda x:float(x)),	# arbitrary precision, limited index precision
	"choice":("s",textconv),			# string from a fixed enumerated list, eg "yes","no","maybe"
	"string":("s",textconv),			# a string indexed as a whole, may have an extensible enumerated list or be arbitrary
	"text":("s",textconv),			# freeform text, fulltext (word) indexing, UNICODE PERMITTED
	"time":("s",timeconv),			# HH:MM:SS
	"date":("s",dateconv),			# yyyy/mm/dd
	"datetime":("s",datetimeconv),		# yyyy/mm/dd HH:MM:SS
	"intlist":(None,lambda y:[int(x) for x in y]),		# list of integers
	"floatlist":(None,lambda y:[float(x) for x in y]), # list of floats
	"stringlist":(None,lambda y:[str(x) for x in y]),	# list of enumerated strings
	"url":("s",lambda x:str(x)),			# link to a generic url
	"hdf":("s",lambda x:str(x)),			# url points to an HDF file
	"image":("s",lambda x:str(x)),			# url points to a browser-compatible image
	"binary":("s",lambda y:[str(x) for x in y]),				# url points to an arbitrary binary... ['bdo:....','bdo:....','bdo:....']
	"binaryimage":("s",lambda x:str(x)),		# non browser-compatible image requiring extra 'help' to display... 'bdo:....'
	"child":("child",lambda y:[int(x) for x in y]),	# link to dbid/recid of a child record
	"link":("link",lambda y:[int(x) for x in y]),		# lateral link to related record dbid/recid
	"boolean":("d",boolconv),
	"dict":(None, lambda x:x), 
	"user":("s",lambda x:str(x)), # ian 09.06.07
	"userlist":(None,lambda y:[str(x) for x in y])
}

# Valid physical property names
# The first item in the value tuple is ostensibly a default, but this
# will generally be provided by the ParamDef. It may be that
# synonyms should be combined in a better way
# valid_properties = { 
# "count":(None,{"k":1000, "K":1000, "pixels":1}),
# "unitless":(None,{"n/a": None}),
# "length":("meter",{"m":1.,"meters":1,"km":1000.,"kilometer":1000.,"cm":0.01,"centimeter":0.01,"mm":0.001,
#		"millimeter":0.001, "um":1.0e-6, "micron":1.0e-6,"nm":1.0e-9,"nanometer":1.0e-9,"angstrom":1.0e-10,
#		"A":1.0e-10}),
# "area":("m^2",{"m^2":1.,"cm^2":1.0e-4}),
# "volume":("m^3",{"m^3":1,"cm^3":1.0e-6,"ml":1.0e-6,"milliliter":1.0e-6,"l":1.0e-3, "ul":1.0e-9, "uL":1.0e-9}),
# "mass":("gram",{"g":1.,"gram":1.,"mg":.001,"milligram":.001,"Da":1.6605387e-24,"KDa":1.6605387e-21, "dalton":1.6605387e-24}),
# "temperature":("K",{"K":1.,"kelvin":1.,"C":lambda x:x+273.15,"F":lambda x:(x+459.67)*5./9.,
#		"degrees C":lambda x:x+273.15,"degrees F":lambda x:(x+459.67)*5./9.}),
# "pH":("pH",{"pH":1.0}),
# "voltage":("volt",{"V":1.0,"volt":1.0,"kv":1000.0,"kilovolt":1000.0,"mv":.001,"millivolt":.001}),
# "current":("amp",{"A":1.0,"amp":1.0,"ampere":1.0}),
# "resistance":("ohm",{"ohm":1.0}),
# "inductance":("henry",{"H":1.0,"henry":1.0}),
# "transmittance":("%T",{"%T":1.0}),
# "relative_humidity":("%RH",{"%RH":1.0}),
# "velocity":("m/s",{"m/s":1.0}),
# "momentum":("kg m/s",{"kg m/s":1.0}),
# "force":("N",{"N":1.0,"newton":1.0}),
# "energy":("J",{"J":1.0,"joule":1.0}),
# "angle":("degree",{"degree":1.0,"deg":1.0,"radian":180.0/pi, "mrad":0.18/pi}),
# "concentration":("mg/ml", {"mg/ml":1.0, "p/ml":1.0, "pfu":1.0}),
# "resolution":('A/pix', {'A/pix':1.0}),
# "bfactor":('A^2', {"A^2":1.0, "A2":1.0}),
# "dose":('e/A2/sec', {'e/A2/sec':1.0}),
# "currentdensity":('Pi Amp/cm2', {'Pi Amp/cm2':1.0}),
# "filesize": ('bytes', {'bytes':1.0, 'kb':1.0e3, 'Mb':1.0e6, 'GB':1.0e9}),
# "percentage":('%', {'%':1.0}),
# "currency":("dollars",{"dollars":1.0}),
# "pressure":("Pa",{"Pa":1.0,"pascal":1.0,"bar":1.0e-5,"atm":9.8692327e-6,"torr":7.500617e-6,"mmHg":7.500617e-6,"psi":1.450377e-4}),
# "unitless":("unitless",{"unitless":1})
# }

# ian: added better synonym support
valid_properties = {

"time": ("s", 
	{"s":1., "min":60, "hour":3600, "day":86400, "ms":.001, "us":1.0e-6, "ns": 1.0e-9}, 
	{"hours":"hour", "mins":"min", "seconds":"s", "secs":"s", "sec": "s", "days":"day", "nanosecond":"ns", "nanoseconds":"ns", "microseconds":"us", "microsecond":"us", "milliseconds":"ms", "millisecond":"ms"}
	),

"length":("m",
	{"m":1., "km":1000., "cm":0.01, "mm":0.001, "um":1.0e-6, "nm":1.0e-9, "A":1.0e-10},
	{"meters":"m", "meter":"m", "kilometer":"km", "kilometers":"km", "centimeter":"cm", "centimeters":"cm", "millimeter":"mm", "millimeters":"mm", "micron":"um", "microns": "um", "nanometer":"nm", "nanometers":"nm", "angstrom":"A", "Angstroms":"A", "angstroms":"A", "Angstrom":"A"}
	),
	
"count":("count",
	{"K":1000, "pixels":1, "count":1},
	{"k":"K"}
	),
	
"unitless":("unitless", 
	{"unitless":1},
	{}
	),
	
"area":("m^2",
	{"m^2":1.,"cm^2":1.0e-4},
	{}
	),
	
"volume":("m^3",
	{"m^3":1,"ml":1.0e-6,"l":1.0e-3,"ul":1.0e-9,"ul":1.0e-9},
	{"cm^3":"ml", "milliliter":"ml", "milliliters":"ml", "uL":"ul"}
	),	
	
"mass":("gram",
	{"g":1.,"mg":.001,"Da":1.6605387e-24,"KDa":1.6605387e-21, "MDa":1.6605387e-18},
	{"gram":"g", "grams":"g", "milligram":"mg", "milligrams":"mg", "dalton":"Da", "daltons":"Da", "kilodaltons":"KDa", "kilodalton":"KDa", "megadaltons":"MDa", "megadalton":"MDa"}
	),	

"temperature":("K",
	{"K":1.,"C":lambda x:x+273.15,"F":lambda x:(x+459.67)*5./9.},
	{"kelvin":"K","degrees C":"C", "degrees F":"F"}
	),
	
"pH":("pH",
	{"pH":1.0},
	{}
	),
	
"voltage":("volt",
	{"V":1.0,"kv":1000.0,"mv":.001},
	{"volt":"V", "volts":"V", "kilovolt":"kv", "kilovolts":"kv", "millivolt":"mv", "millivolts":"mv"}
	),	
	
"current":("amp",
	{"amp":1.0},
	{"ampere":"amp"}
	),
	
"resistance":("ohm",
	{"ohm":1.0},
	{}
	),
	
"inductance":("henry",
	{"H":1.0},
	{"henry":"H"}
	),
	
"transmittance":("%T",
	{"%T":1.0},
	{}
	),
	
"relative_humidity":("%RH",
	{"%RH":1.0},
	{}
	),
	
"velocity":("m/s",
	{"m/s":1.0},
	{}
	),
	
"momentum":("kg m/s",
	{"kg m/s":1.0},
	{}
	),
	
"force":("N",
	{"N":1.0},
	{"newton":"N"}
	),
	
"energy":("J",
	{"J":1.0},
	{"joule":"J"}
	),
	
"angle":("degree",
	{"degree":1.0,"radian":old_div(180.0,pi), "mrad":old_div(0.18,pi)},
	{"deg":"degree", "degrees":"degree"}
	),
	
"concentration":("mg/ml", 
	{"mg/ml":1.0, "p/ml":1.0, "pfu":1.0},
	{}
	),

"resolution":('A/pix', 
	{'A/pix':1.0},
	{}
	),
	
"bfactor":('A^2',
	{"A^2":1.0},
	{"A2":"A^2"}
	),
	
"dose":('e/A2/sec',
	{'e/A2/sec':1.0},
	{'e/A^2/sec':'e/A2/sec'}
	),
	
"exposure":('e/A2',
	{'e/A2':1.0},
	{'e/A^2':'e/A2'}
	),
	
"currentdensity":('Pi Amp/cm2',
	{'Pi Amp/cm2':1.0},
	{}
	),
	
"filesize": ('bytes',
	{'bytes':1.0, 'kB':1.0e3, 'MB':1.0e6, 'GB':1.0e9, 'KiB':1024, 'MiB': 1048576, 'GiB': 1073741824},
	{'B':'bytes'}
	),
	
"percentage":('%',
	{'%':1.0},
	{}
	),
	
"currency":("dollars",
	{"dollars":1.0},
	{}
	),
	
"pressure":("Pa",
	{"Pa":1.0,"bar":1.0e-5,"atm":9.8692327e-6,"torr":7.500617e-6,"psi":1.450377e-4},
	{"pascal":"Pa", "mmHg":"torr"}
	),
	
}


class ParamDef(DictMixin) :
	"""This class defines an individual data Field that may be stored in a Record.
	Field definitions are related in a tree, with arbitrary lateral linkages for
	conceptual relationships. The relationships are handled externally by the
	Database object. Fields may only be modified by the administrator once
	created, and then, they should only be modified for clarification
	
	Note that if the vartype is boolean then sometimes the attribute "dependents" will be assigned, for example in emsprworkflow.py.
	This is just so dependent widgets can be automatically disabled/enabled in emform.py
	
	""" 
	
	# non-admin users may only update descs and choices
	attr_user = set(["desc_long","desc_short","choices"])
	attr_admin = set(["name","vartype","defaultunits","property","creator","creationtime","creationdb"])
	attr_all = attr_user | attr_admin
	
	# name may be a dict; this allows backwards compat dictionary initialization
	def __init__(self,name=None,vartype=None,desc_short=None,desc_long=None,
						property=None,defaultunits=None,choices=None):
		self.name=name					# This is the name of the paramdef, also used as index
		self.vartype=vartype			# Variable data type. List of valid types in the module global 'vartypes'
		self.desc_short=desc_short		# This is a very short description for use in forms
		self.desc_long=desc_long		# A complete description of the meaning of this variable
		self.property=property			# Physical property represented by this field, List in 'properties'
		self.defaultunits=defaultunits	# Default units (optional)
		self.choices=choices			# choices for choice and string vartypes, a tuple
		self.creator=None				# original creator of the record
		self.creationtime=time.strftime("%Y/%m/%d %H:%M:%S")	# creation date
		self.creationdb=None			# dbid where paramdef originated
		
		if isinstance(name,dict):
			self.update(name)


	#################################		
	# repr methods
	#################################			

	def __str__(self):
		return format_string_obj(self.__dict__,["name","vartype","desc_short","desc_long","property","defaultunits","","creator","creationtime","creationdb"])


	#################################		
	# mapping methods
	#################################			

	def __getitem__(self,key):
		try:
			return self.__dict__[key]
		except:
			pass
			#print "no key %s"%key
		
	def __setitem__(self,key,value):
		if key in self.attr_all:
			self.__dict__[key]=value
		else:
			raise KeyError("Invalid key: %s"%key)
			
	def __delitem__(self,key):
		raise KeyError("Key deletion not allowed")
		
	def keys(self):
		return tuple(self.attr_all)
		
		
	#################################		
	# ParamDef methods
	#################################				

	def items_dict(self):
		ret = {}
		for k in self.attr_all:
			ret[k]=self.__dict__[k]
		return ret		


	#################################		
	# validation methods
	#################################				
	
	def validate(self):
		
		if set(self.__dict__.keys())-self.attr_all:
			raise AttributeError("Invalid attributes: %s"%",".join(set(self.__dict__.keys())-self.attr_all))
		
		if str(self.name) == "":
			raise ValueError("name required")

		if self.vartype != None and not str(self.vartype) in valid_vartypes:
			raise ValueError("Invalid vartype; not in valid_vartypes")
			
		if str(self.desc_short) == "":
			raise ValueError("Short description (desc_short) required")

		if str(self.desc_long) == "":
			raise ValueError("Long description (desc_long) required")

		if self.property != None and str(self.property) not in valid_properties:
			#raise ValueError,"Invalid property; not in valid_properties"
			print("Warning: Invalid property; not in valid_properties")

		if self.defaultunits != None:
			a=[]
			for q in list(valid_properties.values()):
				# ian
				a.extend(list(q[1].keys())) 
				a.extend(list(q[2].keys()))
			if not str(self.defaultunits) in a and str(self.defaultunits) != '':
				#raise ValueError,"Invalid defaultunits; not in valid_properties"
				print("Warning: Invalid defaultunits; not in valid_properties")
			
		if self.choices != None:
			try:
				list(self.choices)
				for i in self.choices:
					str(i)
			except:
				raise ValueError("choices must be a list of strings")
			#if isinstance(self.choices,basestring):
			# raise ValueError,"choices must be strings"
			#for i in self.choices:
