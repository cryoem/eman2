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

from UserDict import DictMixin
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
		raise Exception,"Invalid boolean %s"%str(x)

# validation/conversion for text, permitting unicode
def textconv(x):
	try: return str(x)
	except: return unicode(x)

def datetimeconv(x):
	return str(x)
	
def timeconv(x):
	return str(x)

def dateconv(x):
	return str(x)

DEBUG=0

def tojson(o):
	if isinstance(o,(Database.Record,Database.ParamDef,Database.User,Database.RecordDef)):
		return dict(o)
	return o

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
	"intlist":(None,lambda y:map(lambda x:int(x),y)),		# list of integers
	"floatlist":(None,lambda y:map(lambda x:float(x),y)), # list of floats
	"stringlist":(None,lambda y:map(lambda x:str(x),y)),	# list of enumerated strings
	"url":("s",lambda x:str(x)),			# link to a generic url
	"hdf":("s",lambda x:str(x)),			# url points to an HDF file
	"image":("s",lambda x:str(x)),			# url points to a browser-compatible image
	"binary":("s",lambda y:map(lambda x:str(x),y)),				# url points to an arbitrary binary... ['bdo:....','bdo:....','bdo:....']
	"binaryimage":("s",lambda x:str(x)),		# non browser-compatible image requiring extra 'help' to display... 'bdo:....'
	"child":("child",lambda y:map(lambda x:int(x),y)),	# link to dbid/recid of a child record
	"link":("link",lambda y:map(lambda x:int(x),y)),		# lateral link to related record dbid/recid
	"boolean":("d",boolconv),
	"dict":(None, lambda x:x), 
	"user":("s",lambda x:str(x)), # ian 09.06.07
	"userlist":(None,lambda y:map(lambda x:str(x),y))
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
	{"degree":1.0,"radian":180.0/pi, "mrad":0.18/pi},
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
			raise KeyError,"Invalid key: %s"%key
			
	def __delitem__(self,key):
		raise KeyError,"Key deletion not allowed"
		
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
			raise AttributeError,"Invalid attributes: %s"%",".join(set(self.__dict__.keys())-self.attr_all)
		
		if str(self.name) == "":
			raise ValueError,"name required"

		if self.vartype != None and not str(self.vartype) in valid_vartypes:
			raise ValueError,"Invalid vartype; not in valid_vartypes"
			
		if unicode(self.desc_short) == "":
			raise ValueError,"Short description (desc_short) required"

		if unicode(self.desc_long) == "":
			raise ValueError,"Long description (desc_long) required"

		if self.property != None and str(self.property) not in valid_properties:
			#raise ValueError,"Invalid property; not in valid_properties"
			print "Warning: Invalid property; not in valid_properties"

		if self.defaultunits != None:
			a=[]
			for q in valid_properties.values():
				# ian
				a.extend(q[1].keys()) 
				a.extend(q[2].keys())
			if not str(self.defaultunits) in a and str(self.defaultunits) != '':
				#raise ValueError,"Invalid defaultunits; not in valid_properties"
				print "Warning: Invalid defaultunits; not in valid_properties"
			
		if self.choices != None:
			try:
				list(self.choices)
				for i in self.choices:
					unicode(i)
			except:
				raise ValueError,"choices must be a list of strings"
			#if isinstance(self.choices,basestring):
			# raise ValueError,"choices must be strings"
			#for i in self.choices:

def parseparmvalues(text,noempty=0):
	"""This will extract parameter names $param or $param=value """
	# Ed 05/17/2008 -- cleaned up
	#srch=re.findall('<([^> ]*) ([^=]*)="([^"]*)" *>([^<]*)</([^>]*)>' ,text)
	#srch=re.findall('\$\$([^\$\d\s<>=,;-]*)(?:(?:=)(?:(?:"([^"]*)")|([^ <>"]*)))?',text)
	srch=re.finditer('\$\$([a-zA-Z0-9_\-]*)(?:(?:=)(?:(?:"([^"]*)")|([^ <>"]*)))?',text)
	params, vals = ret=[[],{}]
	
	for name, a, b in (x.groups() for x in srch):
		if name is '': continue
		else:
			params.append(name)
			if a is None: val=b
			else: val=a
			vals[name] = val
	return ret


class RecordDef(DictMixin) :
	"""This class defines a prototype for Database Records. Each Record is a member of
	a RecordClass. This class contains the information giving meaning to the data Fields
	contained by the Record"""
	
	attr_user = set(["mainview","views","private","typicalchld"])
	attr_admin = set(["name","params","paramsK","owner","creator","creationtime","creationdb"])
	attr_all = attr_user | attr_admin
	
	def __init__(self,d=None):
		self.name=None				# the name of the current RecordDef, somewhat redundant, since also stored as key for index in Database
		self.views={"recname":"$$rectype $$creator $$creationtime"}				# Dictionary of additional (named) views for the record
		self.mainview="$$rectype $$creator $$creationtime"			# a string defining the experiment with embedded params
									# this is the primary definition of the contents of the record
		self.private=0				# if this is 1, this RecordDef may only be retrieved by its owner (which may be a group)
									# or by someone with read access to a record of this type
		self.typicalchld=[]			# A list of RecordDef names of typical child records for this RecordDef
									# implicitly includes subclasses of the referenced types

		self.params={}				# A dictionary keyed by the names of all params used in any of the views
									# values are the default value for the field.
									# this represents all params that must be defined to have a complete
									# representation of the record. Note, however, that such completeness
									# is NOT REQUIRED to have a valid Record
		self.paramsK = []			# ordered keys from params()
		self.owner=None				# The owner of this record
		self.creator=0				# original creator of the record
		self.creationtime=None		# creation date
		self.creationdb=None		# dbid where recorddef originated
		
		if (d):
			self.update(d)

		self.findparams()
					

	def __setattr__(self,key,value):
		self.__dict__[key] = value
		if key == "mainview": self.findparams()			


	#################################		
	# pickle methods
	#################################
	
	def __setstate__(self,dict):
		"""restore unpickled values to defaults after unpickling"""
		self.__dict__.update(dict)
		if not dict.has_key("typicalchld") : self.typicalchld=[]		


	#################################		
	# repr methods
	#################################

	def __str__(self):
		return "{ name: %s\nmainview:\n%s\nviews: %s\nparams: %s\nprivate: %s\ntypicalchld: %s\nowner: %s\ncreator: %s\ncreationtime: %s\ncreationdb: %s}\n"%(
			self.name,self.mainview,self.views,self.stringparams(),str(self.private),str(self.typicalchld),self.owner,self.creator,self.creationtime,self.creationdb)

	def stringparams(self):
		"""returns the params for this recorddef as an indented printable string"""
		r=["{"]
		for k,v in self.params.items():
			r.append("\n\t%s: %s"%(k,str(v)))
		return "".join(r)+" }\n"
	

	#################################		
	# mapping methods
	#################################
			
	def __getitem__(self,key):
		return self.__dict__[key]
		
	def __setitem__(self,key,value):
		if key in self.attr_all:
			self.__dict__[key]=value
		else:
			raise AttributeError,"Invalid key: %s"%key
			
	def __delitem__(self,key):
		raise AttributeError,"Key deletion not allowed"
		
	def keys(self):
		return tuple(self.attr_all)


	#################################		
	# RecordDef methods
	#################################

	def findparams(self):
		"""This will update the list of params by parsing the views"""
		t,d=parseparmvalues(self.mainview)
		for i in self.views.values():
			t2,d2=parseparmvalues(i)
			for j in t2:
				# ian: fix for: empty default value in a view unsets default value specified in mainview
				if not d.has_key(j):
					t.append(j)
					d[j] = d2[j]
#			d.update(d2)

		self.params=d
		self.paramsK=tuple(t)


	def items_dict(self):
		ret = {}
		for k in self.attr_all:
			ret[k]=self.__dict__[k]
		return ret
		
		
	def fromdict(self,d):
		for k,v in d.items():
			self.__setattr__(k,v)
		self.validate()

		
	#################################		
	# validate methods
	#################################		
		
	def validate(self): 
		
		if set(self.__dict__.keys())-self.attr_all:
			raise AttributeError,"Invalid attributes: %s"%",".join(set(self.__dict__.keys())-self.attr_all)
		
		try:
			if str(self.name) == "" or self.name==None:
				raise Exception
		except: 
			raise ValueError,"name required; must be str or unicode"

		try:
			dict(self.views)
		except:
			raise ValueError,"views must be dict"

		try:
			list(self.typicalchld)
		except:
			raise ValueError,"Invalid value for typicalchld; list of recorddefs required."
						
		try: 
			if unicode(self.mainview) == "": raise Exception
		except:
			raise ValueError,"mainview required; must be str or unicode"

		if not dict(self.views).has_key("recname"):
			print "Warning: recname view strongly suggested"

		for k,v in self.views.items():
			if not isinstance(k,str) or not isinstance(v,basestring):
				raise ValueError,"Views names must be strings; view defs may be unicode"

		if self.private not in [0,1]:
			raise ValueError,"Invalid value for private; must be 0 or 1"


class Record(DictMixin):
	"""This class encapsulates a single database record. In a sense this is an instance
	of a particular RecordDef, however, note that it is not required to have a value for
	every field described in the RecordDef, though this will usually be the case.
	
	To modify the params in a record use the normal obj[key]= or update() approaches. 
	Changes are not stored in the database until commit() is called. To examine params, 
	use obj[key]. There are a few special keys, handled differently:
	creator,creationtime,permissions,comments

	Record instances must ONLY be created by the Database class through retrieval or
	creation operations. self.context will store information about security and
	storage for the record.
	
	Mechanisms for changing existing params are a bit complicated. In a sense, as in a 
	physical lab notebook, an original value can never be changed, only superceded. 
	All records have a 'magic' field called 'comments', which is an extensible array
	of text blocks with immutable entries. 'comments' entries can contain new field
	definitions, which will supercede the original definition as well as any previous
	comments. Changing a field will result in a new comment being automatically generated
	describing and logging the value change.
	
	From a database standpoint, this is rather odd behavior. Such tasks would generally be
	handled with an audit log of some sort. However, in this case, as an electronic
	representation of a Scientific lab notebook, it is absolutely necessary
	that all historical values are permanently preserved for any field, and there is no
	particular reason to store this information in a separate file. Generally speaking,
	such changes should be infrequent.
	
	Naturally, as with anything in Python, anyone with code-level access to the database
	can override this behavior by changing 'params' directly rather than using
	the supplied access methods. There may be appropriate uses for this when constructing
	a new Record before committing changes back to the database.
	"""
	
	attr_user = set([])
	attr_admin = set(["recid","dbid","rectype"])
	attr_private = set(["_Record__params","_Record__comments","_Record__oparams",
								 "_Record__creator","_Record__creationtime","_Record__permissions",
								 "_Record__ptest","_Record__context"])
	attr_restricted = attr_private | attr_admin
	attr_all = attr_user | attr_admin | attr_private
	
	param_special = set(["recid","rectype","comments","creator","creationtime",
						 "permissions"]) # dbid # "modifyuser","modifytime",
	
	#publicparams = property(lambda self: set(self.__params) - set(self.restrictedparams))
		
	def __init__(self,d=None,ctx=None, **kwargs):
		"""Normally the record is created with no parameters, then setContext is called by the
		Database object. However, for initializing from a dictionary (ie - XMLRPC call, this
		may be done at initiailization time."""
		kwargs.update(d or {})
		## recognized keys
		# recid -- 32 bit integer recordid (within the current database)
		# rectype -- name of the RecordDef represented by this Record
		# comments -- a List of comments records
		# creator -- original creator of the record
		# creationtime -- creation date
		# dbid -- dbid where this record resides (any other dbs have clones)
		# params -- a Dictionary containing field names associated with their data
		# permissions -- permissions for 
						# read access, comment write access, full write access, 
						# and administrative access. Each element is a tuple of 
						# user names or group id's,
						# Group -3 includes any logged in user,
						# Group -4 includes any user (anonymous)
		
		self.recid=kwargs.get('recid')
		self.rectype=kwargs.get('rectype', '')				
		self.__comments=kwargs.get('comments',[])			
		self.__creator=kwargs.get('creator',0)
		self.__creationtime=kwargs.get('creationtime')
		self.__permissions=kwargs.get('permissions',((),(),(),()))

		self.dbid=kwargs.get('dbid',None)
		self.__params={}#kwargs.get('params',{}).copy()
		for key in set(kwargs.keys())-self.param_special:
			key=str(key).strip().lower()
			self.__params[key]=kwargs[key]

		self.__oparams={}

		self.__ptest=[1,1,1,1] 
			# ian: changed to [1,1,1,1] so you can create instances (recid=None) directly
			# ed: NOTE: should this be passed in the init dictionary? 
			# Results of security test performed when the context is set
			# correspond to, read,comment,write and owner permissions, return from setContext
			
		self.__context=None # Validated access context 
		if ctx != None: ctx = ctx	 
		else: ctx = kwargs.get('context')
		if (ctx!=None): self.setContext(ctx)

			
			
	def validate(self):
		self.validate_recid()
		self.validate_rectype()
		self.validate_permissions()

		if self.__context.db:
			self.validate_permissions_users()
			self.validate_params()
			self.validate_security()

		

	#def validate_security(self):
		#cp = self.changedparams()

		#if "recid" in cp:
			#raise ValueError,"Cannot change recid"			
		#if "creator" in cp or "creationtime" in cp:
			#raise ValueError,"Cannot change creation info directly"
		#if "rectype" in cp:
			#raise ValueError,"Cannot change rectype"
		##if "permissions" in cp:
		##	self['permissions'] = self.__oparams['permissions']
			
		#if not self.__ptest[2]:
			#if cp != ["comments"]:
				#raise SecurityError,"Insufficient permission to change field values (%d)"%record.recid
			#if not p[1]:
				#raise SecurityError,"Insufficient permission to add comments to record (%d)"%record.recid
				
					
	def validate_recid(self):
		try: 
			if self.recid != None: int(self.recid)
		except:
			raise ValueError,"recid must be positive integer"


	def validate_rectype(self):			
		if str(self.rectype) == "":
			raise ValueError,"rectype must not be empty"


	def validate_permissions(self):		
		try:
			self.__permissions = tuple((tuple(i) for i in self.__permissions))
		except:
			raise ValueError,"permissions must be 4-tuple of tuples"


	def validate_permissions_users(self):
		u=set()
		users=set(self.__context.db.getusernames(self.__context.ctxid,host=self.__context.host))		
		for j in self.__permissions: u|=set(j)
		u -= set([0,-1,-2,-3])
		if u-users:
			print "Warning: undefined users: %s"%",".join(map(str, u-users))
			
			

	def validate_param(self, value, pd):
		print "___ validate %s ___"%(pd.name)
		#print "vartype: %s"%pd.vartype
		#print "property: %s"%pd.property
		#print "defaultunits: %s"%pd.defaultunits

		if pd.property and isinstance(value,basestring):
			print "______ checking units __________"
			value,units=re.compile("([0-9.,]+)?(.*)").search(value).groups()
			value=float(value)
			units=units.strip()
			defaultunits=valid_properties[pd.property][0]

			#print self[pd.name]
			#print value
			#print units
			#print defaultunits

			if pd.defaultunits != None:
				defaultunits=pd.defaultunits

			if units in valid_properties[pd.property][1].keys():
				pass

			elif units in valid_properties[pd.property][2].keys():
				units=valid_properties[pd.property][2][units]

			elif units == "":
				units = defaultunits

			else:
				raise ValueError,"Unknown units: %s"%units

			try:
				# convert units
				value = value * ( valid_properties[pd.property][1][units] / valid_properties[pd.property][1][defaultunits] )
				print "newval: %s"%value
			except:
				raise ValueError,"Unable to convert %s = %s; skipping value"%(pd.name,value)

		try:
			if value != None:
				value=valid_vartypes[pd.vartype][1](value)
		except:
			raise ValueError,"Error converting datatype: %s, %s"%(pd.name,pd.vartype)
	
		if pd.vartype=="choice" and value!=None:
			if value.lower() not in [i.lower() for i in pd.choices]:
				raise ValueError,"%s not in %s"%(value,pd.choices)
		
		
		# Save validated value
		if self[pd.name]!=value:
			self[pd.name]=value

		return
		
		
	def validate_params(self):
		if self.__params.keys():
			pds=self.__context.db.getparamdefs(self.__params.keys())
			for i,pd in pds.items():
				self.validate_param(self[i],pd)


	def changedparams(self):

		print "===== changed params ====="

		cp = set()
		for k in self.keys():
			if k not in (self.param_special-set(["comments"])) and self[k] != self.__oparams.get(k,None):

				if k == "comments":
					print "%s : ..."%(k)
				else:
					print "%s : %s --> %s"%(k,self.__oparams.get(k,None),self[k])
					
				cp.add(k)

		print "changedparams result: %s"%cp
		return cp		
		
	#################################		
	# pickle methods
	#################################
	
	def __getstate__(self):
		"""the context and other session-specific information should not be pickled"""
		odict = self.__dict__.copy() # copy the dict since we change it
		odict["_Record__oparams"]={}
		
		if not odict.has_key("localcpy") :
			try: del odict['_Record__ptest']
			except: pass
		
		try: del odict['_Record__context']
		except: pass

		return odict
	
	def __setstate__(self,dict):
		"""restore unpickled values to defaults after unpickling"""	
		#print "unpickle: _Record__oparams"
		#if dict["_Record__oparams"]	!= {}:
		#	print dict["recid"],dict["_Record__oparams"]
			
		# this is properly handled by putrecord	
		# try:
		# 	p=dict["_Record__params"]
		#		dict["_Record__params"]={}
		#		for i,j in p.items(): 
		#			if j!=None and j!="None" : dict["_Record__params"][i.lower()]=j
		#except:
		#	traceback.print_exc(file=sys.stdout)
		#dict["rectype"]=dict["rectype"].lower()
		
		#if dict.has_key("localcpy") :
		#	del dict["localcpy"]
		self.__dict__.update(dict)	
		#	self.__ptest=[1,1,1,1]
		#else:
		#	self.__dict__.update(dict)	
		#	self.__ptest=[0,0,0,0]
		self.__ptest=[0,0,0,0]
		self.__context=None

	
	
	#################################		
	# repr methods
	#################################
	
	def __unicode__(self):
		"A string representation of the record"
		ret=["%s (%s)\n"%(str(self.recid),self.rectype)]
		for i,j in self.items(): 
			ret.append(u"%12s:	%s\n"%(str(i),unicode(j)))
		return u"".join(ret)
	
	def __str__(self):
		return self.__unicode__().encode('utf-8')

	def __repr__(self):
		return "<Record id: %s recdef: %s at %x>" % (self.recid, self.rectype, id(self))
		
		
	def json_equivalent(self):
		return self.items_dict()

	#################################		
	# mapping methods; 
	#		DictMixin provides the remainder
	#################################
	
	def __getitem__(self,key):
		"""Behavior is to return None for undefined params, None is also
		the default value for existant, but undefined params, which will be
		treated identically"""
		#if not self.__ptest[0] : raise SecurityError,"Permission Denied (%d)"%self.recid
		# key = key.encode('utf-8', 'replace').strip().lower()
		key = str(key).strip().lower()
		result = None
		if   key=="comments" : result = self.__comments
		elif key=="recid" : result = self.recid
		elif key=="rectype" : result = self.rectype
		elif key=="creator" : result = self.__creator
		elif key=="creationtime" : result = self.__creationtime
		elif key=="permissions" : result = self.__permissions
		else: result = self.__params.get(key)
		return result

	def __setitem__(self,key,value):
		"""This and 'update' are the primary mechanisms for modifying the params in a record
		Changes are not written to the database until the commit() method is called!"""
		# comments may include embedded field values if the user has full write access

		key = str(key).strip().lower()

		#try: valuelower = value.lower()
		#except:	valuelower = ""
		
		if value == "None": 
			value = None
		try:
			if len(value) == 0:
				value = None
		except:
			pass

		if key != "comments" and not self.__oparams.has_key(key):
			self.__oparams[key]=self[key]

		if key not in self.param_special:
			self.__params[key]=value

		elif key=="comments":
			self.addcomment(value)
		
		elif key == 'permissions':
			self.__permissions = tuple(tuple(x) for x in value)
			#self.__permissions = tuple([ tuple(set(x) | set(y))	 for (x,y) in zip(value, self.__permissions)])


	def __delitem__(self,key):
		
		#if not self.__ptest[1] : raise SecurityError,"Permission Denied (%d)"%self.recid
		key = str(key).strip().lower()

		if key not in self.param_special and self.__params.has_key(key):
			self.__setitem__(key,None)
		else:
			raise KeyError,"Cannot delete key %s"%key 

	def keys(self):
		"""All retrievable keys for this record"""
		#if not self.__ptest[0] : raise SecurityError,"Permission Denied (%d)"%self.recid		
		return tuple(self.__params.keys())+tuple(self.param_special)
		
	def has_key(self,key):
		if str(key).strip().lower() in self.keys(): return True
		return False
	
	def get(self, key, default=None):
		return DictMixin.get(self, key) or default

	#################################		
	# record methods
	#################################

	def addcomment(self, value):
		if not isinstance(value,basestring):
			print "Warning: invalid comment"
			return

		#print "Updating comments: %s"%value
		d=parseparmvalues(value,noempty=1)[1]

		if d.has_key("comments"):
			print "Warning: cannot set comments inside a comment"			
			return
			
		self.__comments.append((self.__context.user,time.strftime("%Y/%m/%d %H:%M:%S"),value))	# store the comment string itself		

		# now update the values of any embedded params
		for i,j in d.items():
			self.__setitem__(i,j)


	def getparamkeys(self):
		"""Returns parameter keys without special values like owner, creator, etc."""
		return self.__params.keys()		
		

	# ian: hard coded groups here need to be in sync with the various check*ctx methods.
	def setContext(self,ctx):
		"""This method may ONLY be used directly by the Database class. Constructing your
		own context will not work to see if a ctx(a user context) has the permission to access/write to this record
		"""
		#self.__context__=ctx
		self.__context = ctx
		
		if self.__creator==0:
			self.__creator=ctx.user
			self.__creationtime=time.strftime("%Y/%m/%d %H:%M:%S")
			self.__permissions=((),(),(),(ctx.user,))
		
		# test for owner access in this context.
		if (-1 in ctx.groups) : self.__ptest=[1,1,1,1]
		else:
			# we use the sets module to do intersections in group membership
			# note that an empty set tests false, so u1&p1 will be false if
			# there is no intersection between the 2 sets
			p1=set(self.__permissions[0]+self.__permissions[1]+self.__permissions[2]+self.__permissions[3])
			p2=set(self.__permissions[1]+self.__permissions[2]+self.__permissions[3])
			p3=set(self.__permissions[2]+self.__permissions[3])
			p4=set(self.__permissions[3])
			u1=set(ctx.groups+[-4])				# all users are permitted group -4 access
			
			if ctx.user!=None : u1.add(-3)		# all logged in users are permitted group -3 access
			
			# test for read permission in this context
			if (-2 in u1 or ctx.user in p1 or u1&p1) : self.__ptest[0]=1
	
			# test for comment write permission in this context
			if (ctx.user in p2 or u1&p2): self.__ptest[1]=1
						
			# test for general write permission in this context
			if (ctx.user in p3 or u1&p3) : self.__ptest[2]=1
			
			# test for administrative permission in this context
			if (ctx.user in p4 or u1&p4) : self.__ptest[3]=1
		
		return self.__ptest
		
		
	# from items_dict 
	
	def fromdict(self,d):
		if d.has_key("comments"):
			self.__comments=d["comments"]
			del d["comments"]
		for k,v in d.items():
			#print "setting %s = %s"%(k,v)
			self[k]=v
		# generally there will be no context set at this point; validation is short form
		# ed: no validation here!!! 
		# self.validate()		
	
	def items_dict(self):
		"""Returns a dictionary of current values, __dict__ wouldn't return the correct information"""
		#if not self.__ptest[0] : raise SecurityError,"Permission Denied (%d)"%self.recid		
		ret={}
		ret.update(self.__params)
		for i in self.param_special:
			try:
				ret[i]=self[i]
			except:
				pass
		return ret

	def commit(self,host=None):
		"""This will commit any changes back to permanent storage in the database, until
		this is called, all changes are temporary. host must match the context host or the
		putrecord will fail"""
		return self.__context.db.putrecord(self,self.__context.ctxid,host=host)

	def setoparams(self,d):
		self.__oparams=d.copy()

	def isowner(self):
		return self.__ptest[3]	

	def writable(self):
		"""Returns whether this record can be written using the given context"""
		return self.__ptest[2]
		

	def commentable(self):
		"""Does user have level 1 permissions? Required to comment or link."""
		return self.__ptest[1]	

