#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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
import traceback
#traceback.print_stack()

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
import sys
from math import *
from sys import exit
import os
import time
import shelve
import re
import pickle
import zlib
import socket
import subprocess
from EMAN2_cppwrap import *
from EMAN2_meta import *
#import EMAN2db
import EMAN2jsondb
import argparse, copy
import glob
import random
from struct import pack,unpack
import json
from collections import OrderedDict
import traceback
from pathlib import Path

import threading
#from Sparx import *

### If we ever need to add 'cleanup' exit code, this is how to do it. Drawn from the old BDB code.
## if the program exits nicely, close all of the databases
#atexit.register(DB_cleanup)

## if we are killed 'nicely', also clean up (assuming someone else doesn't grab this signal)
#signal.signal(2, DB_cleanup)
#signal.signal(15, DB_cleanup)

os.environ['QT_MAC_WANTS_LAYER'] = '1'


def e2gethome():
	"""platform independent path with '/'"""
	if (sys.platform != 'win32'):
		url = os.getenv("HOME")
	else:
		if (os.getenv("HOMEPATH") == '\\'):
			# Contributed by Alexander Heyne <AHEYNE@fmp-berlin.de>
			url = os.getenv("USERPROFILE")
			url = url.lstrip('CDEFG:')  # could also use substr
			url = url.replace("\\", "/")
		else:
			url = os.getenv("HOMEPATH")
			url = url.replace("\\", "/")
	return url


def e2getcwd():
	"""platform independent path with '/'"""
	url = os.getcwd()
	url = url.replace("\\", "/")
	return url

HOMEDB=None

# This next line is to initialize the Transform object for threadsafety. Utterly stupid approach, but a functional hack
T=Transform({"type":"2d","alpha":0})

# When generating bispectral invariants, we need 2 parameters, which must be used consistently throughout the system
bispec_invar_parm=(32,10)

# These are processors which don't support in-place operation
outplaceprocs=["math.bispectrum.slice","math.harmonic","misc.directional_sum"]

# Without this, in many countries Qt will set things so "," is used as a decimal
# separator by sscanf and other functions, which breaks CTF reading and some other things
try:
	os.putenv("LC_CTYPE","en_US.UTF-8")
	os.putenv("LC_ALL","en_US.UTF-8")
except: pass

# This block attempts to open the standard EMAN2 database interface
# if it fails, it sets db to None. Applications can then alter their
# behavior appropriately
#try:
#import EMAN2db
#from EMAN2db import EMAN2DB,db_open_dict,db_close_dict,db_remove_dict,db_list_dicts,db_check_dict,db_parse_path,db_convert_path,db_get_image_info,e2gethome, e2getcwd
from EMAN2jsondb import JSDict,js_open_dict,js_close_dict,js_remove_dict,js_list_dicts,js_check_dict,js_one_key
#except:
#	HOMEDB=None

XYData.__len__=XYData.get_size

# Who is using this? Transform3D is deprecated use the Transform instead
#Transform3D.__str__=lambda x:"Transform3D(\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g\n\t\t%7.4g\t%7.4g\t%7.4g)\nPretrans:%s\nPosttrans:%s"%(x.at(0,0),x.at(0,1),x.at(0,2),x.at(1,0),x.at(1,1),x.at(1,2),x.at(2,0),x.at(2,1),x.at(2,2),str(x.get_pretrans()),str(x.get_posttrans()))

try:
	if __IPYTHON__ : GUIMode=True
	from PyQt5 import QtGui, QtWidgets
	app=QtWidgets.qApp
except:
	GUIMode=False
	app = 0

GUIbeingdragged=None
originalstdout = sys.stdout

# Aliases
EMData.get=EMData.get_value_at
EMData.set=EMData.set_value_at
EMData.numpy=EMNumPy.em2numpy
from_numpy=EMNumPy.numpy2em
to_numpy=EMNumPy.em2numpy


def emdata_to_string(self):
	"""This returns a compressed string representation of the EMData object, suitable for storage
	or network communication. The EMData object is pickled, then compressed wth zlib. Restore with
	static method from_string()."""
	return zlib.compress(pickle.dumps(self,-1),3)	# we use a lower compression mode for speed

def emdata_from_string(s):
	"""This will restore a serialized compressed EMData object as prepared by as_string()"""
	return pickle.loads(zlib.decompress(s))

EMData.from_string=emdata_from_string
EMData.to_string=emdata_to_string

def list_to_emdata(l):
	"""Converts a 1-D list into a 1-D EMData object (slow)"""
	r=EMData(len(l),1,1)
	for i,j in enumerate(l): r[i]=j
	return r

def timer(fn,n=1):
	a=time.time()
	for i in range(n): fn()
	print(time.time()-a)

# This is to remove stdio buffering, only line buffering is done. This is what is done for the terminal, but this extends terminal behaviour to redirected stdio
# try/except is to prevent errors with systems that already redirect stdio
#try: sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
#except: pass

def stopautoflush():
	""" Return to buffered stdout """
	sys.stdout = originalstdout

# These are very widely used and hard to find, so some shortcuts
# Image file types
IMAGE_MRC = EMUtil.ImageType.IMAGE_MRC
IMAGE_EER = EMUtil.ImageType.IMAGE_EER
IMAGE_EER2X = EMUtil.ImageType.IMAGE_EER2X
IMAGE_EER4X = EMUtil.ImageType.IMAGE_EER4X
IMAGE_SPIDER = EMUtil.ImageType.IMAGE_SPIDER
IMAGE_SINGLE_SPIDER = EMUtil.ImageType.IMAGE_SINGLE_SPIDER
IMAGE_IMAGIC = EMUtil.ImageType.IMAGE_IMAGIC
IMAGE_HDF = EMUtil.ImageType.IMAGE_HDF
IMAGE_DM3 = EMUtil.ImageType.IMAGE_DM3
IMAGE_DM4 = EMUtil.ImageType.IMAGE_DM4
IMAGE_TIFF = EMUtil.ImageType.IMAGE_TIFF
IMAGE_PGM = EMUtil.ImageType.IMAGE_PGM
IMAGE_LST = EMUtil.ImageType.IMAGE_LST
IMAGE_PIF = EMUtil.ImageType.IMAGE_PIF
IMAGE_VTK = EMUtil.ImageType.IMAGE_VTK
IMAGE_PNG = EMUtil.ImageType.IMAGE_PNG
IMAGE_SAL = EMUtil.ImageType.IMAGE_SAL
IMAGE_ICOS = EMUtil.ImageType.IMAGE_ICOS
IMAGE_EMIM = EMUtil.ImageType.IMAGE_EMIM
IMAGE_GATAN2 = EMUtil.ImageType.IMAGE_GATAN2
IMAGE_AMIRA = EMUtil.ImageType.IMAGE_AMIRA
IMAGE_XPLOR = EMUtil.ImageType.IMAGE_XPLOR
IMAGE_EM = EMUtil.ImageType.IMAGE_EM
IMAGE_V4L = EMUtil.ImageType.IMAGE_V4L
IMAGE_UNKNOWN = EMUtil.ImageType.IMAGE_UNKNOWN

# image data storage modes
EM_UNKNOWN = EMUtil.EMDataType.EM_UNKNOWN
EM_CHAR = EMUtil.EMDataType.EM_CHAR
EM_UCHAR = EMUtil.EMDataType.EM_UCHAR
EM_SHORT = EMUtil.EMDataType.EM_SHORT
EM_USHORT = EMUtil.EMDataType.EM_USHORT
EM_INT = EMUtil.EMDataType.EM_INT
EM_UINT = EMUtil.EMDataType.EM_UINT
EM_FLOAT = EMUtil.EMDataType.EM_FLOAT
EM_DOUBLE = EMUtil.EMDataType.EM_DOUBLE
EM_SHORT_COMPLEX = EMUtil.EMDataType.EM_SHORT_COMPLEX
EM_USHORT_COMPLEX = EMUtil.EMDataType.EM_USHORT_COMPLEX
EM_FLOAT_COMPLEX = EMUtil.EMDataType.EM_FLOAT_COMPLEX
EM_COMPRESSED = EMUtil.EMDataType.EM_COMPRESSED


# These map standard names for data types to internal representation, and provide a minimum and maximum value for each type
file_mode_map={
	"int8"  :EMUtil.EMDataType.EM_CHAR,
	"uint8" :EMUtil.EMDataType.EM_UCHAR,
	"int16" :EMUtil.EMDataType.EM_SHORT,
	"uint16":EMUtil.EMDataType.EM_USHORT,
	"int32" :EMUtil.EMDataType.EM_INT,
	"uint32":EMUtil.EMDataType.EM_UINT,
	"float" :EMUtil.EMDataType.EM_FLOAT,
	"compressed": EMUtil.EMDataType.EM_COMPRESSED }

# inverse dictionary for getting printable names
file_mode_imap=dict([[int(v),k] for k,v in list(file_mode_map.items())])

file_mode_intmap={
	1 :EMUtil.EMDataType.EM_CHAR,
	2 :EMUtil.EMDataType.EM_UCHAR,
	3 :EMUtil.EMDataType.EM_SHORT,
	4 :EMUtil.EMDataType.EM_USHORT,
	5 :EMUtil.EMDataType.EM_INT,
	6 :EMUtil.EMDataType.EM_UINT,
	7 :EMUtil.EMDataType.EM_FLOAT,
	8 :EMUtil.EMDataType.EM_COMPRESSED }


#keyed both by type and by the integer version for flexibility
file_mode_range={
	EMUtil.EMDataType.EM_CHAR:(-128,127),
	EMUtil.EMDataType.EM_UCHAR:(0,255),
	EMUtil.EMDataType.EM_SHORT:(-32768,32767 ),
	EMUtil.EMDataType.EM_USHORT:(0,65535 ),
	EMUtil.EMDataType.EM_INT:(-2147483648,2147483647 ),
	EMUtil.EMDataType.EM_UINT:(0,4294967295),
	EMUtil.EMDataType.EM_FLOAT:(-3.40282347e+38,3.40282347e+38 ),
	int(EMUtil.EMDataType.EM_CHAR):(-128,127),
	int(EMUtil.EMDataType.EM_UCHAR):(0,255),
	int(EMUtil.EMDataType.EM_SHORT):(-32768,32767 ),
	int(EMUtil.EMDataType.EM_USHORT):(0,65535 ),
	int(EMUtil.EMDataType.EM_INT):(-2147483648,2147483647 ),
	int(EMUtil.EMDataType.EM_UINT):(0,4294967295),
	int(EMUtil.EMDataType.EM_FLOAT):(-3.40282347e+38,3.40282347e+38 ),
	int(EMUtil.EMDataType.EM_COMPRESSED):(-3.40282347e+38,3.40282347e+38 )
	}

class NotImplementedException(Exception):
	def __init__(self,val=None): pass

def E2init(argv, ppid=-1) :
	"""E2init(argv)
This function is called to log information about the current job to the local logfile. The flags stored for each process
are pid, start, args, progress and end. progress is from 0.0-1.0 and may or may not be updated. end is not set until the process
is complete. If the process is killed, 'end' may never be set."""

	# We go to the end of the file. Record the location, then write a fixed length string
	try:
		hist=open(".eman2log.txt","r+")
		hist.seek(0,os.SEEK_END)
	except:
		try: hist=open(".eman2log.txt","w")
		except: return -1
	n=hist.tell()
	hist.write("%s\tincomplete         \t%6d/%6d\t%s\t%s\n"%(local_datetime(),os.getpid(),ppid,socket.gethostname()," ".join(argv)))

#	hist.flush()
	hist.close()

	#if EMAN2db.BDB_CACHE_DISABLE :
		#print "Note: Cache disabled"

	return n

def E2progress(n,progress):
	"""Updates the progress fraction (0.0-1.0) for a running job. Negative values may optionally be
set to indicate an error exit."""
#	if EMAN2db.BDB_CACHE_DISABLE : return		# THIS MUST REMAIN DISABLED NOW THAT THE CACHE IS DISABLED PERMANENTLY !!!

	try:
		hist=open(".eman2log.txt","r+")
		hist.seek(n+20)
	except:
		return -1
	hist.write("%3d%%         "%(int(progress*100.0)))
#	hist.flush()
	hist.close()

	return n

def E2end(n):
	"""E2end(n)
This function is called to log the end of the current job. n is returned by E2init"""
#	if EMAN2db.BDB_CACHE_DISABLE : return		# THIS MUST REMAIN DISABLED NOW THAT THE CACHE IS DISABLED PERMANENTLY !!!

	try:
		hist=open(".eman2log.txt","r+")
		hist.seek(n+20)
	except:
		return -1
	hist.write("%s"%(local_datetime()))
#	hist.flush()
	hist.close()

	return n

def E2saveappwin(app,key,win):
	"""stores the window geometry using the application default mechanism for later restoration. Note that
	this will only work with Qt windows"""
	try:
		pos=win.pos()
		sz=win.size()
		geom=(pos.x(),pos.y(),win.width(),win.height())

		E2setappval(app,key,geom)
#		print(app,key,geom)
	except:
		print("Error saving window location ",key)

def E2loadappwin(app,key,win):
	"""restores a geometry saved with E2saveappwin"""
	try:
		geom=list(E2getappval(app,key))
		if geom==None : raise Exception
		win.resize(geom[2],geom[3])
		geom[0]=max(32,geom[0])
		geom[1]=max(60,geom[1])
		win.move(geom[0],geom[1])
#		print(app,key,geom)
	except: return

def E2saveprojtype(app,key,win):
	"""stores the project type using the application default mechanism for later restoration. Note that
	this will only work with Qt windows"""
	try: E2setappval(app,key,win.modeCB.currentIndex())
	except: print("Error saving project type")

def E2loadprojtype(app,key,win):
	"""restores project type saved with E2saveappwin"""
	try:
		idx=E2getappval(app,key)
		win.modeCB.setCurrentIndex(idx)
		win._onModeChange(idx)
	except: return

def E2setappval(app,key,value):
	"""E2setappval
This function will set an application default value both in the local directory and ~/.eman2
When settings are read, the local value is checked first, then if necessary, the global value."""
	try:
		app.replace(".","_")
		key.replace(".","_")
	except:
		print("Error with E2setappval, app and key must be strings")
		return

	try:
		db=js_open_dict(".eman2settings.json")
		db[app+"."+key]=value
#		db.close()
	except:
		pass

	try:
		dir=e2gethome()
		dir+="/.eman2"
	except:
		return
	
	try: os.mkdir(dir)
	except: pass

	try:
		db=js_open_dict(dir+"/eman2settings.json")
		db[app+"."+key]=value
#		db.close()
	except:
		return


def E2getappval(app,key,dfl=None):
	"""E2getappval
This function will get an application default by first checking the local directory, followed by
~/.eman2"""
	try:
		app.replace(".","_")
		key.replace(".","_")
	except:
		print("Error with E2getappval, app and key must be strings")
		return None

	try:
		db=js_open_dict(".eman2settings.json")
		ret=db[app+"."+key]
		return ret
	except:
		pass

	try:
		dir=e2gethome()
		dir+="/.eman2"
		db=js_open_dict(dir+"/eman2settings.json")
		ret=db[app+"."+key]
		db.close()

		return ret
	except: pass

	if dfl!=None: E2setappval(app,key,dfl)		# so the user will know what may be available

	return dfl

def E2getappvals():
	"""E2getappvals
This function will return a list of lists containing all currently set application defaults as [program,option,value,global|local]"""

	ret=[]
	ret2=[]

	try:
		db=js_open_dict(".eman2settings.json")
		keys=db.keys()
		ret=[(k.split(".")[0],k.split(".")[1],db[k],"local") for k in keys]
	except:
		pass
	
	try:
		dir=e2gethome()
		dir+="/.eman2"
		dbu=js_open_dict(dir+"/eman2settings.json")
		ret2=[(k.split(".")[0],k.split(".")[1],dbu[k],"user") for k in dbu.keys() if k not in keys]  # only show global when local doesn't exist
		
	except: pass

	return ret2+ret

def e2getinstalldir() :
	"""Final path needs to be computed relative to a path within the installation.
	 An alternative could be to get the installation directory from cmake,
	 but cmake is not run during binary installations."""
	
	this_file_dirname = os.path.dirname(__file__)
	if get_platform() != "Windows":
		rel_path = '../../../'
	else:
		rel_path = '../../Library/'
	
	return os.path.abspath(os.path.join(this_file_dirname, rel_path))

def get_temp_name():
	"""Returns a suitable name for a temporary HDF file in the current directory. Does not create or delete the file."""
	fsp=f"tmp_{random.randint(0,9999999):07d}.hdf"
	if os.path.exists(fsp) : return get_temp_name()		# risky? shouldn't really ever recurse forever...
	return fsp

def num_path_new(prefix):
	"""make a new prefix_xx (or prefix_xxx) folder and return the folder name, underscore added to prefix if not present"""
	if prefix[-1]!="_" : prefix+="_"
	
	pthns=[int(i.rsplit("_",1)[-1]) for i in os.listdir(".") if i[:len(prefix)]==prefix and i.rsplit("_",1)[-1].isdigit()]
	try: newn=max(pthns)+1
	except: newn=0
	path=f"{prefix}{newn:02d}"
	try: os.mkdir(path)
	except:
		raise Exception("Error: could not create "+path)
	
	return path

def num_path_last(prefix,create=False):
	"""find the highest numbered path starting with prefix and return it. If create is set, will create a new one if none exists. underscore added to prefix if not present"""
	if prefix[-1]!="_" : prefix+="_"
	
	pthns=[int(i.rsplit("_",1)[-1]) for i in os.listdir(".") if i[:len(prefix)]==prefix and i.rsplit("_",1)[-1].isdigit()]
	try: n=max(pthns)
	except:
		if create: 
			path=f"{prefix}00"
			try: os.mkdir(path)
			except: raise Exception("Error: could not create "+path)
			return path
		else: raise Exeception("Error: could not find paths beginning with "+prefix)
	
	path=f"{prefix}{n:02d}"
	
	return path

#def numbered_path(prefix,makenew=False):
	#"""Finds or creates folders of the form prefix_NN. If makenew is set, will create a new folder with NN one
	#larger than the largest existing path. If makenew is not set, returns the highest existing numbered path of that form."""
	#if prefix[-1]=="_" : prefix=prefix[:-1]
	#cur=[int(p.split("_")[-1]) for p in os.listdir(".") if p.rsplit("_",1)[0]==prefix and p.split("_")[-1].isdigit()]
	
	#if makenew:
		#cur.append(0)		# in case of no matches
		#path=f"{prefix}_{max(cur)+1:02d}"
		#try: os.mkdir(path)
		#except: 
			#raise Exception(f"ERROR: numbered_path() could not create {path}")
		#return path
	#if len(cur)==0 : raise Exception(f"ERROR: no paths of the form {prefix}_NN found")
	#return f"{prefix}_{max(cur):02d}"

#def get_numbered_directories(prefix,wd=e2getcwd()):
	#'''
	#Gets the numbered directories starting with prefix in the given working directory (wd)
	#A prefix example would be "refine_" or "r2d_" etc. Used originally form within the workflow context
	#'''
	#dirs, files = get_files_and_directories(wd)
	#dirs.sort()
	#l = len(prefix)
	#for i in range(len(dirs)-1,-1,-1):
		#if len(dirs[i]) != l+2:# plus two because we only check for two numbered directories
			#dirs.pop(i)
		#elif dirs[i][:l] != prefix:
			#dirs.pop(i)
		#else:
			#try: int(dirs[i][l:])
			#except: dirs.pop(i)

	## allright everything left in dirs is "refine_??" where the ?? is castable to an int, so we should be safe now

	#return dirs

def get_prefixed_directories(prefix,wd=e2getcwd()):
	'''
	gets directories starting with prefix and without any '.'
	'''
	dirs, files = get_files_and_directories(wd)
	dirs.sort()
	dirs=[i for i in dirs if i.startswith(prefix) and not "." in i]

	return dirs

def get_image_directory():
	dtag = get_dtag()
	
	return e2getinstalldir()+ dtag + "images" + dtag

def get_dtag():
#	pfrm = get_platform()
#	if pfrm == "Windows": return "\\"
#	else: return "/"
	return "/"

def get_files_and_directories(path=".",include_hidden=False):
	if path == ".": l_path = "./"
	else: l_path = path
	if len(l_path) == 0: path = get_dtag()
	elif l_path[-1] not in ["/","\\"]: l_path += get_dtag()

	dirs = []
	files = []
	try:
		entries = os.listdir(l_path)
	except: # something is wrong with the path
		#print "path failed",l_path
		return dirs,files

	for name in entries:
		if len(name) == 0: continue
		if not include_hidden:
			if name[0] == ".":
				continue

		try:
			if os.path.isdir(l_path+name):
				dirs.append(name)
			else:
				files.append(name)
		except:
			pass # something was wrong with the directory
	return dirs,files


def numbered_bdb(bdb_url):
	'''
	give something like "bdb:refine_01#class_indices (which means bdb:refine_01#class_indices_??.bdb files exist)

	will return the next available name bdb:refine_01#class_indices_?? (minus the bdb)
	'''
	print("ERROR: numbered_bdb no longer functions. Please report to developers.")
	sys.exit(1)

def compress_hdf(fsp,bits,nooutliers=False,level=1):
	"""This will take an existing HDF file and rewrite it with the specified compression
	nooutliers can be specified, but minval/maxval are determined automatically. Returns
	immediately if non-HDF filename is provided."""
	if fsp[-4:].lower()!=".hdf" : return
	nm=get_temp_name()
	os.rename(fsp,nm)
	n=EMUtil.get_image_count(nm)
	for i in range(n): EMData(nm,i).write_compressed(fsp,i,bits,nooutliers=nooutliers,level=level)
	os.unlink(nm)

def get_header(filename,i):
	return EMData(filename,i,True).get_attr_dict()

def remove_image(fsp):
	"""This will remove the image file pointed to by fsp. The reason for this function
	to exist is formats like IMAGIC which store data in two files. This insures that
	both files are removed."""

	try:
		os.unlink(fsp)
		if fsp[-4:]==".hed" : os.unlink(fsp[:-3]+"img")
		elif fsp[-4:]==".img" : os.unlink(fsp[:-3]+"hed")
	except: pass

# since there are only 3 odd numbers in the entire list, we remove them, and just limit ourselves to even numbers
# Values thru 512 are carefully calculated as shown on the wiki. Larger values have prime factors 7 or lower.
# 9/2/17 removing low numbers not divisible by 4 for convenience
good_box_sizes=[16, 24, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 84, 96, 100, 104, 112, 120, 128, 132, 140, 168, 180, 192, 196, 208, 216, 220, 224, 240, 256, 260, 288, 300, 320, 352, 360, 384, 416, 440, 448, 480, 512, 540, 560, 576, 588, 600, 630, 640, 648, 672, 686, 700, 720, 750, 756, 768, 784, 800, 810, 840, 864, 882, 896, 900, 960, 972, 980, 1000, 1008, 1024, 1050, 1080, 1120, 1134, 1152, 1176, 1200, 1250, 1260, 1280, 1296, 1344, 1350, 1372, 1400, 1440, 1458, 1470, 1500, 1512, 1536, 1568, 1600, 1620, 1680, 1728, 1750, 1764, 1792, 1800, 1890, 1920, 1944, 1960, 2000, 2016, 2048, 2058, 2100, 2160, 2240, 2250, 2268, 2304, 2352, 2400, 2430, 2450, 2500, 2520, 2560, 2592, 2646, 2688, 2700, 2744, 2800, 2880, 2916, 2940, 3000, 3024, 3072, 3136, 3150, 3200, 3240, 3360, 3402, 3430, 3456, 3500, 3528, 3584, 3600, 3750, 3780, 3840, 3888, 3920, 4000, 4032, 4050, 4096, 4116, 4200, 4320, 4374, 4410, 4480, 4500, 4536, 4608, 4704, 4800, 4802, 4860, 4900, 5000, 5040, 5120, 5184, 5250, 5292, 5376, 5400, 5488, 5600, 5670, 5760, 5832, 5880, 6000, 6048, 6144, 6174, 6250, 6272, 6300, 6400, 6480, 6720, 6750, 6804, 6860, 6912, 7000, 7056, 7168, 7200, 7290, 7350, 7500, 7560, 7680, 7776, 7840, 7938, 8000, 8064, 8100, 8192, 8232, 8400, 8640, 8748, 8750, 8820, 8960, 9000, 9072, 9216, 9408, 9450, 9600, 9604, 9720, 9800, 10000, 10080, 10206, 10240, 10290, 10368, 10500, 10584, 10752, 10800, 10976, 11200, 11250, 11340, 11520, 11664, 11760, 12000, 12096, 12150, 12250, 12288, 12348, 12500, 12544, 12600, 12800, 12960, 13122, 13230, 13440, 13500, 13608, 13720, 13824, 14000, 14112, 14336, 14400, 14406, 14580, 14700, 15000, 15120, 15360, 15552, 15680, 15750, 15876, 16000, 16128, 16200, 16384]
#good_box_sizes=[16,24,32, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 84, 96, 98, 100, 104, 112, 120, 128,130, 132, 140, 150, 154, 168, 180, 182, 192, 196, 208, 210, 220, 224, 240, 250, 256,260, 288, 300, 320, 330, 352, 360, 384, 416, 440, 448, 450, 480, 512, 540, 560, 576, 588, 600, 630, 640, 648, 672, 686, 700, 720, 750, 756, 768, 784, 800, 810, 840, 864, 882, 896, 900, 960, 972, 980, 1000, 1008, 1024, 1050, 1080, 1120, 1134, 1152, 1176, 1200, 1250, 1260, 1280, 1296, 1344, 1350, 1372, 1400, 1440, 1458, 1470, 1500, 1512, 1536, 1568, 1600, 1620, 1680, 1728, 1750, 1764, 1792, 1800, 1890, 1920, 1944, 1960, 2000, 2016, 2048, 2058, 2100, 2160, 2240, 2250, 2268, 2304, 2352, 2400, 2430, 2450, 2500, 2520, 2560, 2592, 2646, 2688, 2700, 2744, 2800, 2880, 2916, 2940, 3000, 3024, 3072, 3136, 3150, 3200, 3240, 3360, 3402, 3430, 3456, 3500, 3528, 3584, 3600, 3750, 3780, 3840, 3888, 3920, 4000, 4032, 4050, 4096, 4116, 4200, 4320, 4374, 4410, 4480, 4500, 4536, 4608, 4704, 4800, 4802, 4860, 4900, 5000, 5040, 5120, 5184, 5250, 5292, 5376, 5400, 5488, 5600, 5670, 5760, 5832, 5880, 6000, 6048, 6144, 6174, 6250, 6272, 6300, 6400, 6480, 6720, 6750, 6804, 6860, 6912, 7000, 7056, 7168, 7200, 7290, 7350, 7500, 7560, 7680, 7776, 7840, 7938, 8000, 8064, 8100, 8192, 8232, 8400, 8640, 8748, 8750, 8820, 8960, 9000, 9072, 9216, 9408, 9450, 9600, 9604, 9720, 9800, 10000, 10080, 10206, 10240, 10290, 10368, 10500, 10584, 10752, 10800, 10976, 11200, 11250, 11340, 11520, 11664, 11760, 12000, 12096, 12150, 12250, 12288, 12348, 12500, 12544, 12600, 12800, 12960, 13122, 13230, 13440, 13500, 13608, 13720, 13824, 14000, 14112, 14336, 14400, 14406, 14580, 14700, 15000, 15120, 15360, 15552, 15680, 15750, 15876, 16000, 16128, 16200, 16384]

def good_size(size):
	"""Will return the next larger 'good' box size (with good refinement performance)"""
#	sizes=[32, 33, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 81, 84, 96, 98, 100, 104, 105, 112, 120, 128,130, 132, 140, 150, 154, 168, 180, 182, 192, 196, 208, 210, 220, 224, 240, 250, 256,260, 288, 300, 320, 330, 352, 360, 384, 416, 440, 448, 450, 480, 512]
	# since there are only 3 odd numbers in the entire list, we remove them, and just limit ourselves to even numbers
	for i in good_box_sizes :
		if i>=size: return i

	return Util.calc_best_fft_size(int(size))

def good_size_small(size):
	"""Will return the next larger 'good' box size (with good refinement performance)"""
#	sizes=[32, 33, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 81, 84, 96, 98, 100, 104, 105, 112, 120, 128,130, 132, 140, 150, 154, 168, 180, 182, 192, 196, 208, 210, 220, 224, 240, 250, 256,260, 288, 300, 320, 330, 352, 360, 384, 416, 440, 448, 450, 480, 512]

	for i in reversed(good_box_sizes) :
		if i<=size: return i

def re_filter_list(listtofilter, regex, invert=False):
	"""
	Filter a list by a regular expression
	"""
	r1 = re.compile(regex,flags=re.I)
	returndict = {}
	for key in listtofilter:
		if bool(r1.search(key)) ^ invert: returndict[key] = listtofilter[key]
	return returndict

def get_optionlist(argv):
	optionlist = []
	for arg1 in argv:
		if arg1[0] == "-":
			argname = arg1.split("=")
			optionlist.append(argname[0].lstrip("-"))
	return optionlist

def intvararg_callback(option, opt_str, value, parser):
	v = [int(i) for i in value.split(',')]
	setattr(parser.values, option.dest, v)
	return

def floatvararg_callback(option, opt_str, value, parser):
	v = [float(i) for i in value.split(',')]
	setattr(parser.values, option.dest, v)
	return

def commandoptions(options,exclude=[]):
	"""This will reconstruct command-line options, excluding any options in exclude"""
	opts=[]
	for opt,val in vars(options).items():
		if opt in exclude or opt=="positionalargs" or val is False or val==None: continue
		if val==True and isinstance(val,bool) : opts.append("--"+opt)
		else: opts.append("--{}={}".format(opt,val))

	return " ".join(opts)


class EMArgumentParser(argparse.ArgumentParser):
	""" subclass of argparser to masquerade as optparser and run the GUI """
	def __init__(self, prog=None,usage=None,description=None,epilog=None,version=None,parents=[],formatter_class=argparse.HelpFormatter,prefix_chars='-',fromfile_prefix_chars=None,argument_default=None,conflict_handler='error',add_help=True,allow_abbrev=True):
		argparse.ArgumentParser.__init__(self,prog=prog,usage=usage,description=description,epilog=epilog,parents=parents,formatter_class=formatter_class,prefix_chars=prefix_chars,fromfile_prefix_chars=fromfile_prefix_chars,argument_default=argument_default,conflict_handler=conflict_handler,add_help=add_help,allow_abbrev=allow_abbrev)

		# A list of options to add to the GUI
		self.optionslist = []
		self.tablist = []

		if version:
			self.add_argument('--version', action='version', version=version)
		self.add_argument("positionalargs", nargs="*")
		self.add_argument('--help-to-html', action='store_true',
		                  help='print this help message in html format')

	def parse_args(self):
		""" Masquerade as optparser parse options """
		if "--help-to-html" in sys.argv[1:]:
			import pandas as pd
			from contextlib import redirect_stdout
			from io import StringIO

			actions = self._get_optional_actions()

			df = pd.DataFrame(columns=['Option', 'Type', 'Description'])

			for i in actions:
				i.type = "None" if not i.type else str(i.type).split("'")[1]
				i.option_strings = ', '.join(i.option_strings)

				if "--help-to-html" in i.option_strings or "--help" in i.option_strings:
					continue

				df = pd.concat([df, pd.DataFrame({'Option': [i.option_strings],
				                                  'Type': [i.type],
				                                  'Description': [i.help]})],
				               ignore_index=True)

			print('<pre>')

			stdout = StringIO()
			with redirect_stdout(stdout):
				self.print_usage()
			print(stdout.getvalue().replace('<', '&lt').replace('>', '&gt'))

			print('</pre>\n')

			print(df.to_html(index=False, justify="center"))

			self.exit()

		parsedargs = argparse.ArgumentParser.parse_args(self)

		return (parsedargs, parsedargs.positionalargs)

	def add_pos_argument(self, **kwargs):
		""" Add a position argument, needed only for the GUI """
		kwargs["positional"]=True
		self.optionslist.append(copy.deepcopy(kwargs))

	def add_header(self, **kwargs):
		""" for the header, you need, title, row, col"""
		kwargs["guitype"]="header"
		self.optionslist.append(copy.deepcopy(kwargs))

	def add_argument(self, *args, **kwargs):

		if "guitype" in kwargs:
			if args[0][:2] == "--":
				kwargs["name"] = args[0][2:]
			else:
				kwargs["name"] = args[1][2:]
			self.optionslist.append(copy.deepcopy(kwargs))
			del kwargs["guitype"]
			del kwargs["row"]
			del kwargs["col"]
			del kwargs["name"]
			if "rowspan" in kwargs: del kwargs["rowspan"]
			if "colspan" in kwargs: del kwargs["colspan"]
			if "expert" in kwargs: del kwargs["expert"]
			if "lrange" in kwargs: del kwargs["lrange"]
			if "urange" in kwargs: del kwargs["urange"]
			if "choicelist" in kwargs: del kwargs["choicelist"]
			if "filecheck" in kwargs: del kwargs["filecheck"]
			if "infolabel" in kwargs: del kwargs["infolabel"]
			if "mode" in kwargs: del kwargs["mode"]
			if "browser" in kwargs: del kwargs["browser"]
			if "dirbasename" in kwargs: del kwargs["dirbasename"]
			if "nosharedb" in kwargs: del kwargs["nosharedb"]
			if "returnNone" in kwargs: del kwargs["returnNone"]

		argparse.ArgumentParser.add_argument(self, *args, **kwargs)

	def getGUIOptions(self):
		return self.optionslist

def parsesym(optstr):
	return Symmetries.get(optstr)

#	# FIXME - this function is no longer necessary since I overwrite the Symmetry3D::get function (on the c side). d.woolford
#	[sym, dict] = parsemodopt(optstr)
#	if sym[0] in ['c','d','h']:
#		dict["nsym"] = int(sym[1:])
#		sym = sym[0]
#
#	return Symmetries.get(sym, dict)

parseparmobj1=re.compile("([^\(]*)\(([^\)]*)\)")	# This parses test(n=v,n2=v2) into ("test","n=v,n2=v2")
parseparmobj2=re.compile("([^=,]*)=([^,]*)")		# This parses "n=v,n2=v2" into [("n","v"),("n2","v2")]
parseparmobj3=re.compile("[^:]\w*=*[-\w.]*") # This parses ("a:v1=2:v2=3") into ("a", "v1=2", "v2=3")
parseparmobj4=re.compile("\w*[^=][\w.]*") # This parses ("v1=2") into ("v1", "2")

def parse_transform(optstr):
	"""This is used so the user can provide the rotation information and get the transform matrix in a convenient form.
	It will parse "dot0,dot1,dot2" and type:name=val:name=val:... then return the tranform matrix t=transform(az,alt,phi)"""

	# We first check for the 'EMAN-style' special case
	tpl=optstr.split(",")
	if len(tpl)==3 :
		try: tpl=[float(i) for i in tpl]
		except:
			raise Exception("Invalid EMAN transform: %s"%optstr)
		return Transform({"type":"eman","az":tpl[0],"alt":tpl[1],"phi":tpl[2]})

	# Now we must assume that we have a type:name=val:... specification
	tpl=optstr.split(":")
	parms={"type":tpl[0]}

	# loop over parameters and add them to the dictionary
	for parm in tpl[1:]:
		s=parm.split("=")
		try : parms[s[0]]=float(s[1])
		except :
			raise Exception("Invalid transform parameter: %s"%parm)

	try: ret=Transform(parms)
	except:
		raise Exception("Invalid transform: %s"%optstr)

	return ret

def unparsemodopt(tupl):
	"""This takes a 2-tuple of the form returned by parsemodopt and returns a corresponding string representation"""

	try:
		if tupl[0]==None : return ""
		if tupl[1]==None or len(tupl[1])==0 : return str(tupl[0])
		parm=["{}={}".format(k,v) for k,v in list(tupl[1].items())]
		return str(tupl[0])+":"+":".join(parm)
	except:
		return ""

def parsemodopt(optstr=None):
	"""This is used so the user can provide the name of a comparator, processor, etc. with options
	in a convenient form. It will parse "dot:normalize=1:negative=0" and return
	("dot",{"normalize":1,"negative":0})"""

	if optstr is None or len(optstr)==0 : return (None,{})
	if optstr.lower()=="none" : return None					# special case doesn't return a tuple

	op2=optstr.split(":")
	if len(op2)==1 or op2[1]=="" : return (op2[0],{})		# name with no options

	r2={}
	for p in op2[1:]:
		try: k,v=p.split("=")
		except:
			print("ERROR: Command line parameter parsing failed on ",optstr)
			print("must have the form name:key=value:key=value")
			return(None,None)

#		v=v.replace("bdb%","bdb:")
		if v.lower()=="true" : v=1
		elif v.lower()=="false" : v=0
		else:
			try: v=int(v)
			except:
				try: v=float(v)
				except:
					if len(v)>2 and v[0]=='"' and v[-1]=='"' : v=v[1:-1]
		r2[k]=v

	return (op2[0],r2)

def parsedict(dictstr):
	cmpparms = {}
	for d in dictstr:
		keyval = d.split(":")
		try:
			cmpparms[keyval[0]] = float(keyval[1])
		except:
			cmpparms[keyval[0]] = keyval[1]
	return cmpparms

parseparmobj_op = re.compile("\+=|-=|\*=|\/=|%=")
parseparmobj_logical = re.compile(">=|<=|==|~=|!=|<|>") 	# finds the logical operators <=, >=, ==, ~=, !=, <, >
parseparmobj_op_words = re.compile("\w*[^=+-\/\*%][\w.]*") # splits ("v1?2") into ("v1","2") where ? can be any combination of the characters "=!<>~"
parseparmobj_logical_words = re.compile("\w*[^=!<>~][\w.]*") # splits ("v1?2") into ("v1","2") where ? can be any combination of the characters "=!<>~"
def parsemodopt_logical(optstr):

	if not optstr or len(optstr)==0 : return (None)

	p_1 = re.findall( parseparmobj_logical_words, optstr )

	if len(p_1)==0: return (optstr,{})
	if len(p_1) != 2:
		print("ERROR: parsemodopt_logical currently only supports single logical expressions")
		print("Could not handle %s" %optstr)
		return (None,None,None)

	p_2 = re.findall( parseparmobj_logical, optstr )

	if len(p_2) != 1:
		print("ERROR: could not find logical expression in %s" %optstr)
		return (None,None,None)


	if ( p_2[0] not in ["==", "<=", ">=", "!=", "~=", "<", ">"] ):
		print("ERROR: parsemodopt_logical %s could not extract logical expression" %(p_2[0]))
		print("Must be one of \"==\", \"<=\", \">=\", \"<\", \">\" \"!=\" or \~=\" ")
		return (None,None,None)

	return p_1[0], p_2[0], p_1[1]


def parsemodopt_operation(optstr):

	if not optstr or len(optstr)==0 : return (None)

	p_1 = re.findall( parseparmobj_op_words, optstr )
	if len(p_1)==0: return (optstr,{})

	if len(p_1) != 2:
		print("ERROR: parsemodopt_logical currently only supports single logical expressions")
		print("Could not handle %s" %optstr)
		return (None,None,None)

	p_2 = re.findall( parseparmobj_op, optstr )
	if len(p_2) != 1:
		print("ERROR: could not find logical expression in %s" %optstr)
		return (None,None,None)


	if ( p_2[0] not in ["+=", "-=", "*=", "/=", "%="]):
		print("ERROR: parsemodopt_logical %s could not extract logical expression" %(p_2[0]))
		print("Must be one of", "+=", "-=", "*=", "/=", "%=")
		return (None,None,None)

	return (p_1[0], p_2[0], p_1[1])

def read_number_file(path):
	"""This will read a text file containing a list of integers. The integers may be separated by any character(s). Any contiguous
	sequence of (0-9) in the file will be treated as a number. '#' is NOT respected as a comment character."""

	try:
		regex = re.compile("[0-9]+")
		return [int(i) for i in regex.findall(open(path,"r").read())]
	except:
		return []


def parse_list_arg(*possible_types):
	"""Return a function that can be passed to argparse.add_argument() as an argument to 'type'.
	The returned function's return type is a list of instances in possible_types, i.e., the parsed
	arguments are converted to the types given in possible_types.

	>>> parse_list_arg(int,float)("3, 4")
	[3, 4.0]
	>>> parse_list_arg(int,int)("3, 4")
	[3, 4]
	>>> parse_list_arg([int,int],[int,int,int])("3,4,5")
	[3, 4, 5]
	>>> parse_list_arg([int,int],[int,int,int])("3,4")
	[3, 4]
	>>> parse_list_arg([int,str],[int,int,int])("3,4")
	[3, '4']
	>>> parse_list_arg([int,str],[int,int,int])("3,4,5")
	[3, 4, 5]
	>>> parse_list_arg([int],[int,int],[int,int,int])("3,4")
	[3, 4]
	>>> parse_list_arg([int],[int,int],[int,int,int])("3")
	[3]
	"""

	types_dict = {}

	# If a single possible set of types is given, they don't have to be wrapped in a list/tuple
	# To determine that, the first item in possible_types is checked if it is a list/tuple instance
	if not isinstance(possible_types[0], (tuple, list)):
		types_dict[len(possible_types)] = possible_types
	else:
		for i in possible_types:
			types_dict[len(i)] = i

	# This is the function that will be passed to argparse.add_argument()'s 'type' argument
	def arg_to_list(s):
		user_input_str = s.split(',')

		if not any(len(user_input_str) == k for k in types_dict.keys()):
			raise argparse.ArgumentTypeError(f"provide {' or '.join(str(i) for i in types_dict.keys())} arguments! See --help for details.")
		else:
			types = types_dict[len(user_input_str)]

		return [type_to_convert(usr_input_val)
				for type_to_convert, usr_input_val
				in zip(types, user_input_str)]

	return arg_to_list


def parse_string_to_slices(seq):
	"""
	>>> parse_string_to_slices("")
	[slice(None, None, None)]
	>>> parse_string_to_slices("3,4")
	[3, 4]
	>>> parse_string_to_slices("3,4,2:5")
	[3, 4, slice(2, 5, None)]
	>>> parse_string_to_slices("3,4,2:5:")
	[3, 4, slice(2, 5, None)]
	>>> parse_string_to_slices("3,4,2:5:2")
	[3, 4, slice(2, 5, 2)]
	>>> parse_string_to_slices("3,4,2:5:2,8,14")
	[3, 4, slice(2, 5, 2), 8, 14]
	>>> parse_string_to_slices("3,4,2::2,8,14")
	[3, 4, slice(2, None, 2), 8, 14]
	>>> parse_string_to_slices("3,4,::2,8,14")
	[3, 4, slice(None, None, 2), 8, 14]
	>>> parse_string_to_slices("3,4,:11:2,8,14")
	[3, 4, slice(None, 11, 2), 8, 14]
	"""

	if not seq:
		return [slice(None, None, None)]

	seq = filter(None, seq.split(','))
	slices = []

	for s in seq:
		if s.isnumeric():
			slices.append(int(s))
		elif ":" in s:
			sl = [None] * 3
			for i, k in enumerate(s.split(":")):
				sl[i] = int(k) if len(k) > 0 else None
			slices.append(slice(*sl))
		else:
			with open(s) as fin:
				file_content = fin.read().replace('\n', ',').replace(' ', ',')

			slices.extend(parse_string_to_slices(file_content))

	return slices


def parse_infile_arg(arg):
	"""
	Parses a string of file name with inclusion and exclusion lists
	Returns the file name and a list of image indices

	Format of the input string:

	filename:inclusion_list^exclusion_list

	inclusion_list/exclusion_list may contain comma-separated
	1. integers
	2. slices (Python's built-in slice representation)
	3. file names

	No spaces are allowed in the passed string.
	Negative integers are allowed as in Python slices.
	However, listed files may contain multi-lines and may contain spaces.

	Slice examples for an image file with N images where N>20:
	2:    -> 2,3,...,N-1
	:3:   -> 0,1,2
	::5   -> 0,5,10,15,...,<N
	:10:2 -> 0,2,4,6,8
	::-1  -> N-1,N-2,N-3,...,2,1,0
	2:9   -> 2,3,4,5,6,7,8

	If a file named test_images.hdf contains 100 images,
	"test_images.hdf:2,4:11:,9,6^25:28,30" will be interpreted as follows.
	filename: test_images.hdf
	inclusion_list: 2,4:11:,19,26:32:2 -> 2,3,4,5,6,7,8,9,10,19,26,28,30
	exclusion_list: 25:28,30      -> 25,26,27,30

	Returned index list: [2,3,4,5,6,7,8,9,10,19,26,28,30] - [25,26,27,30] ->
					   : [2,3,4,5,6,7,8,9,10,19,28]
	"""

	fname, _, seq = arg.partition(':')

	if not (fname and os.path.isfile(fname)):
		raise Exception(f"'{fname}' is not an existing regular file!")

	seq_inc, _, seq_exc = seq.partition('^')

	slices_inc = parse_string_to_slices(seq_inc)
	slices_exc = parse_string_to_slices(seq_exc) if seq_exc else []

	nimg = EMUtil.get_image_count(fname)

	idxs = OrderedDict()
	for i in slices_inc:
		if isinstance(i, int):
			idxs.update({i: None})
		else:
			for k in range(*(i.indices(nimg))):
				idxs.update({k: None})

	ids_exc = set()
	for i in slices_exc:
		if isinstance(i, int):
			ids_exc.add(i)
		else:
			ids_exc.update(range(*(i.indices(nimg))))

	for i in ids_exc:
		if i in idxs:
			idxs.pop(i)

	return fname, tuple(idxs.keys())

def parse_range(rangestr,maxval=None):
	"""parses strings like "1,3,4,6-9,11" and returns (1,3,4,6,7,8,9,11), maxval will support n- without upper values
	if not provided will still work as long as n- isn't used"""

	ret=[]
	for s in rangestr.split(","):
		try: ret.append(int(s))
		except:
			try:
				v1,v2=s.split("-")
				ret.extend(range(int(v1),int(v2)+1))
			except:
				v1=int(s.split("-")[0])
				ret.extend(range(v1,maxval+1))

	return ret


def parse_outfile_arg(arg):
	"""
	Input:
	<filename>:<outbits>[:f|o|<min>|<min>s:<max>|<max>s]
	f - full range
	o - full range with extreme outlier removal
	<min> is absolute min value
	<min>s mean-min*sigma
	if only outbits specified, "o" is implied
	if bits==0 -> floating point with lossless compression
	if bits<0 -> floating point with no compression
	returns:
	(filename,outbits,f|o|None,<min>,<max>,<mins>,<maxs>)
	filename:outbits:rendermin:rendermax

	out.hdf:6:20:1000
	out.hdf:6:-3s:5s
	out.hdf:6:f
	out.hdf:6

	"""

	parm=arg.split(":")
	if len(parm)==0 or len(parm[0])==0:
		raise Exception("Output filename required")

	try: parm[1]=int(parm[1])
	except:	raise Exception("Output filename <filename>:<bits>")

	# only bits specified, default behavior is "o"
	if len(parm)==2: return(parm[0],parm[1],"o",None,None,None,None)

	# special modes
	if   parm[2].lower()=="f" : return(parm[0],parm[1],"f",0,0,0,0)
	elif parm[2].lower()=="o" : return(parm[0],parm[1],"o",0,0,0,0)
	
	# no special mode, either values or sigma coefficients
	vals=[None,None,None,None]
	if "s" in parm[2]:
		try: vals[2]=float(parm[2][:-1])
		except: raise Exception("Output filename: <filename>:<outbits>[:f|o|<min>|<min>s[:<max>|<max>s]]")
	else:
		try: vals[0]=float(parm[2])
		except: raise Exception("Output filename: <filename>:<outbits>[:f|o|<min>|<min>s[:<max>|<max>s]]")

	if "s" in parm[3]:
		try: vals[3]=int(parm[3][:-1])
		except: raise Exception("Output filename: <filename>:<outbits>[:f|o|<min>|<min>s[:<max>|<max>s]]")
	else:
		try: vals[0]=int(parm[2])
		except: raise Exception("Output filename: <filename>:<outbits>[:f|o|<min>|<min>s[:<max>|<max>s]]")

	return(parm[0],parm[1],None,*vals)

def angle_ab_sym(sym,a,b,c=None,d=None):
	"""Computes the angle of the rotation required to go from Transform A to Transform B under symmetry,
	such that the smallest symmetry-related angle is returned. sym may be either a list of Transforms
	or a symmetry specifier, eg "c4". For the two orientations, specify either
	two Transform objects, or the symmetry followed by four floats in the order 
	AltA,AzA,AltB,AzB. Return in degrees."""
	
	if c!=None :
		A=Transform({"type":"eman","alt":a,"az":b})
		B=Transform({"type":"eman","alt":c,"az":d})
	else :
		A=a
		B=b
	
	# easier to do it here
	Bi=B.inverse()
		
	# needs to be a list of Transforms
	if isinstance(sym,str):
		sym=parsesym(sym).get_syms()
		
	#if not (isinstance(A,Transform) and isinstance(B,Transform)):
		#raise Exception,"angle_ab_sym requries two transforms or 4 angles"

	return min([(A*s*Bi).get_rotation("spin")["omega"] for s in sym])

def sock_sendobj(sock,obj):
	"""Sends an object as a (binary) size then a binary pickled object to a socket file object"""
	if obj==None :
		sock.sendall(pack("<Q",0))
		return
	strobj=pickle.dumps(obj,-1)
	sock.sendall(pack("<Q",len(strobj)))
	sock.sendall(strobj)

	# uses socket 'file'
	#if obj==None :
		#sock.write(pack("<I",0))
		#return
	#strobj=pickle.dumps(obj,-1)
	#sock.write(pack("<I",len(strobj)))
	#sock.write(strobj)

	
def sock_recvobj(sock):
	"""receives a packed length followed by a binary (pickled) object from a socket file object and returns"""
	l=sock.recv(8,socket.MSG_WAITALL)
	
	try :
		datlen=unpack("<Q",l)[0]
	except:
		print("Format error in unpacking (%d) '%s'"%(len(l),l))
		raise Exception("Network error receiving object")
	if datlen<=0 :return None
	return pickle.loads(sock.recv(datlen,socket.MSG_WAITALL))

	# uses socket 'file'
	#l=sock.read(4)
	#try :
		#datlen=unpack("<I",l)[0]
	#except:
		#print("Format error in unpacking (%d) '%s'"%(len(l),l))
		#raise Exception("Network error receiving object")
	#if datlen<=0 :return None
	#return pickle.loads(sock.read(datlen))

display_magic=None		# VERY basic security for e3display

def e3display(data,vtype="auto",vname=None,dname="Unknown",settings={},port=31980):
	"""Server-based display function, mainly designed for use with Jupyter Lab sessions, but availble for
	background monitoring of running jobs as well. 
	data - data object of any (appropriate) type, if None will simply update the widget and return a PNG, 
	vtype - "image","imagemx","volume","plot2d","plot3d","histogram"
	vname - name of a new or existing (type specific) display widget to use, if None 'default' will be used 
	dname - name for the data set, may be used as window title, or to support multiple object display in the same widget
	settings - a dictionary of widget-specific settings

	image - 2-D image display, one at a time. data may be a single 2-D or 3-D image, list/tuple of images or
		a list of 4 or 5 np.arrays which will be passed to add a vector overlay to the image display

	imagemx - multiple tiled 2-D image display, pass a list of 2-D images or a 3-D image

	volume - 3-D display, pass a 3-D EMData object

	plot2d - X/Y scatter/line plot. Pass a list of N numpy arrays

	plot3d - X/Y/Z 3-D scatter plot. Pass a list of N numpy arrays

	histogram - 1-D histogram. Pass a single numpy array or a list of N numpy arrays
	"""
	import ipywidgets
	global display_magic

	if display_magic==None:
		magicpath=f"{e2gethome()}/.eman2/server_magic"
		with open(magicpath,"rb") as fin:
			display_magic=fin.read(8)

	try: sock=socket.create_connection(("localhost",port))
	except:
		print(sys.exc_info())
		print("\nFailed to connect. Do you have e2display.py --server running?")
		return False

	sock.sendall(display_magic)

	sock_sendobj(sock,(data,vtype,vname,dname,settings))
	ret=sock_recvobj(sock)
	if ret[0]=="error": raise Exception(f"e3display error: {ret[1]}")
	elif ret[0]=="png": 
		nx,ny,png=ret[1]
		
	return ipywidgets.Image(value=png,format="png",width=nx,height=ny)



def display(img,force_2d=False,force_plot=False):
	"""Generic display function for images or sets of images. You can force images to be displayed in 2-D or as a plot with
	the optional flags"""
	if GUIMode:
		from eman2_gui import emimage
		if isinstance(img,tuple) : img=list(img)
		image = emimage.EMImageWidget(img,None,app,force_2d,force_plot)
		image.show()

		try:
			image.optimally_resize()
		except: pass

		try: image.raise_()
		except: pass

		return image
	else:
		# In non interactive GUI mode, this will display an image or list of images with e2display
		path="%s/tmp"%e2gethome()
		fsp=path+"/tmp.hdf"

		if not os.path.exists(path) : os.mkdir(path)

		try: os.unlink(fsp)
		except: pass

		if isinstance(img,list) or isinstance(img,tuple) :
			for i,j in enumerate(img): j.write_image(fsp,i)
		else:
			img.write_image(fsp,0)

		os.system("e2display.py %s %s"%(("","--single","--plot","--single --plot")[int(force_2d)+2*int(force_plot)],fsp))

def euler_display(emdata_list):
	if len(emdata_list) == 0: return
	if GUIMode:
		from e2eulerxplor import EMEulerWidget
		widget=EMEulerWidget(auto=False,sparse_mode=True)
		module = widget.model
		if isinstance(emdata_list[0],EMData): module.set_emdata_list_as_data(emdata_list)
		elif isinstance(emdata_list[0],Transform):
			module.specify_eulers(emdata_list)
			module.regen_dl()
		widget.show()
	else:
		print("gui mode is disabled")

def browse():
	if GUIMode:
		from eman2_gui.emselector import EMBrowser
		browser = EMBrowser()
		browser.show()
		#app.attach_child(browser)
		#app.show_specific(browser)
	else:
		os.system("e2display.py")

class EMImage(object):
	"""This is basically a factory class that will return an instance of the appropriate EMImage* class """
	def __new__(cls,data=None,old=None,parent=1):
		"""This will create a new EMImage* object depending on the type of 'data'. If
		old= is provided, and of the appropriate type, it will be used rather than creating
		a new instance."""
		if GUIMode:
			from eman2_gui import emimage
			image = emimage.EMImageWidget(data,old,app)
			image.show()
			#app.show_specific(image)
			try: image.optimally_resize()
			except: pass
			return image
		else: print("can not instantiate EMImage in non gui mode")

def plot_image_similarity(im1,im2,skipzero=True,skipnearzero=False):
	"""Will plot pixels in the first image on x vs the same pixel in the second image on y
	with or without zero pixels"""
	n=im1["nx"]*im1["ny"]*im1["nz"]	# the x and y arrays will be this long (- zeros)
	x=[]
	y=[]
	s1=im1["sigma"]
	s2=im2["sigma"]
	for i in range(n):
		if skipzero and (im1[i]==0 or im2[i]==0) : continue
		if skipnearzero and (fabs(im1[i])<old_div(s1,10.0) or fabs(im2[i])<old_div(s2,10.0)) : continue
		x.append(im1[i])
		y.append(im2[i])

	plot((x,y))
	return (x,y)

def plot(data,data2=None,data3=None,show=1,size=(800,600),path="plot.png"):
	"""plots an image or an array using the matplotlib library"""
	if GUIMode:
		from eman2_gui.emplot2d import EMPlot2DWidget
		plotw=EMPlot2DWidget(application=app)
		plotw.set_data(data,"interactive")
		if data2!=None : plotw.set_data(data2,"interactive2")
		if data3!=None : plotw.set_data(data3,"interactive3")
		plotw.setWindowTitle("2D Plot")
		plotw.show()
		try: plotw.raise_()
		except: pass
#		app.show_specific(plotw)
		return plotw
	else :
		import matplotlib
		matplotlib.use('Agg')
		import pylab
		pylab.figure(figsize=(old_div(size[0],72.0),old_div(size[1],72.0)),dpi=72)
		if isinstance(data,EMData) :
			a=[]
			for i in range(data.get_xsize()):
				a.append(data[i])
			pylab.plot(a)
			if data2!=None:
				a=[]
				for i in range(data2.get_xsize()):
					a.append(data2[i])
				pylab.plot(a)
			if data3!=None:
				a=[]
				for i in range(data.get_xsize()):
					a.append(data3[i])
				pylab.plot(a)
		elif isinstance(data,XYData) :
			pylab.plot(data.get_xlist(),data.get_ylist())
		elif isinstance(data,list) or isinstance(data,tuple):
			if isinstance(data[0],list) or isinstance(data[0],tuple) :
				if len(data)>2 :
					for d in data: pylab.plot(d)
				else:
					pylab.plot(data[0],data[1])
					if data2!=None: pylab.plot(data2[0],data2[1])
					if data3!=None: pylab.plot(data3[0],data3[1])
			else:
				try:
					a=float(data[0])
					pylab.plot(data)
					if data2!=None: pylab.plot(data2)
					if data3!=None: pylab.plot(data3)
				except:
					print("List, but data isn't floats")
					return
		else :
			print("I don't know how to plot that type (%s)"%(str(type(data))))
			return

		pylab.savefig(path)
		if show:
			try: os.system("display "+path)
			except: pass

def kill_process(pid):
	'''
	platform independent way of killing a process
	'''
	import os
	import platform
	platform_string = get_platform()
	if platform_string == "Windows":
		# taken from http://www.python.org/doc/faq/windows/#how-do-i-emulate-os-kill-in-windows
		import win32api
		try:
			# I _think_ this will work but could definitely be wrong
			handle = win32api.OpenProcess(1, 0, pid)
			win32api.TerminateProcess(handle,-1)
			win32api.CloseHandle(handle)
			return 1
		except:
			return 0
	else:
		try:
			os.kill(pid,1)
			return 1
		except:
			return 0

def run(cmd,quiet=False):
	"""shortcut redefined all over the place. Might as well put it in one location."""
	if not quiet: print(cmd)
	return launch_childprocess(cmd)

def launch_childprocess(cmd,handle_err=0):
	'''
	Convenience function to launch child processes
	'''
	p = subprocess.Popen(str(cmd)+" --ppid=%d"%os.getpid(), shell=True)

	if get_platform() == 'Windows':
		error = p.wait()	#os.waitpid does not work on windows
	else:
		error = os.waitpid(p.pid, 0)[1]

	if error and handle_err: 
		print("Error {} running: {}".format(error,cmd))
		sys.exit(1)
		
	return error

def process_running(pid):
	'''
	Platform independent way of checking if a process is running, based on the pid
	'''
	import platform
	import os
	platform_string = get_platform()
	if platform_string == "Windows":
		# taken from http://www.python.org/doc/faq/windows/#how-do-i-emulate-os-kill-in-windows
		import win32api
		try:
			# I _think_ this will work but could definitely be wrong
			handle = win32api.OpenProcess(1, 0, pid)
			#win32api.CloseHandle(handle) # is this necessary?
			return 1
		except:
			return 0
	else:
		try:
			os.kill(pid,0)
			return 1
		except:
			return 0

def memory_stats():
	'''
	Returns [total memory in GB,available memory in GB]
	if any errors occur while trying to retrieve either of these values their retun value is -1
	'''
	import platform
	platform_string = get_platform()
	mem_total = -1
	mem_avail = -1
	if platform_string == "Linux":
		try:
			f = open("/proc/meminfo")
			a = f.readlines()
			mt = a[0].split()
			if mt[0] == "MemTotal:":
				mem_total = old_div(float(mt[1]),1000000.0)
			ma = a[1].split()
			if ma[0] == "MemFree:":
				mem_avail = old_div(float(ma[1]),1000000.0)
		except:
			pass

	elif platform_string == "Darwin":
		import subprocess
		status_total, output_total = subprocess.getstatusoutput("sysctl hw.memsize")
		status_used, output_used = subprocess.getstatusoutput("sysctl hw.usermem")
		total_strings = output_total.split()
		if len(total_strings) >= 2: # try to make it future proof, the output of sysctl will have to be double checked. Let's put it in a unit test in
			total_len = len("hw.memsize")
			if total_strings[0][:total_len] == "hw.memsize":
				try:
					mem_total = old_div(float(total_strings[1]),1000000000.0) # on Mac the output value is in bytes, not kilobytes (as in Linux)
				except:pass # mem_total is just -1

		used_strings = output_used.split()
		mem_used = -1
		if len(used_strings) >=2:
			used_len = len("hw.usermem")
			if used_strings[0][:used_len] == "hw.usermem":
				try:
					mem_used = old_div(float(used_strings[1]),1000000000.0) # on Mac the output value is in bytes, not kilobytes (as in Linux)
				except:pass # mem_used is just -1

		if mem_used != -1 and mem_total != -1:
			mem_avail = mem_total - mem_used

	return [mem_total,mem_avail]

def free_space(p=os.getcwd()):
	'''
	Get the free space on the drive that "p" is on
	return value is in bytes
	Works on Linux - needs checking on other platforms
	'''
	s = os.statvfs(p)
	return s.f_bsize*s.f_bavail

def num_cpus():
	'''
	Returns an estimate of the number of cpus available on the current platform
	'''
	import platform
	platform_string = get_platform()
	if platform_string == "Linux":
		try:
			maxphys=0
			cores=1
			for l in open("/proc/cpuinfo","r"):
				if "physical id" in l: maxphys=max(maxphys,int(l.split(":")[-1]))
				if "cpu cores" in l: cores=int(l.split(":")[-1])
			return cores*(maxphys+1)
		except:
			return 2
	elif platform_string == "Windows":
		try:
			cores = os.getenv("NUMBER_OF_PROCESSORS")
			if cores < 1: return 1 # just for safety
			else: return int(cores)
		except:
			return 2
	elif platform_string == "Darwin":
		import subprocess
		status, output = subprocess.getstatusoutput("sysctl hw.logicalcpu")
		strings = output.split()
		cores = 1 # this has to be true or else it's a really special computer ;)
		if len(strings) >=2:
			used_len = len("hw.logicalcpu")
			if strings[0][:used_len] == "hw.logicalcpu": # this essentially means the system call worked
				try:
					cores = int(strings[1])
				except:pass # mem_used is just -1

		if cores < 1:
			print("warning, the number of cpus was negative (%i), this means the MAC system command (sysctl) has been updated and EMAN2 has not accommodated for this. Returning 1 for the number of cores." %cores)
			cores = 2# just for safety, something could have gone wrong. Maybe we should raise instead
		return cores

	else:
		print("error, in num_cpus - unknown platform string:",platform_string," - returning 2")
		return 2

def gimme_image_dimensions2D( imagefilename ):
	"""returns the dimensions of the first image in a file (2-D)"""

	e = EMData()
	e.read_image(imagefilename,0,True)
	return (e.get_xsize(),e.get_ysize())

# get the three dimensions of a an image
def gimme_image_dimensions3D( imagefilename ):

	#pdb.set_trace()

	read_header_only = True
	e = EMData()
	e.read_image(imagefilename,0, read_header_only)
	d = e.get_attr_dict()
	xsize = d["nx"]
	ysize = d["ny"]
	zsize = d["nz"]

	return (xsize, ysize,zsize)

def get_supported_3d_stack_formats():
	'''
	@return a list of the IMAGE formats in EMAN2 that support 3D stacks
	Return is formatted like ['hdf','spi']
	Note using http://blake.bcm.tmc.edu/emanwiki/EMAN2ImageFormats as of April 23 2009
	Note that "bdb:" is not returned though, represented an aberration. Currently the
	calling function must be aware of this
	@Note - it would be nice if this was automatic, not hard coded.
	'''
	return ["hdf","spi","pif","emim"]


def get_supported_2d_stack_formats():
	'''
	@return a list of the IMAGE formats in EMAN2 that support 2D stacks
	Return is formatted like ['hdf','spi']
	Note using http://blake.bcm.tmc.edu/emanwiki/EMAN2ImageFormats as of April 23 2009
	Note that "bdb:" is not returned though, represented an aberration. Currently the
	calling function must be aware of this
	@Note - it would be nice if this was automatic, not hard coded.
	'''
	return ["hdf","spi","pif","emim","img"]


def get_supported_2d_write_formats():
	'''
	@return a list of the IMAGE formats in EMAN2 that can write to disk as a 2D image
	Note using http://blake.bcm.tmc.edu/emanwiki/EMAN2ImageFormats as of July 1 2009
	Note that "bdb:" is not returned though, represented an aberration. Currently the
	calling function must be aware of this
	@Note - it would be nice if this was automatic, not hard coded.
	'''
	return ["mrc","spi","img","hdf"]

def get_supported_3d_formats():
	'''
	@return a list of the IMAGE formats in EMAN2 that 3D images
	Return is formatted like ['hdf','spi']
	Note using http://blake.bcm.tmc.edu/emanwiki/EMAN2ImageFormats as of April 23 2009
	Note that "bdb:" is not returned though, represented an aberration. Currently the
	calling function must be aware of this
	@Note - it would be nice if this was automatic, not hard coded.
	'''
	return ["hdf","spi","mrc","pif","emim","img","vtk","icos","xplor","em","fits"]

def remove_file( file_name, img_couples_too=True ):
	'''
	A function for removing a file from disk. Works for database style file names.
	Adapted so that if you want to remove an img/hed couple it will automatically remove both.
	You can disable the img coupling feature with the the image_couples_too flag
	'''

	if os.path.exists(file_name):

		parts = file_name.split('.')

		file_tag = parts[len(parts)-1]

		if file_tag == 'hed' or file_tag == 'img':
			# get everything that isn't the tag
			name = ''
			for i in range(0,len(parts)-1):
				name = name + parts[i] + '.'

			if img_couples_too:
				os.remove(name+'hed')
				os.remove(name+'img')
			else: os.remove(name+file_tag)
		else:
			if os.path.isfile(file_name): os.remove(file_name)
			elif os.path.isdir(file_name):
				import shutil
				shutil.rmtree(file_name)
			else:
				raise RuntimeError("Unknown url %s" %url)

		return True
	elif db_check_dict(file_name):
		db_remove_dict(file_name)
	else:
		print("Warning, attempt to remove file (%s) that does not exist. No action taken." %file_name)
		return False

# returns the local date and time as a string
def local_datetime(secs=-1):
	"""Returns a timestamp as yyyy/mm/dd hh:mm:ss"""
	from time import localtime,strftime
	if secs<0 : secs=time.time()
	t=localtime(secs)
	return strftime("%Y/%m/%d %H:%M:%S",t)
#	return "%04d/%02d/%02d %02d:%02d:%02d"%t[:6]

def timestamp_diff(t1,t2):
	"""Takes two timestamps in local_datetime() format and subtracts them, result in seconds.
	difftime() available for convenient display."""
	from time import strptime,mktime

	try:
		tt1=mktime(strptime(t1,"%Y/%m/%d %H:%M:%S"))
		tt2=mktime(strptime(t2,"%Y/%m/%d %H:%M:%S"))
	except:
#		print "time error ",t1,t2
		return 0

	return tt2-tt1

def difftime(secs):
	"""Returns a string representation of a time difference in seconds as a Dd hh:mm:ss style string"""

	d=int(floor(old_div(secs,86400)))
	secs-=d*86400
	h=int(floor(old_div(secs,3600)))
	secs-=h*3600
	m=int(floor(old_div(secs,60)))
	secs=int(secs-m*60)

	if d>0 : return "%dd %2d:%02d:%02d"%(d,h,m,secs)
	if h>0 : return "%2d:%02d:%02d"%(h,m,secs)
	return "   %2d:%02d"%(m,secs)

# returns gm time as a string. For example if it's 11:13 pm on the 18th of June 2008 this will return something like
# '23:13:25.14 18/6/2008'
def gm_time_string():
	'''
	Returns gm time as a string. For example if it's 11:13 pm on the 18th of June 2008 this will return something like '23:13:25.14 18-6-2008'
	The use of '/' is intentionally avoided
	'''

	from time import gmtime,time
	a = time()
	b = gmtime(a)
	astr = str(a)
	idx = str.find(astr,'.')
	decimalseconds = astr[idx:len(astr)]

	val = str(b[3])+':'+str(b[4])+':'+str(b[5])+decimalseconds +' '+str(b[2])+'-'+str(b[1])+'-'+str(b[0])
	return val

def is_2d_image_mx(filename):
	'''
	Returns true if the filename exists, is an EM image type, is 2D, and has more than zero images in it
	Note that a single 2D image is by definition a in "image matrix" that has only one member
	'''
	if not file_exists(filename): return False, "File doesn't exist : "+filename

	a = EMData()
	# this is guaranteed not to raise unless the image has been removed from disk since the last call to check_files_are_em_images
	a.read_image(filename,0,True)
	if a.get_ndim() != 2:
		return False, "Image is not 2D :", filename
	elif EMUtil.get_image_count(filename) < 1:
		return False, "Image has not particles in it :", filename
	else:
		return True, "Image is a 2D stack"

def check_files_are_2d_images(filenames):
	'''
	Checks that the files exist, are valid EM types, and are non matrix 2D images
	'''
	fine, message = check_files_are_em_images(filenames)
	if not fine:
		return fine, message
	else:
		for name in filenames:
			if EMUtil.get_image_count(name) > 1:
				return False, "Image contains more than one image :", name

			else:
				read_header_only = True
				a = EMData()
				# this is guaranteed not to raise unless the image has been removed from disk since the last call to check_files_are_em_images
				a.read_image(name,0,read_header_only)
				if a.get_ndim() != 2:
					return False, "Image is not 2D :", name

		return True, "Images are all valid and 2D"


def check_files_are_em_images(filenames):
	'''
	Checks a list of filenames to first verify that they exist, and secondly to verify that they are valid EM images
	Returns bool, string where bool means the filenames are good, and the string is an error message or similar
	'''
	for file in filenames:
		if not os.path.exists(file):
			try:
				is_db = db_check_dict(file)
				if not is_db: raise
			except: return False, "File doesn't exist:"+file

		read_header_only = True
		a = EMData()
		try:
			a.read_image(file,0,read_header_only)
		except:
			return False, "File is not a valid EM image:"+file

	return True,"images are fine"

abs_path=os.path.abspath
#def abs_path(name):
	#'''
	#wraps os.path.absname but detects bdb naming
	#'''
	#return os.path.abspath(name)

def base_name( file_name,extension=False,bdb_keep_dir=False,nodir=False ):
	'''
	wraps os.path.basename but returns something sensible for bdb syntax
	if nodir is set, then the last path element will never be included, otherwise it is included following a set of standard rules.
	'''
	if extension : print("base_name() with extension. please check")
	if bdb_keep_dir : print("base_name() with bdb_keep_dir. please check")

	file_name=str(file_name)

	apath=os.path.relpath(file_name).replace("\\","/").split("/")
	# for specific directories, we want any references to the same micrograph to share an id
	if nodir or (len(apath)>1 and apath[-2] in ("sets","particles","micrographs","movies","movieparticles","ddd","raw","info", "tiltseries", "tomograms", "particles3d", "segmentations")) :
		if extension :
			return os.path.basename(file_name)
		else :
			return os.path.splitext(os.path.basename(file_name))[0].split("__")[0].replace("_ptcls","").replace("_info","")		# double underscore is used to mark tags added to micrograph names

	# but for other files, like classes_xx which users might make selection lists on, we want to include the
	# subdirectory name, to prevent mixing between different refinement directories
	if extension : return "-".join(apath[-2:])
	else : return "-".join(apath[-2:]).rsplit(".",1)[0]

def info_name(file_name,nodir=False):
	"""This will return the name of the info file associated with a given image file, in the form info/basename_info.js"""
	return "info/{}_info.json".format(base_name(file_name,nodir=nodir))

def file_exists( file_name ):
	'''
	A function for checking if a file exists
	basically wraps os.path.exists, but when an img or hed file is the argument, it
	checks for the existence of both images
	Also checks if the argument is a valid dictionary
	'''

	if ( os.path.exists(file_name) ):
		parts = file_name.split('.')
		file_tag = parts[len(parts)-1]

		# get everything that isn't the tag
		name = ''
		for i in range(0,len(parts)-1):
			name = name + parts[i] + '.'

		if ( file_tag == 'hed' ):
			if ( not os.path.exists(name+'img') ):
				print("Warning - %s does not exist" %(name+'img'))
				return False
			else: return True;
		elif (file_tag == 'img'):
			if (not os.path.exists(name+'hed')):
				print("Warning - %s does not exist" %(name+'hed'))
				return False
			else: return True;
		else:
			return True
	else:
		try:
			if db_check_dict(file_name) and EMUtil.get_image_count(file_name) != 0: # a database can exist but have no images in it, in which case we consider it to not exist
				return True
			else: return False
		except: return False


def strip_after_dot(file_name):
	""" a function for stripping the contents of a filename so that all
	 that remains is up to the first '.'
	 eg if given image.sh4. mrc this functions strips the 'sh4.mrc' and returns 'image'
	 FIXME it's probably easiest to do this with regular expressions... """
	idx = str.find(file_name,'.')
	return file_name[0:idx]


def get_platform():
	'''
	wraps platform.system but accommodates for its internal bad programming.
	Returns Darwin, Windows, Linux, or potentially something else
	'''
	import platform
	while True:
		try:
			return platform.system()
		except:
			pass

if get_platform() == "Darwin":
	glut_inited = True # this is a hack for the time being

def display_path(path):
	"""Will generate a suitable reduced path for use in title-bars on windows, etc."""
	
	try: full=os.path.abspath(path)
	except: full=path
	
	full=full.replace("\\","/").split("/")
	if len(full)==1: return full
	return "/".join(full[-3:])

def remove_directories_from_name(file_name,ntk=0):
	'''
	Removes the directories from a file name.
	ntk indicates how many trailing path elements to keep from the normalized path
	'''
	base=os.path.basename(file_name)
	if ntk>0 :
		full=os.path.abspath(file_name)
		if "\\" in full : pre="\\".join(full.split("\\")[-ntk-1:-1])
		else : pre="/".join(full.split("/")[-ntk-1:-1])
		return pre+"/"+base
	return base

def name_has_no_tag(file_name):
	'''
	A convenient way of asking if the file name in has no tag. i.e.
	/home/tmp.jpg would have a tag but /home/tmp would not. Of course
	this function will return true if the argument is the name of a
	folder, but that was not the original intended use
	'''
	idx1 = file_name.rfind("/")
	idx2 = file_name.rfind(".")

	if idx1 >= idx2: return True # -1 == -1
	else: return False


def check_eman2_type(modoptstring, object, objectname, verbose=True):
	''' a function for testing whether a type of Averager, Aligner, Comparitor, Projector, Reconstructor etc
	 can be created from the command line string. Returns false if there are problems
	 examples
	 if not check_eman2_type(options.simaligncmp,Cmps,"Comparitor") : exit(1)
	 if not check_eman2_type(options.simalign,Aligners,"Aligner"): exit(1)
	'''
	if modoptstring == None:
		if verbose:
			print("Error: expecting a string but got python None, was looking for a type of %s" %objectname)
		return False

	if modoptstring == "":
		if verbose:
			print("Error: expecting a string was not empty, was looking for a type of %s" %objectname)
		return False

	try:
		p = parsemodopt(modoptstring)
		if p[0] == None:
			if verbose:
				print("Error: Can't interpret the construction string %s" %(modoptstring))
			return False
		object.get(p[0], p[1])
	except RuntimeError:
		if (verbose):
			print("Error: the specified %s (%s) does not exist or cannot be constructed" %(objectname, modoptstring))
		return False

	return True

def check_eman2_type_string(modoptstring, object, objectname):
	''' a function for testing whether a type of Averager, Aligner, Comparitor, Projector, Reconstructor etc
	can be created from the command line string. Returns false if there are problems
	examples
	error =  check_eman2_type_string(options.simaligncmp,Cmps,"Comparitor")
	if error != None:
		 print error
		 exit(1)
	Another example:
	error = check_eman2_type_string(options.simalign,Aligners,"Aligner") ETC
	'''
	if modoptstring == None:
		return "Error: expecting a string but got python None, was looking for a type of %s" %objectname

	if modoptstring == "":
		return "Error: expecting a string was not empty, was looking for a type of %s" %objectname

	try:
		p = parsemodopt(modoptstring)
		if p[0] == None:
			return "Error: Can't interpret the construction string %s" %(modoptstring)
		object.get(p[0], p[1])
	except RuntimeError:
		return "Error: the specified %s (%s) does not exist or cannot be constructed" %(objectname, modoptstring)


	return None

def qplot(img):
	"""This will plot a 1D image using qplot
	Note that display(img) will automatically plot 1D images.
	"""
	out=open("/tmp/plt.txt","w")
	for i in range(img.get_xsize()):
		out.write("%d\t%f\n"%(i,img.get_value_at(i,0)))
	out.close()
	os.system("qplot /tmp/plt.txt")

def error_exit(s) :
	"""A quick hack until I can figure out the logging stuff. This function
	should still remain as a shortcut"""
	print(s)
	exit(1)

def write_test_boxing_images(name="test_box",num_im=10,type=0,n=100):
	'''

	'''
	if type == 0:
		window_size=(128,128)
		image_size=(4096,4096)
	elif type == 1:
		window_size=(128,128)
		image_size=(4482,6678)
	elif type == 2:
		window_size=(128,128)
		image_size=(8964,13356)

	for i in range(num_im):
		im = test_boxing_image(window_size,image_size,n)
		im.write_image(name+"_"+str(i)+".mrc",0,EMUtil.ImageType.IMAGE_UNKNOWN, False,None,EMUtil.EMDataType.EM_SHORT)

def test_boxing_image(window_size=(128,128),image_size=(4096,4096),n=100):
	'''
	Returns an image useful for testing boxing procedures
	A randomly oriented and randomly flipped scurve test image is inserted into a larger image n times.
	The large image has gaussian noise added to it. Function arguments give control of the window size
	(of the scurve image), the size of the large returned image, and the number of inserted scurve images.
	'''
	ret = EMData(*image_size)

	x_limit = window_size[0]
	y_limit = window_size[0]
	scurve = test_image(0,window_size)
	for i in range(n):
		window = scurve.copy()
		da = Util.get_frand(0,360)
		flip = Util.get_irand(0,1)

		t = Transform({"type":"2d","alpha":da})
		t.set_mirror(flip)
		window.transform(t)
		#if flip:
			#window.process_inplace("xform.flip",{"axis":"x"})

		p = (Util.get_irand(x_limit,image_size[0]-x_limit),Util.get_irand(y_limit,image_size[1]-y_limit))

		ret.insert_clip(window,p)


	noise = EMData(*image_size)
	noise.process_inplace("testimage.noise.gauss")
	ret.add(noise)

	return ret

def test_image_array(subim,size=(0,0),n_array=(1,1)):
	"""This will build a packed array of 2-D images (subim) into a new
image of the specified size. The number of subimages in X and y can be
independently defined. This is useful for making simulated 2-D crystals
(with an orthogonal geometry). Size will default to just fit the specified
number of subimages"""
	if size==(0,0): size=(subim.get_xsize()*n_array[0],subim.get_ysize()*n_array[1])

	out=EMData(*size)

	for x in range(n_array[0]):
		for y in range(n_array[1]):
			out.insert_clip(subim,(old_div(size[0],2)+int((x-old_div(n_array[0],2.0))*subim.get_xsize()),old_div(size[1],2)+int((y-old_div(n_array[1],2.0))*subim.get_ysize())))

	return out

def test_image(type=0,size=(128,128)):
	"""Returns a simple standard test image
	type=0  scurve
	type=1	gaussian noise, 0 mean, sigma 1
	type=2  square
	type=3  hollow square
	type=4  circular sinewave
	type=5  axes
	type=6  linewave
	type=7  scurve plus x,y gradient
	type=8  scurve translated
	type=9  scurve with gaussian noise(mean 0, sigma 1)
	size=(128,128) """
	ret=EMData()
	ret.set_size(*size)
	if type==0 :
		ret.process_inplace("testimage.scurve")
	elif type==1 :
		ret.process_inplace("testimage.noise.gauss")
	elif type==2:
		ret.process_inplace("testimage.squarecube",{"edge_length":old_div(size[0],2)})
	elif type==3:
		ret.process_inplace("testimage.squarecube",{"fill":1,"edge_length":old_div(size[0],2)})
	elif type==4:
		ret.process_inplace("testimage.scurve")
		t = EMData()
		t.set_size(*size)
		t.process_inplace("testimage.linewave",{"period":Util.get_irand(43,143)})
		ret.add(t)
	elif type==5:
		ret.process_inplace("testimage.axes")
	elif type==6:
		ret.process_inplace("testimage.linewave",{"period":Util.get_irand(43,143)})
	elif type==7:
		ret.process_inplace("testimage.scurve")
		ret.mult(10)
		tmp = EMData(*size)
		tmp.process_inplace("testimage.gradient")
		ret.add(tmp)
		tmp.process_inplace("testimage.gradient",{"axis":"y"})
		ret.add(tmp)
		ret.process_inplace("normalize.edgemean")
	elif type==8:
		ret.process_inplace("testimage.scurve")
		t = Transform({"type":"2d","alpha":Util.get_frand(0,360)})
		s = int(old_div(size[0],10))
		t.set_trans(Util.get_irand(-s,s),Util.get_irand(-s,s))
		t.set_mirror(Util.get_irand(0,1))
		ret.transform(t)
	elif type==9:
		ret.process_inplace("testimage.scurve")
		tmp = EMData(*size)
		tmp.process_inplace("testimage.noise.gauss")
		ret.add(tmp)
	elif type==10:
		ret.process_inplace("testimage.scurve")
		t = Transform({"type":"2d","alpha":Util.get_frand(0,360)})
		s = int(old_div(size[0],10))
		t.set_trans(Util.get_irand(-s,s),Util.get_irand(-s,s))
		#t.set_mirror(Util.get_irand(0,1))
		ret.transform(t)
	else:
		raise

	return ret
test_image.broken = True

def test_image_3d(type=0,size=(128,128,128)):
	"""Returns a simple standard test image
	type=0  axes
	type=1  spherical waves
	type=2  test tomo image, size parameter ignored
	type=3  square
	type=4  sphere
	type=5  ellipse with holes
	type=6  random gaussian noise
	type=7  gaussian ellipsoid 3 axes different
	type=8	gaussian ellipsoid Z long x,y same
	type=9	gaussian ellipsoid Z short x,y same
	size=(128,128,128) """
	ret=EMData()
	if len(size) != 3:
		print("error, you can't create a 3d test image if there are not 3 dimensions in the size parameter")
		return None
	if type != 2: ret.set_size(*size)
	else:
		if size!=(256,256,64) : print("Warning, size set to 256x256x64")
		ret.set_size(256,256,64)
	if type==0 :
		ret.process_inplace("testimage.axes")
	elif type==1:
		tmp = EMData()
		tmp.set_size(*size)
		tmp.process_inplace("testimage.sphericalwave",{"wavelength":old_div(size[0],7.0),"phase":0})

		a = tmp.copy()
		a.translate(old_div(size[0],7.0),old_div(size[1],7.0),0)
		ret.process_inplace("testimage.sphericalwave",{"wavelength":old_div(size[0],11.0),"phase":0})
		ret.add(a)

		a = tmp.copy()
		a.translate(old_div(-size[0],7.0),old_div(size[1],7.0),0)
		ret.add(a)

		a = tmp.copy()
		a.translate(old_div(-size[0],7.0),old_div(-size[1],7.0),0)
		ret.add(a)

		a = tmp.copy()
		a.translate(old_div(size[0],7.0),old_div(-size[1],7.0),0)
		ret.add(a)

		ret.process_inplace("normalize")

	elif type==2:
		ret.process_inplace("testimage.tomo.objects")
	elif type==3:
		ret.process_inplace("testimage.squarecube",{"fill":1,"edge_length":old_div(size[0],2)})
	elif type==4:
		ret.process_inplace("testimage.circlesphere",{"radius":int(size[0]*.375)})
	elif type==5:

		t = Transform({"type":"eman","az":60,"alt":30})
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],3),"b":old_div(size[1],5),"c":old_div(size[2],4),"transform":t})

		t = Transform({"type":"eman","az":-45})
		t.set_trans(0,0,0)
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],2),"b":old_div(size[1],16),"c":old_div(size[2],16),"transform":t,"fill":0})

		t.set_trans(0,0,old_div(size[2],6))
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],2),"b":old_div(size[1],16),"c":old_div(size[2],16),"transform":t,"fill":0})

		t.set_trans(0,0,old_div(-size[2],6))
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],2),"b":old_div(size[1],16),"c":old_div(size[2],16),"transform":t,"fill":0})

		t = Transform({"type":"eman","alt":-45})
		t.set_trans(old_div(-size[0],8),old_div(size[1],4),0)
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],16),"b":old_div(size[1],2),"c":old_div(size[2],16),"transform":t,"fill":0})

		t = Transform({"type":"eman","alt":-45})
		t.set_trans(old_div(size[0],8),old_div(-size[1],3.5))
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],16),"b":old_div(size[1],2),"c":old_div(size[2],16),"transform":t,"fill":0})

	elif type==6 :
		ret.process_inplace("testimage.noise.gauss")

	elif type==7:
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],6),"b":old_div(size[0],5),"c":old_div(size[0],3)})

	elif type==8:
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],6),"b":old_div(size[0],6),"c":old_div(size[0],3)})

	elif type==9:
		ret.process_inplace("testimage.ellipsoid",{"a":old_div(size[0],3),"b":old_div(size[0],3),"c":old_div(size[0],6)})

	return ret

# get a font renderer
def get_3d_font_renderer():
	try:
		from libpyGLUtils2 import EMFTGL
		font_renderer = EMFTGL()
		font_renderer.set_face_size(32)
		font_renderer.set_using_display_lists(True)
		font_renderer.set_depth(2)
		pfm = get_platform()
		if pfm in ["Linux","Darwin"]:
			font_renderer.set_font_file_name(e2getinstalldir()+"/fonts/DejaVuSerif.ttf")
			#font_renderer.set_font_file_name(e2getinstalldir()+"/fonts/SourceCodePro-Light.ttf")
		elif pfm == "Windows":
			font_renderer.set_font_file_name("C:\\WINDOWS\\Fonts\\arial.ttf")
		else:
			print("unknown platform:",pfm)
		return font_renderer
	except ImportError:
		#print "Unable to import EMFTGL. The FTGL library may not be installed. Text on 3D and some 2D viewers may not work."
		return None

class EMAbstractFactory(object):
	'''
	see http://blake.bcm.edu/emanwiki/Eman2FactoriesInPython
	'''

	def register(self, methodName, constructor, *args, **kargs):
		"""register a constructor"""
		_args = [constructor]
		_args.extend(args)
#		setattr(self, methodName,Functor(_args, kargs))
		setattr(self, methodName, EMFunctor(*_args, **kargs))

	def unregister(self, methodName):
		"""unregister a constructor"""
		delattr(self, methodName)

class EMFunctor(object):
	'''
	Taken from http://code.activestate.com/recipes/86900/
	'''
	def __init__(self, function, *args, **kargs):
		assert callable(function), "function should be a callable obj"
		self._function = function
		self._args = args
		self._kargs = kargs

	def __call__(self, *args, **kargs):
		"""call function"""
		_args = list(self._args)
		_args.extend(args)
		_kargs = self._kargs.copy()
		_kargs.update(kargs)
		return self._function(*_args, **_kargs)

def isosurface(marchingcubes, threshhold, smooth=False):
	"""Return the Isosurface points, triangles, normals(smooth=True), normalsSm(smooth=False)"""
	marchingcubes.set_surface_value(threshhold)
	d = marchingcubes.get_isosurface(smooth)
	return d['points'], d['faces'], d['normals']

# determines if all of the data dimensions are a power of val
# e.g. if ( data_dims_power_of(emdata,2) print "the dimensions of the image are all a power of 2"
def data_dims_power_of(data,val):

	test_cases = [data.get_xsize()]

	if ( data.get_ndim() >= 2 ):
		test_cases.append(data.get_ysize())
	if ( data.get_ndim == 3):
		test_cases.append(data.get_zsize())

	for i in test_cases:
		x = i
		while ( x > 1 and x % val == 0 ): x /= val
		if x != 1:
			return False

	return True

def pixelprint(self) : return 'Pixel(%d,%d,%d,%g)'%(self.x,self.y,self.z,self.value)
Pixel.__repr__=pixelprint
Pixel.__str__=pixelprint

def ctfprint(self) : return 'EMAN2Ctf().from_dict({"defocus":%1.5f,"bfactor":%1.1f,"ampcont":%1.2f,"apix":%1.3f,"voltage":%1.2f,"cs":%1.2f}) ...'%(self.defocus,self.bfactor,self.ampcont,self.apix,self.voltage,self.cs)
EMAN2Ctf.__repr__=ctfprint
EMAN2Ctf.__str__=ctfprint

def vec2fprint(self) : return "Vec2f(%1.2f,%1.2f)"%(self[0],self[1])
Vec2f.__repr__=vec2fprint
Vec2f.__str__=vec2fprint

def vec3fprint(self) : return "Vec3f(%1.2f,%1.2f,%1.2f)"%(self[0],self[1],self[2])
Vec3f.__repr__=vec3fprint
Vec3f.__str__=vec3fprint

def xformprint(self) :
	try :
		p=self.get_params("2d")
		return "Transform({'tx':%1.2f,'ty':%1.2f,'alpha':%1.3f,'mirror':%1d,'scale':%1.4f,'type':'2d'})"%(p["tx"],p["ty"],p["alpha"],int(p["mirror"]),p["scale"])
	except:
		p=self.get_params("eman")
		return "Transform({'az':%1.3f,'alt':%1.3f,'phi':%1.3f,'tx':%1.2f,'ty':%1.2f,'tz':%1.2f,'mirror':%1d,'scale':%1.4f,'type':'eman'})"%(p["az"],p["alt"],p["phi"],p["tx"],p["ty"],p["tz"],int(p["mirror"]),p["scale"])

Transform.__repr__=xformprint
Transform.__str__=xformprint

def regionprint(self) :
	try : return "Region(%s)"%self.get_string().replace(";",",")
	except: return "Region(Error)"
Region.__repr__=regionprint
Region.__str__=regionprint

def clear_dead_cudajobs():
	# obviously this only works for POSIX system, but then again we don't support CUDA on windows
	locks = glob.glob("%s*"%EMData.getcudalock())
	for lock in locks:
		f = open(lock,"r")
		cpid = int(f.readline())
		f.close()
		try:
			os.kill(cpid, 0)
		except OSError:
			print("removing deadfile ", lock)
			os.unlink(lock)

### Very odd function, I can't find it used anywhere, so I'm commenting it out.
#def set_emdata_array(img, array):
	#"""
	#Return a new EMData object, set its data with a list of float, all all attribute are the same as input image.
	#The array's size must be the same as the img's size, i.e. nx*ny or nx*ny*nz.
	#img - a EMData object
	#array - a list of float data
	#"""
	#import numpy
	#dct = img.get_attr_dict()

	#if len(array) != dict['nz']*dict['ny']*dict['nx']:
		#print "Error: Array's size does not match nx*ny*nz"
		#return

	#numpy_array = numpy.reshape(numpy.array(array, numpy.float32), (dct['nz'], dct['ny'], dct['nx']))
	#new_img = EMNumPy.numpy2em(numpy_array)
	#new_img.set_attr_dict(dct)
	#return new_img

def write_FSC_file(fsc,filename) :
	"""Convenience function takes a standard FSC/FRC resulting from EMData.calc_fourier_shell_correlation and write it to
	disk as a s,FSC text file."""
	sz=len(fsc)//3
	out=open(filename,"w")
	for i in range(sz):
		out.write("{}\t{}\n".format(fsc[i],fsc[i+sz]))

def initializeCUDAdevice():
	# Initialize CUDA upon EMAN2 import. If cuda is not compiled an error will be thrown an nothing will happen
	try:
		clear_dead_cudajobs()
		EMData.cuda_initialize()
	except:
		pass

class LSXFile(object):
	"""This class will manage writing entries to LSX files, which are text files with a defined record length for
rapid access. Each line contains an image number, a filename, and an optional comment, referencing a particle
in another actual image file. Files MUST use the Unix /n convention, not the Windows (/r/n) convention.

The file begins with
#LSX
# Whole file comment, one line of any length
# Line length (including \n)
number<\t>filename<\t>comment
...
"""
	def __init__(self,path,ifexists=False, comments=""):
		"""Initialize the object using the .lst file in 'path'. If 'ifexists' is set, an exception will be raised
if the lst file does not exist."""

		self.path=path
		if len(comments)==0:
			comments="# This file is in fast LST format. All lines after the next line have exactly the number of characters shown on the next line. This MUST be preserved if editing."

		if os.path.isfile(path):
			self.ptr=open(path,"r+")		# file exists
		else:
			if ifexists: raise Exception("Error: lst file {} does not exist".format(path))

			try: os.makedirs(os.path.dirname(path))
			except: pass
			self.ptr=open(path,"w+")	# file doesn't exist
			self.ptr.write("#LSX\n{}\n# 20\n".format(comments))
			self.ptr.flush()

		self.ptr.seek(0)
		l=self.ptr.readline()
		if l==0 or l!="#LSX\n" :
			if l=="#LST\n" :
				#### This is very similar to rewrite(), but is used to convert LST files to LSX files
				self.seekbase=self.ptr.tell()
				tmpfile=open(self.path+".tmp","w")
				tmpfile.write("#LSX\n{}\n".format(comments))

				# we read the entire file, checking the length of each line
				maxlen=0
				while 1:
					ln=self.ptr.readline().strip()
					if len(ln)==0 : break
					maxlen=max(maxlen,len(ln))

				self.linelen=maxlen+1+4						# we make the lines 4 characters longer than necessary to reduce rewrite calls as "n" gets bigger
				tmpfile.write("# {}\n".format(self.linelen))	# the new line length
				newseekbase=tmpfile.tell()

				fmtstr="{{:<{}}}\n".format(self.linelen-1)	# string for formatting

				self.ptr.seek(self.seekbase)
				while 1:
					ln=self.ptr.readline().strip()
					if len(ln)==0 : break
					tmpfile.write(fmtstr.format(ln))

				# close both files
				tmpfile=None
				self.ptr=None
				self.seekbase=newseekbase

				# rename the temporary file over the original
				os.unlink(self.path)
				os.rename(self.path+".tmp",self.path)
				self.ptr=open(self.path,"r+")
				self.ptr.readline()

			else: raise Exception("ERROR: The file {} is not in #LSX format".format(self.path))
		self.filecomment=self.ptr.readline().strip()
		try: self.linelen=int(self.ptr.readline()[1:])
		except:
			print("ERROR: invalid line length in #LSX file {}".format(self.path))
			raise Exception
		self.seekbase=self.ptr.tell()

		# legacy LST file support
		if self.filecomment.startswith("#keys: "):
			print("WARNING: legacy .lst file containing old-style parameters. Support may be removed in future. Consider rewriting file with e2proclst.py")
			self.filekeys=self.filecomment[7:].split(';')
		else: self.filekeys=None

		# potentially time consuming, but this also gives us self.n
		self.normalize()
		self.lock=threading.Lock()

	def __del__(self):
		self.close()

	def __getitem__(self,n):
		return self.read(n)

	def __setitem__(self,n,tupl):
		if len(tupl)==3 : self.write(n,tupl[0],tupl[1],tupl[2])
		else : self.write(n,tupl[0],tupl[1])

	def close(self):
		"""Once you call this, you should not try to access this object any more"""
		if self.ptr!=None :
			self.normalize()
			self.ptr=None

	def write(self,n,nextfile,extfile,jsondict=None):
		"""Writes a record to any location in a valid #LSX file.
n : image number in #LSX file, -1 appends, as does n>= current file len
nextfile : the image number in the referenced image file
extfile : the path to the referenced image file (can be relative or absolute, depending on purpose)
jsondict : optional string in JSON format or a JSON compatible dictionary. values will override header values when an image is read.
"""

		self.lock.acquire()
		if jsondict==None : 
			outln="{}\t{}".format(nextfile,extfile)
		elif isinstance(jsondict,str) and jsondict[0]=="{" and jsondict[-1]=='}' : 
			outln="{}\t{}\t{}".format(nextfile,extfile,jsondict)
		else:
			if not isinstance(jsondict,dict) and jsondict!=None: 
				jsondict={"__default__":jsondict}
			if jsondict!=None:
				jss=json.dumps(jsondict,indent=None,sort_keys=True,separators=(',',':'),default=EMAN2jsondb.obj_to_json)			
				outln="{}\t{}\t{}".format(nextfile,extfile,jss)
			else: outln="{}\t{}".format(nextfile,extfile)
			
		# We can't write in the middle of the file if the existing linelength is too short
		if len(outln)+1>self.linelen : self.rewrite(len(outln))


		fmtstr="{{:<{}}}\n".format(self.linelen-1)	# string for formatting
		outln=fmtstr.format(outln)					# padded output line

		if n<0 or n>=self.n :
			self.ptr.seek(0,os.SEEK_END)		# append
			self.n+=1
		else : self.ptr.seek(self.seekbase+self.linelen*n)		# otherwise find the correct location

		self.ptr.write(outln)
		self.lock.release()

	def read(self,n):
		"""Reads the nth record in the file. Note that this does not read the referenced image, which can be
performed with read_image either here or in the EMData class. Returns a tuple (n extfile,extfile,dict). dict
contains decoded information from the stored JSON dictionary. Will also read certain other legacy comments
and translate them into a dictionary."""
		if n>=self.n : raise IndexError("Attempt to read record {} from #LSX {} with {} records".format(n,self.path,self.n))
		self.lock.acquire()
		n=int(n)
		self.ptr.seek(self.seekbase+self.linelen*n)
		ln=self.ptr.readline().strip().split("\t")
		if len(ln)==2 : ln.append("")
		try: ln[0]=int(ln[0])
		except:
			print(f"Error LSXFile.read({n}). {self.seekbase},{self.linelen},{ln}")
			self.lock.release()
			raise(Exception)
		if len(ln[2])<2: ln[2]={}
		else:
			try: ln[2]=json.loads(ln[2],object_hook=EMAN2jsondb.json_to_obj)
			except:
				if self.filekeys!=None:
					vals=ln[2].split(";")
					ln[2]={self.filekeys[i]:eval(vals[i]) for i in range(len(self.filekeys))}
				elif ';Transform' in ln[2]:
					score=float(ln[2].split(";")[0])
					xf=eval(ln[2].split(";")[1])
					ln[2]={"score_align":score,"xform.projection":xf}
				elif ln[2][:9]=="Transform":
					ln[2]={"xform.projection":eval(ln[2])}
				else:
					ln[2]={"lst_comment":ln[2]}
		self.lock.release()
		return ln

	def read_image(self,N,hdronly=False,region=None):
		"""This reads the image referenced by the nth record in the #LSX file. The same task can be accomplished with EMData.read_image,
but this method prevents multiple open/close operations on the #LSX file."""

		n,fsp,jsondict=self.read(N)
#		print(self.path,n,fsp,jsondict,hdronly,region)
		ret=EMData()
		ret.read_image_c(fsp,n,hdronly,region)
		ret["source_path"]=self.path
		ret["source_n"]=N
		if len(jsondict)>0 :
			for k in jsondict: ret[k]=jsondict[k]

		return ret

	def read_into_image(self,ret=None,N=0,hdronly=False,region=None,is_3d=False,imgtype=IMAGE_UNKNOWN):
		"""This reads the image referenced by the nth record in the #LSX file. The same task can be accomplished with EMData.read_image,
but this method prevents multiple open/close operations on the #LSX file."""

		n,fsp,jsondict=self.read(N)
#		print(self.path,n,fsp,jsondict,hdronly,region)
		ret.read_image_c(fsp,n,hdronly,region,is_3d,imgtype)
		ret["data_n"]=ret["source_n"]
		ret["data_source"]=ret["source_path"]
		ret["source_path"]=self.path
		ret["source_n"]=N
		if len(jsondict)>0 :
			for k in jsondict: ret[k]=jsondict[k]

	def read_images(self,nlst=None,hdronly=False):
		"""This reads a set of images referenced by the nth record in the #LSX file. This is used by read_images in Python when the file is a LST file
		if nlst is None, the entire file is read."""

		# organize the images to read by path to take advantage of read_images performance
		# d2r contains tuples (image number in returned array,image number in file (key),extra data dictionary)
		d2r={}
		if nlst is None or len(nlst)==0:
			for i in range(self.n):
				j,p,d=self.read(i)
				try: d2r[p].append((i,j,d,i))
				except: d2r[p]=[(i,j,d,i)]
			ii=self.n
		else:
			# ii is the index of the image in the array we will eventually return, i is the index of the image
			# in the lst file. j is the index in the referenced image file, p. d is the dictionary of 
			# override values from the lst comment field
			for ii,i in enumerate(nlst):
				j,p,d=self.read(int(i))
				try: d2r[p].append((ii,j,d,i))
				except: d2r[p]=[(ii,j,d,i)]
			ii+=1

#		out=open("dbug.txt","w")
		# we read the actual images with calls to read_images for speed
		# then overlay the metadata overrides from the LST file, and put them in the requested read order
		ret=[None]*ii	# ii is left with the total number of images to be returned
		for fsp in d2r:
			tpls=d2r[fsp]
			imgs=EMData.read_images_c(fsp,[i[1] for i in tpls],IMAGE_UNKNOWN,hdronly)
			for i,tpl in enumerate(tpls):
				imgs[i]["source_path"]=self.path
				try: imgs[i]["source_n"]=int(tpl[3])
				except:
					traceback.print_exc()
					raise Exception(f"Error in read_images: {i},{tpl}")
				for k in tpl[2]: 
					imgs[i][k]=tpl[2][k]
				ret[tpl[0]]=imgs[i]
#				out.write(f"{tpl[0]}\t{i}\t{fsp}\n")

		if None in ret: raise(Exception(f"Error reading {nlst} from {self.path}, {ret.index(None)} is None"))

		return ret


	def __len__(self): return self.n

	def normalize(self):
		"""This will read the entire file and insure that the line-length parameter is valid. If it is not,
it will rewrite the file with a valid line-length. """

		self.ptr.seek(self.seekbase)
		self.n=0
		while 1:
			ln=self.ptr.readline()
			if len(ln)==0 :break
			if len(ln)!=self.linelen :
				self.rewrite()
				break
			self.n+=1

	def rewrite(self,minlen=0):
		"""This will reprocess the entire file to make sure the line length is correct. It does take 2 passes,
but should only be used infrequently, so it should be fine. If specified, minlen is the shortest permissible
line length. Used when a line must be added in the middle of the file."""

		self.ptr.seek(0)

		tmpfile=open(self.path+".tmp","w")
		# copy the header lines
		tmpfile.write(self.ptr.readline())
		tmpfile.write(self.ptr.readline())
		self.ptr.readline()

		# we read the entire file, checking the length of each line
		maxlen=minlen
		while 1:
			ln=self.ptr.readline().strip()
			if len(ln)==0 : break
			maxlen=max(maxlen,len(ln))

		self.linelen=maxlen+1+4						# we make the lines 4 characters longer than necessary to reduce rewrite calls as "n" gets bigger
		tmpfile.write("# {}\n".format(self.linelen))	# the new line length
		newseekbase=tmpfile.tell()

		fmtstr="{{:<{}}}\n".format(self.linelen-1)	# string for formatting

		self.ptr.seek(self.seekbase)
		while 1:
			ln=self.ptr.readline().strip()
			if len(ln)==0 : break
			tmpfile.write(fmtstr.format(ln))

		# close both files
		tmpfile=None
		self.ptr=None
		self.seekbase=newseekbase

		# rename the temporary file over the original
		os.unlink(self.path)
		os.rename(self.path+".tmp",self.path)
		self.ptr=open(self.path,"r+")

#		print "rewrite ",self.linelen

def image_eosplit(filename):
	"""This will take an input image stack in LSX or normal image format and produce output .lst (LSX)
files corresponding to even and odd numbered particles. It will return a tuple with two filenames
corresponding to each 1/2 of the data."""

	# create the even and odd data sets
	print("### Creating virtual stacks for even/odd data")
	if open(filename,"r").read(4)=="#LSX" :				# LSX (fast LST) format
		eset=filename.rsplit(".",1)[0]
		oset=eset+"_odd.lst"
		eset+="_even.lst"
		oute=open(eset,"w")
		outo=open(oset,"w")
		inf=open(filename,"r")

		# copy the first 3 lines to both files
		for i in range(3):
			l=inf.readline()
			oute.write(l)
			outo.write(l)

		# then split the rest
		i=0
		while True:
			l=inf.readline()
			if l=="" : break
			if i%2==0 : oute.write(l)
			else: outo.write(l)
			i+=1

		oute=None
		outo=None
	else :								# This means we have a regular image file as an input
		n=EMUtil.get_image_count(filename)
		eset=filename.rsplit(".",1)[0]+"_even.lst"
		oset=filename.rsplit(".",1)[0]+"_odd.lst"

		try : os.unlink(eset)
		except: pass

		oute=LSXFile(eset)
		for i in range(0,n,2): oute.write(-1,i,filename)
		oute=None

		try : os.unlink(oset)
		except: pass

		oute=LSXFile(oset)
		for i in range(1,n,2): oute.write(-1,i,filename)
		oute=None

	return (eset,oset)

def save_lst_params(lst,fsp, overwrite=True):
	"""Saves a LSX file (fsp) with metadata represented by a list of dictionaries (lst).
	each dictionary must contain 'src', the image file containing the actual image and
	'idx' the index in that file. Additional keys will be stored in the LSX metadata
	region. Overwrite existing file by default."""
	if len(lst)==0: raise(Exception,"ERROR: save_lst_params with empty list")
	
	if overwrite:
		if os.path.isfile(fsp):
			os.remove(fsp)
	
	lsx=LSXFile(fsp)
	for d in lst:
		dct=d.copy()
		p=dct.pop("src")
		n=dct.pop("idx")
		lsx.write(-1,n,p,dct)

def load_lst_params(fsp , imgns=None):
	"""Reads the metadata for all of the images in an LSX file (fsp) with an optional list of
	image numbers (imgns, iterable or None)"""
	lsx=LSXFile(fsp,True)
	if imgns==None or len(imgns)==0: imgns=range(lsx.n)
	
	ret=[]
	for i in imgns:
		n,p,d=lsx.read(i)
		ret.append({"idx":n,"src":p,**d})
		
	return ret
	
##########
#### replace a few EMData methods with python versions to intercept 'bdb:' filenames
##########
def db_emd_init(self, *parms):
	"""
	__init__( (object)arg1) -> object :

		C++ signature :
			void* __init__(_object*)

	__init__( (object)arg1, (object)that) -> object :
		Construct from an EMData (copy constructor).
		Performs a deep copy.

		that - the EMData to copy

		C++ signature :
			void* __init__(_object*,EMAN::EMData)

	__init__( (object)arg1, (str) filespec,(int)image_idx,[(int) header_only],[(Region) region],[(int) is_3d]) -> object :
		Construct from an image file.

		filename - the image file name
		image_index the image index for stack image file(default = 0)

		C++ signature :
			void* __init__(_object*,std::string [,int])

	__init__( (object)arg1, (object)nx, (object)ny [, (object)nz [, (object)is_real]]) -> object :
		makes an image of the specified size, either real or complex.
		For complex image, the user would specify the real-space dimensions.

		nx - size for x dimension
		ny - size for y dimension
		nz size for z dimension(default=1)
		is_real - boolean to specify real(true) or complex(false) image(default=True)

		C++ signature :
			void* __init__(_object*,int,int [,int [,bool]])
"""
	if len(parms) > 0 and len(parms) < 5 and isinstance(parms[0], str):
		self.__initc()
		self.read_image(*parms)
		return

	self.__initc(*parms)
	return


EMData.__initc = EMData.__init__
EMData.__init__ = db_emd_init


# def transform_to_str(self):
## I added the initial return so that it looks good in the eman2 file browser
# s = "\n%.3f %.3f %.3f %.3f\n" %(self.at(0,0),self.at(0,1),self.at(0,2),self.at(0,3))
# s += "%.3f %.3f %.3f %.3f\n" %(self.at(1,0),self.at(1,1),self.at(1,2),self.at(1,3))
# s += "%.3f %.3f %.3f %.3f\n" %(self.at(2,0),self.at(2,1),self.at(2,2),self.at(2,3))
# s += "0.000 0.000 0.000 1.000"
# return s

# Transform.__str__ = transform_to_str

#lsxcache=None


def compressible_formats():
	return ('.hdf', '.jpeg', '.jpg', '.mrc', '.mrcs', '.png', '.tiff', '.tif', '.df3', '.pgm')


def is_file_compressible(fsp):
	return Path(fsp).suffix.lower() in compressible_formats()


def db_read_image(self, fsp, *parms, **kparms):
	"""read_image(filespec,image #,[header only],[region],[is_3d],[imgtype])

	This function can be used to read a set of images from a file or bdb: database. Pass in
	the filename, or bdb: specification, optionally a list of image numbers to read (or None),
	and a flag indicating that only the image headers should be read in. If only the headers
	are read, accesses to the image data in the resulting EMData objects will be invalid."""
	#	print "RI ",fsp,str(parms)

	if fsp is None or fsp=="" : raise Exception("read_image without filename")

	if fsp[:4].lower() == "bdb:":
		print("ERROR: BDB is not supported in this version of EMAN2. You must use EMAN2.91 or earlier to access legacy data.")
		return None
	if fsp[-4:].lower()==".lst":
		return LSXFile(fsp).read_into_image(self,*parms)
		#global lsxcache
		#if lsxcache==None or lsxcache.path!=fsp: lsxcache=LSXFile(fsp,True)
		#return lsxcache.read_image(parms[0])

	fsp, idxs = parse_infile_arg(fsp)

	if len(parms) > 0 and parms[0]:
		idx = idxs[parms[0]] 
	else:
		idx = idxs[0]

	parms = idx, *parms[1:]

	if len(kparms) != 0:
		if 'img_index' not in kparms:
			kparms['img_index'] = 0
		if 'header_only' not in kparms:
			kparms['header_only'] = False
		if 'region' not in kparms:
			kparms['region'] = None
		if 'is_3d' not in kparms:
			kparms['is_3d'] = False
		if 'imgtype' not in kparms:
			kparms['imgtype'] = IMAGE_UNKNOWN

	return self.read_image_c(fsp, *parms, **kparms)

EMData.read_image_c = EMData.read_image
EMData.read_image = db_read_image


def db_read_images(fsp, *parms):
	"""EMData.read_images(filespec,[image # list],[header only])

	This function can be used to read a set of images from a file or bdb: database. Pass in
	the filename, or bdb: specification, optionally a list of image numbers to read (or None),
	and a flag indicating that only the image headers should be read in. If only the headers
	are read, accesses to the image data in the resulting EMData objects will be invalid"""
	if fsp[:4].lower() == "bdb:":
		print("ERROR: BDB is not supported in this version of EMAN2. You must use EMAN2.91 or earlier to access legacy data.")
		return []

	if fsp[-4:].lower()==".lst":
		return LSXFile(fsp).read_images(*parms)
		#global lsxcache
		#if lsxcache==None or lsxcache.path!=fsp: lsxcache=LSXFile(fsp,True)
		#return lsxcache.read_images(*parms)

	fsp, idxs = parse_infile_arg(fsp)

	if len(parms) > 0 and parms[0]:
		idxs = [idxs[i] for i in parms[0]]

	parms = idxs, *parms[1:]

	return EMData.read_images_c(fsp, *parms)


EMData.read_images_c = staticmethod(EMData.read_images)
EMData.read_images = staticmethod(db_read_images)


def db_write_image(self, fsp, *parms):
	"""write_image(fsp,image #,[image type],[header only],[region],[storage type],[use host endian])

	help(EMUtil.ImageType) for a list of valid image types, eg EMUtil.ImageType.IMAGE_MRC
	help(EMUtil.EMDataType) for a list of valid storage types

	Writes an image to a file or a bdb: entry. Note that for bdb: specifications, only image # is supported.
	the remaining options are ignored
	"""
	#	print "In db_write_image, WI ",fsp,str(parms)

	if fsp[:4].lower() == "bdb:":
		print("ERROR: BDB is not supported in this version of EMAN2. You must use EMAN2.91 or earlier to access legacy data.")
		return

	elif ":" in fsp:
		idx = parms[0] if parms else 0

		return self.write_compressed(fsp, idx)

	return self.write_image_c(fsp, *parms)


EMData.write_image_c = EMData.write_image
EMData.write_image = db_write_image

def im_write_compressed(self,fsp,n,bits=8,minval=0,maxval=0,nooutliers=True,level=1,erase=False):
	"""
This flexible image writing routine will write compressed HDF images (or a single image). It may be called
on an instance:

img.write_compressed(fsp,n,bits,...)

or a list of images:

EMData.write_compressed(list,fsp,n_first_img,bits,...)

Compression itself is lossless, but uses variable bit reduction which is (usually) lossy. 

Ignores any render_min/render_max already in the image header. To use those values, pass them as arguments
and do not specify nooutliers.

If nooutliers is set, this will override any specified minval/maxval. If maxval<=minval then the full range of 
image values will be included in the bit reduced image. If large outliers are present, this may effectively remove 
almost all of the information in the image. Nooutliers will eliminate up to ~0.01% of the extreme pixel values from the
histogram by clamping the values at the extrema. For integrating mode raw micrograph data, nooutliers is 
highly recommended, though it is probably unnecessary with counting mode images, unless there is an extremely
bad normalization value.

Compression level (0-9) is passed to the gzip routine. While all levels are lossless, higher levels will
achieve progressively less additional compression and progressively more time. The default level of 1
provides good compression without having a significant performance impact.

Somewhat counterintuitively, the noisier the data, the fewer bits that are required to fully represent the
image. That is, raw micrographs can safely be stored as 3-4 bits, whereas a reconstructed, filtered volume
may require 8 or more bits.

If erase is set, the file will be deleted if it exists, before writing. Obviously this must be used with caution,
but a problem with overwriting existing compressed HDF images is that the original images are not truly deleted
and the file size will increase.
"""
	is_compression_syntax = (":" in fsp)

	if is_compression_syntax:
		fsp, bits, mode, rendermin_abs, rendermax_abs, rendermin_s, rendermax_s = parse_outfile_arg(fsp)

		if not is_file_compressible(fsp):
			raise Exception(f"Only {[i.strip('.') for i in compressible_formats()]} "
			                f"formats are supported by write_compressed()")
		if mode=="f" or mode is None: 
			nooutliers=False
			minval,maxval=0,0
		elif mode=="o" : nooutliers=True
		#print(fsp, bits, mode, rendermin_abs, rendermax_abs, rendermin_s, rendermax_s)

	if isinstance(self,EMData):
		self=[self]
	
	if erase:
		try: os.unlink(fsp)
		except: pass
	
	if n==-1: 
		try: n=EMUtil.get_image_count(fsp)
		except: n=0
	
	for i,im in enumerate(self):
		if not isinstance(im,EMData) : raise(Exception,"write_compressed() requires a list of EMData objects")

		# no compression
		if bits<0 :
			im.del_attr("render_min")
			im.del_attr("render_max")
			im.del_attr("render_bits")
			im.del_attr("render_compress_level")
			im.write_image_c(fsp,i+n)	# bits<0 implies no compression
			continue

		im["render_bits"]=bits
		im["render_compress_level"]=level
		# float compression, no range change
		if bits==0 :
			im.del_attr("render_min")
			im.del_attr("render_max")
			im.write_image_c(fsp,i+n,EMUtil.ImageType.IMAGE_UNKNOWN,0,None,EMUtil.EMDataType.EM_COMPRESSED)
			continue

		### This is an important option, as it will be the default in many cases. It makes an effort to intelligently
		### eliminate extreme outliers, while not modifying cases with highly nongaussian distributions
		if nooutliers:		# nooutliers may be specified either via function call or in filename, same code
			nxyz=im.get_size()
			maxout=max(nxyz//20000,8)		# at most 1 in 20000 pixels should be considered an outlier on each end
			h0=im["minimum"]
			h1=im["maximum"]
			hs=(h1-h0)/4095
			hist=im.calc_hist(4096,h0,h1)
			
			#ok, we're doing this approximately
			for sl in range(2048):
				if sum(hist[:sl+1])>maxout: break
			for sh in range(4095,2048,-1):
				if sum(hist[sh:])>maxout: break

			im["render_min"]=sl*hs+h0
			im["render_max"]=sh*hs+h0
		elif is_compression_syntax and not mode=="f":	# specified range in filename
			if rendermin_abs is not None : im["render_min"]=max(rendermin_abs,im["minimum"])
			else: im["render_min"]=max(im["mean"]-rendermin_s*im["sigma"],im["minimum"])
			if rendermax_abs is not None : im["render_max"]=min(rendermax_abs,im["maximum"])
			else: im["render_max"]=min(im["mean"]+rendermax_s*im["sigma"],im["maximum"])
		elif maxval<=minval:		# Full range
			im["render_min"]=im["minimum"]
			im["render_max"]=im["maximum"]
		else:				# specified range in function call
			im["render_min"]=max(float(minval),im["minimum"])
			im["render_max"]=min(float(maxval),im["maximum"])

		# integer-only images are handled differently during compression, so if value truncation occurs in
		# an integer only image, truncation must be to an integer
		if im["all_int"]:
			im["render_min"]=round(im["render_min"])
			im["render_max"]=round(im["render_max"])
		
		#print(f"PY {im['render_min']} - {im['render_max']} {im['minimum']} - {im['maximum']}   {im['render_bits']}")
		# would like to use the new write_images, but it isn't quite ready yet.
		im.write_image_c(fsp,i+n,EMUtil.ImageType.IMAGE_UNKNOWN,0,None,EMUtil.EMDataType.EM_COMPRESSED)
	
EMData.write_compressed=im_write_compressed


def db_get_image_count(fsp):
	if ":" in fsp:
		fsp, idxs = parse_infile_arg(fsp)
		return len(idxs)
	else:
		return EMUtil.get_image_count_c(fsp)

EMUtil.get_image_count_c = staticmethod(EMUtil.get_image_count)
EMUtil.get_image_count = db_get_image_count


__doc__ = \
"EMAN classes and routines for image/volume processing in \n\
single particle reconstructions.\n\
\n\
The following classes are defined: \n\
  EMData - the primary class to process electronic microscopy images. \n\
 \n\
  Quaternion - implements mathematical quaternion. \n\
  Region - defines a rectangular 2D/3D region. \n\
  Transform - defines a transformation including rotation, translation, and different Euler angles. \n\
  Vec3i - a 3-element integer vector. \n\
  Vec3f - a 3-element floating number vector. \n\
\n\
  EMObject - A wrapper class for int, float, double, string, EMData, XYData, list etc. \n\
  Pixel - defines a pixel's 3D coordinates and its value. \n\
  SimpleCtf - defines EMAN CTF parameter model. \n\
  XYData - implements operations on a series of (x y) data pair. \n\
\n\
  Aligners - Aligner factory. Each Aligner alignes 2D/3D images. \n\
  Averagers - Averager factory. Each Averager averages a set of images. \n\
  Cmps  - Cmp factory. Each Cmp defines a way to do image comparison. \n\
  Processors - Processor factory. Each processor implements an image-processing algorithm. \n\
  Projectors - Projector factory. Each Projector implements an algorithm for 3D image projection. \n\
  Reconstructors - Reconstructor factory. Each Reconstructor implements an algorithm for image reconstruction. \n\
\n\
  EMUtil - Defines utility functions for EMData-related operations. \n\
  TestUtil - Unit test utility functions. \n\
  Util - Generic utility functions. \n\
\n\
  EMNumPy - Defines functions for conversions between EMData and numpy array. \n\
  Log - Log output at different verbose level. \n\
  PointArray - Point array. \n\
\n\
  dump_aligners() - Print out all Aligners and their parameters. \n\
  dump_averagers() - Print out all Averagers and their sparameters. \n\
  dump_cmps() - Print out all Cmps and their parameters. \n\
  dump_processors() - Print out all Processor`s and their parameters. \n\
  dump_projectors() - Print out all Projectors and their parameters. \n\
  dump_reconstructors() - Print out all Reconstructors and their parameters. \n\
  dump_analyzers() - Print out all Analyzers anf their parameters. \n\
"
