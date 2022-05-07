#!/usr/bin/env python
#
# Author: Steven Ludtke, 06/16/2008 (sludtke@bcm.edu)
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

from future import standard_library

standard_library.install_aliases()
from builtins import range
from builtins import object
import atexit
import weakref
from pickle import loads, dumps
from zlib import compress, decompress
import os
import os.path
import signal
import sys
import time
import fnmatch
import random
import threading
import traceback
from EMAN2jsondb import js_open_dict

print("ERROR: If you are seeing this message, you may be using a legacy program which makes use of the BDB: storage system, which was deprecated a decade ago.\
 We are no longer able to support this system in EMAN2 due to issues with the underlying library. If you rely on this program, please contact the developers\
 to let them know it needs to be updated.")

# I started impementing a basic compatibility module using JSONDB below, but didn't have time. Let's see how things go...  --steve

#from libpyEMData2 import EMData
#from libpyUtils2 import EMUtil

#def bdb2new(name):
	#"Converts a bdb: url to the compatibility filename"
	#if name.lower()[:4]=="bdb:" :
		#if "?" in name:
			#print("ERROR: BDB name with '?' not supported in compatibility module. Exiting to avoid data loss.")
			#sys.exit(1)
		#return "bdb_"+name[4:]+".json"
	#return name

#def e2filemodtime(path):
	#"""This will determine the last modification time for a file or bdb: database in seconds"""
	#path=bdb2new(path)

	#return int(os.stat(path)[8])


## replicated in EMAN2.py
#def e2gethome():
	#"""platform independent path with '/'"""
	#if (sys.platform != 'win32'):
		#url = os.getenv("HOME")
	#else:
		#if (os.getenv("HOMEPATH") == '\\'):
			## Contributed by Alexander Heyne <AHEYNE@fmp-berlin.de>
			#url = os.getenv("USERPROFILE")
			#url = url.lstrip('CDEFG:')  # could also use substr
			#url = url.replace("\\", "/")
		#else:
			#url = os.getenv("HOMEPATH")
			#url = url.replace("\\", "/")
	#return url

## replicated in EMAN2.py
#def e2getcwd():
	#"""platform independent path with '/'"""
	#url = os.getcwd()
	#url = url.replace("\\", "/")
	#return url


#def db_convert_path(path):
	#'''
	#Converts a path to bdb syntax. The path pass must contain "EMAN2DB".
	#For instance, if the path is particles/EMAN2DB/1000_ptcls,
	#will return bdb:particles#1000_ptcls
	#This function is in infancy. It may not foresee all circumstances

	#'''
	#if not isinstance(path, str): raise RuntimeError("Path has to be a string")
	#d = path.split("/")

	#if len(d) < 2: raise ValueError("The path is invalid, try something like 'EMAN2DB/ptcls'")

	#if d[-2] != "EMAN2DB": raise ValueError("The path is invalid, try something like 'EMAN2DB/ptcls'")
	## all right we should be in business

	#dir = ""
	#if len(d) > 2:
		#for i in range(0, len(d) - 2):
			#if i != 0:
				#dir += "/"
			#dir += d[i]  # accrue the directory

	#if path[0] == "/" and dir[0] != "/":  # need to double check that this is needed
		#dir = "/" + dir

	#ret = "bdb:"
	#if dir != "":
		#ret += dir + "#"
	#ret += d[-1]
	#return ret


#def db_parse_path(url):
	#print("ERROR: db_parse_path() unsupported in compatibility")
	#sys.exit(1)

##def db_parse_path(url):
	##"""Takes a bdb: url and splits out (path,dict name,keylist). If keylist is None,
	##implies all of the images in the file. Acceptable urls are of the form:
	##bdb:dict
	##bdb:relative/path#dict
	##bdb:/path/to#dict?key,key,key
	##bdb:/path/to/dict   (also works, but # preferred)
	##"""

	##if url[:4].lower() != "bdb:": raise Exception("Invalid URL, bdb: only (%s)" % url)
	##url = url.replace("~", e2gethome())
	##url = url[4:].rsplit('#', 1)
	##if len(url) == 1: url = url[0].rsplit("/", 1)
	##if len(url) == 1:
		##url = url[0].split("?", 1)
		##if len(url) == 1: return (".", url[0], None)  # bdb:dictname
		##url.insert(0, ".")  # bdb:dictname?keys
	##else:
		##u1 = url[1].split("?", 1)
		##if os.name == 'nt':
			##if url[0] == '.':
				##url[0] = e2getcwd()
			##elif url[0][0] != '/' and url[0][1] != ':':
				##url[0] = e2getcwd() + "/" + url[0]  # relative path
		##else:
			##if url[0][0] != '/': url[0] = e2getcwd() + "/" + url[0]  # relative path
		##if url[0][:3] == "../": url[0] = e2getcwd().rsplit("/", 1)[0] + "/" + url[0]  # find cwd if we start with ../
		##if len(u1) == 1: return (url[0], url[1], None)  # bdb:path/to#dictname
		##url = [url[0]] + url[1].split("?")  # bdb:path/to#dictname?keys

	### if we're here we need to deal with ?keys
	##u2 = url[2].split(",")
	##if len(u2) == 1:
		### u2[0]=u2[0].lower()
		##if u2[0][:7].lower() == "select.":
			##ddb = EMAN2DB.open_db(".")
			##ddb.open_dict("select")
			##if u2[0][7:] not in ddb.select: raise Exception("Unknown selection list %s" % u2[0][7:])
			##return (url[0], url[1], ddb.select[u2[0][7:]])  # bdb:path/to#dict?select/name
		##elif u2[0][:8].lower() == "exclude.":
			##ddb = EMAN2DB.open_db(".")
			##ddb.open_dict("select")
			##all_set = set(range(0, EMUtil.get_image_count("bdb:" + url[0] + "#" + url[1])))
			##if u2[0][8:] not in ddb.select:
				##return (url[0], url[1], list(all_set))  # if the exclusion list is missing, we exclude nothing
			###				raise Exception,"Unknown selection list %s"%u2[0][8:]
			##exc_set = set(ddb.select[u2[0][8:]])
			##rem_set = all_set - exc_set
			##return (url[0], url[1], list(rem_set))  # bdb:path/to#dict?select/name
	###		return url											# bdb:path/to#dict?onekey
	### make sure integer string keys are integers
	##for i in range(len(u2)):
		##try:
			##u2[i] = int(u2[i])
		##except:
			##pass
	##return (url[0], url[1], u2)  # bdb:path/to#dict?k1,k2,k3,...

#def db_open_dict(url, ro=False, with_keys=False):
	#ret=js_open_dict(bdb2new(url))
	#if with_keys: return(ret,ret.keys())
	#return ret

#def db_close_dict(url):
	#js_close_dict(bdb2new(url))


#def db_remove_dict(url):
	#js_close_dict(bdb2new(url))
	#try: os.unlink(bdb2new(url))
	#except: pass

#def db_check_dict(url, readonly=True):
	#"""Checks for the existence of the named dictionary, and whether it can be opened
	#read/write (or just read). Deals only with bdb: urls. Returns false for other specifiers"""

	#if os.access(bdb2new(url), os.W_OK | os.R_OK): return True
	#if os.access(bdb2new(url),os.R_OK) and readonly: return True

	#return False


#def db_list_dicts(url):
	#"""Gives a list of available databases (dicts) at a given path. No '#' specification should be given."""
	#path, dictname, keys = db_parse_path(url)

	#if not "#" in url: path = path + "/" + dictname

	##	if path=="." and dictname!="" : path=dictname
	#try:
		#ld = os.listdir(path + "/EMAN2DB")
	#except:
		#return []

	#ld = [i[:-4] for i in ld if i[-4:] == ".bdb" and i != "00image_counts.bdb"]

	#return ld



#def db_get_image_count(fsp):
	#"""get_image_count(path)

#Takes a path or bdb: specifier and returns the number of images in the referenced stack."""
	#if fsp[:4].lower() == "bdb:":
		#path, dictname, keys = db_parse_path(fsp)
		#if keys == None:
			## This dictionary contains image counts for the others in this directory
			#db2 = db_open_dict("bdb:%s#%s" % (path, "00image_counts"))
			## If this is true, we need to update the dictionary
			#if (dictname in db2 and os.path.getmtime("%s/EMAN2DB/%s.bdb" % (path, dictname)) > db2[dictname][
				#0]) or dictname not in db2:
				#db = db_open_dict(fsp, True)
				#try:
					#im = db[0]
					#sz = (im["nx"], im["ny"], im["nz"])
				#except:
					#sz = (0, 0, 0)
				#db2[dictname] = (time.time(), len(db), sz)

			#return db2[dictname][1]

		#else:  # if the user specifies the key in fsp, we ignore parms
			#db = db_open_dict(fsp, True)
			#n = 0
			#for i in keys:
				#if i in db: n += 1
			#return n
	#try:
		#ret = EMUtil.get_image_count_c(fsp)
	#except:
		##		print"Error with get_image_count on : ",fsp
		#raise Exception(fsp)
	#return ret


## NOTE: Toshio Moriya 2018/08/02
#def db_fix_image_count(fsp):
	#"""fix_image_count(path)

#Takes a path or bdb: specifier and returns the number of images in the referenced stack."""
	#if fsp[:4].lower() == "bdb:":
		#path, dictname, keys = db_parse_path(fsp)
		#if keys == None:
			## This dictionary contains image counts for the others in this directory
			#db2 = db_open_dict("bdb:%s#%s" % (path, "00image_counts"))

			## If this is true, we need to update the dictionary
			## if (dictname in db2 and os.path.getmtime("%s/EMAN2DB/%s.bdb"%(path,dictname))>db2[dictname][0]) or dictname not in db2 :
			##
			## NOTE: Toshio Moriya 2018/08/02
			## Always force to update values of 00image_counts.bdb
			#db = db_open_dict(fsp, True)
			#try:
				#im = db[0]
				#sz = (im["nx"], im["ny"], im["nz"])
			#except:
				#sz = (0, 0, 0)
			#db2[dictname] = (time.time(), len(db), sz)

			#return db2[dictname][1]

		#else:  # if the user specifies the key in fsp, we ignore parms
			#db = db_open_dict(fsp, True)
			#n = 0
			#for i in keys:
				#if i in db: n += 1
			#return n
	#try:
		#ret = EMUtil.get_image_count_c(fsp)
	#except:
		##		print"Error with get_image_count on : ",fsp
		#raise Exception(fsp)
	#return ret


#EMUtil.get_image_count_c = staticmethod(EMUtil.get_image_count)
#EMUtil.get_image_count = staticmethod(db_get_image_count)
#EMUtil.fix_image_count = staticmethod(db_fix_image_count)


#def db_get_image_info(fsp):
	#"""get_image_info(path)

#Takes a bdb: specifier and returns the number of images and image dimensions."""
	#if fsp[:4].lower() == "bdb:":
		#path, dictname, keys = db_parse_path(fsp)
		#if keys == None:
			## This dictionary contains image counts for the others in this directory
			#db2 = db_open_dict("bdb:%s#%s" % (path, "00image_counts"))

			## If this is true, we need to update the dictionary
			#if (dictname in db2 and os.path.getmtime("%s/EMAN2DB/%s.bdb" % (path, dictname)) > db2[dictname][
				#0]) or dictname not in db2 or db2[dictname][2][0] == 0:
				##				print "update ",dictname,os.path.getmtime("%s/EMAN2DB/%s.bdb"%(path,dictname)),db2[dictname][0],db2.has_key(dictname)
				#db = db_open_dict(fsp, True)
				#try:
					#im = db[0]
					#sz = (im["nx"], im["ny"], im["nz"])
				#except:
					#sz = (0, 0, 0)
				#db2[dictname] = (time.time(), len(db), sz)

			#return db2[dictname][1:]

		#else:  # if the user specifies the key in fsp, we ignore parms
			#db = db_open_dict(fsp, True)
			#n = 0
			#for i in keys:
				#if i in db: n += 1
			#sz = db[keys[0]]
			#sz = (sz["nx"], sz["ny"], sz["nz"])
			#return (n, sz)
	#img = EMData(fsp, 0, True)
	#return (EMUtil.get_image_count_c(fsp), (img["nx"], img["ny"], img["nz"]))


#def db_get_all_attributes(fsp, *parms):
	#if fsp[:4].lower() == "bdb:":
		#db = db_open_dict(fsp, True)
		#attr_name = parms[-1]
		#if "?" in fsp:
			#keys = fsp[fsp.rfind("?") + 1:].split(",")
			#for i in range(len(keys)):
				#try:
					#keys[i] = int(keys[i])
				#except:
					#pass
			#return [db.get_attr(i, attr_name) for i in keys]
		#else:
			#if len(parms) == 1:
				#keys = list(range(0, len(db)))
			#else:
				#keys = parms[0]
		#return [db.get_attr(i, attr_name) for i in keys]
	#return EMUtil.get_all_attributes_c(fsp, *parms)


#EMUtil.get_all_attributes_c = staticmethod(EMUtil.get_all_attributes)
#EMUtil.get_all_attributes = staticmethod(db_get_all_attributes)


__doc__ = \
	"""This module supports the concept of a local database for storing data and
	metadata associated with a particular EMAN2 refinement. Data stored in this
	database may be extracted into standard flat-files, but use of a database
	with standard naming conventions, etc. helps provide the capability to log
	the entire refinement process."""




