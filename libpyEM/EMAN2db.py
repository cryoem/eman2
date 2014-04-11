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

import atexit
import weakref
from cPickle import loads,dumps
from zlib import compress,decompress
import os
import os.path
import signal
import sys
import time
import fnmatch
import random
import threading
import traceback

try:
	from bsddb3 import db
except ImportError, e:
#	print "WARNING: Could not import bsddb3; falling back on the older bsddb. Consider installing BerkeleyDB and bsddb3:", e
	from bsddb import db

from libpyEMData2 import EMData
from libpyUtils2 import EMUtil

try:
	a=frozenset()
except:
	from sets import Set
	frozenset=Set

# If set, fairly verbose debugging information will be written to the console
# larger numbers will increase the amount of output
DBDEBUG=0

# Flags used to open the database environment
envopenflags=db.DB_CREATE|db.DB_INIT_MPOOL|db.DB_INIT_LOCK|db.DB_INIT_LOG|db.DB_THREAD
#envopenflags=db.DB_CREATE|db.DB_INIT_MPOOL|db.DB_INIT_LOCK|db.DB_INIT_LOG|db.DB_THREAD|db.DB_REGISTER|db.DB_RECOVER

# If set, databases will be opened without any caching, locking, etc.
# DANGER, unless DB is single-threaded and local, don't do any writing
BDB_CACHE_DISABLE=True			# This is the default for EMAN2.1, since SPARX prefers it and EMAN2 no longer uses BDB in any significant way

# This hardcoded value is the maximum number of DBDict objects which will be simultaneously open
# when more than this are open, we begin closing the oldest ones. They will be reopened on-demand,
# but this will prevent resource exhaustion
MAXOPEN=64

# maximum number of times a task will be restarted before we assume there is something fundamentally
# flawed about the task itself
MAXTASKFAIL=10

cachesize=80000000


def DB_cleanup(signum=None,stack=None):
	if signum in (2,15) :
		if len(DBDict.alldicts)>0 :
			try: nopen=len([i for i in DBDict.alldicts if i.bdb!=None])
			except: nopen=0
			print "Program interrupted (%d), closing %d databases, please wait (%d)"%(signum,nopen,os.getpid())
			for i in DBDict.alldicts.keys(): print i.name
		if stack!=None : traceback.print_stack(stack)
	for d in DBDict.alldicts.keys():
		d.forceclose()
		d.lock.acquire(False)		# prevents reopening
	for e in EMAN2DB.opendbs.values(): e.close()
	for d in DBDict.alldicts.keys():
		try : d.lock.release()	# This will (intentionally) allow the threads to fail, since the environment is now closed
		except : pass
	if signum in (2,15) :
		print "Shutdown complete, exiting"
		sys.stderr.flush()
		sys.stdout.flush()
		#parallel_process_exit()
		#cuda_exit()
		os._exit(1)
	"""
	else:
		parallel_process_exit()
		cuda_exit()
	"""

def cuda_exit():
	try:	# will fail if cuda is not installed
		EMData.cuda_cleanup()
	except:
		pass

def parallel_process_exit():
	# Compete HACK to prevent EMAN2DB creation if one deosn't already exisit. Need to do this b/c when anything in EMAN2PAR gets improted, and EMAN2DB is created!!!
	if os.access('EMAN2DB',os.R_OK):
		# Kill any running process from e2paralle.py running on localhost. If none are running nothing happens
		from EMAN2PAR import EMLocalTaskHandler
		for proc in EMLocalTaskHandler.allrunning.values():
			proc.terminate()
			os.kill(proc.pid,signal.SIGKILL)

# if the program exits nicely, close all of the databases
atexit.register(DB_cleanup)

# if we are killed 'nicely', also clean up (assuming someone else doesn't grab this signal)
signal.signal(2,DB_cleanup)
signal.signal(15,DB_cleanup)

def e2filemodtime(path):
	"""This will determine the last modification time for a file or bdb: database in seconds"""
	if path[:4].lower()=="bdb:" :
		p=db_parse_path(path)
		path=p[0]+"/EMAN2DB/"+p[1]+".bdb"

	return int(os.stat(path)[8])

def e2gethome() :
	"""platform independent path with '/'"""
	if(sys.platform != 'win32'):
		url=os.getenv("HOME")
	else:
		if(os.getenv("HOMEPATH") == '\\'):
			#Contributed by Alexander Heyne <AHEYNE@fmp-berlin.de>
			url=os.getenv("USERPROFILE")
			url=url.lstrip('CDEFG:') #could also use substr
			url=url.replace("\\","/")
		else:
			url=os.getenv("HOMEPATH")
			url=url.replace("\\","/")
	return url

def e2getcwd() :
	"""platform independent path with '/'"""
	url=os.getcwd()
	url=url.replace("\\","/")
	return url

def db_convert_path(path):
	'''
	Converts a path to bdb syntax. The path pass must contain "EMAN2DB".
	For instance, if the path is particles/EMAN2DB/1000_ptcls,
	will return bdb:particles#1000_ptcls
	This function is in infancy. It may not foresee all circumstances

	'''
	if not isinstance(path,str): raise RuntimeError("Path has to be a string")
	d = path.split("/")

	if len(d) < 2: raise ValueError("The path is invalid, try something like 'EMAN2DB/ptcls'")

	if d[-2] != "EMAN2DB": raise ValueError("The path is invalid, try something like 'EMAN2DB/ptcls'")
	# all right we should be in business

	dir = ""
	if len(d) > 2:
		for i in range(0,len(d)-2):
			if i != 0:
				dir += "/"
			dir += d[i] # accrue the directory

	if path[0] == "/" and dir[0] != "/": # need to double check that this is needed
		dir = "/" + dir

	ret = "bdb:"
	if dir != "":
		ret += dir+"#"
	ret += d[-1]
	return ret


def db_parse_path(url):
	"""Takes a bdb: url and splits out (path,dict name,keylist). If keylist is None,
	implies all of the images in the file. Acceptable urls are of the form:
	bdb:dict
	bdb:relative/path#dict
	bdb:/path/to#dict?key,key,key
	bdb:/path/to/dict   (also works, but # preferred)
	"""

	if url[:4].lower()!="bdb:": raise Exception,"Invalid URL, bdb: only (%s)"%url
	url=url.replace("~",e2gethome())
	url=url[4:].rsplit('#',1)
	if len(url)==1 : url=url[0].rsplit("/",1)
	if len(url)==1 :
		url=url[0].split("?",1)
		if len(url)==1 : return(".",url[0],None)			# bdb:dictname
		url.insert(0,".")									# bdb:dictname?keys
	else :
		u1=url[1].split("?",1)
		if os.name=='nt':
			if url[0] == '.' :
				url[0]=e2getcwd()
			elif url[0][0]!='/' and url[0][1]!=':' :
				url[0]=e2getcwd()+"/"+url[0] 	# relative path
		else:
			if url[0][0]!='/' : url[0]=e2getcwd()+"/"+url[0]	# relative path
		if url[0][:3]=="../": url[0]=e2getcwd().rsplit("/",1)[0]+"/"+url[0]	# find cwd if we start with ../
		if len(u1)==1 : return(url[0],url[1],None)			# bdb:path/to#dictname
		url=[url[0]]+url[1].split("?")						# bdb:path/to#dictname?keys

	# if we're here we need to deal with ?keys
	u2=url[2].split(",")
	if len(u2)==1 :
		#u2[0]=u2[0].lower()
		if u2[0][:7].lower()=="select." :
			ddb=EMAN2DB.open_db(".")
			ddb.open_dict("select")
			if not ddb.select.has_key(u2[0][7:]) : raise Exception,"Unknown selection list %s"%u2[0][7:]
			return (url[0],url[1],ddb.select[u2[0][7:]])		# bdb:path/to#dict?select/name
		elif u2[0][:8].lower()=="exclude." :
			ddb=EMAN2DB.open_db(".")
			ddb.open_dict("select")
			all_set = set(range(0,EMUtil.get_image_count("bdb:"+url[0]+"#"+url[1])))
			if not ddb.select.has_key(u2[0][8:]) :
				return (url[0],url[1],list(all_set))			# if the exclusion list is missing, we exclude nothing
#				raise Exception,"Unknown selection list %s"%u2[0][8:]
			exc_set = set(ddb.select[u2[0][8:]])
			rem_set = all_set - exc_set
			return (url[0],url[1],list(rem_set))		# bdb:path/to#dict?select/name
#		return url											# bdb:path/to#dict?onekey
	# make sure integer string keys are integers
	for i in range(len(u2)) :
		try: u2[i]=int(u2[i])
		except: pass
	return (url[0],url[1],u2)								# bdb:path/to#dict?k1,k2,k3,...


def db_open_dict(url,ro=False,with_keys=False):
	"""opens a DB through an environment from a db:/path/to/db#dbname string. If you want to specify a specific image by key,
	you can specify the key as:  db:/path/to/db#dbname?key
	If key is an integer, it will be converted to an integer before lookup. For commands like read_images, key may also be a
	comma-separated list of keys. Thus it is impossible to access data items
	with keys like '1' instead of (int)1 using this mechanism. ro is a read only flag, which will disable caching as well."""

	path,dictname,keys=db_parse_path(url)

	ddb=EMAN2DB.open_db(path)
	ddb.open_dict(dictname,ro=ro)

	if with_keys :
		if ro : return (ddb[dictname+"__ro"],keys)
		return (ddb[dictname],keys)

	if ro : return ddb[dictname+"__ro"]
	return ddb[dictname]

def db_close_dict(url):
	"""Closes a named dictionary to free resources. Otherwise open dictionaries are cached.
	After closing, a dict CAN be reopened, but references to existing dict objects should not be used
	after closing. Ignores non bdb: urls"""

	try: path,dictname,keys=db_parse_path(url)
	except: return

	ddb=EMAN2DB.open_db(path)
	ddb.close_dict(dictname)

	return

def db_remove_dict(url):
	"""closes and deletes a database using the same specification as db_open_dict"""
	path,dictname,keys=db_parse_path(url)

	ddb=EMAN2DB.open_db(path)
	ddb.remove_dict(dictname)

	return

def db_check_dict(url,readonly=True):
	"""Checks for the existence of the named dictionary, and whether it can be opened
	read/write (or just read). Deals only with bdb: urls. Returns false for other specifiers"""

	if len(url) < 4 or url[:4] != "bdb:": return False
  	path,dictname,keys=db_parse_path(url)

	path=path+"/EMAN2DB/"+dictname+".bdb"
#	print path
	if os.access(path,os.W_OK|os.R_OK) : return True
	if os.access(path,os.R_OK) and readonly : return True

	return False

def db_list_dicts(url):
	"""Gives a list of available databases (dicts) at a given path. No '#' specification should be given."""
	path,dictname,keys=db_parse_path(url)

	if not "#" in url : path=path+"/"+dictname

#	if path=="." and dictname!="" : path=dictname
	try:
		ld=os.listdir(path+"/EMAN2DB")
	except: return []

	ld=[i[:-4] for i in ld if i[-4:]==".bdb" and i!="00image_counts.bdb"]

	return ld

##########
#### replace a few EMData methods with python versions to intercept 'bdb:' filenames
##########
def db_emd_init(self,*parms):
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

    __init__( (object)arg1, (object)filename [[, (object)image_index],header_only]) -> object :
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
	if len(parms)<5 and len(parms)>0 and isinstance(parms[0],str) and parms[0][:4].lower()=="bdb:":
		self.__initc()
		self.read_image(*parms)
		return
	else :
#		print "toC:", parms
		if len(parms)>2 and isinstance(parms[0],str) :
			self.__initc()
			self.read_image_c(*parms)			# this handles Region reading, which isn't supported in the C++ constructor
		else : self.__initc(*parms)
		return

EMData.__initc=EMData.__init__
EMData.__init__=db_emd_init

#def transform_to_str(self):
	## I added the initial return so that it looks good in the eman2 file browser
	#s = "\n%.3f %.3f %.3f %.3f\n" %(self.at(0,0),self.at(0,1),self.at(0,2),self.at(0,3))
	#s += "%.3f %.3f %.3f %.3f\n" %(self.at(1,0),self.at(1,1),self.at(1,2),self.at(1,3))
	#s += "%.3f %.3f %.3f %.3f\n" %(self.at(2,0),self.at(2,1),self.at(2,2),self.at(2,3))
	#s += "0.000 0.000 0.000 1.000"
	#return s

#Transform.__str__ = transform_to_str

def db_read_image(self,fsp,*parms,**kparms):
	"""read_image(filespec,image #,[header only],[region],[is_3d])

	This function can be used to read a set of images from a file or bdb: database. Pass in
	the filename, or bdb: specification, optionally a list of image numbers to read (or None),
	and a flag indicating that only the image headers should be read in. If only the headers
	are read, accesses to the image data in the resulting EMData objects will be invalid."""
#	print "RI ",fsp,str(parms)

	if fsp[:4].lower()=="bdb:" :
		db,keys=db_open_dict(fsp,True,True)

		if len(parms)>1 and parms[1] : nodata=1
		else: nodata=0

		if len(parms)>2 and parms[2] : region=parms[2]
		else: region = None

		if keys:
			key=keys[parms[0]]
		else:
			try: key=parms[0]
			except: key=0

		x=db.get(key,target=self,nodata=nodata,region=region)
		if x==None : raise Exception("Could not access "+str(fsp)+" "+str(key))
#		except:
#			raise Exception("Could not access "+str(fsp)+" "+str(key))
		return None
	if len(kparms)!=0:
		if not kparms.has_key('img_index'):
			kparms['img_index'] = 0
		if not kparms.has_key('header_only'):
			kparms['header_only'] = False
		if not kparms.has_key('region'):
			kparms['region'] = None
		if not kparms.has_key('is_3d'):
			kparms['is_3d'] = False
	return self.read_image_c(fsp,*parms, **kparms)

EMData.read_image_c=EMData.read_image
EMData.read_image=db_read_image

def db_read_images(fsp,*parms):
	"""EMData.read_images(filespec,[image # list],[header only])

	This function can be used to read a set of images from a file or bdb: database. Pass in
	the filename, or bdb: specification, optionally a list of image numbers to read (or None),
	and a flag indicating that only the image headers should be read in. If only the headers
	are read, accesses to the image data in the resulting EMData objects will be invalid"""
	if fsp[:4].lower()=="bdb:" :
		db,keys=db_open_dict(fsp,True,True)
		if len(parms)>1: nodata=parms[1]
		else: nodata=0
		if keys:
			if len(parms)>0 :
				if not parms[0] or len(parms[0])==0 : parms[0]=keys
				return [db.get(keys[i]) for i in parms[0]]
			return [db.get(i,nodata=nodata) for i in keys]
		else :
			if len(parms)==0 : keys=range(0,len(db))
			else : keys=parms[0]
			if not keys or len(keys)==0 : keys=range(len(db))
		return [db.get(i,nodata=nodata) for i in keys]
	if len(parms)>0 and (parms[0]==None or len(parms[0])==0) :
		parms=(range(EMUtil.get_image_count(fsp)),)+parms[1:]
	return EMData.read_images_c(fsp,*parms)

EMData.read_images_c=staticmethod(EMData.read_images)
EMData.read_images=staticmethod(db_read_images)

def db_write_image(self,fsp,*parms):
	"""write_image(fsp,image #,[image type],[header only],[region],[storage type],[use host endian])

	help(EMUtil.ImageType) for a list of valid image types, eg EMUtil.ImageType.IMAGE_MRC
	help(EMUtil.EMDataType) for a list of valid storage types

	Writes an image to a file or a bdb: entry. Note that for bdb: specifications, only image # is supported.
	the remaining options are ignored
	"""
#	print "In db_write_image, WI ",fsp,str(parms)
		
	if fsp[:4].lower()=="bdb:" :
		db,keys=db_open_dict(fsp,False,True)
		if keys :			# if the user specifies the key in fsp, we ignore parms
			if len(keys)>1 : raise Exception,"Too many keys provided in write_image %s"%str(keys)
			if isinstance(keys[0],int) and keys[0]<0 : raise Exception,"Negative integer keys not allowed %d"%keys[0]
			db[keys[0]]=self
			return
		if len(parms)==0 : parms=[0]
		if parms[0]<0 : parms=(len(db),)+parms[1:]

		key = 0
		region=None
		if len(parms) > 0: key = parms[0]
		#if len(parms) > 1: read_header = parms[1]
		# parsm[2] is unused : image_type=EMUtil.ImageType.IMAGE_UNKNOWN
		if len(parms) > 3: region = parms[3]

		db.set(self,key=key,region=region,txn=None) # this makes region writing work
#		db[parms[0]] = self # this means region writing does not work
		return
	return self.write_image_c(fsp,*parms)

EMData.write_image_c=EMData.write_image
EMData.write_image=db_write_image

def db_get_image_count(fsp):
	"""get_image_count(path)

Takes a path or bdb: specifier and returns the number of images in the referenced stack."""
	if fsp[:4].lower()=="bdb:" :
		path,dictname,keys=db_parse_path(fsp)
		if keys==None :
			# This dictionary contains image counts for the others in this directory
			db2=db_open_dict("bdb:%s#%s"%(path,"00image_counts"))

			# If this is true, we need to update the dictionary
			if (db2.has_key(dictname) and os.path.getmtime("%s/EMAN2DB/%s.bdb"%(path,dictname))>db2[dictname][0]) or not db2.has_key(dictname) :
				db=db_open_dict(fsp,True)
				try:
					im=db[0]
					sz=(im["nx"],im["ny"],im["nz"])
				except : sz=(0,0,0)
				db2[dictname]=(time.time(),len(db),sz)

			return db2[dictname][1]

		else :			# if the user specifies the key in fsp, we ignore parms
			db=db_open_dict(fsp,True)
			n=0
			for i in keys:
				if i in db : n+=1
			return n
	try:
		ret=EMUtil.get_image_count_c(fsp)
	except:
#		print"Error with get_image_count on : ",fsp
		raise Exception,fsp
	return ret

EMUtil.get_image_count_c=staticmethod(EMUtil.get_image_count)
EMUtil.get_image_count=staticmethod(db_get_image_count)

def db_get_image_info(fsp):
	"""get_image_info(path)

Takes a bdb: specifier and returns the number of images and image dimensions."""
	if fsp[:4].lower()=="bdb:" :
		path,dictname,keys=db_parse_path(fsp)
		if keys==None :
			# This dictionary contains image counts for the others in this directory
			db2=db_open_dict("bdb:%s#%s"%(path,"00image_counts"))

			# If this is true, we need to update the dictionary
			if (db2.has_key(dictname) and os.path.getmtime("%s/EMAN2DB/%s.bdb"%(path,dictname))>db2[dictname][0]) or not db2.has_key(dictname) or db2[dictname][2][0]==0:
#				print "update ",dictname,os.path.getmtime("%s/EMAN2DB/%s.bdb"%(path,dictname)),db2[dictname][0],db2.has_key(dictname)
				db=db_open_dict(fsp,True)
				try:
					im=db[0]
					sz=(im["nx"],im["ny"],im["nz"])
				except :
					sz=(0,0,0)
				db2[dictname]=(time.time(),len(db),sz)

			return db2[dictname][1:]

		else :			# if the user specifies the key in fsp, we ignore parms
			db=db_open_dict(fsp,True)
			n=0
			for i in keys:
				if i in db : n+=1
			sz=db[keys[0]]
			sz=(sz["nx"],sz["ny"],sz["nz"])
			return (n,sz)
	img=EMData(fsp,0,True)
	return (EMUtil.get_image_count_c(fsp),(img["nx"],img["ny"],img["nz"]))


def db_get_all_attributes(fsp,*parms):
	if fsp[:4].lower()=="bdb:" :
		db=db_open_dict(fsp,True)
		attr_name = parms[-1]
		if "?" in fsp:
			keys=fsp[fsp.rfind("?")+1:].split(",")
			for i in range(len(keys)):
				try: keys[i]=int(keys[i])
				except: pass
			return [db.get_attr(i, attr_name) for i in keys]
		else :
			if len(parms)==1 : keys=range(0,len(db))
			else : keys=parms[0]
		return [db.get_attr(i, attr_name) for i in keys]
	return EMUtil.get_all_attributes_c(fsp,*parms)

EMUtil.get_all_attributes_c=staticmethod(EMUtil.get_all_attributes)
EMUtil.get_all_attributes=staticmethod(db_get_all_attributes)

#############
###  Task Management classes
#############

class EMTaskQueue:
	"""This class is used to manage active and completed tasks through an
	EMAN2DB object. Pass it an initialized EMAN2DB instance. Tasks are each
	assigned a number unique in the local database. The active task 'max'
	keeps increasing. When a task is complete, it is shifted to the
	tasks_done list, and the 'max' value is increased if necessary. There is
	no guarantee in either list that all keys less than 'max' will exist."""

	lock=threading.Lock()
	caching=False		# flag to prevent multiple simultaneous caching

	def __init__(self,path=None,ro=False):
		"""path should point to the directory where the disk-based task queue will reside without bdb:"""
		if path==None or len(path)==0 : path="."
		self.path=path
		self.active=db_open_dict("bdb:%s#tasks_active"%path,ro)		# active tasks keyed by id
		self.complete=db_open_dict("bdb:%s#tasks_complete"%path,ro)	# complete tasks
		self.nametodid=db_open_dict("bdb:%s#tasks_name2did"%path,ro)	# map local data filenames to did codes
		self.didtoname=db_open_dict("bdb:%s#tasks_did2name"%path,ro)	# map data id to local filename
		self.precache=db_open_dict("bdb:%s#precache_files"%path,ro)		# files to precache on clients, has one element "files" with a list of paths

		#if not self.active.has_key("max") : self.active["max"]=-1
		#if not self.active.has_key("min") : self.active["min"]=0
		#if not self.complete.has_key("max") : self.complete["max"]=0

	def __len__(self) : return len(self.active)

	def get_task(self,clientid=0):
		"""This will return the next task waiting for execution"""
		EMTaskQueue.lock.acquire()
		for tid in sorted(self.active.keys()):
			task=self.active[tid]
			if isinstance(task,int) : continue
			if task==None :
				print "Missing task ",tid
				continue
			if task.starttime==None:
				task.starttime=time.time()
				task.clientid=clientid
				self.active[tid]=task
				EMTaskQueue.lock.release()
				return task

		EMTaskQueue.lock.release()
		return None

	def add_group(self):
		"""returns a new (unique) group id to be used for related tasks"""
		EMTaskQueue.lock.acquire()
		try: ret=self.active["grpctr"]+1
		except: ret=1
		self.active["grpctr"]=ret
		EMTaskQueue.lock.release()
		return ret

	def todid(self,name):
		"""Returns the did for a path, creating one if not already known"""
		fmt=e2filemodtime(name)
		try :
			did=self.nametodid[name]			# get the existing did from the cache (modtime,int)
			if fmt!=did[0]	:					# if the file has been changed, we need to assign a new did
				del self.didtoname[did]
				EMTaskQueue.lock.release()
				raise Exception
		except:
			did=(fmt,random.randint(0,999999))	# since there may be multiple files with the same timestamp, we also use a random int
			while (self.didtoname.has_key(did)):
				did=(fmt,random.randint(0,999999))

		self.nametodid[name]=did
		self.didtoname[did]=name

		return did


	def add_task(self,task):
		"""Adds a new task to the active queue, scheduling it for execution. If parentid is
		specified, a doubly linked list is established. parentid MUST be the id of a task
		currently in the active queue. parentid and wait_for may be set in the task instead"""
		if not isinstance(task,EMTask) : raise Exception,"Invalid Task"
		#self.active["max"]+=1
		#tid=self.active["max"]

		EMTaskQueue.lock.acquire()
		try: tid=max(self.active["maxrec"],self.complete["maxrec"])+1
		except: tid=1
		task.taskid=tid
		task.queuetime=time.time()

		# map data file specifiers to ids
		for j,k in task.data.items():
			try:
				if k[0]!="cache" : continue
			except: continue
			try:
				did=self.todid(k[1])
			except:
				print "Invalid data item %s: %s"%(str(j),str(k))
				print str(task)
				os._exit(1)
			try: k[1]=did
			except:
				task.data[j]=list(k)
				task.data[j][1]=did

		self.active[tid]=task		# store the task in the queue
		try: EMTaskQueue.lock.release()
		except: print "Warning: lock re-released in add_task. Not serious, but shouldn't happen."
		return tid

	def task_progress(self,tid,percent):
		"""Update task progress, unless task is already complete/aborted. Returns True if progress
		update successful"""
		try :
			task=self.active[tid]
			if task==None : raise Exception
		except :
			try:
				task=self.complete[tid]
				return False
			except:
				print "ERROR: Progress, No such task : ",tid,percent
				return False
		task.progtime=(time.time(),percent)
		self.active[tid]=task
		return True

	def task_check(self,tid):
		"""This will check the status of a task. It will return -1 if a task is queued but not yet running,
		0-99.999 while running (% complete) or exactly 100 when the task is done"""
#		print "task_check ",tid
		try :
			task=self.active[tid]
			if task==None: raise Exception
		except:
			task=self.complete[tid]		# if we succeed in retrieving it from the complete list, it's done (or aborted)
			return 100

		if task.starttime==None or task.starttime<1 : return -1
		if task.progtime==None : return 0
		return task.progtime[1]

	def task_done(self, tid):
		"""Mark a Task as complete, by shifting a task to the tasks_complete queue"""
		EMTaskQueue.lock.acquire()
		try:
			task=self.active[tid]
			if task==None:
				print "*** Warning, task %d was already complete"%tid
				EMTaskQueue.lock.release()
				return
		except:
			EMTaskQueue.lock.release()
			return

		task.progtime=(time.time(),100)
		task.endtime=time.time()
		self.complete[tid]=task
		del self.active[tid]
		EMTaskQueue.lock.release()
		#if self.complete["max"]<tid : self.complete["max"]=tid
		#self.active["min"]=min(self.active.keys())


	def task_rerun(self,taskid):
		"""If a task has been started somewhere, but execution fails due to a problem with the target node, call
		this and the task will be returned to the queue, unless it has failed MAXTASKFAIL times."""
		try:
			task=self.active[taskid]
			if task==None: raise Exception
			cpl=False
		except:
			try :
				task=self.complete[taskid]
				cpl=True
			except:
				return

		if task==None :
			print "Warning: tried to requeue task ",taskid," but couldn't find it"
			return

		if task.failcount==MAXTASKFAIL :
			self.task_aborted(taskid)
			return

		task.failcount+=1
		task.starttime=None
		task.progtime=None
		task.endtime=None
		task.clientid=None
		task.exechost=None

		if cpl :
			print "Completed task %d requeued (%d failures)"%(taskid,task.failcount)
			del self.complete[taskid]

		self.active[taskid]=task

		return

	def task_aborted(self, taskid):
		"""Mark a Task as being aborted, by shifting a task to the tasks_complete queue"""
		EMTaskQueue.lock.acquire()
		try:
			task=self.active[taskid]
		except:
			EMTaskQueue.lock.release()
			return

		#if self.active["min"]==taskid : self.active["min"]=min(self.active.keys())
		self.complete[taskid]=task
		self.active[taskid]=None
		EMTaskQueue.lock.release()

class EMTask:
	"""This class represents a task to be completed. Generally it will be subclassed. This is effectively
	an abstract superclass to define common member variables and methods. Note that the data dictionary,
	which contains a mix of actual data elements, and ["cache",filename,#|min,max|(list)] image references, will
	be transparently remapped on the client to similar specifiers which are valid locally. Over the network
	such data requests are remapped into data identifiers (did), then translated back into valid filenames
	in the remote cache.  When subclassing this class, avoid defining new member variables, as EMTask objects
	get transmitted over the network. Make use of command, data and options instead. """
	def __init__(self,command=None,data=None,options=None,user=None,):
		self.taskid=None		# unique task identifier (in this directory)
		self.queuetime=None		# Time (as returned by time.time()) when task queued
		self.starttime=None		# Time when execution began
		self.progtime=None		# (time,% complete) from client
		self.endtime=None		# Time when task completed
		self.clientid=None		# id number of client where the job was executed
		self.exechost=None		# hostname where task was executed
		self.user=None			# Username from customer
		self.group=None			# group this task is in for task management purposes
		self.command=command	# This is a one word description of the purpose of the task, should be set in __init__
		self.data=data			# dictionary of named data specifiers value may be (before transmission):
								# - actual data object (no caching)
								# - ['cache',filename,#]
								# - ['cache',filename/url,min,max]  (max is exclusive, not inclusive)
								# - ['cache',filename/url,(list)]
								# (after transmission):
								# - actual data item
								# - ['cache',didEMD, did is a (modtime,int) tuple
		self.modtimes={}		# Used by MPI parallelism. A dictionary of the last modification times of files specified for caching in data
		self.options=options	# dictionary of options
		self.wait_for=None		# in the active queue, this identifies an exited class which needs to be rerun when all wait_for jobs are complete
		self.failcount=0		# Number of times this task failed to reach completion after starting
		self.errors=[]			# a list of errors (strings) that occured during task execution. Normally empty !
		self.ppid = os.getpid()		# Rrcords the PID to send to to children as PPID (sent to e2parallel.py)

	def execute(self): return

##########
### This is the 'database' object, representing a BerkeleyDB environment
##########

class EMAN2DB:
	"""This class implements a local database of data and metadata for EMAN2"""

	opendbs={}
	lock=threading.Lock()

	def open_db(path=None):
		"""This is an alternate constructor which may return a cached (already open)
		EMAN2DB instance"""
		# check the cache of opened dbs first
#		if not path : path=e2getcwd()
		EMAN2DB.lock.acquire()

		if not path : path=e2gethome()+"/.eman2"
		if path=="." or path=="./" : path=e2getcwd()
		if EMAN2DB.opendbs.has_key(path) :
			EMAN2DB.lock.release()
			return EMAN2DB.opendbs[path]
		ret=EMAN2DB(path)
		EMAN2DB.lock.release()
		return ret


	open_db=staticmethod(open_db)

	def __init__(self,path=None):
		"""path points to the directory containin the EMAN2DB subdirectory. None implies the current working directory"""
		global BDB_CACHE_DISABLE
		#if recover: xtraflags=db.DB_RECOVER
#		if not path : path=e2getcwd()
		if not path : path=e2gethome()+"/.eman2"
		if path=="." or path=="./" : path=e2getcwd()
		self.path=path

		# Keep a cache of opened database environments
		EMAN2DB.opendbs[self.path]=self

		# Make the database directory
		if not os.access("%s/EMAN2DB"%self.path,os.F_OK) :
			try:
				os.makedirs("%s/EMAN2DB"%self.path)
			except:
				# perhaps there is another process running that just made it?
				if not os.access("%s/EMAN2DB"%self.path,os.F_OK):
					raise RuntimeError("Error - there was a problem creating the EMAN2DB directory")

		# make the shared cache directory in /tmp
		if(sys.platform != 'win32'):
			if (not BDB_CACHE_DISABLE):
				if(not os.access("/tmp/eman2db-%s"%os.getenv("USER","anyone"),os.F_OK)):
					try: os.makedirs("/tmp/eman2db-%s"%os.getenv("USER","anyone"))
					except: pass
		else:
			if (not BDB_CACHE_DISABLE):
				if(not os.access("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"),os.F_OK)):
					try: os.makedirs("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"))
					except: pass

		if BDB_CACHE_DISABLE:
			self.dbenv=None
		else :
			self.dbenv=db.DBEnv()
			self.dbenv.set_cachesize(0,cachesize,4)		# gbytes, bytes, ncache (splits into groups)
#			self.dbenv.set_cachesize(1,0,8)		# gbytes, bytes, ncache (splits into groups)
			self.dbenv.set_data_dir("%s/EMAN2DB"%self.path)
			self.dbenv.set_lk_detect(db.DB_LOCK_DEFAULT)	# internal deadlock detection
			self.dbenv.set_lk_max_locks(20000)		# if we don't do this, we can easily run out when dealing with large numbers of files
			try: self.dbenv.set_lg_regionmax(5000000)
			except: print "Could not alter log region size. Please run e2bdb.py -c"
			self.dbenv.set_lk_max_objects(20000)

			try:
				if(sys.platform != 'win32'):
					self.dbenv.open("/tmp/eman2db-%s"%os.getenv("USER","anyone"),envopenflags)
				else:
					self.dbenv.open("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"),envopenflags)
			except:
				try:
					print """Cache open failed. Retyring one time."""
					# retry once
					time.sleep(2)
					if(sys.platform != 'win32'):
						self.dbenv.open("/tmp/eman2db-%s"%os.getenv("USER","anyone"),envopenflags)
					else:
						self.dbenv.open("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"),envopenflags)
				except:
					print "/tmp/eman2db-%s"%os.getenv("USER","anyone")
					traceback.print_exc()
					print """
========
ERROR OPENING DATABASE CACHE      (This most often occurs if you upgrade EMAN2 without running 'e2bdb.py -c' first. It could
also indicate that a program crashed in a bad way causing potential database corruption. You can try running 'e2bdb.py -c', and
see if that fixes the problem, otherwise you may need to 'rm -rf /tmp/eman2db-*' (or on windows remove the corresponding folder).
While there is a small possibility that this will prevent recovery of image files that were corrupted by the crash, it may be the
only practical option.)
"""
					os._exit(1)

		self.dicts={}
		#if self.__dbenv.DBfailchk(flags=0) :
			#self.LOG(1,"Database recovery required")
			#sys.exit(1)


	def close(self):
		"""close the environment associated with this object. This should ONLY be called if no other
		EMAN2 programs are currently running. Even then it is optional. see e2bdb.py --cleanup"""
		try: self.dbenv.close()
		except: pass
		self.dbenv=None

	def __del__(self):
		if not self.dbenv: return
		for i in self.dicts.keys() : self.close_dict(i)
		self.close()

	def __getitem__(self,key):
#		print "get ",key
		try: return self.dicts[key]
		except:
			# if the database is open read/write, return that instead
			if key[-4:]=="__ro" : return self.dicts[key[:-4]]
			raise KeyError

	def open_dict(self,name,ro=False):
#		print "open ",name,ro
		if self.dicts.has_key(name) : return
		self.dicts[name]=DBDict(name,dbenv=self.dbenv,path=self.path+"/EMAN2DB",parent=self,ro=ro)
		self.__dict__[name]=self.dicts[name]

	def close_dict(self,name):
		"this will close a dictionary"
		try:
			self.__dict__[name].close()
		except: pass

	def remove_dict(self,name):
		self.open_dict(name)
		self.__dict__[name].realopen()
		self.__dict__[name].bdb.truncate()
		self.__dict__[name].close()
#		try: os.unlink(self.path+"/EMAN2DB/"+name+".bdb")		# this is unsafe as far at the DB is concerned, but the method below freezes :^(
#		except: pass

		#if self.dicts.has_key(name) : self.dicts[name].close()
		#print self.path+"/EMAN2DB/"+name+".bdb"
		#self.dbenv.dbremove(self.path+"/EMAN2DB/"+name+".bdb")		# this is the 'correct' way, but using it seems to cause a process to lockup on FUTEX_WAIT  :^(

#if not (name in self.dicts.keys()) :
		#if name in self.dicts.keys():
			#self.__dict__[name].close()
		#try: os.unlink(self.path+"/EMAN2DB/"+name+".bdb")
		#except: pass
		for f in os.listdir(self.path+"/EMAN2DB"):
			if fnmatch.fnmatch(f, name+'_*'):
				try: os.unlink(self.path+"/EMAN2DB/"+f)
				except: pass

##########
### This represents a 'dictionary' within a 'database', in BerkeleyDB parlance this is a B-tree (though it could also be a Hash)
##########

class DBDict:
	"""This class uses BerkeleyDB to create an object much like a persistent Python Dictionary,
	keys and data may be arbitrary pickleable types, however, additional functionality is provided
	for EMData objects. Specifically, if integer keys are used, set_attr and get_attr may be used
	to efficiently get and set attributes for images with reduced i/o requirements (in certain cases)."""

	alldicts=weakref.WeakKeyDictionary()
	nopen=0
	closelock=threading.Lock()
	fixedkeys=frozenset(("nx","ny","nz","minimum","maximum","mean","sigma","square_sum","mean_nonzero","sigma_nonzero"))
	def __init__(self,name,filename=None,dbenv=None,path=None,parent=None,ro=False):
		"""This is a persistent dictionary implemented as a BerkeleyDB Hash
		name is required, and will also be used as a filename if none is
		specified. Note that the database is not actually opened until it's used."""

		from EMAN2 import e2getcwd

		global dbopenflags
		DBDict.alldicts[self]=1		# we keep a running list of all trees so we can close everything properly
		self.name = name
		self.parent=parent
		self.dbenv=dbenv
		self.lock=threading.Lock()	# used to reduce thread conflicts
		self.file=filename
		self.rohint=ro
		self.lasttime=time.time()		# last time the database was accessed
		self.opencount=0				# number of times the database has needed reopening
		if path : self.path = path
		else : self.path=e2getcwd()
		self.txn=None	# current transaction used for all database operations
		self.bdb=None
		self.isro=False


	def __str__(self): return "<EMAN2db DBDict instance: %s>" % self.name

	def __del__(self):
		self.close()

	def updateold(self,lfile,ro=False):
		"""Called to update old 4.2 databases with funny problem"""
		self.bdb=db.DB(self.dbenv)
		print "Old format DB detected (%s). Do not be alarmed, this must be done 1 time for each database file. Please wait."%lfile
		try: os.unlink(self.path+"/"+lfile.replace(".bdb",".old"))
		except: pass
		os.rename(self.path+"/"+lfile,self.path+"/"+lfile.replace(".bdb",".old"))
		try:
			tmpdb=db.DB()
			tmpdb.open(self.path+"/"+lfile.replace(".bdb",".old"),self.name)
		except:
			print "Error updating %s. Please contact sludtke@bcm.edu."%lfile
			os._exit(1)
		try:
			self.bdb.open(self.path+"/"+lfile,self.name,db.DB_BTREE,db.DB_CREATE|db.DB_THREAD)
		except:
			print "Error 2 updating %s. Please contact sludtke@bcm.edu."%lfile
			os._exit(1)

		for k in tmpdb.keys():
			self.bdb[k]=tmpdb[k]

		tmpdb.close()
		print "Conversion complete. (-old file retained as an emergency backup, safe to remove)"

	def realopen(self,ro=False):
		"""This actually opens the database (unless already open), if ro is set and the database is not already
		open read-write, it will be opened read-only"""

#		print "open ",self.name

		global DBDEBUG
		if DBDEBUG:
			while not self.lock.acquire(False) :
				print"DB %s locked. Waiting"%self.name
				time.sleep(1)
		else : self.lock.acquire()
		self.lasttime=time.time()
		if self.bdb!=None :
			if ro==True or self.isro==False :
				if DBDEBUG : print "already open",self.name
				self.lock.release()
				return  		# return if the database is already open and in a compatible read-only mode
			if DBDEBUG : print "reopening R/W ",self.name
			self.lock.release()
			self.close()	# we need to reopen read-write
			self.lock.acquire()

		#if DBDEBUG:
			## look at the locking subsystem stats
			#ls=self.dbenv.lock_stat()
			#print "lock_stat:\t",
			#for i in ("nlocks","maxlocks","nlockers","maxlockers","nobjects","maxobjects","maxnlocks"): print ls[i],"\t",
			#print


		self.bdb=db.DB(self.dbenv)		# we don't check BDB_CACHE_DISABLE here, since self.dbenv will already be None if its set
		if self.file==None : lfile=self.name+".bdb"
		else : lfile=self.file
#		print "open ",self.path+"/"+file,self.name,ro
		if ro :
			try:
				self.bdb.open(self.path+"/"+lfile,self.name,db.DB_BTREE,db.DB_RDONLY|db.DB_THREAD)
			except db.DBInvalidArgError:
				self.updateold(lfile,ro)
			except db.DBNoSuchFileError:
				self.bdb=None
				self.lock.release()
				if DBDEBUG : traceback.print_exc()
				raise Exception, "Cannot open or find %s"%self.name

			#except:
				## try one more time... this shouldn't be necessary...
				#time.sleep(1)
##				try:
				#self.bdb.open(self.path+"/"+file,self.name,db.DB_BTREE,db.DB_RDONLY|db.DB_THREAD)
##				except:
##					raise Exception,"Cannot open database : %s"%self.path+"/"+file
			self.isro=True
		else :
			try:
				self.bdb.open(self.path+"/"+lfile,self.name,db.DB_BTREE,db.DB_CREATE|db.DB_THREAD)
			except db.DBInvalidArgError:
				self.updateold(lfile,ro)
			except:
				try:
					os.makedirs("%s/EMAN2DB"%self.path)
					self.bdb=db.DB(self.dbenv)
					self.bdb.open(self.path+"/"+lfile,self.name,db.DB_BTREE,db.DB_CREATE|db.DB_THREAD)
				except :
					self.bdb=None
					self.lock.release()
					traceback.print_exc()
					print "Unable to open read/write %s (%s/%s)"%(self.name,self.path,lfile)
					return
			#except:
				## try one more time... this shouldn't be necessary...
				#time.sleep(1)
				#try:
					#self.bdb.open(self.path+"/"+file,self.name,db.DB_BTREE,db.DB_CREATE|db.DB_THREAD)
				#except:
					#raise Exception,"Cannot create database : %s"%self.path+"/"+file
			self.isro=False

		self.opencount+=1	# how many times we have had to reopen this database
		DBDict.nopen+=1
		self.lock.release()

		global MAXOPEN
		if DBDict.nopen>MAXOPEN : self.close_one()

		if DBDEBUG : print "Opened ",self.name

#		print "%d open"%DBDict.nopen
#		print "opened ",self.name,ro
#		self.bdb.open(file,name,db.DB_HASH,dbopenflags)
#		print "Init ",name,file,path

	def close_one(self):
		"""Will select and close any excess open databases. Closure is based on the number of times it has been repoened and the
		time it was last used."""
		if not DBDict.closelock.acquire(False) : return

		global MAXOPEN
#		l=[(i.opencount,i.lasttime,i) for i in self.alldicts if i.bdb!=None]		# list of all open databases and usage,time info
		l=[(i.lasttime,i.opencount,i) for i in self.alldicts if i.bdb!=None]		# sort by time
		l.sort()

		global DBDEBUG

#		if len(l)>MAXOPEN :
#			print "%d dbs open, autoclose disabled"%len(l)
#			for j,i in enumerate(l): print j,i[2].name,i[0]-time.time(),i[1]
		if len(l)>MAXOPEN :
			if DBDEBUG: print "DB autoclosing %d/%d "%(len(l)-MAXOPEN,len(l))
			for i in range(len(l)-MAXOPEN):
				if DBDEBUG: print "CLOSE:",l[i][2].name
				if (l[i][2]!=self) :
					l[i][2].close()
				else :
					print "Warning: attempt to autoclose the DB we just opened, please report this"

		DBDict.closelock.release()

	def forceclose(self):
		global DBDEBUG

		# There should be no lock to acquire if the bdb is None
		try :
			if self.bdb == None:
				return
		except: return

		if not self.lock.acquire(False):
			time.sleep(0.1)


#		print "close x ",self.path+"/"+str(self.file),self.name,"XXX"
		try: self.bdb.close()
		except: pass
		self.bdb=None
		DBDict.nopen-=1
		self.lock.release()
		if DBDEBUG : print "Closed ",self.name


	def close(self):
		global DBDEBUG
		n=0
		while not self.lock.acquire(False) and n<3:
			print "Sleep on close ",self.name
			time.sleep(.2)
			n+=1
		if n>=4 : return	# failed too many times, just return and let things fail where they may...

		if self.bdb == None:
			self.lock.release()
			return
#		print "close x ",self.path+"/"+str(self.file),self.name,"XXX"
		self.bdb.close()
		self.bdb=None
		DBDict.nopen-=1
		self.lock.release()
		if DBDEBUG : print "Closed ",self.name

	def sync(self):
		if self.bdb!=None : self.bdb.sync()

	def set_txn(self,txn):
		"""sets the current transaction. Note that other python threads will not be able to use this
		Hash until it is 'released' by setting the txn back to None"""
		if txn==None:
			self.txn=None
			return

		while self.txn :
			time.sleep(.1)
		self.txn=txn

	def get_attr(self,n,attr):
		"""Returns an attribute or set of attributes for an image or set of images. n may be a single key or a list/tuple/set of keys,
		and attr may be a single attribute or a list/tuple/set. Returns the attribute, a dict of attributes or a image keyed dict of dicts keyed by attribute"""
		self.realopen(self.rohint)	# make sure the database is open
		try :
			ret={}
			for i in n:
				d=loads(self.bdb.get(dumps(i,-1),txn=self.txn))
				if getattr(attr, '__iter__', False):
					ret[i]={}
					for a in attr:
						if a in d : ret[i][a]=d[a]
				else:
					try: ret[i]=d[attr]
					except: pass
			return ret
		except:
			if getattr(attr, '__iter__', False):
				d=loads(self.bdb.get(dumps(n,-1),txn=self.txn))
				ret={}
				for a in attr:
					if a in d : ret[a]=d[a]
				return ret
			return loads(self.bdb.get(dumps(n,-1),txn=self.txn))[attr]

	def set_attr(self,n,attr,val=None):
		"""Sets an attribute to val in EMData object 'n'. Alternatively, attr may be a dictionary containing multiple key/value pairs
		to be updated in the EMData object. Unlike with get_attr, n must always refer to a single EMData object in the database."""
		self.realopen()
		a=loads(self.bdb.get(dumps(n,-1),txn=self.txn))
		if isinstance(attr,dict) :
			a.update(attr)
		else: a[attr]=val
		self[n]=a

	def __len__(self):
		self.realopen(self.rohint)
		try: return self["maxrec"]+1
		except: return 0
#		return self.bdb.stat(db.DB_FAST_STAT)["nkeys"]
#		return len(self.bdb)

	def put(self,key,val,txn=None):
		"""This performs the bdb.put function with some error detection and retry capabilites.
Retrying should not be necessary, but we have been completely unable to figure out the cause
of these occasional errors"""
		n=0
		while n<10:
			try:
				self.bdb.put(key,val,txn=txn)
				break
			except:
				if n in (0,9) : traceback.print_exc()
				print "********** Warning: problem writing ",key," to ",self.name,". Retrying (%d/10)"%n
				time.sleep(5)
				n+=1


	def __setitem__(self,key,val):
		self.realopen()
		if (val==None) :
			try:
				self.__delitem__(key)
			except: return
		elif isinstance(val,EMData) :
			# decide where to put the binary data
			ad=val.get_attr_dict()
			pkey="%s/%s_"%(self.path,self.name)
			fkey="%dx%dx%d"%(ad["nx"],ad["ny"],ad["nz"])
#			print "w",fkey
			try :
				n=loads(self.bdb.get(fkey+dumps(key,-1),txn=self.txn))
			except:
				if not self.has_key(fkey) : self[fkey]=0
				else: self[fkey]+=1
				n=self[fkey]
			self.put(fkey+dumps(key,-1),dumps(n,-1),txn=self.txn)		# a special key for the binary location

			# write the metadata
			try: del ad["data_path"]
			except:pass

			t=time.localtime(time.time())
			ad["timestamp"]="%04d/%02d/%02d %02d:%02d:%02d"%t[:6]
			self.put(dumps(key,-1),dumps(ad,-1),txn=self.txn)

			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) :
				self["maxrec"]=key

				# Updates the image count cache
				db2=db_open_dict("bdb:%s#%s"%(self.path[:-8],"00image_counts"))
				try:
					im=self[0]
					sz=(im["nx"],im["ny"],im["nz"])
				except : sz=(0,0,0)
				db2[self.name]=(time.time(),len(self),sz)
			else:
				# Updates the count cache time, to show it's up to date
				try:
					db2=db_open_dict("bdb:%s#%s"%(self.path[:-8],"00image_counts"))
					db2[self.name]=(time.time(),db2[self.name][1],db2[self.name][2])
				except : pass

			# write the binary data
			val.write_data(pkey+fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"])

		else :
			self.put(dumps(key,-1),dumps(val,-1),txn=self.txn)
			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) : self["maxrec"]=key

	def __getitem__(self,key):
		self.realopen(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1),txn=self.txn))
		except: return None
		if isinstance(r,dict) and r.has_key("is_complex_x") :
			pkey="%s/%s_"%(self.path,self.name)
			fkey="%dx%dx%d"%(r["nx"],r["ny"],r["nz"])
#			print "r",fkey
			ret=EMData(r["nx"],r["ny"],r["nz"])
			if r.has_key("data_path"):
				p,l=r["data_path"].split("*")
#				print "read ",os.getcwd(),self.path,p,l
				if p[0]=='/' : ret.read_data(p,int(l))
				else : ret.read_data(self.path+"/"+p,int(l))
			else:
				try: n=loads(self.bdb.get(fkey+dumps(key,-1)))	 # this is the index for this binary data item in the image-dimensions-specific binary data file
				except: raise KeyError,"Undefined data location key for : %s"%key
				ret.read_data(pkey+fkey,n*4*r["nx"]*r["ny"]*r["nz"])
			k=set(r.keys())
			k-=DBDict.fixedkeys

			for i in k:
				if i not in ('nx', 'ny', 'nz'):
					ret.set_attr(i,r[i])

			ret["source_path"]="bdb:"+pkey[:-1].replace("/EMAN2DB/","#")
			ret["source_n"]=key
			return ret
		return r

	def __delitem__(self,key):
		self.realopen()
		self.bdb.delete(dumps(key,-1),txn=self.txn)

	def __contains__(self,key):
		self.realopen(self.rohint)
		return self.bdb.has_key(dumps(key,-1))

	def item_type(self,key):
		self.realopen(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1),txn=self.txn))
		except: return None
		if isinstance(r,dict) and r.has_key("is_complex_x") : return EMData
		return type(r)

	def keys(self):
		self.realopen(self.rohint)
		try: return [loads(x) for x in self.bdb.keys() if x[0]=='\x80']
		except:
			traceback.print_exc()
			print "This is a serious error, which should never occur during normal usage.\n Please report it (with the error text above) to sludtke@bcm.edu. Please also read the database warning page in the Wiki."
			sys.exit(1)

	def values(self):
		self.realopen(self.rohint)
		return [self[k] for k in self.keys()]

	def items(self):
		self.realopen(self.rohint)
		return [(k,self[k]) for k in self.keys()]

	def has_key(self,key):
		self.realopen(self.rohint)
		return self.bdb.has_key(dumps(key,-1))

	def get_data_path(self,key):
		"""returns the path to the binary data as "path*location". Only valid for EMData objects."""
		self.realopen(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1)))
		except: return None
		if isinstance(r,dict) and r.has_key("is_complex_x") :
			pkey="%s/%s_"%(self.path,self.name)
			fkey="%dx%dx%d"%(r["nx"],r["ny"],r["nz"])
			if r.has_key("data_path"):
				if r["data_path"][0]=="/" : return r["data_path"]
				return self.path+"/"+r["data_path"]
			else :
				n=loads(self.bdb.get(fkey+dumps(key,-1)))
				return "%s*%d"%(pkey+fkey,n*4*r["nx"]*r["ny"]*r["nz"])
		return None

	def get(self,key,dfl=None,txn=None,target=None,nodata=0,region=None,idx=0):
		"""Alternate method for retrieving records. Permits specification of an EMData 'target'
		object in which to place the read object"""
		self.realopen(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1),txn=txn))
		except: return dfl
		if isinstance(r,dict) and r.has_key("is_complex_x") :
			pkey="%s/%s_"%(self.path,self.name)
			rnx,rny,rnz = r["nx"],r["ny"],r["nz"]
			fkey="%dx%dx%d"%(rnx,rny,rnz)

#			print "r",fkey
			if region != None:
				size = region.get_size()
				# zeros are annoyingly necessary
				for i in range(len(size)):
					if size[i] == 0: size[i] = 1
				nx,ny,nz = int(size[0]),int(size[1]),int(size[2])
			else:
				nx,ny,nz = rnx,rny,rnz

			if target : ret = target
			else: ret = EMData()

			# metadata
			k=set(r.keys())
#			k-=DBDict.fixedkeys
			for i in k:
				if i not in ('nx', 'ny', 'nz'):
					ret.set_attr(i,r[i])
			ret["source_path"]="bdb:"+pkey[:-1].replace("/EMAN2DB/","#")
			ret["source_n"]=key

			# binary data
			ret.set_size(nx,ny,nz,nodata)
			if not nodata:

				if region != None: ret.to_zero() # this has to occur in situations where the clip region goes outside the image
				if r.has_key("data_path"):
					p,l=r["data_path"].split("*")
					if p[0]=='/' or p[0]=='\\' or p[1]==':': ret.read_data(p,int(l),region,rnx,rny,rnz)		# absolute path
					else :
						ret.read_data(self.path+"/"+p,int(l),region,rnx,rny,rnz) 		# relative path
				else:
					try: n=loads(self.bdb.get(fkey+dumps(key,-1)))	 # this is the index for this binary data item in the image-dimensions-specific binary data file
					except: raise KeyError,"Undefined data location key %s for %s"%(key,pkey+fkey)
					try: ret.read_data(pkey+fkey,n*4*rnx*rny*rnz,region,rnx,rny,rnz)	# note that this uses n, NOT 'key'. Images cannot be located in the binary file based on their numerical key
					except :
						import socket
						print "Data read error (%s) on %s (%d)"%(socket.gethostname(),pkey+fkey,key*4*rnx*rny*rnz)
						traceback.print_exc()
						sys.stderr.flush()
						sys.stdout.flush()
						os._exit(1)
			return ret
		return r

	def get_header(self,key,txn=None,target=None):
		"""Alternate method for retrieving metadata for EMData records."""
		self.realopen(self.rohint)
		try: return loads(self.bdb.get(dumps(key,-1),txn=txn))
		except: return None

	def set(self,val,key,region=None,txn=None):
		'''
		I have to support the image_type and read_header parameters even though they are not used -
		this is so this function can be used interchangeably with EMData.write_image

		Alternative to x[key]=val with transaction set'''
		self.realopen()
		if (val==None) :
			try:
				self.__delitem__(key)
			except: return
		elif isinstance(val,EMData) :
			# decide where to put the binary data
			if region:
				ad = self.get_header(key)
				#except: raise RuntimeError("If you're using a region the file must already exist")
			else:
				ad=val.get_attr_dict()

			pkey="%s/%s_"%(self.path,self.name)
			fkey="%dx%dx%d"%(ad["nx"],ad["ny"],ad["nz"])
#			print "w",fkey
			try :
				n=loads(self.bdb.get(fkey+dumps(key,-1),txn=txn))
			except:
				if not self.has_key(fkey) :
					self[fkey]=0
				else: self[fkey]+=1
				n=self[fkey]
			self.put(fkey+dumps(key,-1),dumps(n,-1),txn=txn)		# a special key for the binary location

			# write the metadata
			try: del ad["data_path"]
			except:pass

			t=time.localtime(time.time())
			ad["timestamp"]="%04d/%02d/%02d %02d:%02d:%02d"%t[:6]
			self.put(dumps(key,-1),dumps(ad,-1),txn=txn)

			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) :
				self["maxrec"]=key

				# update the image count cache
				db2=db_open_dict("bdb:%s#%s"%(self.path[:-8],"00image_counts"))
				try:
					im=self[0]
					sz=(im["nx"],im["ny"],im["nz"])
				except : sz=(0,0,0)
				db2[self.name]=(time.time(),len(self),sz)
			else:
				# Updates the count cache time, to show it's up to date
				try:
					db2=db_open_dict("bdb:%s#%s"%(self.path[:-8],"00image_counts"))
					db2[self.name]=(time.time(),db2[self.name][1],db2[self.name][2])
				except : pass

			# write the binary data
			if region:
				val.write_data(pkey+fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"],region,ad["nx"],ad["ny"],ad["nz"])
			else:
				val.write_data(pkey+fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"])
#			print "WI ",self.path,self.name,key,val,BDB_CACHE_DISABLE,pkey+fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"]

		else :
			self.put(dumps(key,-1),dumps(val,-1),txn=txn)
			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) : self["maxrec"]=key

	def set_header(self,key,val,txn=None):
		"Alternative to x[key]=val with transaction set"
		self.realopen()
		# make sure the object exists and is an EMData object
		try: r=loads(self.bdb.get(dumps(key,-1),txn=txn))
		except: raise Exception,"set_header can only be used to update existing EMData objects"
		if not isinstance(r,dict) or not r.has_key("is_complex_ri") :
			raise Exception,"set_header can only be used to update existing EMData objects"

		if (val==None) : raise Exception,"You cannot delete an EMData object header"
		elif isinstance(val,EMData) :
			# write the metadata
			ad=val.get_attr_dict()
			self.bdb.put(dumps(key,-1),dumps(ad,-1),txn=txn)
		elif isinstance(val,dict) :
			self.bdb.put(dumps(key,-1),dumps(val,-1),txn=txn)
		else : raise Exception,"set_header is only valid for EMData objects or dictionaries"


	def update(self,dict):
		self.realopen()
		for i,j in dict.items(): self[i]=j


#def DB_cleanup():
	#"""This does at_exit cleanup. It would be nice if this were always called, but if python is killed
	#with a signal, it isn't. This tries to nicely close everything in the database so no recovery is
	#necessary at the next restart"""
	#sys.stdout.flush()
	#print >>sys.stderr, "Closing %d BDB databases"%(len(BTree.alltrees)+len(IntBTree.alltrees)+len(FieldBTree.alltrees))
	#if DEBUG>2: print >>sys.stderr, len(BTree.alltrees), 'BTrees'
	#for i in BTree.alltrees.keys():
		#if DEBUG>2: sys.stderr.write('closing %s\n' % str(i))
		#i.close()
		#if DEBUG>2: sys.stderr.write('%s closed\n' % str(i))
	#if DEBUG>2: print >>sys.stderr, '\n', len(IntBTree.alltrees), 'IntBTrees'
	#for i in IntBTree.alltrees.keys():
		#i.close()
		#if DEBUG>2: sys.stderr.write('.')
	#if DEBUG>2: print >>sys.stderr, '\n', len(FieldBTree.alltrees), 'FieldBTrees'
	#for i in FieldBTree.alltrees.keys():
		#i.close()
		#if DEBUG>2: sys.stderr.write('.')
	#if DEBUG>2: sys.stderr.write('\n')
## This rmakes sure the database gets closed properly at exit
#atexit.register(DB_cleanup)


__doc__ = \
"""This module supports the concept of a local database for storing data and
metadata associated with a particular EMAN2 refinement. Data stored in this
database may be extracted into standard flat-files, but use of a database
with standard naming conventions, etc. helps provide the capability to log
the entire refinement process."""




