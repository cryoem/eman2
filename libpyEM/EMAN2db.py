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

from EMAN2 import *
import atexit
import weakref
from cPickle import loads,dumps
from zlib import compress,decompress
import os
import signal
import sys
import fnmatch
import random
import traceback
try:
	from bsddb3 import db
except:
	from bsddb import db

try:
	a=frozenset()
except:
	from sets import Set
	frozenset=Set

# Flags used to open the database environment
envopenflags=db.DB_CREATE|db.DB_INIT_MPOOL|db.DB_INIT_LOCK|db.DB_INIT_LOG|db.DB_THREAD
#envopenflags=db.DB_CREATE|db.DB_INIT_MPOOL|db.DB_INIT_LOCK|db.DB_INIT_LOG|db.DB_THREAD|db.DB_REGISTER|db.DB_RECOVER

# If set, databases will be opened without any caching, locking, etc.
# DANGER, unless DB is single-threaded and local, don't do any writing
MPIMODE=0

# This hardcoded value is the maximum number of DBDict objects which will be simultaneously open
# when more than this are open, we begin closing the oldest ones. They will be reopened on-demand,
# but this will prevent resource exhaustion
MAXOPEN=20

# maximum number of times a task will be restarted before we assume there is something fundamentally
# flawed about the task itself
MAXTASKFAIL=10

cachesize=80000000


def DB_cleanup(a1=None,a2=None):
	if a1 in (2,15) :
		print "Program interrupted, closing databases, please wait (%d)"%os.getpid() 
	for d in DBDict.alldicts.keys(): d.close()
	for e in EMAN2DB.opendbs.values(): e.close()
	if a1 in (2,15) :
		print "Databases closed, exiting" 
		sys.exit(1)

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
		if url[0][0]!='/' : url[0]=e2getcwd()+"/"+url[0]	# relative path
		if url[0][:3]=="../": url[0]=e2getcwd().rsplit("/",1)[0]+"/"+url[0]	# find cwd if we start with ../
		if len(u1)==1 : return(url[0],url[1],None)			# bdb:path/to#dictname
		url=[url[0]]+url[1].split("?")						# bdb:path/to#dictname?keys
	
	# if we're here we need to deal with ?keys
	u2=url[2].split(",")
	if len(u2)==1 :
		if u2[0][:7]=="select." :
			ddb=EMAN2DB.open_db(".")
			ddb.open_dict("select")
			if not ddb.select.has_key(u2[0][7:]) : raise Exception,"Unknown selection list %s"%u2[0][7:]
			return (url[0],url[1],ddb.select[u2[0][7:]])		# bdb:path/to#dict?select/name
		elif u2[0][:8]=="exclude." :
			ddb=EMAN2DB.open_db(".")
			ddb.open_dict("select")
			if not ddb.select.has_key(u2[0][8:]) : raise Exception,"Unknown selection list %s"%u2[0][8:]
			exc_set = set(ddb.select[u2[0][8:]])
			all_set = set(range(0,EMUtil.get_image_count("bdb:"+url[0]+"#"+url[1])))
			rem_set = all_set - exc_set
			return (url[0],url[1],list(rem_set))		# bdb:path/to#dict?select/name
		return url											# bdb:path/to#dict?onekey
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
	after closing."""

	path,dictname,keys=db_parse_path(url)
	
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
	
	ld=[i[:-4] for i in ld if i[-4:]==".bdb"]
	
	return ld
	
##########
#### replace a few EMData methods with python versions to intercept 'bdb:' filenames
##########
def db_emd_init(self,*parms):
	if len(parms)<4 and len(parms)>0 and isinstance(parms[0],str) and parms[0][:4].lower()=="bdb:":
		self.__initc()
		self.read_image(*parms)
		return
	else : 
#		print "toC:", parms
		self.__initc(*parms)	
		return

EMData.__initc=EMData.__init__
EMData.__init__=db_emd_init

def transform_to_str(self):
	# I added the initial return so that it looks good in the eman2 file browser
	s = "\n%.3f %.3f %.3f %.3f\n" %(self.at(0,0),self.at(0,1),self.at(0,2),self.at(0,3))
	s += "%.3f %.3f %.3f %.3f\n" %(self.at(1,0),self.at(1,1),self.at(1,2),self.at(1,3))
	s += "%.3f %.3f %.3f %.3f\n" %(self.at(2,0),self.at(2,1),self.at(2,2),self.at(2,3))
	s += "0.000 0.000 0.000 1.000"
	return s

Transform.__str__ = transform_to_str

def db_read_image(self,fsp,*parms):
	"""read_image(filespec,image #,[header only],[region],[is_3d])
	
	This function can be used to read a set of images from a file or bdb: database. Pass in
	the filename, or bdb: specification, optionally a list of image numbers to read (or None),
	and a flag indicating that only the image headers should be read in. If only the headers
	are read, accesses to the image data in the resulting EMData objects will be invalid.
	Region reading is not supported for bdb:entries yet."""
	if fsp[:4].lower()=="bdb:" :
		db,keys=db_open_dict(fsp,True,True)
		 
		if len(parms)>1 and parms[1] : nodata=1
		else: nodata=0
		
		if len(parms)>2 and parms[2] : region=parms[2]
		else: region = None
		
		if keys:
			key=keys[parms[0]]
		else: key=parms[0]
#		try :
		x=db.get(key,target=self,nodata=nodata,region=region)
#		except: 
#			raise Exception("Could not access "+str(fsp)+" "+str(key))
		return None
	return self.read_image_c(fsp,*parms)

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
	return EMData.read_images_c(fsp,*parms)

EMData.read_images_c=staticmethod(EMData.read_images)
EMData.read_images=staticmethod(db_read_images)


def db_write_image(self,fsp,*parms):
	"""write_image(fsp,image #,[image type],[header only],[region],[storage type],[use host endian])

	Writes and image to a file or a bdb: entry. Note that for bdb: specifications, only image # is supported.
	the remaining options are ignored"
	"""
	if fsp[:4].lower()=="bdb:" :
		db,keys=db_open_dict(fsp,False,True)
		if keys :			# if the user specifies the key in fsp, we ignore parms
			if len(keys)>1 : raise Exception,"Too many keys provided in write_image %s"%str(keys)
			if isinstance(keys[0],int) and keys[0]<0 : raise Exception,"Negative integer keys not allowed %d"%keys[0]
			db[keys[0]]=self
			return
		if parms[0]<0 : parms=(len(db),)+parms[1:]
		db[parms[0]]=self
		return
	return self.write_image_c(fsp,*parms)

EMData.write_image_c=EMData.write_image
EMData.write_image=db_write_image

def db_get_image_count(fsp):
	"""get_image_count(path)

Takes a path or bdb: specifier and returns the number of images in the referenced stack."""
	if fsp[:4].lower()=="bdb:" :
		db,keys=db_open_dict(fsp,True,True)
		if keys :			# if the user specifies the key in fsp, we ignore parms
			n=0
			for i in keys:
				if i in db : n+=1
			return n
		return len(db)
	return EMUtil.get_image_count_c(fsp)


EMUtil.get_image_count_c=staticmethod(EMUtil.get_image_count)
EMUtil.get_image_count=staticmethod(db_get_image_count)


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
	
	def __init__(self,path=None,ro=False):
		"""path should point to the directory where the disk-based task queue will reside without bdb:"""
		if path==None or len(path)==0 : path="."
		self.path=path
		self.active=db_open_dict("bdb:%s#tasks_active"%path,ro)		# active tasks keyed by id
		self.complete=db_open_dict("bdb:%s#tasks_complete"%path,ro)	# complete tasks
		self.nametodid=db_open_dict("bdb:%s#tasks_name2did"%path,ro)	# map local data filenames to did codes
		self.didtoname=db_open_dict("bdb:%s#tasks_did2name"%path,ro)	# map data id to local filename

		#if not self.active.has_key("max") : self.active["max"]=-1
		#if not self.active.has_key("min") : self.active["min"]=0
		#if not self.complete.has_key("max") : self.complete["max"]=0
	
	def __len__(self) : return len(self.active)
	
	def get_task(self,clientid=0):
		"""This will return the next task waiting for execution"""
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
				return task 
			
		return None
	
	def add_group(self):
		"""returns a new (unique) group id to be used for related tasks"""
		try: ret=self.active["grpctr"]+1
		except: ret=1
		self.active["grpctr"]=ret
		
		return ret
	
	def add_task(self,task):
		"""Adds a new task to the active queue, scheduling it for execution. If parentid is
		specified, a doubly linked list is established. parentid MUST be the id of a task
		currently in the active queue. parentid and wait_for may be set in the task instead"""
		if not isinstance(task,EMTask) : raise Exception,"Invalid Task"
		#self.active["max"]+=1
		#tid=self.active["max"]
		try: tid=max(self.active["maxrec"],self.complete["maxrec"])+1
		except: tid=1
		task.taskid=tid
		task.queuetime=time.time()

		# map data file specifiers to ids
		for j,k in task.data.items():
			try: 
				if k[0]!="cache" : continue
			except: continue
			
			fmt=e2filemodtime(k[1])
			try : 
				did=self.nametodid[k[1]]			# get the existing did from the cache (modtime,int)
				if fmt!=did[0]	:					# if the file has been changed, we need to assign a new did
					del self.didtoname[did]
					raise Exception
			except: 
				did=(fmt,random.randint(0,999999))	# since there may be multiple files with the same timestamp, we also use a random int
				while (self.didtoname.has_key(did)):
					did=(fmt,random.randint(0,999999))
			
			self.nametodid[k[1]]=did
			self.didtoname[did]=k[1]
			try: k[1]=did
			except:
				task.data[j]=list(k)
				task.data[j][1]=did

		self.active[tid]=task		# store the task in the queue
		return tid
		
	def task_progress(self,taskid,percent):
		"""Update task progress, unless task is already complete/aborted. Returns True if progress
		update successful"""
		try : task=self.active[tid]
		except :
			task=self.complete[tid]
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
		try:
			task=self.active[tid]
			if task==None:
				print "*** Warning, task %d was already complete"%tid
				return
		except:
			return
		
		task.endtime=time.time()
		self.complete[tid]=task
		del self.active[tid]
		#if self.complete["max"]<tid : self.complete["max"]=tid
		#self.active["min"]=min(self.active.keys())
		

	def task_rerun(self,taskid):
		"""If a task has been started somewhere, but execution fails due to a problem with the target node, call
		this and the task will be returned to the queue, unless it has failed MAXTASKFAIL times."""
		try:
			task=self.active[tid]
		except:
			return

		if task.failcount==MAXTASKFAIL :
			self.task_aborted(taskid)
			return
		
		task.failcount+=1
		task.starttime=None
		return
		
	def task_aborted(self, taskid):
		"""Mark a Task as being aborted, by shifting a task to the tasks_complete queue"""
		try:
			task=self.active[tid]
		except:
			return
		
		#if self.active["min"]==taskid : self.active["min"]=min(self.active.keys())
		self.complete[tid]=task
		self.active[tid]=None
		

		
	
class EMTask:
	"""This class represents a task to be completed. Generally it will be subclassed. This is effectively
	an abstract superclass to define common member variables and methods. Note that the data dictionary,
	which contains a mix of actual data elements, and ["cache",filename,#|min,max|(list)] image references, will
	be transparently remapped on the client to similar specifiers which are valid locally. Over the network
	such data requests are remapped into data identifiers (did), then translated back into valid filenames
	in the remote cache.  When subclassing this class, avoid defining new member variables, as EMTask objects
	get transmitted over the network. Make use of command, data and options instead. """
	def __init__(self,user=None,data=None,options=None):
		self.taskid=None		# unique task identifier (in this directory)
		self.queuetime=None		# Time (as returned by time.time()) when task queued
		self.starttime=None		# Time when execution began
		self.progtime=None		# (time,% complete) from client
		self.endtime=None		# Time when task completed
		self.clientid=None		# id number of client where the job was executed
		self.exechost=None		# hostname where task was executed
		self.user=None			# Username from customer
		self.group=None			# group this task is in for task management purposes
		self.command="Nothing"	# This is a one word description of the purpose of the task, should be set in __init__
		self.data=data			# dictionary of named data specifiers value may be (before transmission):
								# - actual data object (no caching)
								# - ['cache',filename,#]
								# - ['cache',filename/url,min,max]  (max is exclusive, not inclusive)
								# - ['cache',filename/url,(list)]
								# (after transmission):
								# - actual data item
								# - ['cache',didEMD, did is a (modtime,int) tuple
		self.options=options	# dictionary of options
		self.wait_for=None		# in the active queue, this identifies an exited class which needs to be rerun when all wait_for jobs are complete
		self.failcount=0		# Number of times this task failed to reach completion after starting
		self.errors=[]			# a list of errors (strings) that occured during task execution. Normally empty !

	def execute(self): return

##########
### This is the 'database' object, representing a BerkeleyDB environment
##########

class EMAN2DB:
	"""This class implements a local database of data and metadata for EMAN2"""
	
	opendbs={}
	
	def open_db(path=None):
		"""This is an alternate constructor which may return a cached (already open)
		EMAN2DB instance"""
		# check the cache of opened dbs first
#		if not path : path=e2getcwd()
		if not path : path=e2gethome()+"/.eman2"
		if path=="." or path=="./" : path=e2getcwd()
		if EMAN2DB.opendbs.has_key(path) : return EMAN2DB.opendbs[path]
		return EMAN2DB(path)
	
	open_db=staticmethod(open_db)
	
	def __init__(self,path=None):
		"""path points to the directory containin the EMAN2DB subdirectory. None implies the current working directory"""
		global MPIMODE
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
			if (not MPIMODE) and not os.access("/tmp/eman2db-%s"%os.getenv("USER","anyone"),os.F_OK) : os.makedirs("/tmp/eman2db-%s"%os.getenv("USER","anyone"))
		else:
			if (not MPIMODE) and not os.access("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"),os.F_OK) : os.makedirs("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"))

		if MPIMODE:
			self.dbenv=None
		else :
			self.dbenv=db.DBEnv()
			self.dbenv.set_cachesize(0,cachesize,4)		# gbytes, bytes, ncache (splits into groups)
#			self.dbenv.set_cachesize(1,0,8)		# gbytes, bytes, ncache (splits into groups)
			self.dbenv.set_data_dir("%s/EMAN2DB"%self.path)
			self.dbenv.set_lk_detect(db.DB_LOCK_DEFAULT)	# internal deadlock detection

			if(sys.platform != 'win32'):
				self.dbenv.open("/tmp/eman2db-%s"%os.getenv("USER","anyone"),envopenflags)
			else:
				self.dbenv.open("/tmp/eman2db-%s"%os.getenv("USERNAME","anyone"),envopenflags)

		self.dicts={}
		#if self.__dbenv.DBfailchk(flags=0) :
			#self.LOG(1,"Database recovery required")
			#sys.exit(1)
		
		
	def close(self):
		"""close the environment associated with this object"""
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
		if name in self.dicts.keys():
			self.__dict__[name].close()
		try: os.unlink(self.path+"/EMAN2DB/"+name+".bdb")
		except: pass
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
	fixedkeys=frozenset(("nx","ny","nz","minimum","maximum","mean","sigma","square_sum","mean_nonzero","sigma_nonzero"))
	def __init__(self,name,file=None,dbenv=None,path=None,parent=None,ro=False):
		"""This is a persistent dictionary implemented as a BerkeleyDB Hash
		name is required, and will also be used as a filename if none is
		specified. Note that the database is not actually opened until it's used."""
		
		global dbopenflags
		DBDict.alldicts[self]=1		# we keep a running list of all trees so we can close everything properly
		self.name = name
		self.parent=parent
		self.dbenv=dbenv
		self.file=file
		self.rohint=ro
		self.lasttime=time.time()		# last time the database was accessed
		self.opencount=0				# number of times the database has needed reopening
		if path : self.path = path
		else : self.path=e2getcwd()
		self.txn=None	# current transaction used for all database operations
		self.bdb=None
		self.isro=False
		

	def __str__(self): return "<EMAN2db DBHash instance: %s>" % self.name

	def __del__(self):
		self.close()

	def open(self,ro=False):
		"""This actually opens the database (unless already open), if ro is set and the database is not already
		open read-write, it will be opened read-only"""

		self.lasttime=time.time()
		if self.bdb!=None :
			if ro==True or self.isro==False : 
#				print "already open",self.name
				return  		# return if the database is already open and in a compatible read-only mode
			self.close()	# we need to reopen read-write
		
		global MAXOPEN
		if DBDict.nopen>MAXOPEN : self.close_one() 

		self.bdb=db.DB(self.dbenv)		# we don't check MPIMODE here, since self.dbenv will already be None if its set
		if self.file==None : file=self.name+".bdb"
		else : file=self.file
#		print "open ",self.path+"/"+file,self.name,ro
		self.opencount+=1	# how many times we have had to reopen this database
		if ro : 
			try:
				self.bdb.open(self.path+"/"+file,self.name,db.DB_BTREE,db.DB_RDONLY|db.DB_THREAD)
				self.isro=True
			except: raise Exception,"No such database : %s"%self.path+"/"+file
		else : 
			try: 
				self.bdb.open(self.path+"/"+file,self.name,db.DB_BTREE,db.DB_CREATE|db.DB_THREAD)
			except: raise Exception,"Cannot create database : %s"%self.path+"/"+file
			self.isro=False
			
		DBDict.nopen+=1
#		print "%d open"%DBDict.nopen
#		print "opened ",self.name,ro
#		self.bdb.open(file,name,db.DB_HASH,dbopenflags)
#		print "Init ",name,file,path

	def close_one(self):
		"""Will select and close any excess open databases. Closure is based on the number of times it has been repoened and the
		time it was last used."""
		global MAXOPEN
		l=[(i.opencount,i.lasttime,i) for i in self.alldicts if i.bdb!=None]		# list of all open databases and usage,time info
		l.sort()

		if len(l)>MAXOPEN :
#			print "closing ",len(l)-MAXOPEN
			for i in range(len(l)-MAXOPEN): l[i][2].close()

	def close(self):
		if self.bdb == None: return
#		print "close x ",self.path+"/"+str(self.file),self.name,"XXX"
		self.bdb.close()
		self.bdb=None
		DBDict.nopen-=1
	
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
		self.open(self.rohint)	# make sure the database is open
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
		self.open()
		a=loads(self.bdb.get(dumps(n,-1),txn=self.txn))
		if isinstance(attr,dict) :
			a.update(attr)
		else: a[attr]=val
		self[n]=a

	def __len__(self):
		self.open(self.rohint)
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
				if n==9 : traceback.print_exc()
				print "********** Warning: problem writing ",key," to ",self.name,". Retrying (%d/10)"%n
				time.sleep(5)
				n+=1


	def __setitem__(self,key,val):
		self.open()
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
			self.put(fkey+dumps(key,-1),dumps(n,-1))		# a special key for the binary location
			
			# write the metadata
			try: del ad["data_path"]
			except:pass

			self.put(dumps(key,-1),dumps(ad,-1),txn=self.txn)

			if not self.has_key("maxrec") or key>self["maxrec"] : self["maxrec"]=key
			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) : self["maxrec"]=key
			
			# write the binary data
			val.write_data(pkey+fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"])
			
		else :
			self.put(dumps(key,-1),dumps(val,-1),txn=self.txn)
			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) : self["maxrec"]=key
				
	def __getitem__(self,key):
		self.open(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1),txn=self.txn))
		except: return None
		if isinstance(r,dict) and r.has_key("is_complex_x") :
			pkey="%s/%s_"%(self.path,self.name)
			fkey="%dx%dx%d"%(r["nx"],r["ny"],r["nz"])
#			print "r",fkey
			ret=EMData(r["nx"],r["ny"],r["nz"])
			if r.has_key("data_path"):
				p,l=r["data_path"].split("*")
				ret.read_data(p,int(l))
			else:
				n=loads(self.bdb.get(fkey+dumps(key,-1)))
				ret.read_data(pkey+fkey,n*4*r["nx"]*r["ny"]*r["nz"])
			k=set(r.keys())
			k-=DBDict.fixedkeys
			for i in k: ret.set_attr(i,r[i])
			return ret
		return r

	def __delitem__(self,key):
		self.open()
		self.bdb.delete(dumps(key,-1),txn=self.txn)

	def __contains__(self,key):
		self.open(self.rohint)
		return self.bdb.has_key(dumps(key,-1))

	def item_type(self,key):
		self.open(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1),txn=self.txn))
		except: return None
		if isinstance(r,dict) and r.has_key("is_complex_x") : return EMData
		return type(r)

	def keys(self):
		self.open(self.rohint)
		return [loads(x) for x in self.bdb.keys() if x[0]=='\x80']

	def values(self):
		self.open(self.rohint)
		return [self[k] for k in self.keys()]

	def items(self):
		self.open(self.rohint)
		return [(k,self[k]) for k in self.keys()]

	def has_key(self,key):
		self.open(self.rohint)
		return self.bdb.has_key(dumps(key,-1))

	def get_data_path(self,key):
		"""returns the path to the binary data as "path*location". Only valid for EMData objects."""
		self.open(self.rohint)
		try: r=loads(self.bdb.get(dumps(key,-1)))
		except: return None
		if isinstance(r,dict) and r.has_key("is_complex_x") :
			pkey="%s/%s_"%(self.path,self.name)
			fkey="%dx%dx%d"%(r["nx"],r["ny"],r["nz"])
			n=loads(self.bdb.get(fkey+dumps(key,-1)))
			if r.has_key("data_path"): return r["data_path"]
			else : return "%s*%d"%(pkey+fkey,n*4*r["nx"]*r["ny"]*r["nz"])
		return None

	def get(self,key,dfl=None,txn=None,target=None,nodata=0,region=None):
		"""Alternate method for retrieving records. Permits specification of an EMData 'target'
		object in which to place the read object"""
		self.open(self.rohint)
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
			ret = EMData()
			if target : ret = target
			ret.set_size(nx,ny,nz)
			
			# metadata
			k=set(r.keys())
			k-=DBDict.fixedkeys
			for i in k: 
				ret.set_attr(i,r[i])
			ret.set_attr("nx",nx)
			ret.set_attr("ny",ny)
			ret.set_attr("nz",nz)
				
			# binary data
			if not nodata:
				

				if region != None: ret.to_zero() # this has to occur in situations where the clip region goes outside the image
				if r.has_key("data_path"):
					p,l=r["data_path"].split("*")
					ret.read_data(p,int(l),region,rnx,rny,rnz)
				else:
					try: n=loads(self.bdb.get(fkey+dumps(key,-1)))
					except: raise KeyError,"Undefined data location key for : ",key
					ret.read_data(pkey+fkey,n*4*nx*ny*nz,region,rnx,rny,rnz)
			return ret
		return r
		
	def get_header(self,key,txn=None,target=None):
		"""Alternate method for retrieving metadata for EMData records."""
		self.open(self.rohint)
		try: return loads(self.bdb.get(dumps(key,-1),txn=txn))
		except: return None

	def set(self,key,val,txn=None):
		"Alternative to x[key]=val with transaction set"
		self.open()
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
				n=loads(self.bdb.get(fkey+dumps(key,-1),txn=txn))
			except:
				if not self.has_key(fkey) : self[fkey]=0
				else: self[fkey]+=1 
				n=self[fkey]
			self.put(fkey+dumps(key,-1),dumps(n,-1))		# a special key for the binary location
			
			# write the metadata
			try: del ad["data_path"]
			except:pass

			self.put(dumps(key,-1),dumps(ad,-1),txn=txn)

			if not self.has_key("maxrec") or key>self["maxrec"] : self["maxrec"]=key
			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) : self["maxrec"]=key
			
			# write the binary data
			val.write_data(pkey+fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"])
			
		else :
			self.put(dumps(key,-1),dumps(val,-1),txn=txn)
			if isinstance(key,int) and (not self.has_key("maxrec") or key>self["maxrec"]) : self["maxrec"]=key

	def set_header(self,key,val,txn=None):
		"Alternative to x[key]=val with transaction set"
		self.open()
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
		self.open()
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


