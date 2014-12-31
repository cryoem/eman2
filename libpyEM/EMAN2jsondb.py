#!/usr/bin/env python
#
# Author: Steven Ludtke, 05/01/2013 (sludtke@bcm.edu)
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
import json
import cPickle
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
import re

from libpyEMData2 import EMData
from libpyUtils2 import EMUtil


# If set, fairly verbose debugging information will be written to the console
# larger numbers will increase the amount of output
DBDEBUG=0

def js_one_key(url,key):
	"""Opens a JSON file and returns a single key before closing the file. Not really faster, but conserves memory by not leaving the file open"""

	return JSDict.one_key(url,key)

def js_open_dict(url):
	"""Opens a JSON file as a dict-like database object. The interface is almost identical to the BDB db_* functions.
If opened. Writes to JDB dictionaries may be somewhat inefficient due to the lack of a good model (as BDB has) for
multithreaded access. Default behavior is to write the entire dictionary to disk when any element is changed. File
locking is attempted to avoid conflicts, but may not work in all situations. read-only access is a meaningless concept
because file pointers are not held open beyond discrete transations. While it is possible to store images in JSON files
it is not recommended due to inefficiency, and making files which are difficult to read."""

	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"

	return JSDict.open_db(url)

def js_close_dict(url):
	"""This will free some resources associated with the database. Not associated with closing a file pointer at present."""

	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"

	ddb=JSDict.get_db(url)
	if ddb!=None : ddb.close()

	return

def js_remove_dict(url):
	"""closes and deletes a database using the same specification as db_open_dict"""

	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"

	js_close_dict(url)
	try : os.unlink(url)
	except OSError: pass

	return

def js_check_dict(url,readonly=True):
	"""Checks for the existence of the named JSON file and insures that it can be opened for reading [and writing].
It does not check the contents of the file, just for its exsistence and permissions."""

	if url==None : return False
	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"

	if readonly and os.access(url,os.R_OK) : return True
	if os.access(url,os.W_OK|os.R_OK) : return True

	return False

def js_list_dicts(url):
	"""Gives a list of readable json files at a given path."""

	try:
		ld=os.listdir(url)
	except: return []

	ld=[i for i in ld if i[-5:]==".json" and os.access("{}/{}".format(url,i),os.R_OK)]

	return ld

############
### JSON support for specific objects
############
from libpyEMData2 import *
from libpyAligner2 import *
from libpyTransform2 import *
import base64,zlib

def emdata_to_jsondict(obj):
	"""This is tacked on to EMData objects to give them non-pickle JSON support"""
	ret=obj.get_attr_dict()
	ret["__class__"]="EMData"
	ret["~bindata~"]=base64.encodestring(zlib.compress(obj.get_data_string(),1))		# we use ~ here as a delimiter because it's alphabetically after letters
	return ret

EMData.to_jsondict=emdata_to_jsondict		# we hack this into the EMData object

def emdata_from_jsondict(dct):
	"""This returns a new EMData object reconstituted from a JSON file"""
	fixedkeys=frozenset(("nx","ny","nz","minimum","maximum","mean","sigma","square_sum","mean_nonzero","sigma_nonzero","__class__","~bindata~"))
	ret=EMData(dct["nx"],dct["ny"],dct["nz"])
	ret.set_data_string(zlib.decompress(base64.decodestring(dct["~bindata~"])))
	for k in fixedkeys:
		try: del dct[k]
		except: pass

	for k in dct.keys():
		ret[k]=dct[k]

	ret.update()
	return ret

def eman2ctf_to_jsondict(obj):
	ret=obj.to_dict()
	ret["__class__"]="EMAN2Ctf"
	return ret

EMAN2Ctf.to_jsondict=eman2ctf_to_jsondict

def eman2ctf_from_jsondict(dct):
	ret=EMAN2Ctf()
	#ret.from_dict(dct)				# for some reason this crashes under certain situations, so we do it in python
	ret.defocus = dct["defocus"]
	ret.dfdiff = dct["dfdiff"]
	ret.dfang = dct["dfang"]
	ret.bfactor = dct["bfactor"]
	ret.ampcont = dct["ampcont"]
	ret.voltage = dct["voltage"]
	ret.cs = dct["cs"]
	ret.apix = dct["apix"]
	ret.dsbg = dct["dsbg"]
	ret.background = dct["background"]
	ret.snr = dct["snr"]
	return ret

def transform_to_jsondict(obj):
	ret={"__class__":"Transform"}
	ret["matrix"]=obj.get_matrix()
	return ret

Transform.to_jsondict=transform_to_jsondict

def transform_from_jsondict(dct):
	ret=Transform()
	ret.set_matrix(dct["matrix"])
	return ret

############
### File locking routines, hopefully platform independent
############

### First we try Linux/Mac
try:
	import fcntl		# Linux/Unix file locking

	def file_lock(fileobj, readonly=True):
		"""Unix file locking. Not truly enforced, but useful in thread synchronization. Multiple threads can have read-locks, but only one can have writelock"""

		if readonly : fcntl.lockf(fileobj.fileno(),fcntl.LOCK_SH)
		else : fcntl.lockf(fileobj.fileno(),fcntl.LOCK_EX)

	def file_unlock(fileobj):
		fcntl.lockf(fileobj.fileno(),fcntl.LOCK_UN)

### If that fails, we try windows
except ImportError:
	try:
		import msvcrt	# Windows file locking

		def file_lock(fileobj, readonly=True):
			"""Windows file locking (I hope)"""

			if readonly : l=msvcrt.LK_NBRLCK
			else : l=msvcrt.LK_NBLCK

			# We try for 30 seconds before giving up
			for i in xrange(30):
				try :
					msvcrt.locking(fileobj.fileno(), l, 1)
					break
				except:
					time.sleep(1)

			if i==30 :
				print "WARNING: Could not lock %s. Continuing without lock.  Please report this as a bug."%fileobj.name

		def file_unlock(fileobj) :
			try:
				msvcrt.locking(fileobj.fileno(), msvcrt.LK_UNLCK, 1)
			except:
				pass

### If that failed too, we don't lock files
	except:
		print "WARNING: Could not initialize either Linux/Mac or Windows file locking. Disabling locking (risky). Please report this as a bug !"

		def file_lock(fileobj, readonly=True):
			return

		def file_unlock(fileobj):
			return

#############
###  Task Management classes
#############

class JSTaskQueue:
	"""This class is used to manage active and completed tasks through an
	JSDict object. Tasks are each assigned a number unique in the local file.
	The active task 'max' keeps increasing. When a task is complete, it is shifted to the
	tasks_done list, and the 'max' value is increased if necessary. There is
	no guarantee in either list that all keys less than 'max' will exist."""

	lock=threading.Lock()
	caching=False		# flag to prevent multiple simultaneous caching

	def __init__(self,path=None):
		"""path should point to the directory where the disk-based task queue will reside without bdb:"""
		if path==None or len(path)==0 :
			path="tmp"

		if not os.path.isdir(path) : os.makedirs(path)
		self.path=path
		self.active=js_open_dict("%s/tasks_active.json"%path)		# active tasks keyed by id
		self.complete=file("%s/tasks_complete.txt"%path,"a")	# complete task log
		self.nametodid=js_open_dict("%s/tasks_name2did.json"%path)	# map local data filenames to did codes
		self.didtoname=js_open_dict("%s/tasks_did2name.json"%path)	# map data id to local filename
		self.precache=js_open_dict("%s/precache_files.json"%path)		# files to precache on clients, has one element "files" with a list of paths

		#if not self.active.has_key("max") : self.active["max"]=-1
		#if not self.active.has_key("min") : self.active["min"]=0
		#if not self.complete.has_key("max") : self.complete["max"]=0

	def to_jsondict(self):
		"""The path is really the only thing we need to store. Not sure that we ever really need to do this anyway..."""
		dct={"__class__":"JSTaskQueue"}
		dct["path"]=self.path
		return dct

	@classmethod
	def from_jsondict(cls,data):
		del data["__class__"]
		return cls(data["path"])

	def __len__(self) : return len(self.active)

	def get_task(self,clientid=0):
		"""This will return the next task waiting for execution"""
		JSTaskQueue.lock.acquire()
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
				JSTaskQueue.lock.release()
				return task

		JSTaskQueue.lock.release()
		return None

	def add_group(self):
		"""returns a new (unique) group id to be used for related tasks"""
		JSTaskQueue.lock.acquire()
		try: ret=self.active["grpctr"]+1
		except: ret=1
		self.active["grpctr"]=ret
		JSTaskQueue.lock.release()
		return ret

	def todid(self,name):
		"""Returns the did for a path, creating one if not already known"""
		fmt=e2filemodtime(name)
		try :
			did=self.nametodid[name]			# get the existing did from the cache (modtime,int)
			if fmt!=did[0]	:					# if the file has been changed, we need to assign a new did
				del self.didtoname[did]
				JSTaskQueue.lock.release()
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

		JSTaskQueue.lock.acquire()
		try: tid=self.active["taskctr"]+1
		except: tid=1
		self.active["taskctr"]=tid
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
		try: JSTaskQueue.lock.release()
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
			return 100		# if we don't find it, we assume it's done

		if task.starttime==None or task.starttime<1 : return -1
		if task.progtime==None : return 0
		return task.progtime[1]

	def task_done(self, tid):
		"""Mark a Task as complete, by removing it and logging it"""
		JSTaskQueue.lock.acquire()
		try:
			task=self.active[tid]
			if task==None: raise Exception
		except:
			print "*** Warning, task %d was already complete"%tid
			JSTaskQueue.lock.release()
			return

		# log the completed task
		self.complete.write("{tid}\t{cls}\t{runtime}\t{endtime}\t{starttime}\t{queuetime}\t{host}\n".format
			(tid=tid,cls=task.__class__.__name__,runtime=time.time()-task.starttime,endtime=time.time(),starttime=task.starttime,queuetime=task.queuetime,host=task.exechost))
		self.complete.flush()

		del self.active[tid]		# remove from active queue
		JSTaskQueue.lock.release()

	def task_rerun(self,taskid):
		"""If a task has been started somewhere, but execution fails due to a problem with the target node, call
		this and the task will be returned to the queue, unless it has failed MAXTASKFAIL times."""
		try:
			task=self.active[taskid]
			if task==None: raise Exception
		except:
			print "Fatal error: Could not find task {} to rerun.".format(taskid)

		if task.failcount==MAXTASKFAIL :
			self.task_aborted(taskid)
			return

		task.failcount+=1
		task.starttime=None
		task.progtime=None
		task.endtime=None
		task.clientid=None
		task.exechost=None

		self.active[taskid]=task
		return

	def task_aborted(self, taskid):
		"""Mark a Task as being aborted, by shifting a task to the tasks_complete queue"""
		JSTaskQueue.lock.acquire()
		try:
			task=self.active[taskid]
		except:
			JSTaskQueue.lock.release()
			return

		print "Error running task:\n{tid}\t{cls}\t{runtime}\t{endtime}\t{starttime}\t{queuetime}\t{host}".format(tid=tid,cls=task.__class__.__name__,runtime=time.time()-task.starttime,endtime=time.time(),starttime=task.starttime,queuetime=task.queuetime,host=task.exechost)
		#if self.active["min"]==taskid : self.active["min"]=min(self.active.keys())
		del self.active[taskid]
		JSTaskQueue.lock.release()

class JSTask:
	"""This class represents a task to be completed. Generally it will be subclassed. This is effectively
	an abstract superclass to define common member variables and methods. Note that the data dictionary,
	which contains a mix of actual data elements, and ["cache",filename,#|min,max|(list)] image references, will
	be transparently remapped on the client to similar specifiers which are valid locally. Over the network
	such data requests are remapped into data identifiers (did), then translated back into valid filenames
	in the remote cache.  When subclassing this class, avoid defining new member variables, as EMTask objects
	get transmitted over the network. Make use of command, data and options instead.

	If you subclass this class, make SURE you add an appropriate definition in the 'jsonclasses' dict below, and
	make sure the subclass is imported in EMAN2PAR.py"""
	def __init__(self,command=None,data=None,options=None,user=None):
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

	def to_jsondict(self):
		dct=dict(self.__dict__)		# copies the dict of the object
		dct["__class__"]=self.__class__.__name__
		#if self.__class__.__name__!="JSTask" :
			#print "WARNING : class <{}> must have to_jsondict and from_jsondict methods defined to function properly. See EMAN2jsondb.py.".format(self.__class__.__name__)
		return dct

	@classmethod
	def from_jsondict(cls,data):
		ret=cls()
		del data["__class__"]
		self.__dict__.update(data)
		return ret

	def execute(self): return


##########
### This object represents a single .json file in the filesystem as a persistent dictionary
### New JSDicts should be created by calling the static open_db method.
##########

# These two items are used to compact the file representation. Dict keys are separated by newlines,
# but lists are all on one line.
listrex=re.compile("\[[^\]\{]*\]")		# This regex will find all lists that don't contain dicts

def denl(s):
	"This will replace \n with nothing in a search match"
	return s.group(0).replace("\n","")

class JSDict:
	"""This class provides dict-like access to a JSON file on disk. It goes to some lengths to insure thread/process-safety, even if
performance must be sacrificed. The only case where it may not work is when a remote filesystem which doesn't obey file-locking is used.
Note that when opened, the entire JSON file is read into memory. For this reason (and others) it is not a good idea to store (many) images
in JSON files. JSDict objects are cached in RAM, and will not be removed from the cache unless the close() method is called, or too many
JSDicts are open at one time."""

	opendicts={}
	lock=threading.Lock()		# to make this section threadsafe

	@classmethod
	def open_db(cls,path=None):
		"""This should be used to create a JSDict instance. It caches already open dictionaries to avoid redundancy and conflicts."""

		cls.lock.acquire()

		if not isinstance(path,str) :
			cls.lock.release()
			raise Exception,"Must specify path to open JSONDB"
		if path[-5:]!=".json" :
			cls.lock.release()
			raise Exception,"JSON databases must have .json extension ('{}')".format(path)

		try: normpath=os.path.abspath(path)
		except:
			cls.lock.release()
			raise Exception,"Cannot find path for {}".format(path)

		if cls.opendicts.has_key(normpath) :
			cls.lock.release()
			return cls.opendicts[normpath]

		try : ret=JSDict(path)
		except:
			cls.lock.release()
			raise Exception,"Unable to open "+path

		cls.lock.release()

#		print "JSON: {} open {} kb".format(len(cls.opendicts),sum([i.filesize for i in cls.opendicts.values()])/1024)

		return ret

	@classmethod
	def get_db(cls,path=None):
		"""This will return an existing JSDict for 'path' if any, otherwise None"""

		cls.lock.acquire()

		if not isinstance(path,str) : raise Exception,"Must specify path to open JSONDB"
		if path[-5:]!=".json" :
			raise Exception,"JSON databases must have .json extension ('{}')".format(url)

		try: normpath=os.path.abspath(path)
		except: raise Exception,"Cannot find path for {}".format(path)

		ret=None
		if cls.opendicts.has_key(normpath) :
			ret=cls.opendicts[normpath]

		cls.lock.release()
		return ret

	@classmethod
	def one_key(cls,path,key):
		"""Opens a JSDict for path, reads a single key (if present) then closes the JSDict, to reduce memory issues
		when a large set of Dicts is iterated over in a large project"""

		try:
			db=JSDict.open_db(path)
			ret=db[key]
			db.close()
		except:
#			traceback.print_exc()
			return None

		return ret

	def __init__(self,path=None):
		"""This is a dict-like representation of a JSON file on disk. Warning, the entire file contents are parsed and held
in memory for efficient access. File change monitoring and file locking is used to insure self-consistency across processes.
Due to JSON module, there may be some data types which aren't permitted as values. While this module may be used like a traditional
dictionary for the most part, for efficiency, you may consider using the setval() and get() methods which permit deferring
synchronization with the disk file.

There is no name/path separation as existed with BDB objects. 'path' is a full path to the .json file. A normalized version
of the path is stored as self.normpath"""

		from EMAN2 import e2getcwd

		self.path=path
		try: self.normpath=os.path.abspath(path)
		except: self.normpath=path
		self.filesize=0					# stores the size of the text file on disk for approximate memory management

		self.data={}					# a cached copy of the actual data
		self.changes={}					# a set of changes to merge when next committing to disk
		self.delkeys=set()				# a set of keys to delete on next update
		self.lasttime=0					# last time the database was accessed

		self.sync()
		JSDict.opendicts[self.normpath]=self	# add ourselves to the cache

	def __str__(self): return "<JSDict instance: %s>" % self.path

	def __del__(self):
		if len(self.changes)>0 or len(self.delkeys): self.sync()

	def close(self):
		"""This will free effectively all of the memory associated with the object. It doesn't actually eliminate
		the object entirely since there may be multiple copies around. If the dictionary is accessed again, it will
		be automatically reopened."""
		if len(self.changes)>0 or len(self.delkeys): self.sync()
		self.lasttime=0
		self.data={}
#		del JSDict.opendicts[self.normpath]

	def sync(self):
		"""This is where all of the JSON file access occurs. This one routine handles both reading and writing, with file locking"""

		# We check for the _tmp file first
		try:
			mt2=os.stat(self.normpath[:-5]+"_tmp.json").st_mtime
		except:
			mt2=None

		# We check the modification time on the file to see if we need to read an update
		try:
			mt=os.stat(self.normpath).st_mtime
		except:		# file doesn't exist or we can't stat() it
			try:
				# if we find this, then we probably caught the files at just the wrong instant when another thread was doing an update
				if mt2!=None : time.sleep(0.5)
				mt=os.stat(self.normpath).st_mtime		# if it still doesn't exist, then the _tmp file may be an orphan
			except:
				if mt2!=None :
					# recover an orphaned js file (caused by interupt during write)
					os.rename(self.normpath[:-5]+"_tmp.json",self.normpath)
					mt=mt2
					mt2=None
				else :
					# if we get here there are only 2 possibilities, A) the file doesn't exist or B) we don't have read permission on the directory. In either case, this should be the right thing to do
					try : jfile=file(self.normpath,"w")
					except :
						try: os.makedirs(os.path.dirname(self.normpath))	# Can't open the file for writing, so we try to make sure the full path exists. If this fails, we let the actual Exception get raised
						except : pass
						try: jfile=file(self.normpath,"w")
						except: raise Exception,"Error: Unable to open {} for writing".format(self.normpath)
					file_lock(jfile,readonly=False)
					json.dump({},jfile)
					file_unlock(jfile)
					jfile=None
					mt=time.time()


		### Read entire dict from file
		# If we have unprocessed changes, or if the file has changed since last access
		if len(self.changes)>0 or mt>self.lasttime :
			jfile=file(self.normpath,"r")		# open the file
			file_lock(jfile,readonly=True)		# lock it for reading

			try:
				self.data=json.load(jfile,object_hook=json_to_obj)			# parse the whole JSON file, which should be a single dictionary
			except:
				jfile.seek(0)
				a=jfile.read()
				if len(a.strip())==0 : self.data={}		# json.load doesn't like completely empty files
				else :
					file_unlock(jfile)					# unlock the file
					print "Error in file: ",self.path
					traceback.print_exc()
					raise Exception,"Error reading JSON file : {}".format(self.path)
			self.filesize=jfile.tell()			# our location after reading the data from the file
			file_unlock(jfile)					# unlock the file
			jfile=None							# implicit close

		### Write entire dict to file
		# If we have unprocessed changes, we need to apply them and write back to disk
		if len(self.changes)>0:
			try:
				os.rename(self.normpath,self.normpath[:-5]+"_tmp.json")		# we back up the original file, just in case
			except:
				raise Exception,"WARNING: file '{}' cannot be created, conflict in writing JSON files. You may consider reporting this if you don't know why this happened.".format(self.normpath[:-3]+"_tmp.json")

			### We do the updates and prepare the string in-ram. If someone else tries a write while we're doing this, it should raise the above exception
			self.data.update(self.changes)		# update the internal copy of the data
			self.changes={}
			for k in self.delkeys:
				try: del self.data[k]
				except: pass
			self.delkeys=set()
			jss=json.dumps(self.data,indent=0,sort_keys=True,default=obj_to_json,encoding="ascii")			# write the whole dictionary back to disk
			jss=re.sub(listrex,denl,jss)

			### We do the actual write as a rapid sequence to avoid conflicts
			jfile=file(self.normpath,"w")
			file_lock(jfile,readonly=False)
			jfile.write(jss)
			file_unlock(jfile)
			jfile=None
			os.unlink(self.normpath[:-5]+"_tmp.json")

		self.lasttime=os.stat(self.normpath).st_mtime	# make sure we include our recent change, if made

	def __len__(self):
		"""Ignores any pending updates for speed"""
		return len(self.data)



	def __contains__(self,key):
		self.sync()
		if str(key) in self.data : return True
		return False

	def keys(self):
		self.sync()
		return self.data.keys()

	def values(self):
		self.sync()
		return self.data.values()

	def items(self):
		self.sync()
		return self.data.items()

	def has_key(self,key):
		self.sync()
		if str(key) in self.data : return True
		return False

	def update(self,newdict):
		"""Equivalent to dictionary update(). Performs JSON file update all at once, so substantially better
performance than many individual changes."""

		for k in newdict.keys(): self.setval(k,newdict[k],deferupdate=True)
		self.sync()

	def setdefault(self,key,dfl,noupdate=False):
		"""Returns the value for key if it exists. Otherwise, sets the value to be dfl and returns it"""

		key=str(key)

		if noupdate:
			if self.lasttime==0 : self.sync()		# if DB is closed, sync anyway
			if key in self.delkeys and key not in self.changes and key not in self.data :
				del self.delkeys[key]
				self.changes[key]=dfl
				return dfl
			if key in self.changes : return self.changes[key]
			return self.data[key]

		self.sync()
		if key in self.data : return self.data[key]
		self.changes[key]=dfl
		return dfl

	def getdefault(self,key,dfl,noupdate=False):
		"""Returns the value for key if it exists. Otherwise returns dfl"""

		key=str(key)

		if noupdate:
			if self.lasttime==0 : self.sync()		# if DB is closed, sync anyway
			if key in self.delkeys and key not in self.changes and key not in self.data : return dfl
			if key in self.changes : return self.changes[key]
			return self.data[key]

		self.sync()
		if key in self.data : return self.data[key]
		return dfl


	def get(self,key,noupdate=False):
		"""Alternate method for retrieving records. Can avoid expensive update call upon request."""

		key=str(key)

		if noupdate:
			if self.lasttime==0 : self.sync()		# if DB is closed, sync anyway
			if key in self.delkeys and key not in self.changes and key not in self.data : raise KeyError,key
			if key in self.changes : return self.changes[key]
			return self.data[key]

		self.sync()
		if key in self.data : return self.data[key]

		raise KeyError,key


	def setval(self,key,val,deferupdate=False):
		'''Alternative to assignment operator. Permits deferred writing to improve performance.'''
		# ok, decided to permit non-string keys to pass through and get converted to strings
#		if not isinstance(key,str) : raise Exception,"JSONDB keys must be strings"
		key=str(key)
		if key in self.delkeys : self.delkeys.remove(key)
		self.changes[key]=val
		if not deferupdate : self.sync()

	def delete(self,key,deferupdate=False):
		self.delkeys.add(key)
		if key in self.changes : del self.changes[key]
		if not deferupdate : self.sync()


JSDict.__setitem__=JSDict.setval
JSDict.__getitem__=JSDict.get
JSDict.__delitem__=JSDict.delete

### We must explicitly list any classes which are willing to be stored non-pickled in JSON
# more classes may get added to this dict by external modules when they import this one
jsonclasses = {
	"JSTask":JSTask.from_jsondict,
	"JSTaskQueue":JSTaskQueue.from_jsondict,
	"EMData":emdata_from_jsondict,
	"EMAN2Ctf":eman2ctf_from_jsondict,
	"Transform":transform_from_jsondict
}

def json_to_obj(jsdata):
	"""converts a javascript object representation back to the original python object"""

	if jsdata.has_key("__pickle__") :
		try: return cPickle.loads(str(jsdata["__pickle__"]))
		except: return str(jsdata["__pickle__"])				# This shouldn't happen. Means a module hasn't been loaded. This is an emergency stopgap to avoid crashing
	elif jsdata.has_key("__class__") : return jsonclasses[jsdata["__class__"]](jsdata)
	else: return jsdata

def obj_to_json(obj):
	"""converts a python object to a supportable json type"""
	try:
		return obj.to_jsondict()
	except:
		return {"__pickle__":cPickle.dumps(obj,0)}

__doc__ = \
"""This module provides a dict-like wrapper for JSON files on disk, with full support for file locking and other
internal consistency measures to improve thread/process safety. In some cases performance is sacrificed to insure
consistency. Performance will be substantially worse than BDB, particularly with large databases, but the JSON files
are human-readable, and it will not suffer from the issues with caching and database corruption which frustrated so
many users with BDB."""




