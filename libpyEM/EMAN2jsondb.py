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

from libpyEMData2 import EMData
from libpyUtils2 import EMUtil


# If set, fairly verbose debugging information will be written to the console
# larger numbers will increase the amount of output
DBDEBUG=0

def json_open_dict(url):
	"""Opens a JSON file as a dict-like database object. The interface is almost identical to the BDB db_* functions.
If opened. Writes to JDB dictionaries may be somewhat inefficient due to the lack of a good model (as BDB has) for
multithreaded access. Default behavior is to write the entire dictionary to disk when any element is changed. File
locking is attempted to avoid conflicts, but may not work in all situations. read-only access is a meaningless concept
because file pointers are not held open beyond discrete transations. While it is possible to store images in JSON files
it is not recommended due to inefficiency, and making files which are difficult to read."""

	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"
	
	return JSONDict.open_db(url)

def json_close_dict(url):
	"""This will free some resources associated with the database. Not associated with closing a file pointer at present."""
	
	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"
	
	ddb=JSONDict.get_db(url)
	if ddb!=None : ddb.close()
	
	return

def json_remove_dict(url):
	"""closes and deletes a database using the same specification as db_open_dict"""

	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"
	
	json_close_dict(url)
	os.unlink(url)
	
	return

def json_check_dict(url,readonly=True):
	"""Checks for the existence of the named JSON file and insures that it can be opened for reading [and writing].
It does not check the contents of the file, just for its exsistence and permissions."""

	if url[-5:]!=".json" :
		raise Exception,"JSON databases must have .json extension"
	
	if readonly and os.access(url,os.R_OK) : return True
	if os.access(url,os.W_OK|os.R_OK) : return True
	
	return False
	
def json_list_dicts(url):
	"""Gives a list of readable json files at a given path."""
	
	try:
		ld=os.listdir(url)
	except: return []
	
	ld=[i for i in ld if i[-5:]==".json" and os.access("{}.{}".format(url,i),os.R_OK)]
	
	return ld

### File locking routines, hopefully platform independent
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

### FIXME - this class is not yet implemented. Still has the old BDB implementation
class JSONTaskQueue:
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



##########
### This object represents a single .json file in the filesystem as a persistent dictionary
### New JSONDicts should be created by calling the static open_db method.
##########	
	
class JSONDict:
	"""This class provides dict-like access to a JSON file on disk. It goes to some lengths to insure thread/process-safety, even if
performance must be sacrificed. The only case where it may not work is when a remote filesystem which doesn't obey file-locking is used.
Note that when opened, the entire JSON file is read into memory. For this reason (and others) it is not a good idea to store images
in JSON files."""
	
	opendicts={}
	lock=threading.Lock()		# to make this section threadsafe
	
	@classmethod
	def open_db(cls,path=None):
		"""This should be used to create a JSONDict instance. It caches already open dictionaries to avoid redundancy and conflicts."""
		
		cls.lock.acquire()
		
		if not isinstance(path,str) : raise Exception,"Must specify path to open JSONDB"
		if path[-5:]!=".json" :
			raise Exception,"JSON databases must have .json extension ('{}')".format(path)
		
		try: normpath=os.path.abspath(path)
		except: raise Exception,"Cannot find path for {}".format(path)
		
		if cls.opendicts.has_key(normpath) : 
			cls.lock.release()
			return cls.opendicts[normpath]
		
		ret=JSONDict(path)
		cls.lock.release()
		return ret

	@classmethod
	def get_db(cls,path=None):
		"""This will return an existing JSONDict for 'path' if any, otherwise None"""
		
		cls.lock.acquire()
		
		if not isinstance(path,str) : raise Exception,"Must specify path to open JSONDB"
		if url[-5:]!=".json" :
			raise Exception,"JSON databases must have .json extension ('{}')".format(url)
		
		try: normpath=os.path.abspath(path)
		except: raise Exception,"Cannot find path for {}".format(path)
		
		ret=None
		if cls.opendicts.has_key(normpath) : 
			ret=cls.opendicts[normpath]
		
		cls.lock.release()
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
		
		self.data={}					# a cached copy of the actual data
		self.changes={}					# a set of changes to merge when next committing to disk
		self.delkeys=set()				# a set of keys to delete on next update
		self.lasttime=0					# last time the database was accessed
		
		self.sync()

	def __str__(self): return "<JSONDict instance: %s>" % self.path

	def __del__(self):
		if len(self.changes)>0 : self.sync()
				

	def close(self):
		if len(self.changes)>0 : self.sync()

	def sync(self):
		"""This is where all of the JSON file access occurs. This one routine handles both reading and writing, with file locking"""
	
		# We check the modification time on the file to see if we need to read an update
		try:
			mt=os.stat(self.normpath).st_mtime
		except:
			try:
				mt=os.stat(self.normpath[:-5]+"_tmp.json").st_mtime		# if we find this, then we probably caught the files at just the wrong instant when another thread was doing an update
				time.sleep(0.2)
				mt=os.stat(self.normpath).st_mtime		# if it still doesn't exist, then the _tmp file may be an orphan we should ignore
			except:
				# if we get here there are only 2 possibilities, A) the file doesn't exist or B) we don't have read permission on the directory. In either case, this should be the right thing to do
				jfile=file(self.normpath,"w")
				file_lock(jfile,readonly=False)
				json.dump({},jfile)
				file_unlock(jfile)
				jfile=None
				mt=time.time()
		
		# If we have unprocessed changes, or if the file has changed since last access, we reread from disk
		if len(self.changes)>0 or mt>self.lasttime :
			jfile=file(self.normpath,"r")		# open the file
			file_lock(jfile,readonly=True)		# lock it for reading
			try:
				self.data=json.load(jfile)			# parse the whole JSON file, which should be a single dictionary
			except:
				jfile.seek(0)
				a=jfile.read()
				if len(a.strip())==0 : self.data={}		# json.load doesn't like completely empty files
				else :
					file_unlock(jfile)					# unlock the file
					raise Exception,"Error reading JSON file : {}".format(self.path)
			file_unlock(jfile)					# unlock the file
			jfile=None							# implicit close
			
		# If we have unprocessed changes, we need to apply them and write back to disk
		if len(self.changes)>0:
			try:
				os.rename(self.normpath,self.normpath[:-5]+"_tmp.json")		# we back up the original file, just in case
			except:
				raise Exception,"WARNING: file '{}' cannot be created, conflict in writing JSON files. You may consider reporting this if you don't know what happened.".format(self.normpath[:-5]+"_tmp.json")
			jfile=file(self.normpath,"w")
			file_lock(jfile,readonly=False)
			self.data.update(self.changes)		# update the internal copy of the data
			self.changes={}
			for k in self.delkeys: del self.data[k]
			self.delkeys=set()
			json.dump(self.data,jfile,indent=1,sort_keys=True)			# write the whole dictionary back to disk
			file_unlock(jfile)
			jfile=None
			os.unlink(self.normpath[:-5]+"_tmp.json")
			
		self.lasttime=os.stat(self.normpath).st_mtime	# make sure we include our recent change, if made

	def __len__(self):
		"""Ignores any pending updates for speed"""
		return len(self.data)



	def __contains__(self,key):
		self.sync()
		if key in self.data : return True
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
		if key in self.data : return True
		return False

	def update(self,newdict):
		"""Equivalent to dictionary update(). Performs JSON file update all at once, so substantially better
performance than many individual changes."""

		for k in newdict.keys(): self.setval(k,newdict[k],deferupdate=True)
		self.sync()

	def get(self,key,noupdate=False):
		"""Alternate method for retrieving records. Can avoid expensive update call upon request."""
		if noupdate :
			if key in self.delkeys : raise KeyError
			try: return repairval(self.changes[key])
			except: return repairval(self.data[key])
			
		self.sync()
		return repairval(self.data[key])
	
	def setval(self,key,val,deferupdate=False):
		'''Alternative to assignment operator. Permits deferred writing to improve performance.'''
		if not isinstance(key,str) : raise Exception,"JSONDB keys must be strings"
		if key in self.delkeys : self.delkeys.remove(key)
		if isinstance(val,str) :
			if val[0]=="$" : val="$"+value
		elif not isinstance(val,int) and not isinstance(val,float) : 
			try : tmp=json.dumps(val)
			except: val="$ "+cPickle.dumps(val,protocol=0)
		self.changes[key]=val
		if not deferupdate : self.sync()

	def delete(self,key,deferupdate=False):
		self.delkeys.add(key)
		if key in self.changes : del self.changes[key]
		if not deferupdate : self.sync()
		

JSONDict.__setitem__=JSONDict.setval
JSONDict.__getitem__=JSONDict.get
JSONDict.__delitem__=JSONDict.delete

def repairval(v):
	"""Objects with no JSON representation are stored as pickled strings or other specialty types, with a $ prefix
as an identifier"""
	if not isinstance(v,str) and not isinstance(v,unicode): return v
	if v[:2]=="$$" : return v[1:]
	if v[:2]=="$ " : return cPickle.loads(str(v[2:]))
	return v

__doc__ = \
"""This module provides a dict-like wrapper for JSON files on disk, with full support for file locking and other
internal consistency measures to improve thread/process safety. In some cases performance is sacrificed to insure
consistency. Performance will be substantially worse than BDB, particularly with large databases, but the JSON files
are human-readable, and it will not suffer from the issues with caching and database corruption which frustrated so
many users with BDB."""




