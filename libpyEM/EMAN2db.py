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
try:
	from bsddb3 import db
except:
	from bsddb import db
	

envopenflags=db.DB_CREATE|db.DB_INIT_MPOOL|db.DB_INIT_LOCK|db.DB_INIT_LOG|db.DB_THREAD
#dbopenflags=db.DB_CREATE
cachesize=10000000

class EMAN2DB:
	"""This class implements a local database of data and metadata for EMAN2"""
	
	def __init__(self,path=None):
		"""path points to the directory containin the EMAN2DB subdirectory. None implies the current working directory"""
		#if recover: xtraflags=db.DB_RECOVER
		if not path : path=os.getcwd()
		self.path=path
		
		if not os.access("./EMAN2DB/cache",os.F_OK) : os.makedirs("./EMAN2DB/cache")
		self.dbenv=db.DBEnv()
		self.dbenv.set_cachesize(0,cachesize,4)		# gbytes, bytes, ncache (splits into groups)
#		self.dbenv.set_cachesize(1,0,8)		# gbytes, bytes, ncache (splits into groups)
		self.dbenv.set_data_dir("%s/EMAN2DB"%self.path)
		self.dbenv.set_lk_detect(db.DB_LOCK_DEFAULT)	# internal deadlock detection
		self.dicts=[]
		#if self.__dbenv.DBfailchk(flags=0) :
			#self.LOG(1,"Database recovery required")
			#sys.exit(1)
			
		self.dbenv.open("%s/EMAN2DB/cache"%self.path,envopenflags)

	def __del__(self):
		self.dbenv=None
		for i in self.dicts : self.close_dict(i)

	def open_dict(self,name):
		self.__dict__[name]=DBDict(name,dbenv=self.dbenv,path=self.path+"/EMAN2DB")
		self.dicts.append(name)
	
	def close_dict(self,name):
		self.__dict__[name].close()
		del(self.__dict__[name])
		self.dicts.remove(name)

class DBDict:
	"""This class uses BerkeleyDB to create an object much like a persistent Python Dictionary,
	keys and data may be arbitrary pickleable types, however, additional functionality is provided
	for EMData objects. Specifically, if integer keys are used, set_attr and get_attr may be used
	to efficiently get and set attributes for images with reduced i/o requirements (in certain cases)."""
	
	allhashes=weakref.WeakKeyDictionary()
	fixedkeys=frozenset(("nx","ny","nz","minimum","maximum","mean","sigma","square_sum","mean_nonzero","sigma_nonzero"))
	def __init__(self,name,file=None,dbenv=None,path=None):
		"""This is a persistent dictionary implemented as a BerkeleyDB Hash
		name is required, and will also be used as a filename if none is
		specified. """
		
		global dbopenflags
		DBDict.allhashes[self]=1		# we keep a running list of all trees so we can close everything properly
		self.name = name
		if path : self.path = path
		else : self.path=os.getcwd()
		self.txn=None	# current transaction used for all database operations
		self.bdb=db.DB(dbenv)
		if file==None : file=name+".bdb"
		try: self.bdb.open(file,name,db.DB_BTREE,0)
		except: 
			self.bdb.open(file,name,db.DB_BTREE,db.DB_CREATE)
			self.bdb.put(dumps("_freebinary",-1),dumps([],-1),txn=self.txn)
			print "new db created"
#		self.bdb.open(file,name,db.DB_HASH,dbopenflags)

	def __str__(self): return "<EMAN2db DBHash instance: %s>" % self.name

	def __del__(self):
		self.close()

	def close(self):
		if self.bdb == None: return
		self.bdb.close()
		self.bdb=None
	
	def sync(self):
		if self.bdb : self.bdb.sync()
		
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
		return loads(self.bdb.get(dumps(n,-1),txn=self.txn))[attr]
		
	def set_attr(self,n,attr,val):
		a=loads(self.bdb.get(dumps(n,-1),txn=self.txn))
		a[attr]=val
		self[n]=a

	def __len__(self):
		return self.bdb.stat(db.DB_FAST_STAT)["nkeys"]
#		return len(self.bdb)

	def __setitem__(self,key,val):
		if (val==None) :
			try: self.__delitem__(key)
			except: return
		elif isinstance(val,EMData) : 
			# decide where to put the binary data
			ad=val.get_attr_dict()
			fkey="%s/%s%dx%dx%d"%(self.path,self.name,ad["nx"],ad["ny"],ad["nz"])
#			print "w",fkey
			try :
				o=loads(self.bdb.get(dumps(key,-1),txn=self.txn))
				n=o["binary"]
			except:
				if not self.has_key(fkey) : self[fkey]=0
				else: self[fkey]+=1 
				n=self[fkey]
			ad["binary"]=n
			
			# write the metadata
			self.bdb.put(dumps(key,-1),dumps(ad,-1),txn=self.txn)
			if not self.has_key("maxrec") or key>self["maxrec"] : self["maxrec"]=key
			
			# write the binary data
			val.write_data(fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"])
			
		else :
			self.bdb.put(dumps(key,-1),dumps(val,-1),txn=self.txn)
				
	def __getitem__(self,key):
		try: r=loads(self.bdb.get(dumps(key,-1),txn=self.txn))
		except: return None
		if isinstance(r,dict) and r.has_key("binary") :
			fkey="%s/%s%dx%dx%d"%(self.path,self.name,r["nx"],r["ny"],r["nz"])
#			print "r",fkey
			ret=EMData(r["nx"],r["ny"],r["nz"])
			ret.read_data(fkey,r["binary"]*4*r["nx"]*r["ny"]*r["nz"])
			k=set(r.keys())
			k-=DBDict.fixedkeys
			for i in k: ret.set_attr(i,r[i])
			return ret
		return r

	def __delitem__(self,key):
		self.bdb.delete(dumps(key,-1),txn=self.txn)

	def __contains__(self,key):
		return self.bdb.has_key(dumps(key,-1),txn=self.txn)

	def keys(self):
		return map(lambda x:loads(x),self.bdb.keys())

	def values(self):
		return map(lambda x:loads(decompress(x)),self.bdb.values())

	def items(self):
		return map(lambda x:(loads(x[0]),loads(decompress(x[1]))),self.bdb.items())

	def has_key(self,key):
		return self.bdb.has_key(dumps(key,-1))

	def get(self,key,txn=None):
		try: r=loads(self.bdb.get(dumps(key,-1),txn=txn))
		except: return None
		if isinstance(r,dict) and r.has_key("binary") :
			fkey="%s/%s%dx%dx%d"%(self.path,self.name,r["nx"],r["ny"],r["nz"])
#			print "r",fkey
			ret=EMData(r["nx"],r["ny"],r["nz"])
			ret.read_data(fkey,r["binary"]*4*r["nx"]*r["ny"]*r["nz"])
			k=set(r.keys())
			k-=DBDict.fixedkeys
			for i in k: ret.set_attr(i,r[i])
			return ret
		return r

	def set(self,key,val,txn=None):
		"Alternative to x[key]=val with transaction set"
		if (val==None) :
			try: self.__delitem__(key)
			except: return
		elif isinstance(val,EMData) : 
			# decide where to put the binary data
			ad=val.get_attr_dict()
			fkey="%s/%s%dx%dx%d"%(self.path,self.name,ad["nx"],ad["ny"],ad["nz"])
#			print "w",fkey
			try :
				o=loads(self.bdb.get(dumps(key,-1),txn=txn))
				n=o["binary"]
			except:
				if not self.has_key(fkey) : self[fkey]=0
				else: self[fkey]+=1 
				n=self[fkey]
			ad["binary"]=n
			
			# write the metadata
			self.bdb.put(dumps(key,-1),dumps(ad,-1),txn=txn)
			if not self.has_key("maxrec") or key>self["maxrec"] : self["maxrec"]=key
			
			# write the binary data
			val.write_data(fkey,n*4*ad["nx"]*ad["ny"]*ad["nz"])
			
		else :
			self.bdb.put(dumps(key,-1),dumps(val,-1),txn=txn)

	def update(self,dict):
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
