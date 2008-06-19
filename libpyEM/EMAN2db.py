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
try:
	from bsddb3 import db
except:
	from bsddb import db

envopenflags=db.DB_CREATE|db.DB_INIT_MPOOL|db.DB_INIT_LOCK|db.DB_INIT_LOG|db.DB_THREAD
dbopenflags=db.DB_CREATE
cachesize=100000000

class EMAN2DB:
	"""This class implements a local database of data and metadata for EMAN2"""
	
	def __init__(self):
		#if recover: xtraflags=db.DB_RECOVER
		if not os.access("EMAN2DB/dbcache",os.F_OK) : os.makedirs("EMAN2DB/dbcache")
		self.dbenv=db.DBEnv()
		self.dbenv.set_cachesize(0,cachesize,4)		# gbytes, bytes, ncache (splits into groups)
		self.dbenv.set_data_dir("EMAN2DB")
		self.dbenv.set_lk_detect(db.DB_LOCK_DEFAULT)	# internal deadlock detection
		#if self.__dbenv.DBfailchk(flags=0) :
			#self.LOG(1,"Database recovery required")
			#sys.exit(1)
			
		self.__dbenv.open("EMAN2DB",envopenflags|xtraflags)

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
"This module supports the concept of a local database for storing data and
metadata associated with a particular EMAN2 refinement. Data stored in this
database may be extracted into standard flat-files, but use of a database
with standard naming conventions, etc. helps provide the capability to log
the entire refinement process. 
"
