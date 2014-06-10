#!/usr/bin/env python

#
# Author: Steven Ludtke, 03/06/2009 (sludtke@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# e2parallel.py Steven Ludtke
# This program implements, via various options, the parallelism system for EMAN2

from EMAN2db import EMTaskQueue, EMTask
from EMAN2 import *
from EMAN2PAR import *
from math import *
import time
import os
import sys
import socket
import traceback
import subprocess

debug=False
logid=None

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	commandlist=("dcserver","dcclient","realdcclient","dckill","dckillclients","servmon","rerunall","killall","precache","localclient","mpiclient")
	usage = """prog [options] <command> ...
	
This program implements much of EMAN2's coarse-grained parallelism mechanism. There are several flavors available via
different options in this program. Note that the simple threaded parallelism used on a single computer with multiple
CPUs does not require use of this program. See the Wiki for more information.  This tool primarily serves the needs
of the distributed computing (dc) parallelism mechanism.

<command> is one of: dcserver, dcclient, dckillclients, dckill, rerunall, precache, killall, servmon

run e2parallel.py servmon to run a GUI server monitor. This MUST run on the same machine in the same directory as the server.

client-server DC system:
run e2parallel.py dcserver on the machine containing your data and project files
run e2parallel.py dcclient on as many other machines as possible, pointing at the server machine

"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--server",type=str,help="Specifies host of the server to connect to",default="localhost")
	parser.add_argument("--port",type=int,help="Specifies server port, default is automatic assignment",default=-1)
	parser.add_argument("--clientpath",type=str,help="Scratch directory for DC clients. Default is current directory.",default=None)
	parser.add_argument("--clientid",type=int,help="Internal use only",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--taskin", type=str,help="Internal use only. Used when executing local threaded tasks.")
	parser.add_argument("--taskout", type=str,help="Internal use only. Used when executing local threaded tasks.")
	parser.add_argument("--scratchdir", type=str,help="Internal use only. Used by the MPI client",default="/tmp")
	parser.add_argument("--cache", action="store_true",help="Internal use only. Used by the MPI client",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
#	parser.add_argument("--cpus",type=int,help="Number of CPUs/Cores for the clients to use on the local machine")
#	parser.add_argument("--idleonly",action="store_true",help="Will only use CPUs on the local machine when idle",default=False)
	
	(options, args) = parser.parse_args()
	# Initialize CUDA iof needed
	initializeCUDAdevice()
	
	if len(args)<1 or args[0] not in commandlist: parser.error("command required: "+str(commandlist))

	if args[0]=="dcserver" :
		print "Sorry : e2parallel DC mode is not supported in EMAN2.1. We may revisit this if users complain, though. Please email sludtke@bcm.edu."
		sys.exit(1)
#		rundcserver(options.port,options.verbose-1)
		
	elif args[0]=="dcclient" :
		print "Sorry : e2parallel DC mode is not supported in EMAN2.1. We may revisit this if users complain, though. Please email sludtke@bcm.edu."
		sys.exit(1)
		if options.clientpath!=None : 
			try: os.makedirs(options.clientpath)
			except: pass
			os.chdir(options.clientpath)
		rundcclients(options.server,options.port,options.verbose-1)
		
	elif args[0]=="realdcclient" :
		rundcclient(options.server,options.port,options.verbose-1,options.clientid)
		
	elif args[0]=="dckill" :
		killdcserver(options.server,options.port,options.verbose-1)
	
	elif args[0]=="dckillclients" :
		killdcclients(options.server,options.port,options.verbose-1)

	elif args[0]=="precache" :
		precache(args[1:])

	elif args[0]=="localclient" :
		runlocaltask(options.taskin,options.taskout)

	elif args[0]=="mpiclient" :
		runinmpi(options.scratchdir,options.verbose)
		
	elif args[0]=="rerunall":
		rerunall()

	elif args[0]=="killall":
		killall()

	elif args[0]=="servmon" :
		runservmon()
	
	sys.exit(0)
	
def progcb(val):
	return True

def runinmpi(scratchdir,verbose):
	"""This function can only be used when called by mpirun from a valid MPI environment"""
#	from mpi import *
#	mpi_barrier(MPI_COMM_WORLD)
	client=EMMpiClient(scratchdir)
	client.test(verbose)		# don't skip this. It's necessary to identify node names
	client.run(verbose)

	sys.exit(0)

def runlocaltask(taskin,taskout):
	"""Exectues a task on the local machine. Reads the pickled task from 'taskfile'. Returns results to taskout. """
	from cPickle import load,dump

#	print "Executing %s (%s)"%(taskin,taskout)
	from e2classaverage import ClassAvTask
	from e2simmx import EMSimTaskDC
	from e2project3d import EMProject3DTaskDC
#	from e2tomoaverage import EMTomoAlignTaskDC
	
	task=load(file(taskin,"rb"))
	
	try: dump(task.execute(progcb),file(taskout,"wb"),-1)
	except:
		traceback.print_exc()
		sys.exit(1)		# Error !
#	print "Done %s (%s)"%(taskin,taskout)
	
	sys.exit(0)

def rundcserver(port,verbose):
	"""Launches a DCServer. If port is <1 or None, will autodetermine. Does not return."""
	import EMAN2db
	# The following was causing issues with the multithreaded parallelism server. Seems like we need to insure the server and the customer
	# are running on the same physical computer !!!
#	EMAN2db.BDB_CACHE_DISABLE=1	# this diables caching on the server so the customer knows it can freely write to local database files
	server=runEMDCServer(port,verbose)			# never returns

def killdcclients(server,port,verbose):
	import EMAN2db
	server=runEMDCServer(port,verbose,True)			# never returns

def rundcclients(host,port,verbose):
	clientid=random.randint(1,2000000000)
	while 1:
		rc=subprocess.call(["e2parallel.py","realdcclient","--server="+str(host),"--port="+str(port),"--verbose="+str(verbose),"--clientid="+str(clientid)])
		if rc : 
			if rc==1 : print "Client exited at server request"
			else : print "Client exited with status code %s"%str(rc)
			break

def rundcclient(host,port,verbose,clientid):
	"""Starts a DC client running, runs forever"""
#	while (1):
	client=EMDCTaskClient(host,port,verbose,myid=clientid)
	r=client.run(onejob=True)
	sys.exit(r)
#	print "New client (%d alloced)"%EMData.totalalloc

def precache(files):
	"""Adds a list of filenames to the precaching queue. Precaching will occur before jobs are started."""
	q=EMTaskQueue()
	print len(files)," files queued for precaching" 
	q.precache["files"]=files
	

def rerunall():
	"""Requeues all active (incomplete) tasks"""
	q=EMTaskQueue()
	e=q.active.keys()
	e=[i for i in e if isinstance(i,int) and q.active[i]!=None]
	
	for i in e: q.task_rerun(i)
	
	print "Requeued %d tasks"%len(e)

def killall():
	"""Requeues all active (incomplete) tasks"""
	q=EMTaskQueue()
	e=q.active.keys()
	e=[i for i in e if isinstance(i,int) and q.active[i]!=None]
	
	for i in e: q.task_aborted(i)

	print "Killed %d tasks"%len(e)

def killdcserver(server,port,verbose):
	EMDCsendonecom(server,port,"QUIT")


# We import Qt even if we don't need it
try:
	from PyQt4 import QtCore, QtGui
	from PyQt4.QtCore import Qt
except:
	class dummy:
		"A dummy class for use when Qt not installed"
		def __init__(self,quiet=False):
			if not quiet :
				print "ERROR: Qt4 could not be imported, check your PyQt installation"
				import traceback
				traceback.print_exc()

	QtGui=dummy(True)
	QtGui.QWidget=dummy
	QtGui.QMainWindow=dummy
	QtCore=dummy(True)
	QtCore.QAbstractTableModel=dummy
	
def runservmon():
	import EMAN2db
	# we changed the meaning of the variable to disable writing to cache altogether, pap 9-01
	EMAN2db.BDB_CACHE_DISABLE=1

	queue=EMAN2db.EMTaskQueue(".",ro=True)

#	activedata=TaskData(queue.active)
#	completedata=TaskData(queue.complete)

	app = QtGui.QApplication([])
	window = GUIservmon()
	
#	ui.tableView.setModel(data)

	window.show()
	window.set_data(queue)
	
	app.exec_()
	
class GUIservmon(QtGui.QMainWindow):
	"""A DC server monitor GUI"""
	def __init__(self):
		QtGui.QWidget.__init__(self,None)

		self.cw=QtGui.QWidget()
		self.setCentralWidget(self.cw)
		self.vbl=QtGui.QVBoxLayout(self.cw)
		
		self.tabs = QtGui.QTabWidget()
		self.vbl.addWidget(self.tabs)
		self.tabs.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
	
		self.activetab=QtGui.QWidget()
		self.vblat=QtGui.QVBoxLayout(self.activetab)
		self.actview=QtGui.QTableView()
		self.vblat.addWidget(self.actview)
		self.tabs.addTab(self.activetab,"Active")
		
		self.donetab=QtGui.QWidget()
		self.vbldt=QtGui.QVBoxLayout(self.donetab)
		self.doneview=QtGui.QTableView()
		self.vbldt.addWidget(self.doneview)
		self.tabs.addTab(self.donetab,"Complete")
		
		self.clienttab=QtGui.QWidget()
		self.vblct=QtGui.QVBoxLayout(self.clienttab)
		self.clientview=QtGui.QTableView()
		self.vblct.addWidget(self.clientview)
		self.tabs.addTab(self.clienttab,"Clients")
		
		self.startTimer(10000)
		
	def timerEvent(self,event):
		if self.tabs.currentIndex()==0 :
			self.actmodel.load()
			self.actview.setModel(None) 				# suboptimal, but a hack for now
			self.actview.setModel(self.actmodel)
			self.actview.resizeColumnsToContents()
		elif self.tabs.currentIndex()==1 :
			self.donemodel.load()
			self.doneview.setModel(None)
			self.doneview.setModel(self.donemodel)
			self.doneview.resizeColumnsToContents()
		
	def set_data(self,queue):
		"""This takes an EMTaskQueue object and displays it"""
		
		self.actmodel=TaskData(queue.active)
		self.actview.setModel(self.actmodel)
		self.donemodel=TaskData(queue.complete)
		self.doneview.setModel(self.donemodel)
		
		self.actmodel.load()
		self.donemodel.load()
		
#		self.vbl.addWidget(self.tabs)

class TaskData(QtCore.QAbstractTableModel):
	def __init__(self,target):
		QtCore.QAbstractTableModel.__init__(self)
		self.target=target
		self.nrows=0
		self.rows=[]

	def load(self):
		"""Updates the cached display from source"""
		keys=self.target.keys()
		keys.sort()
		keys=keys[:-2]
		self.nrows=len(keys)
		self.rows=[]
		for r in range(self.nrows):
			task=self.target[keys[r]]
			self.rows.append([self.col(task,i) for i in range(8)])
		
	def col(self,task,n):
		"""gets a single table entry"""
		if not isinstance(task,EMTask) : 
			print loc.row(),keys[loc.row()]
			print self.target[keys[loc.row()]]
			return QtCore.QVariant("???")
			
		if n==0 : ret=task.taskid
		elif n==1: 
			try: 
				if task.progtime==None or task.progtime[1]==-1 : ret = "-"
				elif task.progtime[1]==0 : ret= "#"
				elif task.progtime[1]<100 : ret= "#"*(1+task.progtime[1]/5)
				else : ret = "DONE"
				if task.progtime[0]-time.time()>300 : ret+=" ?"
			except: ret="?"
		elif n==2 : ret=task.command
		elif n==3 : ret=local_datetime(task.queuetime)
		elif n==4 : ret=local_datetime(task.starttime)
		elif n==5 : 
			try: ret=difftime(task.endtime-task.starttime)
			except: ret = "incomplete"
		elif n==6 : 
			ret=task.exechost
			if ret==None : ret=task.clientid
		elif n==7 :
			try:
				if   task.command=="e2classaverage.py" : ret= "Class %d"%task.class_idx
				elif task.command=="e2simmx.py" : ret= "Range: %d - %d : %d - %d"%(task.data["references"][2],task.data["references"][3],task.data["particles"][2],task.data["particles"][3])
				elif task.command=="e2project3d.py" : ret="Proj: %d - %d"%(min(task.data["indices"]),max(task.data["indices"]))
				else : ret = task.command
			except : ret = str(task)
			
		return QtCore.QVariant(str(ret))

	def data(self,loc,role):
		if not loc.isValid() or role != QtCore.Qt.DisplayRole : return QtCore.QVariant()
		try : return self.rows[loc.row()][loc.column()]
		except : return QtCore.QVariant("---")
		
	def rowCount(self,parent):
		if parent.isValid() : return 0
		return self.nrows

	def columnCount(self,parent):
		if parent.isValid() : return 0
		return 8


if __name__== "__main__":
	main()
	

	
