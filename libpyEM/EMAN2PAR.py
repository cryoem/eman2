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

# This file contains functions related to running jobs in parallel in EMAN2

DBUG=False		# If set will dump a bunch of debugging output, normally should be False

import os.path
import time
import socket
import sys
import random
import signal
import traceback
import shutil
import subprocess
import thread,threading
import getpass
import select

from EMAN2 import test_image,EMData,abs_path,local_datetime,EMUtil,Util,get_platform
from EMAN2db import e2filemodtime
from EMAN2jsondb import JSTask,JSTaskQueue,js_open_dict
from e2classaverage import ClassAvTask
from e2classifytree import TreeClassifyTask

from e2tiltvalidate import CompareToTiltTask
from e2simmx import EMSimTaskDC
from e2project3d import EMProject3DTaskDC

from e2spt_classaverage import Align3DTask
from e2spt_hac import Align3DTaskAVSA
from e2spt_simulation import SubtomoSimTask
from e2symsearch3d import SymAlignTask
from e2tvrecon import TVReconTask
from e2classifytree import TreeClassifyTask

from e2initialmodel import InitMdlTask
import SocketServer
from cPickle import dumps,loads,dump,load
from struct import pack,unpack

# If we can't import it then we probably won't be trying to use MPI
try :
	from mpi import *
	from mpi_eman import *
except : pass

# if there is a problem with zlib just don't compress
try: from zlib import compress,decompress
except:
	print "Warning: no compression available, please install zlib"
	def compress(s) : return(s)
	def decompress(s) : return(s)

# used to make sure servers and clients are running the same version
EMAN2PARVER=14

# This is the maximum number of active server threads before telling clients to wait
DCMAXTHREADS=7

def DCcustomer_alarm(signum=None,stack=None):
#		if stack!=None : traceback.print_stack(stack)
	print "ALARM - network interrupt"
#	if stack!=None : traceback.print_stack(stack)
	return

class EMTaskCustomer:
	"""This will communicate with the specified task server on behalf of an application needing to
	have tasks completed"""
	def __init__(self,target):
		"""Specify the type and target host of the parallelism server to use.
	dc[:hostname[:port]] - default hostname localhost, default port 9990
	thread:nthreads[:scratch_dir]
	mpi:ncpu[:scratch_dir_on_nodes]
	"""
		target=target.lower()
		self.servtype=target.split(":")[0]

		if self.servtype=="dc" :
			self.addr=target.split(":")[1:]
			if len(self.addr)==0 : self.addr=("",9990)
			elif len(self.addr)==1 : self.addr.append(9990)
			else : self.addr[1]=int(self.addr[1])
			self.addr=tuple(self.addr)
			signal.signal(signal.SIGALRM,DCcustomer_alarm)	# this is used for network timeouts
		elif self.servtype=="thread":
			self.groupn=0
			self.maxthreads=int(target.split(":")[1])
			try: self.scratchdir=target.split(":")[2]
			except: self.scratchdir="/tmp"
			self.handler=EMLocalTaskHandler(self.maxthreads,self.scratchdir)
		elif self.servtype=="mpi":
			self.maxthreads=int(target.split(":")[1])
			try: self.scratchdir=target.split(":")[2]
			except: self.scratchdir="/tmp"
			try:
				# Caching no longer used at all with MPI
				if target.split(":")[3].lower()=="nocache" :
					self.cache=False
				else: self.cache=False
			except: self.cache=False
			self.handler=EMMpiTaskHandler(self.maxthreads,self.scratchdir)
		else : raise Exception,"Only 'dc', 'thread' and 'mpi' servertypes currently supported"

	def __del__(self):
		if self.servtype=="thread" :
			print "Cleaning up thread server. Please wait."
			self.handler.stop()
		elif self.servtype=="mpi" :
			print "Stopping MPI. Please wait."
			try: self.handler.stop()
			except: print "Error: MPI environment was never established, cannot clean it up"

	def wait_for_server(self,delay=10):
		print "%s: Server communication failure, sleeping %d secs"%(local_datetime(),delay)
		time.sleep(delay)
		try:
			x=EMDCsendonecom(self.addr,"TEST",None)
			if (x[0]=="TEST") :
				print "%s: Server is ok now"%local_datetime()
				return
		except: pass

		self.wait_for_server(delay*2)
		return

	def precache(self,filelist):
		"""This will cause a list of files to be precached on the compute nodes before actual computation begins"""
		if self.servtype=="dc" :
			filelist=[abs_path(i) for i in filelist]

			try: EMDCsendonecom(self.addr,"CACH",filelist)
			except:
				self.wait_for_server()
				EMDCsendonecom(self.addr,"CACH",filelist)
		elif self.servtype=="mpi":
			self.handler.precache(filelist)

	def cpu_est(self,wait=True):
		"""Returns an estimate of the number of available CPUs based on the number
		of different nodes we have talked to. Doesn't handle multi-core machines as
		separate entities yet. If wait is set, it will not return until ncpu > 1"""
		if self.servtype =="thread" : return self.maxthreads
		if self.servtype =="mpi" : return self.maxthreads-1

		if self.servtype=="dc" :
			n=0
			while (n==0) :
				try:
					signal.alarm(120)
					n = EMDCsendonecom(self.addr,"NCPU",None)
					signal.alarm(0)
				except:
					signal.alarm(120)
					self.wait_for_server()
					signal.alarm(120)
					n=EMDCsendonecom(self.addr,"NCPU",None)
					signal.alarm(0)
				if not wait : return n
				if n==0 :
					print "Server reports no CPUs available. I will try again in 60 sec"
			return n

	def new_group(self):
		"""request a new group id from the server for use in grouping subtasks"""

		if self.servtype in ("thread","mpi"):
			self.groupn+=1
			return self.groupn

		if self.servtype=="dc" :
			try: return EMDCsendonecom(self.addr,"NGRP",None)
			except:
				self.wait_for_server()
				return EMDCsendonecom(self,addr,"NCPU",None)

	def rerun_task(self,tid):
		"""Trigger an already submitted task to be re-executed"""
		if self.servtype in ("thread","mpi") :
			self.handler.stop()
			raise Exception,"MPI/Threaded parallelism doesn't support respawning tasks"

		if self.servtype=="dc":
			while (1):
				signal.alarm(120)
				try:
					ret=EMDCsendonecom(self.addr,"RQUE",tid)
					signal.alarm(0)
					return ret
				except :
					print "Requeue failure on ",tid
					time.sleep(10)
					continue

	def send_tasks(self,tasks):
		"""Send a group of tasks to the server. Returns a list of taskids."""

		ret=[]

		for task in tasks:
			try: task.user=getpass.getuser()
			except: task.user="anyone"

		if self.servtype in ("thread","mpi"):
			return [self.handler.add_task(t) for t in tasks]


		if self.servtype=="dc" :
			for task in tasks:
				for k in task.data:
					try :
						if task.data[k][0]=="cache" :
							task.data[k]=list(task.data[k])
							task.data[k][1]=abs_path(task.data[k][1])
					except: pass

			# Try to open a socket until we succeed
			try:
				signal.alarm(360)		# sometime there is a lot of congestion here...
				ret=EMDCsendonecom(self.addr,"TSKS",tasks)
				signal.alarm(0)
				return ret
			except:
				traceback.print_exc()
				print "***************************  ERROR SENDING TASK"
				signal.alarm(120)
				self.wait_for_server()
				signal.alarm(0)
				return EMDCsendonecom(self.addr,"TSKS",tasks)

		raise Exception,"Unknown server type"

	def send_task(self,task):
		"""Send a task to the server. Returns a taskid."""

		try: task.user=getpass.getuser()
		except: task.user="anyone"

		if self.servtype in ("thread","mpi"):
			return self.handler.add_task(task)

		if self.servtype=="dc" :
			for k in task.data:
				try :
					if task.data[k][0]=="cache" :
						task.data[k]=list(task.data[k])
						task.data[k][1]=abs_path(task.data[k][1])
				except: pass

			try:
				signal.alarm(120)
				ret=EMDCsendonecom(self.addr,"TASK",task)
				signal.alarm(0)
				return ret
			except:
				traceback.print_exc()
				print "***************************  ERROR SENDING TASK"
				signal.alarm(120)
				self.wait_for_server()
				signal.alarm(0)
				return EMDCsendonecom(self.addr,"TASK",task)

		raise Exception,"Unknown server type"

	def check_task(self,taskid_list):
		"""Check on the status of a list of tasks. Returns a list of ints, -1 to 100. -1 for a task
		that hasn't been started. 0-99 for tasks that have begun, but not completed. 100 for completed tasks."""
		if self.servtype in ("thread","mpi") :
			return self.handler.check_task(taskid_list)

		if self.servtype=="dc":
			try:
				signal.alarm(180)
				ret=EMDCsendonecom(self.addr,"STAT",taskid_list)
				signal.alarm(0)
				return ret
			except:
				signal.alarm(180)
				self.wait_for_server()
				signal.alarm(0)
				return EMDCsendonecom(self.addr,"STAT",taskid_list)
		print self.servtype
		raise Exception,"Unknown server type"

	def get_results(self,taskid,retry=True):
		"""Get the results for a completed task. Returns a tuple (task object,dictionary}."""

		if self.servtype in ("thread","mpi") :
			return self.handler.get_results(taskid)

		if self.servtype=="dc":
			try:
				signal.alarm(180)
				sock,sockf=openEMDCsock(self.addr,retry=10)
				sockf.write("RSLT")
				sendobj(sockf,taskid)
				sockf.flush()

				signal.alarm(180)
				task=recvobj(sockf)

				k=0
				rd={}
				while 1:
					signal.alarm(360)
					k=recvobj(sockf)
					if k==None: break
					v=recvobj(sockf)
					rd[k]=v

#				print "TASK ",task.__dict__
#				print "RESULT ",rd
				sockf.write("ACK ")
				sockf.flush()
				signal.alarm(0)
				return (task,rd)
			except:
				traceback.print_exc()
				if not retry :
					print "************************* Failed to retrieve results, aborting attempt"
					raise Exception,"Unable to retrieve results for %s"%str(taskid)
				print "***************************  ERROR RETRIEVING RESULTS - retrying"
				self.wait_for_server()
				return self.get_results(taskid,False)

class EMTaskHandler:
	"""This is the actual server object which talks to clients and customers. It coordinates task execution
 acts as a data clearinghouse. This parent class doesn't contain any real functionality. Subclasses are always
 used for acutual servers."""
	queue=None

	def __init__(self,path=None):
		if EMTaskHandler.queue==None : EMTaskHandler.queue=JSTaskQueue(path)
		self.queue=EMTaskHandler.queue
#		self.queue=EMTaskQueue(path)

class EMTaskClient:
	"""This class executes tasks on behalf of the server. This parent class implements the actual task functionality.
Communications are handled by subclasses."""
	def __init__(self):
		self.myid=random.randint(1,2000000000)

	def process_task(self,task,callback):
		"""This method implements the actual image processing by calling appropriate module functions"""

		return task.execute(callback)

def image_range(a,b=None):
	"""This is an iterator which handles the (#), (min,max), (1,2,3,...) image number convention for passed data"""
	if b!=None:
		for i in range(a,b): yield i
	elif isinstance(a,int) : yield a
	else:
		for i in a : yield i

class EMTestTask(JSTask):
	"""This is a simple example of a EMTask subclass that actually does something"""

	def __init__(self):
		JSTask.__init__(self)
		self.command="test_task"

	def execute(self):
		# test command. takes a set of images, inverts them and returns the results
		# data should contain one element "input"
		data=self.data["input"]
		cname=data[1]
		cache=db_open_dict(cname)

		ret=[]
		for i in image_range(*data[2:]):		# this allows us to iterate over the specified image numbers
			ret.append(cache[i]*-1)

		return ret


#######################
#  Here are classes for implementing xmlrpc based parallelism

from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler


def runXMLRPCServer(port,verbose):
	"""This will create a ThreadingTCPServer instance and execute it"""
	try: EMDCTaskHandler.verbose=int(verbose)
	except: EMDCTaskHandler.verbose=0

	if port!=None and port>0 :
		server = SimpleXMLRPCServer(("", port),SimpleXMLRPCRequestHandler,False,allow_none=True)
		print "Server started on %s port %d"%(socket.gethostname(),port)
	# EMAN2 will use ports in the range 9900-9999
	else :
		for port in range(9990,10000):
			try:
				server = SimpleXMLRPCServer(("", port),SimpleXMLRPCRequestHandler,False,allow_none=True)
				print "Server started on %s port %d"%(socket.gethostname(),port)
			except:
				if verbose>1 : print "Port %d unavailable"%port
				continue
			break
	server.register_introspection_functions()
	handler=XMLRPCTaskHandler(verbose)
	server.register_instance(handler)

	server.serve_forever()

class XMLRPCTaskHandler(EMTaskHandler):
	def __init__(self, verbose=0):
		# if a port is specified, just try to open it directly. No error detection

		EMTaskHandler.__init__(self)
		self.verbose=verbose

	# Customer routines
	def customer_task_add(self,task):
		"""Enqueue a new task to be executed. Any data passed by reference in the task.data must be
accessible to the server."""
		pass

	def customer_task_check(self,tidlist):
		"""Check a set of tasks for completion. Returns a list of completed task ids from tidlist."""
		pass

	def customer_task_results(self,tid):
		"""Retrieve the return value from a completed task. This should be done only once per task. """
		pass

	# Client routines
	def client_task_get(self):
		"""This is how the client asks the server for work to perform. A task object will be returned."""
		pass

	def client_data_get(self,did):
		pass

	def client_task_status(self,state):
		pass

	# Utility routines
	def test(self,data):
		print "Test message (%d) : %s"%(self.verbose,data)
		return JSTask()

	def quit(self):
		pass

#######################
# Here we define the classes for local threaded parallelism
class EMLocalTaskHandler():
	"""Local threaded Taskserver. This runs as a thread in the 'Customer' and executes tasks. Not a
	subclass of EMTaskHandler for efficient local processing and to avoid data name translation."""
	lock=threading.Lock()
	allrunning = {}	# Static dict of running local tasks. Used for killing thses task upon parent kill
	def __init__(self,nthreads=2,scratchdir="/tmp"):
		self.maxthreads=nthreads
		self.running=[]			# running subprocesses
		self.completed=set()	# completed subprocesses
		self.scratchdir="%s/e2tmp.%d"%(scratchdir,random.randint(1,2000000000))
		self.maxid=0
		self.nextid=0
		self.doexit=0


		os.makedirs(self.scratchdir)
		self.thr=threading.Thread(target=self.run)
		self.thr.start()

	def stop(self):
		"""Called externally (by the Customer) to nicely shut down the task handler"""
		self.doexit=1
		self.thr.join()
		shutil.rmtree(self.scratchdir,True)

	def add_task(self,task):
		EMLocalTaskHandler.lock.acquire()
		if not isinstance(task,JSTask) : raise Exception,"Non-task object passed to EMLocalTaskHandler for execution"
		dump(task,file("%s/%07d"%(self.scratchdir,self.maxid),"wb"),-1)
		ret=self.maxid
		self.maxid+=1
		EMLocalTaskHandler.lock.release()
		return ret

	def check_task(self,id_list):
		"""Checks a list of tasks for completion. Note that progress is not currently
		handled, so results are always -1, 0 or 100 """
		ret=[]
		for i in id_list:
			if i>=self.nextid : ret.append(-1)
			elif i in self.completed : ret.append(100)
			else: ret.append(0)
		return ret

	def get_results(self,taskid):
		"""This returns a (task,dictionary) tuple for a task, and cleans up files"""
#		print "Retrieve ",taskid
		if taskid not in self.completed : raise Exception,"Task %d not complete !!!"%taskid

		task=load(file("%s/%07d"%(self.scratchdir,taskid),"rb"))
		results=load(file("%s/%07d.out"%(self.scratchdir,taskid),"rb"))

		os.unlink("%s/%07d.out"%(self.scratchdir,taskid))
		os.unlink("%s/%07d"%(self.scratchdir,taskid))
		self.completed.remove(taskid)

		return (task,results)

	def run(self):

		while(1):
			time.sleep(1)		# only go active at most 1/second
			if self.doexit==1:
#				shutil.rmtree(self.scratchdir)
				break

			for i,p in enumerate(self.running):
				# Check to see if the task is complete
				if p[0].poll()!=None :

					# This means that the task failed to execute properly
					if p[0].returncode!=0 :
						print "Error running task : ",p[1]
						thread.interrupt_main()
						sys.stderr.flush()
						sys.stdout.flush()
						os._exit(1)

#					print "Task complete ",p[1]
					# if we get here, the task completed
					self.completed.add(p[1])
					del(EMLocalTaskHandler.allrunning[p[1]])

			self.running=[i for i in self.running if i[1] not in self.completed]	# remove completed tasks

			while self.nextid<self.maxid and len(self.running)<self.maxthreads:
#				print "Launch task ",self.nextid
				EMLocalTaskHandler.lock.acquire()
				#There is the issue that when shell=True, popen.pid return shell pid and not process, so we set shell=false (there will be issues on Windows, but we don't support paralellization on windows
				#proc=subprocess.Popen("e2parallel.py" + " localclient" + " --taskin=%s/%07d"%(self.scratchdir,self.nextid) + " --taskout=%s/%07d.out"%(self.scratchdir,self.nextid), shell=True)
				if get_platform() == 'Windows':
					proc=subprocess.Popen(["python", "%s\\bin\\e2parallel.py"%os.getenv('EMAN2DIR'),"localclient","--taskin=%s/%07d"%(self.scratchdir,self.nextid),"--taskout=%s/%07d.out"%(self.scratchdir,self.nextid)])
				else:
					proc=subprocess.Popen(["e2parallel.py","localclient","--taskin=%s/%07d"%(self.scratchdir,self.nextid),"--taskout=%s/%07d.out"%(self.scratchdir,self.nextid)])
				self.running.append((proc,self.nextid))
				EMLocalTaskHandler.allrunning[self.nextid] = proc
				self.nextid+=1
				EMLocalTaskHandler.lock.release()


#######################
#  Here we define the classes for MPI parallelism

class EMMpiClient():
	"""MPI communications are a bit complicated. An instance of EMMpiTaskHandler is created by the
	customer. This object spawns the actual MPI job by executing mpirun. This MPI job as implemented
	in e2parallel.py will make use of the run() method of EMMpiClient for its main loop. MPI rank 0
	is assumed to run on the same host as the EMMpiTaskHandler instance, and will communicate with it
	vi a pair of named pipes (FIFOs). Command requests are queued in the scratch directory and tags
	passed to rank 0 of EMMpiClient, which sends them to appropriate nodes for computation. Rank 0 does
	no computation itself, but is simply a communications hub. The reason for the EMMpiTaskHandler is
	as a wrapper for the mpirun command, which we cannot do at a higher level due to the modular parallelism
	system."""

	def __init__(self,scratchdir="/tmp"):
		mpi_init(0, [])
		self.rank=mpi_comm_rank(MPI_COMM_WORLD)
		self.nrank=mpi_comm_size(MPI_COMM_WORLD)
		self.scratchdir=scratchdir

		self.queuedir=self.scratchdir+"/queue"
		self.cachedir=None

		self.lastupdate=0		# last time we sent an update to rank 0
		if self.rank==0 : self.logfile=file(self.scratchdir+"/rank0.log","w")
		elif self.rank==1 : self.logfile=file(self.scratchdir+"/rank1.log","w")
#		elif self.rank==12 : self.logfile=file(self.scratchdir+"/rank12.log","a")
		else: self.logfile=None

		self.rankmap={}			# key=rank, value=hostname
		self.noderanks={}		# key=hostname, value=rank. Provides one rank/node to be used when precaching
		if DBUG : print "Run EMMpiClient in: ",os.getcwd()
		mpi_barrier(MPI_COMM_WORLD)		# make sure all ranks are up before we start

	def log(self,s):
		if self.logfile!=None:
			self.logfile.write("%s\t%s\n"%(local_datetime(),s))
			self.logfile.flush()

	def test(self,verbose):
		"""This routine will test MPI communications by broadcasting HELO to all of the cpus,
		then waiting for an OK response."""

		if verbose and self.rank==0: print "Testing MPI Communications"
		mpi_barrier(MPI_COMM_WORLD)		# make sure all ranks are up before we start

		# A little test to make sure MPI communications are really established. Also identifies node names
		# for each rank
		if self.rank==0:
			mpi_bcast_send("HELO")

			allsrc=set(range(1,self.nrank))	# we use this to make sure we get a reply from all nodes
			while (1):
				if len(allsrc)==0 : break
				mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD)
				com,b,src=mpi_eman2_recv(MPI_ANY_SOURCE)
				self.log("Rank %d = %s"%(src,b))
				if com!="OK  " :
					print "MPI: Failed receive from node=%d"%src
					mpi_finalize()
					sys.stderr.flush()
					sys.stdout.flush()
					os._exit(1)

				self.rankmap[src]=str(b)[3:]
				self.noderanks[str(b)[3:]]=src		# we just need 1 random rank on each node
				allsrc.remove(src)

			if verbose>1 : print "Successful HELO to all MPI nodes !"
		else:
			a=mpi_bcast_recv(0)
			if a!="HELO" :
				print "MPI: Failed receive on node=%d"%self.rank
				mpi_finalize()
				sys.stderr.flush()
				sys.stdout.flush()
				os._exit(1)
			mpi_eman2_send("OK  ",socket.gethostname(),0)

		mpi_barrier(MPI_COMM_WORLD)		# make sure all ranks are done before we move on


	def run(self,verbose):

		# rank 0 is responsible for communications and i/o, and otherwise does no real work
		if self.rank==0:
			if verbose: print "MPI running on %d processors"%self.nrank
			self.log("MPI on %d processors"%self.nrank)

			# FIFOs to talk to the controlling process, creation order is important !
			#self.fmcon=file("%s/tompi"%self.scratchdir,"rb",0)
			#self.tocon=file("%s/fmmpi"%self.scratchdir,"wb",0)

			# Using a UNIX domain socket due to odd problems with the named FIFO pairs deadlocking
			self.mpisock=socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
			time.sleep(2)
			self.mpisock.connect("%s/mpisock"%self.scratchdir)
			self.mpifile=self.mpisock.makefile()
			self.log("Connected to Controller")

			# Initial handshake to make sure we're both here
			if (self.mpifile.read(4)!="HELO") :
				print "Fatal error establishing MPI controller communications"
				sys.stderr.flush()
				sys.stdout.flush()
				os._exit(1)
			self.log( "Controller said HELO")
			self.mpifile.write("HELO")
			self.mpifile.flush()
			self.log("Said HELO back")

			self.rankjobs=[-1 for i in xrange(self.nrank)]		# Each element is a rank, and indicates which job that rank is currently running (-1 if idle)
			self.rankjobs[0]=-2					# this makes sure we don't try to send a job to ourself
			self.maxjob=-1						# current highest job number waiting for execution
			self.nextjob=1						# next job waiting to run
			self.status={}						# status of each job

			while 1:
				# Look for a command from our controlling process
				if select.select([self.mpifile],[],[],0)[0]:
					com,data=load(self.mpifile)
					if com=="EXIT" :
						dump("OK",self.mpifile,-1)
						self.mpifile.flush()
						self.log("Normal EXIT")
						for i in range(1,self.nrank):
							r=mpi_eman2_send("EXIT","",i)

						self.mpifile.close()
						self.mpisock.close()

						break
					elif com=="NEWJ" :
						dump("OK",self.mpifile,-1)
						self.mpifile.flush()
						if verbose>1 : print "New job %d from customer"%data
						self.log("New job %d from customer"%data)
						self.maxjob=data	# this is the highest number job currently assigned

					elif com=="CHEK" :
						for i in xrange(len(data)):
							if data[i] not in self.status : data[i]=-1
							else: data[i]=self.status[data[i]]

						dump(data,self.mpifile,-1)
						self.mpifile.flush()
					elif com=="CACH" :
						if verbose>1 : print "(ignored) Cache request from customer: "
						dump("OK",self.mpifile,-1)
						self.mpifile.flush()

					else : print "Unknown command from client '%s'"%com
					continue


				# Finally, see if we have any jobs that need to be executed
				if self.nextjob<=self.maxjob :

					if -1 in self.rankjobs :
						rank=self.rankjobs.index(-1)
						if verbose>1 : print "Sending job %d to rank %d (%d idle)"%(self.nextjob,rank,self.rankjobs.count(-1))

						task = file("%s/%07d"%(self.queuedir,self.nextjob),"rb").read()		# we don't unpickle
						self.log("Sending task %d to rank %d (%s)"%(self.nextjob,rank,str(type(task))))
						r=mpi_eman2_send("EXEC",task,rank)

						# if we got here, the task should be running
						self.rankjobs[rank]=self.nextjob
						self.nextjob+=1
						self.log("Sending task rank %d done"%(rank))
						continue

				# Now look for any requests from existing running jobs
				info=mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD)
				if info:
					com,data,src=mpi_eman2_recv(MPI_ANY_SOURCE)
					if com=="DONE" :
						taskid=self.rankjobs[src]
						file("%s/%07d.out"%(self.queuedir,taskid),"wb").write(data)	# assume what we got back is already pickled
						self.status[taskid]=100
						self.rankjobs[src]=-1
						self.log('Task %s complete on rank %d'%(taskid,src))
					elif com=="PROG" :
						if data[1]<0 or data[1]>99 :
							print "Warning: Invalid progress report :",data
						else :
							try : self.status[data[0]]=data[1]
							except: print "Warning: Invalid progress report :",data
					else : print "Warning: unknown task command ",com
					continue

				time.sleep(2)

		# all other ranks handle executing jobs
		else :

			while 1:
				r=mpi_iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD)		# listen for a message from rank 0, this one doesn't block
				#if self.logfile!=None : self.logfile.write( "Listen response (%s)\n"%str(r))
				if r:
					com,data,src=mpi_eman2_recv(0)

					if com=="EXIT":
						if verbose>1 : print "rank %d: I was just told to exit"%self.rank
						break
					
					if com=="EXEC":
						if self.logfile!=None : self.logfile.write( "EXEC\n")
						task=loads(data)		# just for clarity
						if not isinstance(task,JSTask) : raise Exception,"Non-task object passed to MPI for execution ! (%s)"%str(type(task))
						if verbose>1 : print "rank %d: I just got a task to execute (%s):"%(self.rank,socket.gethostname()),task.command,str(task.options)

						self.taskfile="%s/taskexe.%d"%(self.queuedir,os.getpid())
						self.taskout="%s/taskout.%d"%(self.queuedir,os.getpid())

						# Execute the task
						self.task=task	# for the callback
						try: ret=task.execute(self.progress_callback)
						except:
							print "ERROR in executing task"
							traceback.print_exc()
							break

						# return results to rank 0
						if verbose : print "rank %d: Process done :"%self.rank,self.task.taskid
						r=mpi_eman2_send("DONE",dumps(ret,-1),0)

						######## OpenMPI wasn't so hot on the idea of a process that fork()ed, even if it wasn't doing MPI directly
						## Run the job
						#try: os.makedirs("/".join(self.taskfile.split("/")[:-1]))       # this should not be required ?
						#except: pass
						#dump(task,file(self.taskfile,"w"),-1)		# write the task to disk


						#cmd="e2parallel.py localclient --taskin=%s --taskout=%s"%(self.taskfile,self.taskout)

						#self.job=subprocess.Popen(cmd, stdin=None, stdout=None, stderr=None, shell=True)
						#if verbose>2 : print "rank %d: started my job:"%self.rank,cmd

#						if verbose>2 : print "rank %d: running but not done"%self.rank
				else:
					time.sleep(0.1)
#					if self.logfile!=None : self.logfile.write( "Sleep\n")
#					if verbose>2 : print "rank %d: waiting"%self.rank

		mpi_finalize()

	def progress_callback(self,prog):
		""" This gets progress callbacks from the task. We need to make sure we haven't been asked
		to exit if we get this, and we want to update the progress on rank 0 """
		if DBUG : print "progress called",prog
		r=mpi_iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD)		# listen for a message from rank 0
		if r :
			if DBUG: print "probed ",r
			com,data,src=mpi_eman2_recv(0)

			if com=="EXIT":
				print "rank %d: I was just told to exit during processing"%self.rank
				mpi_eman2_send("OK",0,2)
				mpi_finalize()
				sys.exit(0)
			else:
				print "ERROR: Got mysterious command during processing: ",com,data
				mpi_finalize()
				sys.stderr.flush()
				sys.stdout.flush()
				os._exit(1)

		if time.time()-self.lastupdate>120 :
			if DBUG : print "Sending progress"
			mpi_eman2_send("PROG",(self.task.taskid,prog),0)
			self.lastupdate=time.time()

		if DBUG : print "progress done"
		return True

	def pathtocache(self,path):
		"""This will convert a remote filename to a local (unique) cache-file bdb:path"""

		if self.cachedir==None : return path	# caching disabled, return original filename

		# We strip out any punctuation, particularly '/', and take the last 40 characters, which hopefully gives us something unique
		return "%s/%s.hdf"%(self.cachedir,path.translate(None,"#/\\:.!@$%^&*()-_=+")[-40:])

	def pathtocachename(self,path):
		"""This will convert a remote filename to a local (unique) cache-file name"""

		# We strip out any punctuation, particularly '/', and take the last 40 characters, which hopefully gives us something unique
		return str(path.translate(None,"#/\\:.!@$%^&*()-_=+")[-40:]+".hdf")

def filesenum(lst):
	"""This is a generator that takes a list of (name,#), (name,(#,#,#)), or (name,min,max) specifiers
	and yields (name,#) until the list is exhausted"""

	for i in lst:
		for j in fileenum(i): yield j

def fileenum(lst):
	"""This is a generator that takes a single (name,#), (name,(#,#,#)), or (name,min,max) specifier
	and yields (name,#) until the numbers are exhausted"""
	if len(lst)==3 :
		for i in xrange(lst[1],lst[2]) : yield (lst[0],i)
	elif isinstance(lst[1],int) :
		yield lst
	else :
		for i in lst[1]: yield (lst[0],i)

def imgnumenum(lst):
	"""like fileenum, but skips the 'name' references on return"""
	if len(lst)==3 :
		for i in xrange(lst[1],lst[2]) : yield i
	elif isinstance(lst[1],int) :
		yield lst[1]
	else :
		for i in lst[1]: yield i



class EMMpiTaskHandler():
	"""MPI based task handler. This exists as a thread in the customer and handles communications with the actual
	MPI program, which this handler spawns. We do not subclass the EMTaskHandler because we are using our own
	file caching naming scheme here, since the MPI task is not persistent across jobs. If this handler dies,
	all knowledge of running processes dies with it."""
	lock=threading.Lock()
	def __init__(self,ncpus=2,scratchdir="/tmp"):
		try: user=getpass.getuser()
		except: user="anyone"

		self.scratchdir="%s/eman2mpi-%s"%(scratchdir,user)
		self.queuedir=self.scratchdir+"/queue"
		self.cachedir=None
		try: os.makedirs(self.queuedir)
		except:
			#print "makedirs queue failed (%s). May have failed because already exists. "%self.queuedir
			pass

		self.maxid=1			# Current task counter, points to the next open number
		self.completed={}		# set of completed tasks, key is task id, value is completion status

		self.mpiout=file("%s/mpiout.txt"%self.scratchdir,"w")
		self.mpierr=file("%s/mpierr.txt"%self.scratchdir,"w")

		# Using a UNIX domain socket due to odd problems with the named FIFO pairs deadlocking
		try : os.unlink("%s/mpisock"%self.scratchdir)
		except : pass

		self.mpisock=socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
		self.mpisock.bind("%s/mpisock"%self.scratchdir)

		# Launch the MPI subprocess
		mpiopts=os.getenv("EMANMPIOPTS","-n %d"%ncpus)
		cmd="mpirun %s e2parallel.py mpiclient --scratchdir=%s -v 2 </dev/null"%(mpiopts,self.scratchdir)

		self.mpitask=subprocess.Popen(cmd, stdin=None, stdout=self.mpiout, stderr=self.mpierr, shell=True)

		self.mpisock.listen(1)
		self.mpiconn, self.mpiaddr = self.mpisock.accept()
		self.mpifile=self.mpiconn.makefile()

		# order is important here, since opening for reading will block until the other end opens for writing
		#self.tompi=file("%s/tompi"%self.scratchdir,"wb",0)
		#self.fmmpi=file("%s/fmmpi"%self.scratchdir,"rb",0)

		# Send a HELO and wait for a reply. We then know that the MPI system is setup and available
		print "Say HELO to MPI rank 0"
		self.mpifile.write("HELO")
		self.mpifile.flush()
		rd=self.mpifile.read(4)
		if (rd!="HELO") :
			print "Fatal error establishing MPI communications (%s)",rd
			sys.stderr.flush()
			sys.stdout.flush()
			os._exit(1)
		print "Rank 0 said HELO back"

		print "MPI initialized"

	def stop(self):
		"""Called externally (by the Customer) to nicely shut down the task handler"""
		self.sendcom("EXIT")

		self.mpifile.close()
		self.mpiconn.close()
		shutil.rmtree(self.queuedir,True)

	def precache(self,filelist):
		"""Called by the customer to initiate precaching on the nodes, but we no longer support caching for MPI"""
		return

	def sendcom(self,com,data=None):
		"""Transmits a command to MPI rank 0 and waits for a single object in reply"""
		global DBUG
		if DBUG : self.mpiout.write("{} customer sending {}\n".format(local_datetime(),com))
		dump((com,data),self.mpifile,-1)
		self.mpifile.flush()

		if DBUG : self.mpiout.write("{} customer sent\n".format(local_datetime(),com))
		return load(self.mpifile)

	def add_task(self,task):
		if not isinstance(task,JSTask) : raise Exception,"Non-task object passed to EMLocalTaskHandler for execution"

		task.taskid=self.maxid

		dump(task,file("%s/%07d"%(self.queuedir,self.maxid),"wb"),-1)
		ret=self.maxid
		self.sendcom("NEWJ",self.maxid)
		if DBUG : self.mpiout.write("{} customer NEWJ complete {}\n".format(local_datetime(),self.maxid))

		self.maxid+=1

		return ret

	def check_task(self,id_list):
		"""Checks a list of tasks for completion."""

		chk=self.sendcom("CHEK",id_list)
		for i,j in enumerate(id_list):
			self.completed[j]=chk[i]

		return chk

	def get_results(self,taskid):
		"""This returns a (task,dictionary) tuple for a task, and cleans up files"""
#		print "Retrieve ",taskid

		if DBUG : self.mpiout.write("{} customer results {}\n".format(local_datetime(),taskid))

		try :
			task=load(file("%s/%07d"%(self.queuedir,taskid),"rb"))
			results=load(file("%s/%07d.out"%(self.queuedir,taskid),"rb"))
			os.unlink("%s/%07d.out"%(self.queuedir,taskid))
			os.unlink("%s/%07d"%(self.queuedir,taskid))
			del self.completed[taskid]
		except :
			traceback.print_exc()
			print "Error: asked for results on incomplete job"
			return None


		return (task,results)

#######################
#  Here we define the classes for publish and subscribe parallelism

def openEMDCsock(addr,clientid=0, retry=3):
	alrm=signal.alarm(0)
	addr=tuple(addr)
	for i in range(retry):
		try :
			xch="WAIT"
			while xch=="WAIT" :
				if alrm>0 : signal.alarm(alrm+1)
				sock=socket.socket()
				sock.connect(addr)
				sockf=sock.makefile()
				xch=sockf.read(4)
				signal.alarm(0)
				if xch=="WAIT" :
					time.sleep(random.randint(8,30))
					continue
		except:
			time.sleep(8)
			if i>2 : print "Retrying connect to server (%d)"%i
			continue
		if xch!="EMAN" : raise Exception,"Not an EMAN server"
		break
	else: raise Exception,"Exceeded max retries in opening socket to "+str(addr)

	# Introduce ourselves
	#addr=socket.inet_aton(socket.gethostbyname(socket.gethostname()))	# This often returns 127.0.0.1  :^(
	addr=socket.inet_aton(sock.getsockname()[0])
	sockf.write("EMAN")
	sockf.write(pack("<II",EMAN2PARVER,clientid))
	sockf.write(addr)
	sockf.flush()
	if sockf.read(4)=="BADV" :
		print "ERROR: Server version mismatch ",socket.gethostname()
		sys.exit(1)

	if alrm>0 : signal.alarm(alrm+1)

	return(sock,sockf)

global oseq
oseq=1

def broadcast(sock,obj):
	"""This will broadcast an object through a bound datagram socket on port 9989 by serializing the
	the packets into a set of 1k hunks, assums MTU > 1k. Packets contain
	PKLD+uid+totlen+objseq+pktseq+1kdata """
	global oseq
	p=dumps(obj,-1)
	hdr=pack("<4sIII","EMAN",os.getuid(),len(p),oseq)
	for seq in xrange(1+(len(p)-1)/1024):
		r=sock.sendto(hdr+pack("<I",seq)+p[seq*1024:(seq+1)*1024],("<broadcast>",9989))
		if r<0 :
			print "transmit fail %d"%seq
	                r=sock.sendto(hdr+pack("<I",seq)+p[seq*1024:(seq+1)*1024],("<broadcast>",9989))

	for seq in xrange((len(p)-1)/1024,-1,-1):
		sock.sendto(hdr+pack("<I",seq)+p[seq*1024:(seq+1)*1024],("<broadcast>",9989))
#	for seq in xrange(1+(len(p)-1)/1024):
#		sock.sendto(hdr+pack("<I",seq)+p[seq*1024:(seq+1)*1024],("<broadcast>",9989))
	oseq+=1

def recv_broadcast(sock):
	"""This will receive an object sent using broadcast(). If a partial object is received, then a new object starts,
	the first object will be abandoned, and the second (successful) object will be returned"""

	myuid=os.getuid()
	curobjseq=-1
	curpktseq=-1
	while 1:
		pkt=sock.recv(1044)		# 20 byte header + 1024 bytes of data
		mag,uid,totlen,objseq,pktseq=unpack("<4sIIII",pkt[:20])
		if mag !="EMAN" or uid!=myuid: continue			# not for us
		if pktseq!=0 and objseq!=curobjseq  : continue	# for us, but in the middle of a transmission and we missed some
		if pktseq==0:
			curobjseq=objseq
			payload=pkt[20:]
		else :
			if pktseq!=curpktseq+1 : continue			# we missed some, or got the repeat :^(
			payload+=pkt[20:]
		if len(payload)==totlen :
			ret=payload
			payload=""
			return ret
		curpktseq=pktseq

def sendstr(sock,obj):
	"""Sends a string on a socket. obj must be a string or None. A string of zero length is equivalent to None. Use recvstr() to receive"""
	if obj==None or len(obj)==0:
		sock.write(pack("<I",0))
		return
	sock.write(pack("<I",len(obj)))
	sock.write(obj)

def recvstr(sock):
	"""Receives a (compressed) string on a socket. obj must be a string. Use sendstr() to receive"""
	l=sock.read(4)
	try :
		datlen=unpack("<I",l)[0]
	except:
		print "Format error in unpacking (%d) '%s'"%(len(l),l)
		raise Exception,"Network error receiving string"
	if datlen<=0 :return None
	return sock.read(datlen)


def sendobj(sock,obj):
	"""Sends an object as a (binary) size then a binary pickled object to a socket file object"""
	if obj==None :
		sock.write(pack("<I",0))
		return
	strobj=dumps(obj,-1)
	sock.write(pack("<I",len(strobj)))
	sock.write(strobj)

def recvobj(sock):
	"""receives a packed length followed by a binary (pickled) object from a socket file object and returns"""
	l=sock.read(4)
	try :
		datlen=unpack("<I",l)[0]
	except:
		print "Format error in unpacking (%d) '%s'"%(len(l),l)
		raise Exception,"Network error receiving object"
	if datlen<=0 :return None
	return loads(sock.read(datlen))

def EMDCsendonecom(addr,cmd,data,clientid=0):
	"""Connects to an EMAN EMDCServer sends one command, receives a returned object and returns it.
	addr is the standard (host,port) tuple. cmd is a 4 character string, and data is the payload.
	Returns the response from the server."""
	#addr=tuple(addr)
	#sock=socket.socket()
	#sock.connect(addr)
	#sockf=sock.makefile()
	#if sockf.read(4)!="EMAN" : raise Exception,"Not an EMAN server"
	#sockf.write("EMAN")
	#sockf.write(pack("I4s",EMAN2PARVER,cmd))
	sock,sockf=openEMDCsock(addr,clientid=clientid,retry=12)
	sockf.write(cmd)
	sendobj(sockf,data)
	sockf.flush()

	ret=recvobj(sockf)
	sockf.close()
	return ret


import threading
class DCThreadingMixIn:
	"""EMAN2 Mix-in class uses threads, but sets an upper limit on simultaneous threads"""

	# Decides how threads will act upon termination of the
	# main process
	N=1

	def process_request_thread(self, request, client_address):
		"""Same as in BaseServer but as a thread.

		In addition, exception handling is done here.

		"""
		try:
			self.finish_request(request, client_address)
			self.close_request(request)
		except:
			self.handle_error(request, client_address)
			self.close_request(request)

	def process_request(self, request, client_address):
		"""Start a new thread to process the request."""

		N=DCThreadingMixIn.N
		DCThreadingMixIn.N+=1
		count=0

		# Handled in the thread now
		#if threading.active_count()>10:
		#while threading.active_count()>8 and count<10:
			#time.sleep(2)
			#count +=1


		t = threading.Thread(target = self.process_request_thread,
								args = (request, client_address))
		t.setName(str(N))
		t.start()

#class ThreadingTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer): pass
class ThreadingTCPServer(DCThreadingMixIn, SocketServer.TCPServer): pass

def runEMDCServer(port,verbose,killclients=False):
	"""This will create a ThreadingTCPServer instance and execute it"""

	# place to put results
	try:
		os.makedirs("results")
	except:
		pass

	try: EMDCTaskHandler.verbose=int(verbose)
	except: EMDCTaskHandler.verbose=0
	EMDCTaskHandler.killclients=killclients
	EMDCTaskHandler.clients=db_open_dict("bdb:clients")

	if port!=None and port>0 :
		server=None
		while server==None:
			try:
				server = SocketServer.ThreadingTCPServer(("", port), EMDCTaskHandler)	# "" is the hostname and will bind to any IPV4 interface/address
		#		server = SocketServer.TCPServer(("", port), EMDCTaskHandler)	# "" is the hostname and will bind to any IPV4 interface/address
			except :
				print "Port in use, waiting"
				time.sleep(5)
				continue
		if verbose: print "Server started on port ",port
	# EMAN2 will use ports in the range 9900-9999
	else :
		for port in range(9990,10000):
			try:
				server = SocketServer.ThreadingTCPServer(("", port), EMDCTaskHandler)
#				server = SocketServer.TCPServer(("", port), EMDCTaskHandler)
				print "Server started on %s port %d"%(socket.gethostname(),port)
			except:
				if verbose: print "Port %d unavailable"%port
				continue
			break

	if killclients : print "Client killing mode"
	server.serve_forever()

class EMDCTaskHandler(EMTaskHandler,SocketServer.BaseRequestHandler):
	"""Distributed Computing Taskserver. In this system, clients run on hosts with free cycles and request jobs
	from the server, which runs on a host with access to the data to be processed."""
	verbose=0
	rtcount=0
	datacount=0
	lasthk=0
	killclients=False
	tasklock=threading.Lock()
#	dbugfile=file("debug.out","w")
	rotate=('/','-','\\','|')

	def __init__(self, request, client_address, server):
		# if a port is specified, just try to open it directly. No error detection

		EMTaskHandler.__init__(self)
		self.verbose=EMDCTaskHandler.verbose
		if self.verbose>1 : print len(self.queue)
		self.sockf=request.makefile()		# this turns our socket into a buffered file-like object
		SocketServer.BaseRequestHandler.__init__(self,request,client_address,server)
		self.client_address=client_address

	def housekeeping(self):
		# if we haven't heard from a client in 5 minutes, assume it's gone
		for k in EMDCTaskHandler.clients.keys():
			if k=="maxrec" : continue
			c=EMDCTaskHandler.clients[k]

			try:
				if time.time()-c[1]>300 or k==0 :
					if self.verbose : print "Removing client %d"%k
					del EMDCTaskHandler.clients[k]
			except: continue

		for k in self.queue.active.keys():
			j=self.queue.active[k]
			if isinstance(j,int) or j.starttime==None : continue
			try:
				if time.time()-j.starttime>360 and (j.progtime==None or time.time()-j.progtime[0]>360) : raise Exception
			except:
				if self.verbose : print "Task %s doesn't seem to be making progress, restarting"%str(j.taskid)
				self.queue.task_rerun(j.taskid)

	def handle(self):
		"""Process requests from a client. The exchange is:
	send EMAN,version
	recv EMAN,version
	recv command, data len, data (if len not 0)
	send command, data len, data (if len not 0)
	close"""

		if self.verbose>1 : print "Thread %s start"%threading.currentThread().getName()

		# periodic housekeeping, but not if too busy
		if threading.activeCount()>=DCMAXTHREADS :  EMDCTaskHandler.lasthk = time.time()
		if time.time()-EMDCTaskHandler.lasthk>30 :
			EMDCTaskHandler.lasthk=time.time()
			self.housekeeping()

		if self.verbose>1 : print "connection from %s"%(str(self.client_address))

		# the beginning of a message is a struct containing
		# 0-3  : EMAN
		# 4-7  : int, EMAN2PARVER
		# 8-11 : int, client id
		# 12-15 : 4 char command
		# 16-19: count of bytes in pickled data following header

		# if we have too many threads, we trigger the client to sleep a while before even doing a handshake
		# We'll try not to make localhost wait
		if threading.activeCount()>DCMAXTHREADS and self.client_address!="127.0.0.1":
			self.sockf.write("WAIT")
			self.sockf.flush()
			if self.verbose>1 : print "Telling client to wait ",self.client_address
			return

		# initial exchange to make sure we're talking to a client
		self.sockf.write("EMAN")
		self.sockf.flush()
		msg = self.sockf.read(4)
		if msg!="EMAN" : raise Exception,"Non EMAN client"

		ver=unpack("<I",self.sockf.read(4))[0]
		if ver!=EMAN2PARVER :
			self.sockf.write("BADV")
			self.sockf.flush()
			raise Exception,"Version mismatch in parallelism (%d!=%d)"%(ver,EMAN2PARVER)
		self.sockf.write("ACK ")
		self.sockf.flush()
		client_id=unpack("<I",self.sockf.read(4))[0]
		client_addr=socket.inet_ntoa(self.sockf.read(4))

		while (1):
			cmd = self.sockf.read(4)
			if len(cmd)<4 :
				if self.verbose>1 : print "connection closed %s"%(str(client_addr))
				break

			if cmd=="ACK " :
				print "Warning : Out of band ACK"
				continue		# some sort of protocol error, just ignore

			data = recvobj(self.sockf)

			if self.verbose :
				if cmd=="RDYT" :
					EMDCTaskHandler.rtcount+=1
					print " %s  (%d)   \r"%(EMDCTaskHandler.rotate[EMDCTaskHandler.rtcount%4],len(EMDCTaskHandler.clients.bdb)-1),
					sys.stdout.flush()
				elif cmd=="DATA" :
					EMDCTaskHandler.datacount+=1
					if EMDCTaskHandler.datacount%100==0 :
						print"*** %d data   \r"%EMDCTaskHandler.datacount,
						sys.stdout.flush()
				elif cmd=="STAT" :
					EMDCTaskHandler.datacount+=1
					print"? \r",
					sys.stdout.flush()
				elif cmd=="PROG" :
					print "Command %s (%s): %s %s    \r"%(str(client_addr),str(client_id),cmd,str(data)),
					sys.stdout.flush()
				else :
					try: print "Command %s (%s): %s (%d)  "%(str(client_addr),str(client_id),cmd,len(data))
					except: print "Command %s (%s): %s (-)  "%(str(client_addr),str(client_id),cmd)

			######################  These are issued by clients
			# Ready for a task
			if cmd=="RDYT" :
				if EMDCTaskHandler.killclients:
					sendobj(self.sockf,"EXIT")
					self.sockf.flush()
					r=self.sockf.read(4)
					print "Client killed"
					break

				# keep track of clients
#				EMDCTaskHandler.clients[client_id]=(self.client_address[0],time.time(),cmd)
				EMDCTaskHandler.clients[client_id]=(client_addr,time.time(),cmd)

				EMDCTaskHandler.tasklock.acquire()
				if self.queue.caching :
					sendobj(self.sockf,None)			# clients will listen for cache data while idle
					self.sockf.flush()
					r=recvobj(self.sockf)
					if self.verbose>1 : print "Telling client to wait ",client_addr
					EMDCTaskHandler.tasklock.release()
					return

				#### This implements precaching of large data files
				try :
					files=self.queue.precache["files"]
					if files!=None and len(files)>0:
						self.queue.caching=True
						EMDCTaskHandler.tasklock.release()
						lst=["CACHE",[]]
						for i in files:
							lst[1].append(self.queue.todid(i))

						# send the list of objects to the client and get back a list of which ones the client needs
						sendobj(self.sockf,lst)
						self.sockf.flush()
						needed=recvobj(self.sockf)


						# if none are needed, clear out the list of needed files
						if len(needed)==0:
							sendobj(self.sockf,[]) 		# empty chain list
							self.sockf.flush()
							self.queue.precache["files"]=[]
							self.queue.caching=False
						else :
							if self.verbose : print "Precaching ",needed

							# send a list of all clients to try to talk to
							allclients=set()
							for i in EMDCTaskHandler.clients.keys():
								if i=="maxrec" : continue
								allclients.add(EMDCTaskHandler.clients[i][0])
							allclients=list(allclients)
							if self.verbose : print "Clients: ", allclients
							sendobj(self.sockf,allclients)
							self.sockf.flush()

							# now send data for each needed file
							for i in needed:
#								if self.verbose : print "Cache ",i
								name=self.queue.didtoname[i]			# get the filename back from the did
								n=nimg = EMUtil.get_image_count(name)	# how many images to cache in this file
								a=EMData()
								for j in xrange(n):				# loop over images
									a.read_image(name,j)
									xmit=compress(dumps((i[0],i[1],j,a),-1),3)		# compressed pickled string for efficient transfer
									sendstr(self.sockf,xmit)
									self.sockf.flush()
									rsp=self.sockf.read(4)
									if rsp!="ACK " : print "Odd, non-ACK during caching"
									if self.verbose and j%100==0 :
										print "\r Caching %s: %d / %d        "%(name,j+1,n),
										sys.stdout.flush()

						if self.verbose : print "\nDone caching\n"
						sendstr(self.sockf,"DONE")
						self.sockf.flush()
						self.queue.caching=False
						ack=recvobj(self.sockf)
						if ack != "ACK " :
							print "No ack after caching (%s)"%str(ack)
							return
					else : EMDCTaskHandler.tasklock.release()

				except:
					self.queue.caching=False
					traceback.print_exc()

				#### Now we get back to sending an actual task
				EMDCTaskHandler.tasklock.acquire()

				# Get the first task and send it (pickled)
				task=self.queue.get_task(client_id)		# get the next task

				if task==None :
					sendobj(self.sockf,None)			# no tasks available
				else:
					sendobj(self.sockf,task)
					if self.verbose>1: print "Send task: ",task.taskid
				self.sockf.flush()

				# check for an ACK, if not, requeue
				r=["HUH",0]
				try:
					r=recvobj(self.sockf)
					if r[0]!="ACK " :
						EMDCTaskHandler.tasklock.release()
						raise Exception
				except:
					if self.verbose: print "Task sent, no ACK"
					if task!=None : self.queue.task_rerun(task.taskid)
				if task!=None :
					task.exechost=r[1]
					try: self.queue.active[task.taskid]=task
					except: pass
#					EMDCTaskHandler.dbugfile.write("Task %5d sent to %s (%08X)  (%s)  [%s]\n"%(task.taskid,client_addr,client_id,r,local_datetime()))
#					EMDCTaskHandler.dbugfile.flush()
				EMDCTaskHandler.tasklock.release()


			# Job is completed, results returned
			elif cmd=="DONE" :
				# the first object we get is the task identifier
				tid=data

				if self.verbose>1 : print "DONE task (%08X) "%client_id,tid

				# keep track of clients
				EMDCTaskHandler.clients[client_id]=(client_addr,time.time(),cmd)

				# then we get a sequence of key,value objects, ending with a final None key
#				result=db_open_dict("bdb:%s#result_%d"%(self.queue.path,tid))
				result=file("%s/results/%07d"%(self.queue.path,tid),"w")

				cnt=0
				while (1) :
					try:
						key=recvobj(self.sockf)
						if key==None : break
						val=recvobj(self.sockf)
					except:
						traceback.print_exc()
						print "Error in network communications"
						cnt=-1
						break
#					result[key]=val				# store each key in the database
					dump(key,result,-1)
					dump(val,result,-1)
					cnt+=1
					if self.verbose>3: print key,val

				result.close()
				if self.verbose>2 : print "Task %d: %d data elements"%(tid,cnt)

				if cnt>=0 : self.queue.task_done(tid)		# don't mark the task as done until we've stored the results
#				EMDCTaskHandler.dbugfile.write("Task %5d complete from %08X with %d data elements  [%s]\n"%(tid,client_id,cnt,local_datetime()))
#				EMDCTaskHandler.dbugfile.flush()

			# Request data from the server
			# the request should be of the form ["queue",did,image#]
			# Returns requested object
			elif cmd=="DATA" :
				# if the server is congested, we slow down the data flow to the clients a bit to produce better threading
				#if threading.active_count()>DCMAXTHREADS :
					#if random.randint(1,10)==1 : time.sleep(1)
				try:
					fsp=self.queue.didtoname[data[1]]
					obj=EMData(fsp,data[2])
					sendobj(self.sockf,obj)
					if self.verbose>2 : print "Data sent %s(%d)"%(fsp,data[2])
				except:
					sendobj(self.sockf,None)
					if self.verbose : print "Error sending %s(%d)"%(fsp,data[2])
				self.sockf.flush()
				EMDCTaskHandler.clients[client_id]=(client_addr,time.time(),cmd)

			# Notify that a task has been aborted
			# request should be taskid
			# no return
			elif cmd=="ABOR" :
				self.queue.task_rerun(data)
				if self.verbose : print "Task execution abort : ",data
				del EMDCTaskHandler.clients[client_id]		# assume this client is going away unless we hear from it again

			# Progress message indicating that processing continues
			# data should be (taskid,percent complete)
			# returns "OK  " if processing should continue
			# retunrs "ABOR" if processing should be aborted
			elif cmd=="PROG" :
				if EMDCTaskHandler.killclients:
					sendobj(self.sockf,"ABOR")
					self.sockf.flush()
					print "Client killed"
					break

				ret=self.queue.task_progress(data[0],data[1])
				if self.verbose>2 : print "Task progress report : ",data
				if ret : sendobj(self.sockf,"OK  ")
				else : sendobj(self.sockf,"ABOR")
				self.sockf.flush()

				# keep track of clients
				EMDCTaskHandler.clients[client_id]=(client_addr,time.time(),cmd)

			###################### These are utility commands
			# Returns whatever is sent as data
			elif cmd=="TEST":
				sendobj(self.sockf,("TEST",data))
				self.sockf.flush()

			# Cause this server to exit (cleanly)
			# actually this may be partially broken...
			elif cmd=="QUIT" :
				sendobj(self.sockf,None)
				self.sockf.flush()
				self.server.server_close()
				if self.verbose : print "Server exited cleanly"
				break
#				sys.exit(0)

			################### These are issued by customers
			# A new task to enqueue
			# the request contains the EMTask object
			# return is the task id or None upon error
			elif cmd=="TASK":
				tid=self.queue.add_task(data)
				if self.verbose>1 : print "new TASK %s.%s"%(data.command,str(data.data))
				try:
					sendobj(self.sockf,tid)
					self.sockf.flush()
				except:
					sendobj(self.sockf,None)
					self.sockf.flush()

#				EMDCTaskHandler.dbugfile.write("Task %5d queued [%s]\n"%(tid,local_datetime()))

			################### These are issued by customers
			# A set of new tasks to enqueue
			# the request contains a list of EMTask objects
			# return is a list of task ids or None upon error
			elif cmd=="TSKS":
				grp=self.queue.add_group()
				tids=[]
				for task in data:
					task.group=grp
					tids.append(self.queue.add_task(task))
					if self.verbose>1 : print "new TASK %s.%s"%(task.command,str(data.data))
				try:
					sendobj(self.sockf,tids)
					self.sockf.flush()
				except:
					sendobj(self.sockf,None)
					self.sockf.flush()

#				EMDCTaskHandler.dbugfile.write("Tasks %5d - %5d queued [%s]\n"%(tids[0],tids[-1],local_datetime()))

			# This will requeue a task which has already been run (and may be completed)
			# This will generally happen if there was some unexpected error in the returned results
			# and should not ever really be necessary. It is likely indicative of some sort of
			# problem with one of the compute nodes
			elif cmd=="RQUE":
				if data!= None :
					self.queue.task_rerun(data)
					if self.verbose : print "Requeuing ",data
					sendobj(self.sockf,"OK")
					self.sockf.flush()
				else : print "ERROR: tried to requeue None"

			# This will cause files to be precached on the clients before running more tasks
			elif cmd=="CACH" :
				self.queue.precache["files"]=data
				sendobj(self.sockf,None)
				self.sockf.flush()
				if self.verbose : print "Accepted list for precaching ",data

			# Get an estimate of the number of CPUs available to run jobs
			# At the moment, this is the number of hosts that have communicated with us
			# so it doesn't handle multiple cores
			elif cmd=="NCPU" :
				sendobj(self.sockf,len(EMDCTaskHandler.clients.bdb)-1)
				self.sockf.flush()
				if self.verbose : print len(EMDCTaskHandler.clients.bdb)-1," clients reported"

			elif cmd=="NGRP" :
				sendobj(self.sockf,self.queue.add_group())
				self.sockf.flush()

			# Cancel a pending task
			# request contains the taskid to cancel
			elif cmd=="CNCL":
				self.queue.task_aborted(data)

			# status request
			# request contains a list/tuple of taskids
			# return is a list/tuple with the same number of elements containing % complete for each task
			# if percent complete is exactly 100, then results are ready for pickup. if -1, task not yet started
			elif cmd=="STAT":
				ret=[self.queue.task_check(i) for i in data]
				sendobj(self.sockf,ret)
				self.sockf.flush()

			# Retreieve results for completed task
			# request contains taskid
			# return is a task object followed by a series of key/value pairs terminating with a None key
			# initial return None if task incomplete
			elif cmd=="RSLT":
				try: ret=self.queue.complete[data]
				except:
					sendobj(self.sockf,None)
					self.sockf.flush()
					continue
				sendobj(self.sockf,ret)

				if self.verbose>2 : print "RSLT: ",data
#				result=db_open_dict("bdb:%s#result_%d"%(self.queue.path,data))
				#for k in result.keys():
					#if k=="maxrec": continue
					#sendobj(self.sockf,k)
					#sendobj(self.sockf,result[k])
					#if self.verbose>3 : print k,result[k]

				result=file("%s/results/%07d"%(self.queue.path,data),"r")
				while 1:
					try:
						k=load(result)
						v=load(result)
						sendobj(self.sockf,k)
						sendobj(self.sockf,v)
					except:
						break
					if self.verbose>3 : print k,v
				result.close()

				sendobj(self.sockf,None)
				self.sockf.flush()

				try:
					if self.sockf.read(4)!="ACK " : raise Exception
#					db_remove_dict("bdb:%s#result_%d"%(self.queue.path,data))
					os.unlink("%s/results/%07d"%(self.queue.path,data))
#					EMDCTaskHandler.dbugfile.write("Results for task %5d retrieved [%s]\n"%(data,local_datetime()))
				except:
					if self.verbose: print "No ACK on RSLT. Keeping results for retry."
#					EMDCTaskHandler.dbugfile.write("Results for task %5d FAILED [%s]\n"%(data,local_datetime()))


			else :
				sendobj(self.sockf,"ERROR: Unknown command")
				self.sockf.flush()

		if self.verbose>1 : print "Thread %s exit"%threading.currentThread().getName()
# self.request is the TCP socket connected to the client
#		self.data = self.request.recv(1024).strip()
#		print "%s wrote:" % self.client_address[0]
#		print self.data
#		# just send back the same data, but upper-cased
#		self.request.send(self.data.upper())

def DCclient_alarm(signum=None,stack=None):
	"""for normal use"""
#		if stack!=None : traceback.print_stack(stack)
	print "ALARM - network interrupt"
	if stack!=None : traceback.print_stack(stack)
	return

def DCclient_alarm2(signum=None,stack=None):
	"""for the special case of precaching"""
#		if stack!=None : traceback.print_stack(stack)
	raise Exception,"timeout"


class EMDCTaskClient(EMTaskClient):
	"""Distributed Computing Task Client. This client will connect to an EMDCTaskServer, request jobs to run
 and run them ..."""

	def __init__(self,server,port,verbose=0,myid=None):
		EMTaskClient.__init__(self)
		self.addr=(server,port)
		self.verbose=verbose
		self.lastupdate=0
		self.task=None
		if myid!=None : self.myid=myid
		signal.signal(signal.SIGALRM,DCclient_alarm)	# this is used for network timeouts

	def imalive(self,progress):
		"""Executed code should call this periodically to inidicate that they are still running. This must be called at least once every 5 minutes
		with an integer in the range 0-100. 100 means about to exit, 0 means not started yet. A -1 indicates an error has occurred and the request
		is being aborted. If this function returns False, the client should abort the task in progress."""
		ret=True

		signal.alarm(60)
		if time.time()-self.lastupdate>120 :
			if self.task==None :
				signal.alarm(0)
				return True
			try:
				ret=EMDCsendonecom(self.addr,"PROG",(self.task.taskid,progress),clientid=self.myid)
				self.lastupdate=time.time()
			except: pass

		signal.alarm(0)
		if ret=="ABOR": return False

		return True

	def cachewriter(self,cq):
		"""This writes the incoming broadcast cache data in a separate thread so we don't miss any broadcasts.
		Not using it any more with new chain scheme."""

		n=0
		lname=""
		while 1:
			if len(cq)==0 : continue
#			print len(cq), " in cache list"
			img=cq.pop()
			if img==None :
				if len(cq)==0 : break
				cq.insert(0,None)

			# The data item should be a pickled tuple (time,rand,img#,image)
			try : img=loads(decompress(img))
			except : continue			# bad pickle :^(
			try : cname="bdb:cache_%d.%d"%(img[0],img[1])
			except : continue			# data wasn't what we expected

			if cname!=lname :
				if self.verbose : print "Receiving cache data ",cname
				f=db_open_dict(cname)
				lname=cname
			f[img[2]]=img[3]
#			print "> ",img[2],cname
			n+=1

		if self.verbose and n>0: print n," items cached"

	def connectfromlist(self,hostlist):
		"""Given a list of possible hostnames to connect to, try to connect to each in sequence until one
		answers. Transmit the remainder of the list to the host after connection. Returns the connected socket or
		None if no connection was sucessful. """

		fail=1
		while fail :
			try:
				nexthost=hostlist.pop()
			except:
				sockout,sockoutf=None,None
				break
			sockout=socket.socket()
			signal.alarm(5)
			try :
				sockout.connect((nexthost,9989))
				sockoutf=sockout.makefile()
				sendobj(sockoutf,hostlist)				# First thing we do is send the next node in the chain a list of the remaining nodes
				sockoutf.flush()
				fail=0
			except:
#				traceback.print_exc()
				print "connect %s to %s failed"%(socket.gethostname(),nexthost)

#		if sockout!=None : print "connect %s to %s"%(socket.gethostname(),nexthost)
		signal.alarm(0)
		return (sockout,sockoutf)

	def listencache(self):
		"""This will listen for cached data (idle for up to 30 seconds) or sleep for 30 seconds if someone else on this node is listening"""

		try:
			sock=socket.socket()
			sock.bind(("",9989))
#			print "listening"
		except:
#			print "Sleeping (not listening)"
			sock=None
			time.sleep(60)
			return

		signal.signal(signal.SIGALRM,DCclient_alarm2)	# this is used for network timeouts
		try:
#			signal.alarm(30)
			sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVTIMEO, pack('LL', 45, 0))
			sock.listen(1)
			sock2=sock.accept()[0]
			sockf=sock2.makefile()
#			print "connection !"
		except:
			sock=None
			return
#		sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVTIMEO, pack('LL', 30, 0))

		signal.signal(signal.SIGALRM,DCclient_alarm)	# this is used for network timeouts
		#cq=[]	# a list of returned images, gradually written by the thread
		#thr=threading.Thread(target=self.cachewriter,args=(cq,))
		#thr.start()

#		signal.signal(signal.SIGALRM,DCclient_alarm2)	# this is used for network timeouts

		chainlist=recvobj(sockf)		# The first thing we get is a list of other clients to send data to
		#try: chainlist.remove(socket.gethostbyname(socket.gethostname()))
		try: chainlist.remove(sock2.getsockname()[0])
		except:pass
#		print "chainlist: ",chainlist
		sockout,sockoutf=self.connectfromlist(chainlist)		# try to establish a connection to one of them

		try:
			ret=None
			lname=""
			nrecv=0
			while 1 :
				if nrecv==0 : signal.alarm(15)		# if we haven't gotten anything yet, only wait a little while
				else : signal.alarm(30)				# if we've been receiving data, wait longer for more

#				cq.append(recv_broadcast(sock))		# too slow
#				ret=Util.recv_broadcast(sock.fileno())

				ret=None
				try:
					ret=recvstr(sockf)		# this should contain time,rint,img#,image in a compressed pickled string

					if ret=="DONE" or len(ret)==0: break
					sockf.write("ACK ")			# We do the ACK immediately along the chain for performance
					sockf.flush()
#					print "Got object (%d)"%len(ret)
				except:
					traceback.print_exc()
					print "**** error receiving data on ",socket.gethostname()
					break

				signal.alarm(60)				# Longer timeout for retransmitting down the chain
				# Send the image down the chain
				if sockout!=None :
					try:
						sendstr(sockoutf,ret)
						sockoutf.flush()
					except:
						traceback.print_exc()
						print "Chain broken ! (%s)"%socket.gethostname()
						sockout=None

				# The data item should be a pickled tuple (time,rand,img#,image)
				try : img=loads(decompress(ret))
				except :
					print "ERROR (%s): Bad data on chain"%socket.gethostname()
					continue			# bad pickle :^(
				try : cname="bdb:cache_%d.%d"%(img[0],img[1])
				except :
					print "ERROR (%s): Bad data object on chain"%socket.gethostname()
					continue			# data wasn't what we expected

				if cname!=lname :
					if self.verbose : print "Receiving cache data ",cname
					f=db_open_dict(cname)
					lname=cname
				f[img[2]]=img[3]
				nrecv+=1

				if sockout!=None :
					try:
						if sockoutf.read(4)!="ACK " : raise Exception
					except:
						traceback.print_exc()
						print "Chain broken (NACK) ! (%s)"%socket.gethostname()
						sockout=None


		except:
			traceback.print_exc()		# we shouldn't really get an exception here
			print "**** Exception in outer chain loop"

		if self.verbose :
			print nrecv," total items cached"

		# Tell the chain we're done
		if sockout!=None :
			try:
				sendstr(sockoutf,"DONE")
				sockoutf.flush()
#				sockoutf.close()
#				sockout.close()
			except: pass

		signal.signal(signal.SIGALRM,DCclient_alarm)	# this is used for network timeouts
#		print "Done listening"

		signal.alarm(0)
		#thr.join()			# wait for cache writing to complete

	def docache(self,sock,sockf,clist):
		"""This routine will receive data to cache from the server, then transmit it to the next server in line. This used
		to use a network broadcast mechanism, but it proved to be too unreliable"""
#		if self.verbose: print "Caching starting ",clist
		needed=[]
		for i in clist[1] :			# loop over list of files to cache
			cname="bdb:cache_%d.%d"%(i[0],i[1])
			cache=db_open_dict(cname)
			try: z=cache.get_header(0)
			except: needed.append(i)
			if z==None : needed.append(i)


		sendobj(sockf,needed)
		sockf.flush()
		if self.verbose : print "Caching phase, need : ",needed
#		if len(needed)==0 : return

		#bcast=socket.socket(socket.AF_INET,socket.SOCK_DGRAM)	# One client will broadcast the data on its subnet for mass distribution
		#bcast.bind(("",9989))
		#bcast.setsockopt(socket.SOL_SOCKET,socket.SO_BROADCAST,1)
		chainlist=recvobj(sockf)		# The first thing we get back is a list of other clients to send data to
		#try: chainlist.remove(socket.gethostbyname(socket.gethostname()))
		try: chainlist.remove(sock.getsockname()[0])
		except: pass
#		print "chainlist (%s): %s"%(socket.gethostbyname(socket.gethostname()),chainlist)
		sockout,sockoutf=self.connectfromlist(chainlist)		# try to establish a connection to one of them

		lname=""
		n=0
#		t0,t1,t2,t3=0,0,0,0
		# This loop receives the data from the server, then forwards it to the next host in the chain
		while 1:
			signal.alarm(60)
			xmit=recvstr(sockf)		# this should contain time,rint,img#,image

			if xmit=="DONE" : break
			try : img=loads(decompress(xmit))
			except : print "ERROR : Bad cache data from server"

			# immediate ACK so the server can prepare the next packet
			sockf.write("ACK ")
			sockf.flush()

			try : cname="bdb:cache_%d.%d"%(img[0],img[1])
			except :
				print "Invalid cache data '",img,"'"
				break
			if cname!=lname :
				if self.verbose : print "Receiving cache data ",cname
				f=db_open_dict(cname)
				lname=cname

			# Send the image down the chain
			if sockout!=None :
				try:
#					print "Sending down chain"
					sendstr(sockoutf,xmit)
					sockoutf.flush()
					if sockoutf.read(4)!="ACK " : raise Exception
				except:
					print "Chain broken ! (%s)"%socket.gethostname()
					sockout=None

			f[img[2]]=img[3]		# Save the image in the local cache

			n+=1

		# Tell the chain we're done
		if sockout!=None :
			try:
				sendstr(sockoutf,"DONE")
				sockoutf.flush()
#				sockout.close()
			except: pass

		signal.alarm(0)
#		print "Done precaching"


	def run(self,dieifidle=86400,dieifnoserver=86400,onejob=False):
		"""This is the actual client execution block. dieifidle is an integer number of seconds after
		which to terminate the client if we aren't given something to do. dieifnoserver is the same if we can't find a server to
		communicate with. Both default to 24 hours."""
		count=0
		lastjob=time.time()
		lastserv=time.time()
		retcode=2
		while (1):
			count +=1
			# connect to the server
			if self.verbose>1 : print "Connect to (%s,%d)"%self.addr
			try :
				signal.alarm(60)
				sock,sockf=openEMDCsock(self.addr,clientid=self.myid,retry=3)
				sockf.write("RDYT")
				sendobj(sockf,None)
				sockf.flush()

				# Get a task from the server
				task=recvobj(sockf)
				if self.verbose>1 : print "Task: ",task

				# This means the server wants to use us to precache files on all of the clients, we won't
				# get a task until we finish this
				if isinstance(task,list) and task[0]=="CACHE" :
					self.docache(sock,sockf,task)
					sendobj(sockf,"ACK ")
					sockf.flush()
					task=recvobj(sockf)

				# acknowledge the task even before we know what we got
				sendobj(sockf,("ACK ",socket.gethostname()))
				sockf.flush()
				signal.alarm(0)
				if task=="EXIT" :
					retcode=1
					break
				if self.verbose and task!=None: print "%s running task id %d"%(socket.gethostname(),task.taskid)
				self.task=task
			except :
				print "No response from server, sleeping 30 sec"
				if time.time()-lastserv>dieifnoserver :
					print "No server for too long. Terminating"
					retcode=3
					break
				time.sleep(30)
				continue
			lastserv=time.time()

			if task==None:
				if self.verbose :
					if count%1==0 : print " | \r",
					else : print " - \r",
					sys.stdout.flush()
				sockf.close()
				sock.close()
				if time.time()-lastjob>dieifidle :
					print "Idle too long. Terminating"
					retcode=4
					break
				self.listencache()		# We will listen for precached data for 15 seconds (or sleep if another thread is listening)
				continue
			sockf.flush()

			# Translate and retrieve (if necessary) data for task
			for k,i in task.data.items():
				if self.verbose>1 : print "Data translate ",k,i
#				try:
				if isinstance(i,list) and len(i)>0 and i[0]=="cache" :
					cname="bdb:cache_%d.%d"%(i[1][0],i[1][1])
					cache=db_open_dict(cname)
					if self.verbose>2 : print "Open cache : ",cname

					for j in image_range(*i[2:]):
#						if not cache.has_key(j):
						try:
							z=cache.get_header(j)
							if z==None : raise Exception
						except:
							cache[j]=self.get_data(sockf,i[1],j)
#							print j,cache[j]

						z=cache[j]		# this is here to raise an exception if something went funny in the data retrieval

					i[1]=cname
#				except: pass
			sockf.close()
			sock.close()

			# Execute translated task
#			try:
			ret=self.process_task(task,self.imalive)
#			except Exception,err:
#				ret={"error (%d)"%task.taskid:err}

			if ret==None :
				if self.verbose : print "Task aborted %d"%task.taskid
				sock,sockf=openEMDCsock(self.addr,clientid=self.myid,retry=10)
				sockf.write("ABOR")
				sockf.close()
				lastjob=time.time()
				continue

			# Return results
			if self.verbose : print "Task done %d"%task.taskid
			if self.verbose>3 : print self.__dict__

			retry=True
			retrycount=0
			while retry:
				signal.alarm(120)
				retry=False
				retrycount+=1
				if retrycount>10 :
					print "Failed in 10 attempts to send results, aborting (%d)"%task.taskid
					retcode=10
					break

				try:
					sock,sockf=openEMDCsock(self.addr,clientid=self.myid,retry=10)
					sockf.write("DONE")
				except:
					print "Server communication failure, trying again in 1 minute (%d)"%task.taskid
					time.sleep(60)
					retry=True
					continue

				signal.alarm(120)
				try:
					sendobj(sockf,task.taskid)
					sockf.flush()
				except:
					print "Immediate ERROR (retrying ",task.taskid,")"
					retry=True
					continue

				for k,v in ret.items():
					signal.alarm(120)
					try:
						sendobj(sockf,k)
						sendobj(sockf,v)
					except :
						print "ERROR (retrying ",task.taskid,") on : ",k, " in ",ret.items()
						if isinstance(v,EMData) : v.write_image("error.hdf",-1)
						time.sleep(3)
						retry=True
						retcode=11
						break

				try:
					sendobj(sockf,None)
					sockf.flush()
				except:
					print "Error on flush (%d)"%task.taskid
					retry=True

				try:
					sockf.close()
					sock.close()
				except:
					print "Error on close (%d)"%task.taskid
					retry=True

			signal.alarm(0)

			if retrycount>10 : break			# if we completely failed once, this client should die

			if self.verbose : print "Task returned %d"%task.taskid
			if self.verbose>2 : print "Results :",ret

			if onejob :
				retcode=0		# this means it's ok to spawn more jobs, the current job was a success
				break

			lastjob=time.time()
			time.sleep(3)

		if retcode : print "Client on %s exiting. Bye !"%socket.gethostname()
		return retcode

	def get_data(self,sockf,did,imnum):
		"""request data from the server on an open connection"""
		signal.alarm(240)
		if self.verbose>2 : print "Retrieve %s : %d"%(str(did),imnum)
		sockf.write("DATA")
		sendobj(sockf,["cache",did,imnum])
		sockf.flush()
		try: ret=recvobj(sockf)
		except :
			print "ERROR on %s: could not retrieve %s : %d, retry"%(socket.gethostname(),str(did),imnum)
			time.sleep(3)
			sockf.write("DATA")
			sendobj(sockf,["cache",did,imnum])
			sockf.flush()
			ret=recvobj(sockf)
		signal.alarm(0)

		return ret

	#		self.sockf=

