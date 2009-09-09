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

from EMAN2 import test_image,EMData,abs_path,local_datetime
from EMAN2db import EMTask,EMTaskQueue,db_open_dict,db_remove_dict
from e2classaverage import EMClassAveTaskDC
from e2simmx import EMSimTaskDC
from e2project3d import EMProject3DTaskDC
from e2tomoaverage import EMTomoAlignTaskDC
import SocketServer
from cPickle import dumps,loads,dump,load
from struct import pack,unpack
import os.path
import time
import socket
import sys
import random
import signal
import traceback

# used to make sure servers and clients are running the same version
EMAN2PARVER=11

# This is the maximum number of active server threads before telling clients to wait
DCMAXTHREADS=10

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
	"""
		target=target.lower()
		self.servtype=target.split(":")[0]
		if self.servtype!="dc" : raise Exception,"Only 'dc' servertype currently supported"
		
		self.addr=target.split(":")[1:]
		if len(self.addr)==0 : self.addr=("",9990)
		elif len(self.addr)==1 : self.addr.append(9990)
		else : self.addr[1]=int(self.addr[1])
		self.addr=tuple(self.addr)
		signal.signal(signal.SIGALRM,DCcustomer_alarm)	# this is used for network timeouts

	def wait_for_server(self,delay=10):
		print "%s: Server communication failure, sleeping %d secs"%(time.ctime(),delay)
		time.sleep(delay)
		try:
			x=EMDCsendonecom(self.addr,"TEST",None)
			if (x[0]=="TEST") : 
				print "%s: Server is ok now"%time.ctime()
				return
		except: pass
		
		self.wait_for_server(delay*2)
		return

	def cpu_est(self,wait=True):
		"""Returns an estimate of the number of available CPUs based on the number
		of different nodes we have talked to. Doesn't handle multi-core machines as
		separate entities yet. If wait is set, it will not return until ncpu > 1"""
		if self.servtype=="dc" :
			n=0
			while (n==0) :
				try:
					signal.alarm(60)
					n = EMDCsendonecom(self.addr,"NCPU",None)
					signal.alarm(0)
				except:
					signal.alarm(60)
					self.wait_for_server()
					signal.alarm(60)
					n=EMDCsendonecom(self.addr,"NCPU",None)
					signal.alarm(0)
				if not wait : return n
				if n==0 : 
					print "Server reports no CPUs available. I will try again in 60 sec"
			return n

	def new_group(self):
		"""request a new group id from the server for use in grouping subtasks"""
		
		if self.servtype=="dc" :
			try: return EMDCsendonecom(self.addr,"NGRP",None)
			except:
				self.wait_for_server()
				return EMDCsendonecom(self,addr,"NCPU",None)

	def rerun_task(self,tid):
		"""Trigger an already submitted task to be re-executed"""
		
		if self.servtype=="dc":
			signal.alarm(60)
			try: return EMDCsendonecom(self.addr,"RQUE",tid)
			except:
				signal.alarm(60)
				self.wait_for_server()
				ret=EMDCsendonecom(self.addr,"RQUE",tid)
				signal.alarm(0)
				return ret
				
			signal.alarm(0)

	def send_task(self,task):
		"""Send a task to the server. Returns a taskid."""
		
		try: task.user=os.getlogin()
		except: task.user="unknown"
		
		for k in task.data:
			try :
				if task.data[k][0]=="cache" :
					task.data[k]=list(task.data[k])
					task.data[k][1]=abs_path(task.data[k][1])
			except: pass
		
		if self.servtype=="dc" :
			try: 
				signal.alarm(60)
				ret=EMDCsendonecom(self.addr,"TASK",task)
				signal.alarm(0)
				return ret
			except:
				traceback.print_exc()
				print "***************************  ERROR SENDING TASK"
				signal.alarm(60)
				self.wait_for_server()
				signal.alarm(0)
				return EMDCsendonecom(self.addr,"TASK",task)

		raise Exception,"Unknown server type"

	def check_task(self,taskid_list):
		"""Check on the status of a list of tasks. Returns a list of ints, -1 to 100. -1 for a task
		that hasn't been started. 0-99 for tasks that have begun, but not completed. 100 for completed tasks."""
		if self.servtype=="dc":
			try:
				signal.alarm(60)
				ret=EMDCsendonecom(self.addr,"STAT",taskid_list)
				signal.alarm(0)
				return ret
			except:
				signal.alarm(60)
				self.wait_for_server()
				signal.alarm(0)
				return EMDCsendonecom(self.addr,"STAT",taskid_list)
		print self.servtype
		raise Exception,"Unknown server type"
	
	def get_results(self,taskid):
		"""Get the results for a completed task. Returns a tuple with the task object and dictionary."""
		if self.servtype=="dc":
			try:
				signal.alarm(60)
				sock,sockf=openEMDCsock(self.addr,retry=10)
				sockf.write("RSLT")
				sendobj(sockf,taskid)
				sockf.flush()
				
				signal.alarm(60)
				task=recvobj(sockf)
				
				k=0
				rd={}
				while 1:
					signal.alarm(240)
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
				print "***************************  ERROR RETRIEVING RESULTS - retrying"
				self.wait_for_server()
				return self.get_results(taskid)

class EMTaskHandler:
	"""This is the actual server object which talks to clients and customers. It coordinates task execution
 acts as a data clearinghouse. This parent class doesn't contain any real functionality. Subclasses are always
 used for acutual servers."""
	
	def __init__(self,path=None):
		self.queue=EMTaskQueue(path)
		
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

class EMTestTask(EMTask):
	"""This is a simple example of a EMTask subclass that actually does something"""
	
	def __init__(self):
		EMTask.__init__(self)
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
		return EMTask()

	def quit(self):
		pass
	

#######################
#  Here we define the classes for publish and subscribe parallelism

def openEMDCsock(addr,clientid=0, retry=3):
	alrm=signal.alarm(0)
	addr=tuple(addr)
	for i in range(retry):
		try :
			xch="WAIT"
			while xch=="WAIT" :
				if alrm>0 : signal.alarm(alrm)
				sock=socket.socket()
				sock.connect(addr)
				sockf=sock.makefile()
				xch=sockf.read(4)
				signal.alarm(0)
				if xch=="WAIT" :
					time.sleep(random.randint(5,20))
					continue
		except:
			time.sleep(8)
			if i>2 : print "Retrying connect to server (%d)"%i
			continue
		if xch!="EMAN" : raise Exception,"Not an EMAN server"
		break
	else: raise Exception,"Exceeded max retries in opening socket to "+str(addr)

	# Introduce ourselves and ask for a task to execute
	sockf.write("EMAN")
	sockf.write(pack("I4",EMAN2PARVER))
	sockf.write(pack("I4",clientid))
	sockf.flush()
	if sockf.read(4)=="BADV" : 
		print "ERROR: Server version mismatch ",socket.gethostname()
		sys.exit(1)

	if alrm>0 : signal.alarm(alrm)

	return(sock,sockf)

def sendobj(sock,obj):
	"""Sends an object as a (binary) size then a binary pickled object to a socket file object"""
	if obj==None : 
		sock.write(pack("I",0))
		return
	strobj=dumps(obj,-1)
	sock.write(pack("I",len(strobj)))
	sock.write(strobj)

def recvobj(sock):
	"""receives a packed length followed by a binary (pickled) object from a socket file object and returns"""
	l=sock.read(4)
	try : 
		datlen=unpack("I",l)[0]
	except:
		print "Format error in unpacking '%s'"%l
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
		if verbose: print server
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
	dbugfile=file("debug.out","w")
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

	def handle(self):
		"""Process requests from a client. The exchange is:
	send EMAN,version
	recv EMAN,version
	recv command, data len, data (if len not 0)
	send command, data len, data (if len not 0)
	close"""

		if self.verbose>1 : print "Thread %s start"%threading.currentThread().getName()
		
		# periodic housekeeping
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
		if threading.active_count()>DCMAXTHREADS and self.client_address!="127.0.0.1":
			self.sockf.write("WAIT")
			self.sockf.flush()
			if self.verbose>1 : print "Telling client to wait ",self.client_address
			return

		# initial exchange to make sure we're talking to a client
		self.sockf.write("EMAN")
		self.sockf.flush()
		msg = self.sockf.read(4)
		if msg!="EMAN" : raise Exception,"Non EMAN client"

		ver=unpack("I",self.sockf.read(4))[0]
		if ver!=EMAN2PARVER : 
			self.sockf.write("BADV")
			self.sockf.flush()
			raise Exception,"Version mismatch in parallelism (%d!=%d)"%(ver,EMAN2PARVER)
		self.sockf.write("ACK ")
		self.sockf.flush()
		client_id=unpack("I",self.sockf.read(4))[0]

		while (1):
			cmd = self.sockf.read(4)
			if len(cmd)<4 :
				if self.verbose>1 : print "connection closed %s"%(str(self.client_address))
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
					print "Command %s (%s): %s %s    \r"%(str(self.client_address),str(client_id),cmd,str(data)),
					sys.stdout.flush()
				else :
					try: print "Command %s (%s): %s (%d)  "%(str(self.client_address),str(client_id),cmd,len(data))
					except: print "Command %s (%s): %s (-)  "%(str(self.client_address),str(client_id),cmd)
			
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
				EMDCTaskHandler.clients[client_id]=(self.client_address[0],time.time(),cmd)
				
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
				r="HUH"
				try:
					r=self.sockf.read(4)
					if r!="ACK " : 
						EMDCTaskHandler.tasklock.release()
						raise Exception
				except:
					if self.verbose: print "Task sent, no ACK"
					if Task!=None : self.queue.task_rerun(task.taskid)
				if task!=None : 
					EMDCTaskHandler.dbugfile.write("Task %5d sent to %s (%08X)  (%s)  [%s]\n"%(task.taskid,self.client_address[0],client_id,r,local_datetime()))
					EMDCTaskHandler.dbugfile.flush()
				EMDCTaskHandler.tasklock.release()

			
			# Job is completed, results returned
			elif cmd=="DONE" :
				# the first object we get is the task identifier
				tid=data
				
				if self.verbose>1 : print "DONE task (%08X) "%client_id,tid
				
				# keep track of clients
				EMDCTaskHandler.clients[client_id]=(self.client_address[0],time.time(),cmd)
				
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
				EMDCTaskHandler.dbugfile.write("Task %5d complete from %08X with %d data elements  [%s]\n"%(tid,client_id,cnt,local_datetime()))
				EMDCTaskHandler.dbugfile.flush()

			# Request data from the server
			# the request should be of the form ["queue",did,image#]
			# Returns requested object
			elif cmd=="DATA" :
				try:
					fsp=self.queue.didtoname[data[1]]
					obj=EMData(fsp,data[2])
					sendobj(self.sockf,obj)
					if self.verbose>2 : print "Data sent %s(%d)"%(fsp,data[2])
				except:
					sendobj(self.sockf,None)
					if self.verbose : print "Error sending %s(%d)"%(fsp,data[2])
				self.sockf.flush()
				EMDCTaskHandler.clients[client_id]=(self.client_address[0],time.time(),cmd)
				
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
				EMDCTaskHandler.clients[client_id]=(self.client_address[0],time.time(),cmd)

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
				
				EMDCTaskHandler.dbugfile.write("Task %5d queued [%s]\n"%(tid,local_datetime()))

			# This will requeue a task which has already been run (and may be completed)
			# This will generally happen if there was some unexpected error in the returned results
			# and should not ever really be necessary. It is likely indicative of some sort of 
			# problem with one of the compute nodes
			elif cmd=="RQUE":
				self.queue.task_rerun(data)
				if self.verbose : print "Requeuing ",data

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
					EMDCTaskHandler.dbugfile.write("Results for task %5d retrieved [%s]\n"%(data,local_datetime()))
				except:
					if self.verbose: print "No ACK on RSLT. Keeping results for retry."
					EMDCTaskHandler.dbugfile.write("Results for task %5d FAILED [%s]\n"%(data,local_datetime()))
					
				
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
#		if stack!=None : traceback.print_stack(stack)
	print "ALARM - network interrupt"
	if stack!=None : traceback.print_stack(stack)
	return

class EMDCTaskClient(EMTaskClient):
	"""Distributed Computing Task Client. This client will connect to an EMDCTaskServer, request jobs to run
 and run them ..."""
 
	def __init__(self,server,port,verbose=0):
		EMTaskClient.__init__(self)
		self.addr=(server,port)
		self.verbose=verbose
		self.lastupdate=0
		self.task=None
		signal.signal(signal.SIGALRM,DCclient_alarm)	# this is used for network timeouts

	def imalive(self,progress):
		"""Executed code should call this periodically to inidicate that they are still running. This must be called at least once every 5 minutes
		with an integer in the range 0-100. 100 means about to exit, 0 means not started yet. A -1 indicates an error has occurred and the request
		is being aborted. If this function returns False, the client should abort the task in progress."""
		ret=True
		
		signal.alarm(60)
		if time.time()-self.lastupdate>60 : 
			self.lastupdate=time.time()
			if self.task==None : 
				signal.alarm(0)
				return True 
			ret=EMDCsendonecom(self.addr,"PROG",(self.task.taskid,progress),clientid=self.myid)
		
		signal.alarm(0)
		return ret

	def run(self,dieifidle=86400,dieifnoserver=86400):
		"""This is the actual client execution block. dieifidle is an integer number of seconds after
		which to terminate the client if we aren't given something to do. dieifnoserver is the same if we can't find a server to
		communicate with. Both default to 24 hours."""
		count=0
		lastjob=time.time()
		lastserv=time.time()
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
				sockf.write("ACK ")
				sockf.flush()
				signal.alarm(0)
				if task=="EXIT" : break
				if self.verbose and task!=None: print "%s running task id %d"%(socket.gethostname(),task.taskid)
				self.task=task
			except :
				print "No response from server, sleeping 30 sec"
				if time.time()-lastserv>dieifnoserver :
					print "No server for too long. Terminating"
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
					break
				time.sleep(15)
				continue
			sockf.flush()
			
			# Translate and retrieve (if necessary) data for task
			for k,i in task.data.items():
				if self.verbose>1 : print "Data translate ",k,i
#				try:
				if isinstance(i,list) and i[0]=="cache" :
					cname="bdb:cache_%d.%d"%(i[1][0],i[1][1])
					cache=db_open_dict(cname)
					if self.verbose>2 : print "Open cache : ",cname
					
					for j in image_range(*i[2:]):
#						if not cache.has_key(j):
						try: 
							z=cache[j]
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

			# Return results
			if self.verbose : print "Task done %d"%task.taskid
			if self.verbose>3 : print self.__dict__

			retry=True
			retrycount=0
			while retry:
				signal.alarm(60)
				retry=False
				retrycount+=1
				if retrycount>10 :
					print "Failed in 10 attempts to send results, aborting (%d)"%task.taskid
					break
					
				try:
					sock,sockf=openEMDCsock(self.addr,clientid=self.myid,retry=10)
					sockf.write("DONE")
				except:
					print "Server communication failure, trying again in 1 minute (%d)"%task.taskid
					time.sleep(60)
					retry=True
					continue

				signal.alarm(60)
				try:
					sendobj(sockf,task.taskid)
				except:
					print "Immediate ERROR (retrying ",task.taskid,")"
					retry=True
					continue

				for k,v in ret.items():
					signal.alarm(60)
					try:
						sendobj(sockf,k)
						sendobj(sockf,v)
					except :
						print "ERROR (retrying ",task.taskid,") on : ",k, " in ",ret.items()
						if isinstance(v,EMData) : v.write_image("error.hdf",-1)
						time.sleep(3)
						retry=True
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
			
			lastjob=time.time()
			time.sleep(3)
		
		print "Client on %s exiting. Bye !"%socket.gethostname()
		
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
		
