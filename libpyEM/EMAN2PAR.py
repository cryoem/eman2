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

from EMAN2 import test_image
from EMAN2db import EMTask,EMTaskQueue,db_open_dict
import SocketServer
from cPickle import dumps,loads
from struct import pack,unpack
import time
import socket
import sys

# used to make sure servers and clients are running the same version
EMAN2PARVER=1

class EMTaskCustomer:
	"""This will communicate with the specified task server on behalf of an application needing to
	have tasks completed"""

class EMTaskHandler:
	"""This is the actual server object which talks to clients and customers. It coordinates task execution
 acts as a data clearinghouse"""
	
	def __init__(self,path=None):
		self.queue=EMTaskQueue(path)
		
class EMTaskClient:
	"""This class executes tasks on behalf of the server. This parent class implements the actual task functionality.
Communications are handled by subclasses."""
	def __init__(self):
		pass
	
	def process_task(self,task):
		"""This method implements the actual image processing by calling appropriate module functions"""

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
		print "Server started on port %d"%port
	# EMAN2 will use ports in the range 9900-9999
	else :
		for port in range(9990,10000):
			try: 
				server = SimpleXMLRPCServer(("", port),SimpleXMLRPCRequestHandler,False,allow_none=True)
				if verbose: print "Server started on port %d"%port
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

def openEMDCsock(addr):
	sock=socket.socket()
	sock.connect(addr)
	sockf=sock.makefile()
			
	# Introduce ourselves and ask for a task to execute
	if sockf.read(4)!="EMAN" : raise Exception,"Not an EMAN server"
	sockf.write("EMAN")
	sockf.write(pack("I4",EMAN2PARVER))
	return(sock,sockf)

def sendobj(sock,obj):
	"""Sends an object as a (binary) size then a binary pickled object to a socket file object"""
	if obj==None : 
		sock.write("\0\0\0\0")
		return
	strobj=dumps(obj,-1)
	sock.write(pack("I",len(strobj)))
	sock.write(strobj)

def recvobj(sock):
	"""receives a packed length followed by a binary (pickled) object from a socket file object and returns"""
	datlen=unpack("I",sock.read(4))[0]
	if datlen<=0 :return None
	return loads(sock.read(datlen))

def EMDCsendonecom(host,port,cmd,data):
	"""Connects to an EMAN EMDCServer sends one command then disconnects"""
	# the beginning of a message is a struct containing
	# 0-3  : EMAN
	# 4-7  : int, EMAN2PARVER
	# 8-11 : 4 char command
	# 12-15: count of bytes in pickled data following header
	sock=socket.socket()
	sock.connect((host,port))
	sockf=sock.makefile()
	if sockf.read(4)!="EMAN" : raise Exception,"Not an EMAN server"
	sockf.write("EMAN")
	sockf.write(pack("I4s",EMAN2PARVER,cmd))
	sendobj(sockf,data)
	sockf.flush()
	
	ret=recvobj(sockf)
	sockf.close()
	return ret

class ThreadingTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer): pass

def runEMDCServer(port,verbose):
	"""This will create a ThreadingTCPServer instance and execute it"""
	try: EMDCTaskHandler.verbose=int(verbose)
	except: EMDCTaskHandler.verbose=0

	if port!=None and port>0 : 
		server = SocketServer.ThreadingTCPServer(("", port), EMDCTaskHandler)	# "" is the hostname and will bind to any IPV4 interface/address
		print server
	# EMAN2 will use ports in the range 9900-9999
	else :
		for port in range(9990,10000):
			try: 
				server = SocketServer.ThreadingTCPServer(("", port), EMDCTaskHandler)
				print "Server started on port %d"%port
			except:
				print "Port %d unavailable"%port
				continue
			break

	server.serve_forever()

class EMDCTaskHandler(EMTaskHandler,SocketServer.BaseRequestHandler):
	"""Distributed Computing Taskserver. In this system, clients run on hosts with free cycles and request jobs
	from the server, which runs on a host with access to the data to be processed."""
	verbose=0
	
	def __init__(self, request, client_address, server):
		# if a port is specified, just try to open it directly. No error detection

		EMTaskHandler.__init__(self)
		self.verbose=EMDCTaskHandler.verbose
		if self.verbose>1 : print len(self.queue)
		self.sockf=request.makefile()		# this turns our socket into a buffered file-like object
		SocketServer.BaseRequestHandler.__init__(self,request,client_address,server)

	def handle(self):
		"""Process requests from a client. The exchange is:
	send EMAN,version
	recv EMAN,version
	recv command, data len, data (if len not 0)
	send command, data len, data (if len not 0)
	close"""

		print self.sockf
		if self.verbose>1 : print "connection from %s"%(str(self.client_address))

		# the beginning of a message is a struct containing
		# 0-3  : EMAN
		# 4-7  : int, EMAN2PARVER
		# 8-11 : 4 char command
		# 12-15: count of bytes in pickled data following header

		# initial exchange to make sure we're talking to a client
		self.sockf.write("EMAN")
		self.sockf.flush()
		msg = self.sockf.read(4)
		if msg!="EMAN" : raise Exception,"Non EMAN client"

		ver=unpack("I",self.sockf.read(4))[0]
		if ver!=EMAN2PARVER : raise Exception,"Version mismatch in parallelism (%d!=%d)"%(ver,EMAN2PARVER)

		while (1):
			cmd = self.sockf.read(4)
			if len(cmd)<4 :
				if self.verbose>1 : print "connection closed %s"%(str(self.client_address))
				return
				
			if cmd=="ACK " :
				print "Warning : Out of band ACK"
				continue		# some sort of protocol error, just ignore

			data = recvobj(self.sockf)
				
			if self.verbose : 
				try: print "Command %s: %s (%d)"%(str(self.client_address),cmd,len(data))
				except: print "Command %s: %s (-)"%(str(self.client_address),cmd)
			

			######################  These are issued by clients
			# Ready for a task
			if cmd=="RDYT" :
				# Get the first task and send it (pickled)
				task=self.queue.get_task()
				if task==None :
					sendobj(self.sockf,None)			# no tasks available
				else:
					task.starttime=time.time()			# send the task
					sendobj(self.sockf,task)
					if self.verbose>1: print "Send task: ",task.taskid
				self.sockf.flush()
				
				# check for an ACK, if not, requeue
				try:
					if self.sockf.read(4)!="ACK " : raise Exception
				except:
					if self.verbose: print "Task sent, no ACK"
					if task: task.starttime=None		# requeue
			
			
			# Job is completed, results returned
			elif cmd=="DONE" :
				# the first object we get is the task identifier
				tid=data
				
				if self.verbose>1 : print "DONE task ",tid
				
				# then we get a sequence of key,value objects, ending with a final None key
				result=db_open_dict("bdb:%s#result_%d"%(self.queue.path,tid))
				
				cnt=0
				while (1) :
					key=recvobj(self.sockf)
					if key==None : break
					val=recvobj(self.sockf)
					result[key]=val				# store each key in the database
					cnt+=1

				result.close()
				if self.verbose>2 : print "Task %d: %d data elements"%(tid,cnt)

				self.queue.task_done(tid)		# don't mark the task as done until we've stored the results

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
				
			# Notify that a task has been aborted
			# request should be taskid
			# no return
			elif cmd=="ABOR" : 
				self.queue.task_rerun(data)
				if self.verbose>1 : print "Task execution abort : ",data

			# Progress message indicating that processing continues
			# request should be (taskid,percent complete)
			# returns "OK  " if processing should continue
			# retunrs "ABOR" if processing should be aborted
			elif cmd=="PROG" :
				ret=self.queue.task_progress(data[0],data[1])
				if self.verbose>2 : print "Task progress report : ",data
				if ret : self.sockf.send("OK  ")
				else : self.sockf.send("ABOR")

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
					self.sendobj.flush()
					continue
				sendobj(self.sockf,ret)

				result=db_open_dict("bdb:%s#result_%d"%(self.queue.path,tid))
				for k in result.keys():
					if k=="maxrec": continue
					sendobj(self.sockf,k)
					sendobj(self.sockf,result[k])
					
				sendobj(self.sockf,None)
				self.sockf.flush()
				
		# self.request is the TCP socket connected to the client
#		self.data = self.request.recv(1024).strip()
#		print "%s wrote:" % self.client_address[0]
#		print self.data
#		# just send back the same data, but upper-cased
#		self.request.send(self.data.upper())

class EMDCTaskClient(EMTaskClient):
	"""Distributed Computing Task Client. This client will connect to an EMDCTaskServer, request jobs to run
 and run them ..."""
 
	def __init__(self,server,port,verbose=0):
		EMTaskClient.__init__(self)
		self.addr=(server,port)
		self.verbose=verbose

	def run(self):
		while (1):
			# connect to the server
			if self.verbose>1 : print "Connect to (%s,%d)"%self.addr
			sock,sockf=openEMDCsock(self.addr)
			sockf.write("RDYT")
			sendobj(sockf,None)
			sockf.flush()
			
			# Get a task from the server
			task=recvobj(sockf)
			sockf.write("ACK ")
			if task==None:
				if self.verbose : print "No tasks to run :^("
				sockf.close()
				sock.close()
				time.sleep(10)
				continue
			sockf.flush()
			
			# Translate and retrieve (if necessary) data for task
			for k,i in task.data.items():
				if self.verbose>1 : print "Data translate ",k,i
				try:
					if i[0]=="cache" :
						cname="bdb:cache_%d.%d"%(i[1][0],i[1][1])
						cache=db_open_dict(cname)
						
						# form min,max
						if len(i) == 4:
							for j in range(i[2],i[3]):
								if not cache.has_key(j):
									cache[j]=self.get_data(sockf,i[1],j)
						# form image number
						elif isinstance(i[2],int):
							if not cache.has_key(i[2]) :
								cache[i[2]]=self.get_data(sockf,i[1],i[2])
						# form list of numbers
						else:
							for j in i[2]:
								if cache.has_key(j) : continue
								cache[j]=self.get_data(sockf,i[1],j)
								
						i[1]=cname
				except: pass
			sockf.close()
			sock.close()
			
			# Execute translated task
			
			# Return results
			if self.verbose>1 : print "Task done"
			sock,sockf=openEMDCsock(self.addr)
			sockf.write("DONE")
			sendobj(sockf,task.taskid)
			sendobj(sockf,"result")
			sendobj(sockf,"TEST")
			sendobj(sockf,"img")
			sendobj(sockf,test_image())
			sendobj(sockf,None)
			sockf.flush()
			sockf.close()
			sock.close()
			
			time.sleep(5)
			
	def get_data(self,sockf,did,imnum):
		"""request data from the server on an open connection"""
		if self.verbose>2 : print "Retrieve %s : %d"%(str(did),imnum)
		sockf.write("DATA")
		sendobj(sockf,["cache",did,imnum])
		sockf.flush()
		return recvobj(sockf)
		
	#		self.sockf=
		
