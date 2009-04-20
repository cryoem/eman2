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

from EMAN2db import EMTask,EMTaskQueue
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
#  Here we define the classes for publish and subscribe parallelism

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

def runEMDCServer(port,verbose):
	"""This will create a ThreadingTCPServer instance and execute it"""
	try: EMDCTaskHandler.verbose=int(verbose)
	except: EMDCTaskHandler.verbose=0

	if port!=None and port>0 : 
		server = SocketServer.TCPServer(("", port), EMDCTaskHandler)	# "" is the hostname and will bind to any IPV4 interface/address
		print server
	# EMAN2 will use ports in the range 9900-9999
	else :
		for port in range(9990,10000):
			try: 
				server = SocketServer.TCPServer(("", port), EMDCTaskHandler)
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

		# initial exchange to make sure we're talking to a client
		self.sockf.write("EMAN")
		self.sockf.flush()
		msg = self.sockf.read(4)
		if msg!="EMAN" : raise Exception,"Non EMAN client"
		
		while (1):
			# the beginning of a message is a struct containing
			# 0-3  : EMAN
			# 4-7  : int, EMAN2PARVER
			# 8-11 : 4 char command
			# 12-15: count of bytes in pickled data following header
			msg = self.sockf.read(8)
			if len(msg)<8 :
				if self.verbose>1 : print "connection closed %s"%(str(self.client_address))
				return
				
			ver,cmd=unpack("I4s",msg)
			if ver!=EMAN2PARVER : raise Exception,"Version mismatch in parallelism"
			
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
					sendobj(self.sockf,None)
				else:
					task.starttime=time.time()
					sendobj(self.sockf,task)
				self.sockf.flush()
				
				# check for an ACK, if not, requeue
				try:
					if sel.sockf.read(4)!="ACK " : raise Exception
				except: 
					if task: task.starttime=None		# requeue
			
			
			# Job is completed, results returned
			elif cmd=="DONE" :
				pass
			
			# Request data from the server
			elif cmd=="DATA" :
				datac=
				sendobj(self.sockf,datac)
			
			# Notify that a task has been aborted
			elif cmd=="ABOR" : 
				pass
		
			# Progress message indicating that processing continues
			elif cmd=="PROG" :
				pass

			###################### These are utility commands
			# Returns whatever is sent as data
			elif cmd=="TEST":
				sendobj(self.sockf,("TEST",data))
				self.sockf.flush()
				
			
			# Cause this server to exit (cleanly)
			elif cmd=="QUIT" :
				sendobj(self.sockf,None)
				self.sockf.flush()
				self.server.server_close()
				if self.verbose : print "Server exited cleanly"
#				sys.exit(0)
			
			################### These are issued by customers
			# A new task to enqueue
			elif cmd=="TASK":
				try: 
					self.queue.add_task(data)
					if self.verbose>1 : print "TASK %s.%s"%(data.module,data.command)
				

			# status request
			elif cmd=="STAT":
				pass
		
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
			sock=socket.socket()
			sock.connect(self.addr)
			sockf=sock.makefile()
			
			# Introduce ourselves and ask for a task to execute
			if sockf.read(4)!="EMAN" : raise Exception,"Not an EMAN server"
			sockf.write("EMAN")
			sockf.write(pack("I4s",EMAN2PARVER,"RDYT"))
			sendobj(sockf,None)
			sockf.flush()
			
			# Get a task from the server
			task=recvobj(sockf)
			if task==None:
				if self.verbose : print "No tasks to run :^("
				sockf.close()
				sock.close()
				time.sleep(10)
				continue
			
			# Process the task
			for k,i in task.data.items():
				try:
					if i[0]=="cache" :
						sockf.write("DATA")
						task.data[k]=sockf.sendobj(sockf,i)
				except: pass
			
	#		self.sockf=
		
