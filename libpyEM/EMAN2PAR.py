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
import socket

# used to make sure servers and clients are running the same version
EMAN2PARVER=1

class EMTaskProxy:
	"""This will communicate with the specified task server on behalf of an application needing to
	have tasks completed"""

class EMTaskHandler:
	"""This will respond to client requests"""
	
	def __init__(self,path=None):
		self.queue=EMTaskQueue(path)
		
	
class EMTaskClient:
	"""This class executes tasks on behalf of the server. This parent class implements the actual task functionality.
Communications are handled by subclasses."""
	def process_task(self,task):
		"""This method implements the actual image processing by calling appropriate module functions"""
		
#######################
#  Here we define the classes for publish and subscribe parallelism

def sendobj(sock,obj):
	"""Sends an object as a (binary) size then a binary pickled object"""
	strobj=dumps(obj,-1)
	sock.send(pack("I",len(strobj)))
	sock.send(strobj)

def recvobj(sock):
	"""receives a packed length followed by a binary (pickled) object and returns"""
	datlen=unpack("I",sock.recv(4))
	if datlen<=0 :return None
	return loads(sock.recv(datlen))

def EMDCsendonecom(host,port,cmd,data):
	"""Connects to an EMAN EMDCServer sends one command then disconnects"""
	# the beginning of a message is a struct containing
	# 0-3  : EMAN
	# 4-7  : int, EMAN2PARVER
	# 8-11 : 4 char command
	# 12-15: count of bytes in pickled data following header
	sock=socket.socket()
	sock.connect((host,port))
	if sock.recv(4)!="EMAN" : raise Exception,"Not an EMAN server"
	sock.send("EMAN")
	sock.send(pack("I4s",EMAN2PARVER,cmd))
	sendobj(sock,data)
	
	ret=recvobj(sock)
	sock.close()
	return ret

def runEMDCServer(port,verbose):
	"""This will create a ThreadingTCPServer instance and execute it"""
	try: EMDCTaskHandler.verbose=int(verbose)
	except: EMDCTaskHandler.verbose=0

	if port!=None and port>0 : 
		server = SocketServer.TCPServer(("", port), EMDCTaskHandler)	# "" is the hostname and will bind to any IPV4 interface/address

	# EMAN2 will use ports in the range 9900-9999
	for port in range(9900,10000):
		try: server = SocketServer.TCPServer(("", port), self)
		except: continue
		break

	server.serve_forever()


class EMDCTaskHandler(EMTaskHandler,SocketServer.BaseRequestHandler):
	"""Distributed Computing Taskserver. In this system, clients run on hosts with free cycles and request jobs
	from the server, which runs on a host with access to the data to be processed."""
	verbose=0
	
	def __init__(self, request, client_address, server):
		# if a port is specified, just try to open it directly. No error detection

		EMTaskHandler.__init__(self)
		SocketServer.BaseRequestHandler(request,client_address,server)
		self.verbose=EMDCTaskHandler.verbose


	
	def handle(self):
		"""Process requests from a client. The exchange is:
	send EMAN,version
	recv EMAN,version
	recv command, data len, data (if len not 0)
	send command, data len, data (if len not 0)
	close"""

		while (1):
			# the beginning of a message is a struct containing
			# 0-3  : EMAN
			# 4-7  : int, EMAN2PARVER
			# 8-11 : 4 char command
			# 12-15: count of bytes in pickled data following header
			self.request.send("EMAN")
			msg = sel.request.recv(4)
			if msg!="EMAN" : raise Exception,"Non EMAN client"
			msg = sel.request.recv(12)
			ver,cmd,datlen=unpack("I4sI",msg)
			if ver!=EMAN2PARVER : raise Exception,"Version mismatch in parallelism"
			
			data = loads(sel.request.recv(datlen))
			
			if self.verbose : print "Command %s: %s (%d)"%(str(self.client_address),cmd,datlen)

			######################  These are issued by clients
			# Ready for a task
			if cmd=="RDYT" :
				# Get the first task and send it (pickled)
				task=self.queue.get_task()
				task.starttime=time.time()
				sendobj(self.request,task)
				
				# check for an ACK, if not, requeue
				try:
					if sel.request.recv(4)!="ACK " : raise Exception
				except: 
					task.starttime=None		# requeue
			
			
			# Job is completed, results returned
			elif cmd=="DONE" :
				pass
			
			# Request data from the server
			elif cmd=="DATA" :
				pass
			
			# Notify that a task has been aborted
			elif cmd=="ABOR" : 
				pass
		
			# Progress message indicating that processing continues
			elif cmd=="PROG" :
				pass

			###################### These are utility commands
			# Returns whatever is sent as data
			elif cmd=="TEST":
				ret=dumps(data)
				if data : sel.request.send(dumps(data))
				
			
			# Cause this server to exit (cleanly)
			elif cmd=="QUIT" :
				pass
			
			################### These are issued by customers
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
	"""Distributed Computing Task Client. This client will connect to an EMDCTaskServer and request jobs to run
 and run them ..."""

		