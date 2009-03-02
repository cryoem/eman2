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

# used to make sure servers and clients are running the same version
EMAN2PARVER=1

class EMTaskProxy:
	"""This will communicate with the specified task server on behalf of an application needing to
	have tasks completed"""

class EMTaskServer:
	"""This will serve tasks to clients who execute them"""
	
	def __init__(self,path=None):
		self.queue=EMTaskQueue(path)
		
	
class EMTaskClient:
	"""This class executes tasks on behalf of the server"""

#######################
#  Here we define the classes for publish and subscribe parallelism
class EMDCTaskServer(EMTaskServer,SocketServer.BaseRequestHandler):
	"""Distributed Computing Taskserver. In this system, clients run on hosts with free cycles and request jobs
	from the server, which runs on a host with access to the data to be processed."""
	
	def __init__(self):
		EMTaskServer.__init__(self)
		server = SocketServer.TCPServer((HOST, PORT), self)

		server.serve_forever()

	
	def handle(self):
		"""Process requests from the client"""

		# the beginning of a message is a struct containing
		# 0-3  : EMAN
		# 4-7  : int, EMAN2PARVER
		# 8-11 : 4 char command
		# 12-15: bytes in pickled data following header
		msg = sel.request.recv(4)
		if msg!="EMAN" : raise Exception,"Non EMAN client"
		msg = sel.request.recv(12)
		ver,cmd,datlen=unpack("I4sI",msg)
		if ver!=EMAN2PARVER : raise Exception,"Version mispatch in parallelism"
		
		data = loads(sel.request.recv(datlen))
		
		# Ready for a task
		if cmd=="RDYT" :
			pass
			
		# Notify that a task has been aborted
		elif cmd=="ABOR" : 
			pass
	
		# Progress message indicating that processing continues
		elif cmd=="PROG" :
			pass
		
		
		# self.request is the TCP socket connected to the client
#		self.data = self.request.recv(1024).strip()
#		print "%s wrote:" % self.client_address[0]
#		print self.data
#		# just send back the same data, but upper-cased
#		self.request.send(self.data.upper())

	