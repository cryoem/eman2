#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/25/2009 (sludtke@bcm.edu)
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

# tasktest.py Steven Ludtke
# This program implements, via various options, the parallelism system for EMAN2

from EMAN2 import *
from EMAN2PAR import *
from optparse import OptionParser
from math import *
import time
import os
import sys
import socket

sock=socket.socket()
sock.connect(("localhost",9990))
sockf=sock.makefile()

task=EMTask({1:"Testing"},"Test",{})

if sockf.read(4)!="EMAN" : raise Exception,"Not an EMAN server"
sockf.write("EMAN")
sockf.write(pack("I4s",EMAN2PARVER,"TASK"))
sendobj(sockf,task)
sockf.flush()
tid=recvobj(sockf)

print "Taskid = ",tid

for i in range(10):
	time.sleep(15.0)
	sockf.write("STAT")
	sendobj(sockf,[tid])
	sockf.flush()
	result=recvobj(sockf)[0]
	print result
	if result==100 :
		print "Task completed"
		sockf.write("RSLT")
		sendobj(sockf,tid)
		sockf.flush()
		
		task=recvobj(sockf)
		print task.__dict__
		
		k=0
		while 1:
			k=recvobj(sockf)
			if k==None: break
			v=recvobj(sockf)
			print k,": ",v
		
		break

