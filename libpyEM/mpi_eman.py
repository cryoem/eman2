#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
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

from future import standard_library
standard_library.install_aliases()
import sys
from pickle import dumps,loads
from zlib import compress,decompress
from struct import pack,unpack
import numpy as np

# no longer used
def mpi_dout(data):
	"""Convert data for transmission efficiently"""

	if isinstance(data,str):
		if len(data)>256 : return "C"+compress(data,1)
		else : return "S"+data
	else :
		d2x=dumps(data,-1)
		if len(d2x)>256 : return "Z"+compress(d2x,1)
		else : return "O"+d2x

	
# no longer used
def mpi_din(msg):
	"""Unpack a data payload prepared by mpi_cout"""
	
	if msg[0]=="C" : return decompress((msg[1:]).tostring())
	elif msg[0]=="S" : return (msg[1:]).tostring()
	elif msg[0]=="Z" : return loads(decompress((msg[1:]).tostring()))
	elif msg[0]=="O" : return loads((msg[1:]).tostring())
	else :
		print("ERROR: Invalid MPI message. Please contact developers. (%s)"%str(msg[:20]))
		sys.exit(1)

def bytestoint(data):
	d=[int.from_bytes(data[i:i+4], 'little') for i in range(0,len(data),4)]
	if len(d)==1: d=d[0]
	return d

def mpi_eman2_send(com,data,dest):
	"""Synchronously send 'data' to 'dest' with 4 char string command flag com. data may be any pickleable type.
	Compression and pickling will be performed on the fly when deemed useful. Note that data duplication
	in memory may occur."""
	from mpi import mpi_send, MPI_CHAR, MPI_COMM_WORLD,mpi_comm_rank, MPI_INT
	rank = mpi_comm_rank(MPI_COMM_WORLD)	
	comint=bytestoint(com.encode('utf-8'))
	#mpi_send(com,4,MPI_CHAR,dest,1,MPI_COMM_WORLD)
	#print('eman_send', com, comint, data,dest, type(data))
	if isinstance(data,str) :
		tag=3
		mpi_send([len(data), len(data), rank, tag, comint],5,MPI_INT,dest,1,MPI_COMM_WORLD)
		mpi_send(data, len(data), MPI_CHAR, dest, tag, MPI_COMM_WORLD)
		
	else:
		tag=2
		data = dumps(data,-1)
		l0=len(data)
		data=bytestoint(data)
		mpi_send([l0, len(data), rank, tag, comint],5,MPI_INT,dest,1,MPI_COMM_WORLD)
		mpi_send(data, len(data), MPI_INT, dest, tag, MPI_COMM_WORLD)
	

def mpi_eman2_recv(src):
	"""Synchronously receive a message from 'src' with 'tag'."""
	from mpi import mpi_probe, mpi_get_count, mpi_recv, MPI_CHAR, MPI_COMM_WORLD, MPI_INT
	#com=mpi_recv(4,MPI_CHAR, src,1,MPI_COMM_WORLD)
	#com=b''.join(com).decode('utf-8')
	
	lmsg=mpi_recv(5,MPI_INT, src,1,MPI_COMM_WORLD)		# first get the message length
	#print('recv', type(lmsg), len(lmsg),lmsg)
	l0, l,srank, tag, comint=lmsg
	
	#print('eman_recv',  lmsg, src, comint)
	com=np.uint32(comint).item().to_bytes(4,'little')
	com=com.decode('utf-8')
	
	if tag==3:
		msg=mpi_recv(l,MPI_CHAR, srank,tag,MPI_COMM_WORLD)
		msg=b''.join(msg).decode('utf-8')
		return (com,msg,srank)
	elif tag==2:
		msg=mpi_recv(l,MPI_INT, srank,tag,MPI_COMM_WORLD)	# then the message
		msg=msg.astype(np.uint32).tolist()
		msg=b''.join([i.to_bytes(4,'little') for i in msg])
		msg=msg[:l0]
		return (com,loads(msg),srank)
	
	
def mpi_bcast_send(data):
	"""Unlike the C routine, in this python module, mpi_bcast is split into a send and a receive method. Send must be 
	called on exactly one core, and receive called on all of the others. This routine also coordinates transfer of
	variable length objects."""
	from mpi import mpi_comm_rank, mpi_bcast, MPI_CHAR, MPI_COMM_WORLD, MPI_INT
	rank = mpi_comm_rank(MPI_COMM_WORLD)
	#print('mpi bcast send ', data)
	if isinstance(data,str) :
		mpi_bcast([-1, len(data)],2,MPI_INT, rank, MPI_COMM_WORLD)
		mpi_bcast(data,len(data),MPI_CHAR, rank, MPI_COMM_WORLD)
	else:
		data=dumps(data,-1)
		l0=len(data)
		data=bytestoint(data)
		#l=pack("I",len(data))
		
		mpi_bcast([l0, len(data)],2,MPI_INT, rank, MPI_COMM_WORLD)
		#mpi_bcast(l, len(l), MPI_CHAR, rank, MPI_COMM_WORLD)
		mpi_bcast(data, len(data), MPI_INT, rank, MPI_COMM_WORLD)
	
def mpi_bcast_recv(src):
	"""Unlike the C routine, in this python module, mpi_bcast is split into a send and a receive method. Send must be 
	called on exactly one core, and receive called on all of the others. This routine also coordinates transfer of
	variable length objects. src is the rank of the sender"""
	from mpi import mpi_bcast, MPI_CHAR, MPI_COMM_WORLD, MPI_INT
	
	l = mpi_bcast(None, 2, MPI_INT, src, MPI_COMM_WORLD)
	l0, l=l
	
	if l0<0:
		### string data
		data = mpi_bcast(None, l, MPI_CHAR, src, MPI_COMM_WORLD)
		data=b''.join(data).decode('utf-8')
		#print('mpi bcast recv ', data)
		return data
	else:
		data = mpi_bcast(None, l, MPI_INT, src, MPI_COMM_WORLD)
		data=data.astype(np.uint32).tolist()
		data=b''.join([i.to_bytes(4,'little') for i in data])
		data=data[:l0]
		return loads(data)
