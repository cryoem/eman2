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

import sys
from cPickle import dumps,loads
from zlib import compress,decompress
from struct import pack,unpack

def mpi_dout(data):
	"""Convert data for transmission efficiently"""

	if isinstance(data,str):
		if len(data)>256 : return "C"+compress(data,1)
		else : return "S"+data
	else :
		d2x=dumps(data,-1)
		if len(d2x)>256 : return "Z"+compress(d2x,1)
		else : return "O"+d2x

	
def mpi_din(msg):
	"""Unpack a data payload prepared by mpi_cout"""
	
	if msg[0]=="C" : return decompress((msg[1:]).tostring())
	elif msg[0]=="S" : return (msg[1:]).tostring()
	elif msg[0]=="Z" : return loads(decompress((msg[1:]).tostring()))
	elif msg[0]=="O" : return loads((msg[1:]).tostring())
	else :
		print "ERROR: Invalid MPI message. Please contact developers. (%s)"%str(msg[:20])
		sys.exit(1)


def mpi_eman2_send(data,dest,tag):
	"""Synchronously send 'data' to 'dest' with tag 'tag'. data may be any pickleable type.
	Compression and pickling will be performed on the fly when deemed useful. Note that data duplication
	in memory may occur."""
	from mpi import mpi_send, MPI_CHAR, MPI_COMM_WORLD

	data = mpi_dout(data)
	mpi_send(data, len(data), MPI_CHAR, dest, tag, MPI_COMM_WORLD)

def mpi_eman2_recv(src,tag):
	"""Synchronously receive a message from 'src' with 'tag'."""
	from mpi import mpi_probe, mpi_get_count, mpi_recv, MPI_CHAR, MPI_COMM_WORLD
	
	mpi_probe(src, tag, MPI_COMM_WORLD)
	count = mpi_get_count(MPI_CHAR)
	msg = mpi_recv(count, MPI_CHAR, src, tag, MPI_COMM_WORLD)
	return (mpi_din(msg),src,tag)
	
def mpi_bcast_send(data):
	"""Unlike the C routine, in this python module, mpi_bcast is split into a send and a receive method. Send must be 
	called on exactly one core, and receive called on all of the others. This routine also coordinates transfer of
	variable length objects."""
	from mpi import mpi_comm_rank, mpi_bcast, MPI_CHAR, MPI_COMM_WORLD

	data=mpi_dout(data)
	l=pack("I",len(data))
	rank = mpi_comm_rank(MPI_COMM_WORLD)	
	mpi_bcast(l, len(l), MPI_CHAR, rank, MPI_COMM_WORLD)
	mpi_bcast(data, len(data), MPI_CHAR, rank, MPI_COMM_WORLD)
	
def mpi_bcast_recv(src):
	"""Unlike the C routine, in this python module, mpi_bcast is split into a send and a receive method. Send must be 
	called on exactly one core, and receive called on all of the others. This routine also coordinates transfer of
	variable length objects. src is the rank of the sender"""
	from mpi import mpi_bcast, MPI_CHAR, MPI_COMM_WORLD
	
	l = mpi_bcast(None, 4, MPI_CHAR, src, MPI_COMM_WORLD)
	l=unpack("I",l)[0]
	data = mpi_bcast(None, l, MPI_CHAR, src, MPI_COMM_WORLD)
	return mpi_din(data)
