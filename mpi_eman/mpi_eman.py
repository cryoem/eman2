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


import eman_mpi_c		# these are the actual C functions. We map these selectively to the python namespace so we can enhance their functionality
from cPickle import dumps,loads

# These functions don't require mapping
mpi_init=eman_mpi_c.mpi_init
mpi_comm_rank=eman_mpi_c.mpi_comm_rank
mpi_comm_size=eman_mpi_c.mpi_comm_size
mpi_barrier=eman_mpi_c.mpi_barrier
mpi_finalize=eman_mpi_c.mpi_finalize

def mpi_send(data,dest,tag):
	"""Synchronously send 'data' to 'dest' with tag 'tag'. data may be any pickleable type"""
	if isinstance(data,str):
		
	else :
		d2x=dumps(data,-1)
		if len(d2x>64) 

def mpi_recv

def mpi_bcast
