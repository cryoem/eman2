#!/usr/bin/env python

# recoverctf.py
# Author: Steven Ludtke, 09/16/2011 (sludtke@bcm.edu)
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

import sys
from EMAN2 import *
from EMAN2db import *

# This program will extract CTF parameters from the headers of phase-flipped particles to restore a corrupted
# bdb:e2ctf.parms database.

# To use it:
# e2bdb.py -c
# rm -rf EMAN2DB/e2ctf.parms
# rerun automatic fitting for your project
# run this script. It will replace the automatic fit parameters with the parameters from the phase flipped particles

ptcls=[i for i in db_list_dicts("bdb:particles") if i[-8:]=="ctf_flip"]
db=db_open_dict("bdb:e2ctf.parms")

for i in ptcls:
	print i[:-9]," recovered"
	img=EMData("bdb:particles#%s"%i,0)
	cur=db[i[:-9]]
	cur[0]=img["ctf"].to_string()
	db[i[:-9]]=cur


