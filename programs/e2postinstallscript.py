#!/usr/bin/env python

# Muyuan Chen 2016-01
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
import os
import time

### replacing the timestamp in e2version.py
now_time = time.ctime()
emandir= os.getenv("EMAN2DIR")
e2ver=open("{}/bin/e2version.py".format(emandir),'r')
lines=e2ver.readlines()
e2ver.close()

e2ver=open("{}/bin/e2version.py".format(emandir),'w')
for l in lines:
	e2ver.write(l.replace("BUILD_DATE",now_time))
e2ver.close()

### remove this file since we do not want user to run this script after installation
try: os.remove("{}/bin/e2postinstallscript.py".format(emandir))
except: pass
