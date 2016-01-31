#!/usr/bin/env python

#
# Author: Grant Tang
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

import os
import sys
import platform
from subprocess import *


EMANVERSION="EMAN 2.12"
DATESTAMP="BUILD_DATE"

def main():
	print(EMANVERSION + ' (GITHUB: ' + DATESTAMP +')')

	if sys.platform=='linux2':
		print('Your EMAN2 is running on: {} {}'.format(platform.platform(), os.uname()[2], os.uname()[-1]))

	elif sys.platform=='darwin':
		print('Your EMAN2 is running on: Mac OS {} {}'.format(platform.mac_ver()[0], platform.mac_ver()[2]))

	elif sys.platform=='win32':
		ver = sys.getwindowsversion()
		ver_format = ver[3], ver[0], ver[1]
		win_version = {
					(1, 4, 0): '95',
					(1, 4, 10): '98',
					(1, 4, 90): 'ME',
					(2, 4, 0): 'NT',
					(2, 5, 0): '2000',
					(2, 5, 1): 'XP',
					(2, 5, 2): '2003',
					(2, 6, 0): '2008',
					(2, 6, 1): '7'
				}

		if win_version.has_key(ver_format):
			winsysver = 'Windows' + ' ' + win_version[ver_format]
		else:
			winsysver = 'Windows'

		if 'PROGRAMFILES(X86)' in os.environ:
			winsystype = '64bit'
		else:
			winsystype = '32bit'

		print('Your EMAN2 is running on: {} {}'.format(winsysver, winsystype))

	print('Your Python version is: {}'.format(os.sys.version.split()[0]))

if __name__== "__main__":
	main()
