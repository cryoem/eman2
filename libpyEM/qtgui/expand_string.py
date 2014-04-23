#!/usr/bin/env python

#
# Author: Vernon Williams (vewillia@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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
#-------------------------------------------------------------------------
#
# Python Function  expand_string     Level     5     Date 04/23/14
#
# Purpose:  expand_string is a Python string function to expand any
#           environment variables (beginning with '$' and ended
#           by a blank or special character other than '_') in a
#           string by substituting in the values of the environment
#           variables found.
#
# Sample Call Statement:
#
# expanded_string = expand_string (string)
#
# where:
#
# string is an input character string giving a string in which
# environment variables should be expanded into their values.
#
# The string returned is a static local string of up to 2048
# characters, so the next call to this routine will overwrite
# the string.
#
# External Variables Used:  None.
#
# Called Routines:  None.
#
# Written  by Vernon Williams, February  19, 1996.
# Modified by Dean Knott,      December   2, 1999,
#    to change include file extensions from .INC to .fi, for PC
#    conversion.
# Modified by Sarah Steuber,   February  27, 2001,
#    to rename called routine GETLOG to GTLENV, to avoid a name
#    conflict (apparently, Sun has a system routine called GETLOG).
# Modified by Vernon Williams, May       16, 2006,
#    to fix a small bug using backslashes (to replace the one
#    backslash in string constant envend with two if the compiler
#    translates \\ to \, then it will be the same as before, which
#    works if it does not do this, there will be two backslashes,
#    in the string, which will not hurt so all cases are handled).
# Modified by Vernon Williams, September  8, 2009,
#    to convert from FORTRAN to K, to rename the routine from EXPSTR
#    to w_expand_string.
# Modified by Vernon Williams, April     23, 2014,
#    to make Python version expand_string.py from K version
#    w_expand_string.c.
#
#--------------------------------------------------------------------------*/

import os

def expand_string (string) :
	# Constants:

	blank  = " "

	# Character to signal start of an environment variable:

	envbeg = "$"

	# Characters to signal end of an environment variable:

	envend = ",.<>/?:\'[]{}\\|`~!@#$%^&*[-+= ']\""

	gotenv = False
	inplen = len (string)
	outlen = 0
	outstr = ""

	for ii in xrange (0, inplen + 1) :
		if ii < inplen :
			ch = string [ii]
		else :
			ch = envbeg

		if gotenv :
			if ch in envend :
				gotenv = (ch == envbeg)

				if envnam == "" :
					envval = ""

					error  = False
				else :
					if envnam in os.environ :
						envval = os.getenv (envnam)
					else :
						envval = ""

					error  = False

				if not error :
					vallen = len (envval)

					for iv in xrange (0, vallen) :
						outstr += envval [iv]
						outlen += 1
				if gotenv :
					envnam = ""
					envlen = 0
				else :
					outstr += ch
					outlen += 1
			else :
				envnam += ch
				envlen += 1
		else :
			if ch == envbeg :
				gotenv = True
				envlen = 0
				envnam = ""
			else :
				outstr += ch
				outlen += 1

	return outstr
