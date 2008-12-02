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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# e2history.py   08/08/2005  Steven Ludtke
# This program will dump the local logfile of all EMAN2 programs run in the current
# directory. The file is stored as a python 'shelve' file

import shelve
import sys,os,time

# also defined in EMAN2, but we don't want to have to import it
def local_datetime(secs):
	from time import localtime
	t=localtime(secs)
	return "%04d/%02d/%02d %02d:%02d:%02d"%t[:6]

def time_diff(secs):
	if secs<3600 : return "%d:%02d"%(secs/60,secs%60)
	return "%d:%02d:%02d"%(secs/3600,(secs%3600)/60,secs%60)

if len(sys.argv)>1 and sys.argv[1]=="--help" :
	print "Usage:\ne2history [--all]\n"
	sys.exit(1)

try:
	import EMAN2db
	db=EMAN2db.EMAN2DB.open_db()
	db.open_dict("history")
except:
	db=None

if db:
	try:
		n=int(db.history["count"])
	except:
		print "no logfile"
		sys.exit(0)
	
	if len(sys.argv)>1 and (sys.argv[1]=="--all" or sys.argv[1]=="-A") :
		ah={}
		for i in range(n):
			try: h=db.history[i+1]
			except: print "Entry ",i," missing"
			ah.setdefault(h["path"],[]).append(h)
		for k in ah.keys():
			print "---------- ",k
			for i in ah[k]:
				if i.has_key("end") : print local_datetime(i["start"]),"\t   ",time_diff(i["end"]-i["start"]),"\t"," ".join(i["args"])
				else: print local_datetime(i["start"]),"\tincomplete\t"," ".join(i["args"])
	else:
		for i in range(n):
			try: h=db.history[i+1]
			except: print "Entry ",i," missing"
			if h != None and h.has_key("path") and h["path"]==os.getcwd():
				if h.has_key("end") :print local_datetime(h["start"]),"\t   ",time_diff(h["end"]-h["start"]),"\t"," ".join(h["args"])
				else: print local_datetime(h["start"]),"\tincomplete\t"," ".join(h["args"])
else:
	db=shelve.open(".eman2log")
	try:
		n=int(db["count"])
	except:
		print "no logfile"
		sys.exit(0)
	
	for i in range(n-1):
		print " ".join(db[str(i+1)]["args"])
