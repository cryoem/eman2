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
from EMAN2 import base_name, EMArgumentParser
import EMAN2db

# also defined in EMAN2, but we don't want to have to import it

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	A tool for displaying EMAN2 command history
	"""
	
	parser = EMArgumentParser(usage=usage)
	parser.add_argument("--gui", "-g",default=False, action="store_true",help="Open history in an interface with a sortable table.")
	parser.add_argument("--all", "-a",default=False, action="store_true",help="Show for all directories.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if options.gui:
		from emapplication import EMApp
		app = EMApp()
		hist = HistoryForm(app,os.getcwd())
		app.show()
		app.execute()
		
	else: print_to_std_out(options.all)

class HistoryForm:
	def __init__(self,application,wd):
		'''
		wd is the working directory
		'''
		self.wd = wd
		
		from emform import EMFormWidget
		self.form = EMFormWidget(params=self.get_history_table())
		self.form.setWindowTitle("EMAN2 history")
		
		from PyQt4 import QtGui,QtCore
		self.form.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/feather.png"))
		self.form.resize(640,480)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_cancel)
		
		
 	def get_history_table(self):
 		from emdatastorage import ParamDef
 		try:
			import EMAN2db
			db=EMAN2db.EMAN2DB.open_db()
			db.open_dict("history")
		except:
			db=None
		
		params = []
		
		try:
			n=int(db.history["count"])
		except:
			n = 0
		
		if db == None or n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits="There appears to be no history in this directory",choices=None))
		else:
			from emform import EMParamTable
			start = []
			duration = []
			prgargs = []
			cmd = []
			full_cmd = []
			
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits="Use this form to examine the EMAN2 commands that have occurred in this directory.",choices=None))
			p = EMParamTable(name="commands",desc_short="Historical table of EMAN2 commands",desc_long="") 
			
			for i in range(n):
				try: h=db.history[i+1]
				except: continue
				
				if h != None and h.has_key("path") and h["path"]==self.wd:
					start.append(local_datetime(h["start"]))
					if h.has_key("end") :
						duration.append(time_diff(h["end"]-h["start"]))
					else:
						if h.has_key("progress"):
							try:
								duration.append(str(int(h["progress"]*100))+"%")
							except:
								duration.append("incomplete")
						else: duration.append("incomplete")
					
					
					args = h["args"]
					if len(args) > 0:
						cmd.append(base_name(args[0]))
						full_cmd.append(args[0])
						if len(args) > 1:
							prgargs.append(" ".join(args[1:]))
						else:prgargs.append("")
					else:
						cmd.append("")
						full_cmd.append("")
						prgargs.append("")
						
				
			pcmd = ParamDef(name="cmd",vartype="stringlist",desc_short="Program",desc_long="A shortened version of the name of a Python program that was executed",property=None,defaultunits=None,choices=cmd)
			pstart = ParamDef(name="start",vartype="stringlist",desc_short="Start time",desc_long="The time when the command was first executed",property=None,defaultunits=None,choices=start)
			pduration = ParamDef(name="duration",vartype="stringlist",desc_short="Duration",desc_long="The time taken to execute this command",property=None,defaultunits=None,choices=duration)
			pfull_cmd = ParamDef(name="fullcmd",vartype="stringlist",desc_short="File path",desc_long="The location of program on disk",property=None,defaultunits=None,choices=full_cmd)
			pargs = ParamDef(name="args",vartype="stringlist",desc_short="Arguments",desc_long="The arguments that were given to the program",property=None,defaultunits=None,choices=prgargs)
			
	
			p.append(pcmd)
			p.append(pstart)
			p.append(pduration)
			p.append(pargs)
			p.append(pfull_cmd)
			params.append(p)
	
	
 		return params
 	
 	def on_ok(self):
 		self.form.close()
 		
 	def on_cancel(self):
 		self.form.close()
 			
	
def local_datetime(secs):
	from time import localtime
	t=localtime(secs)
	return "%04d/%02d/%02d %02d:%02d:%02d"%t[:6]

def local_date(secs):
	from time import localtime
	t=localtime(secs)
	return "%04d/%02d/%02d"%t[:3]

def local_time(secs):
	from time import localtime
	t=localtime(secs)
	return "%02d:%02d:%02d"%t[3:6]

def time_diff(secs):
	if secs<3600 : return "%d:%02d"%(secs/60,secs%60)
	return "%d:%02d:%02d"%(secs/3600,(secs%3600)/60,secs%60)

#if len(sys.argv)>1 and sys.argv[1]=="--help" :
#	print "Usage:\ne2history [--all]\n"
#	sys.exit(1)



def print_to_std_out(all):

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
		
		if all  :
			ah={}
			for i in range(n):
				try: h=db.history[i+1]
				except: print "Entry ",i," missing"
				try : ah.setdefault(h["path"],[]).append(h)
				except: continue
			for k in ah.keys():
				print "---------- ",k
				for i in ah[k]:
					if i.has_key("end") : print local_datetime(i["start"]),"\t   ",time_diff(i["end"]-i["start"]),"\t"," ".join(i["args"])
					elif i.has_key("progress") : print local_datetime(i["start"]),"\t ",int(i["progress"]*100)," % done\t"," ".join(i["args"])
					else: print local_datetime(i["start"]),"\tincomplete\t"," ".join(i["args"])
		else:
			for i in range(n):
				try: h=db.history[i+1]
				except: print "Entry ",i," missing"
				if h != None and h.has_key("path") and h["path"]==os.getcwd():
					if h.has_key("end") :print local_datetime(h["start"]),"\t   ",time_diff(h["end"]-h["start"]),"\t"," ".join(h["args"])
					elif h.has_key("progress") :print local_datetime(h["start"]),"\t   ",int(h["progress"]*100)," % done\t"," ".join(h["args"])
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
		print "done"
		
	try:
		fin=file(".eman2log.txt","r")
		print "----- .eman2log.txt"
		for l in fin: print l
		
	except: pass
	
	print "NOTICE: While e2history.py continues to function, log files are now stored in plain text in .eman2log.txt in each directory, and can just be examined directly. e2history.py is primarily useful for looking at historical logfiles."
		
if __name__ == '__main__':
	main()
