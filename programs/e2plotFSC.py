#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2012- Baylor College of Medicine
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
#

from EMAN2 import *
from emplot2d import EMPlot2DWidget,colortypes
from emapplication import EMApp

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog Refinement directory [options]
	Plot FSC curvers produced by e2refine, eotest, etc.
	>"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser.add_pos_argument(name="plot_files",help="List the directories to plot here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)", positional=True, row=0, col=0,rowspan=1, colspan=2)
	parser.add_pos_argument(name="fscdir",help="The refinement directory to use for FSC plotting.", default="", guitype='dirbox', dirbasename='refine|frealign', positional=True, row=0, col=0,rowspan=1, colspan=2)
	parser.add_header(name="filterheader", help='Options below this label are specific to e2plotFSC', title="### e2plotFSC options ###", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--plote2res",action="store_false",help="Plot curves from e2resoltion",default=True,guitype='boolbox',row=2,col=0,rowspan=1,colspan=1)
	parser.add_argument("--plote2eotest",action="store_false",help="Plot curves from e2eotest",default=True,guitype='boolbox',row=2,col=1,rowspan=1,colspan=1)
	parser.add_argument("--plotconvergence",action="store_false",help="Plot curves from refinement convergence",default=True,guitype='boolbox',row=3,col=0,rowspan=1,colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()

	# Make the QT app
	app = EMApp()
	module = EMPlot2DWidget()
	
	# Get data from data base
	db_name = "bdb:"+args[0]+"#convergence.results"
	if not db_check_dict(db_name):
		print "Rubbish!!!, no FSC curves found!!!"
		return
	db = db_open_dict(db_name,ro=True)
	keys = db.keys()
	
	# Load desired FSC curves
	res = []
	eo = []
	conv = []
	# Method to the maddness, I use not here because I need to only plot when open is not presented AND I need to keep presentation on in the GUI
	if not options.plote2res: res= get_e2resolution_results_list(keys)
	if not options.plote2eotest: eo = get_e2eotest_results_list(keys)
	if not options.plotconvergence: conv = get_convergence_results_list(keys)
	
	# Plot FSC curves
	i = 0
	max = len(colortypes)		
	for k in conv:
		module.set_data(db[k],k,color=(i%max),linewidth=1) # there are only a ceratin number of  colors
		i += 1
			
	for plot in [eo,res]:
		for k in plot:
			module.set_data(db[k],k,color=(i%max),linewidth=3) # there are only a ceratin number of  colors
			i += 1
					
	module.show()
	app.exec_()

def get_e2resolution_results_list(keys):
		'''
		Extract the names from the keys that match the e2resolution.py output naming convention
		(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
		'''
		solns = []
		for k in keys:
			if len(k) > 6 and k[-7:] == "res_fsc":
				solns.append(k)
		if not solns: print "Rubbish!!!, no e2resolution FSC curves found!"
		solns.sort()
		return solns
	
def get_e2eotest_results_list(keys):
	'''
	Extract the names from the keys that match the e2eotest.py output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	for k in keys:
		if len(k) > 7 and k[0:8] == "even_odd":
			solns.append(k)
	if not solns: print "Rubbish!!!, no e2eotest FSC curves found!"
	solns.sort()
	return solns

def get_convergence_results_list(keys):
	'''
	Extract the names from the keys that match the e2refine.py convergence plot output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	if "init_00_fsc" in keys:
		solns.append("init_00_fsc")
		
	i = 0
	while True:
		s1 = str(i)
		s2 = str(i+1)
		if len(s1) == 1: s1 = "0"+s1
		if len(s2) == 1: s2 = "0"+s2
		k = s1+"_"+s2+"_fsc"
		if k in keys:
			solns.append(k)
		else:
			break

		i += 1

	return solns
	
if __name__ == "__main__":
	main()