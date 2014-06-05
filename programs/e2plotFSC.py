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
from EMAN2fsc import *
from emplot2d import EMPlot2DWidget,colortypes
from empmwidgets import PMFSCTableWidget
from emapplication import EMApp

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog Refinement directory [options]
	Plot FSC curvers produced by e2refine, eotest, etc.
	>"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#dirbasename='refine|frealign',
	#parser.add_pos_argument(name="plot_files",help="List the directories to plot here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2)
	parser.add_header(name="filterheader", help='There is no help', title="This program is currently not functional. The table below is still useful,\n but for actual plots, suggest using e2evalrefine.py for now.", row=0, col=0, rowspan=1, colspan=2)
	parser.add_pos_argument(name="fscdir",help="The refinement directory to use for FSC plotting.", default="", guitype='fsctable', row=1, col=0,rowspan=1, colspan=2)
	parser.add_header(name="filterheader", help='Options below this label are specific to e2plotFSC', title="### e2plotFSC options ###", row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--plote2res",action="store_false",help="Plot curves from e2resoltion",default=True,guitype='boolbox',row=3,col=0,rowspan=1,colspan=1)
	parser.add_argument("--plote2eotest",action="store_false",help="Plot curves from e2eotest",default=True,guitype='boolbox',row=3,col=1,rowspan=1,colspan=1)
	parser.add_argument("--plotconvergence",action="store_false",help="Plot curves from refinement convergence",default=True,guitype='boolbox',row=4,col=0,rowspan=1,colspan=1)
	parser.add_argument("--ploteoconvergence",action="store_false",help="Plot curves from refine_even_odd convergence",default=True,guitype='boolbox',row=4,col=1,rowspan=1,colspan=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()

	# Make the QT app
	app = EMApp()
	
	# display table
	if len(args) == 0:
		fsctable = PMFSCTableWidget("fsc","",None,resize=True)
		fsctable.show()
		
	# or let user choose FSC plotting pars 
	else:
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
		eoconv = []
		# Method to the maddness, I use not here because I need to only plot when open is not presented AND I need to keep presentation on in the GUI
		if not options.plote2res: res= get_e2resolution_results_list(keys)
		if not options.plote2eotest: eo = get_e2eotest_results_list(keys)
		if not options.plotconvergence: conv = get_convergence_results_list(keys)
		if not options.ploteoconvergence: eoconv = get_e2refine_even_odd_results_list(keys)
	
		# Plot FSC curves
		i = 0
		max = len(colortypes)		
		for k in conv:
			module.set_data(db[k],k,color=(i%max),linewidth=1) # there are only a ceratin number of  colors
			i += 1
		
		# plot e2refine_even_odd curves
		for k in eoconv:
			module.set_data(db[k],k,color=(i%max),linewidth=2) # there are only a ceratin number of  colors
			i += 1
		
		#plot eo test and res
		for plot in [eo,res]:
			for k in plot:
				module.set_data(db[k],k,color=(i%max),linewidth=3) # there are only a ceratin number of  colors
				i += 1
					
		module.show()
	app.exec_()
	
if __name__ == "__main__":
	main()