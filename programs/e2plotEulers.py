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
from emplot2d import EMPolarPlot2DWidget,colortypes
from emapplication import EMApp
import math

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog Refinement directory [options]
	Plot Euler angle distributions for refinement results. Poiont size is proportional to Euler bin count.
	>"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser.add_pos_argument(name="plot_files",help="List the directories to plot here.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2)
	parser.add_pos_argument(name="refinedir",help="The refinement directory to use for FSC plotting.", default="", guitype='dirbox', dirbasename='refine|TiltValidate',  row=0, col=0,rowspan=1, colspan=2)
	parser.add_header(name="eulerheader", help='Options below this label are specific to e2plotEuler', title="### e2plotEuler options ###", row=1, col=0, rowspan=1, colspan=1)
	parser.add_argument("--iteration",type=int,help="Refinement iteration to plot", default=0, guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--pointwidth",type=float,help="The relative scale of the points plotted. The absoule size is dpenedent on particle count", default=1.0, guitype='floatbox', row=2, col=1, rowspan=1, colspan=1)
	parser.add_argument("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.", guitype='symbox', row=3, col=0, rowspan=1, colspan=2)
	parser.add_argument("--norticklabels",action="store_true",help="Disable radius tick labels", guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, default=False)
	parser.add_argument("--nothetaticklabels",action="store_true",help="Disable Theta tick labels", guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
		
	(options, args) = parser.parse_args()

	# Grab the data
	iteration = 0
	data = EMData.read_images("bdb:%s#classify_%02d"%(args[0],iteration))
	projections = EMData.read_images("bdb:%s#projections_%02d"%(args[0],iteration))
	
	# We use a hash data structure to count same Eulers(could also use an array an sort it)
	eulerhash = {}
	for i in xrange(data[0].get_ysize()):
		# Get Eulers
		projnum = int(data[0][i])
		euler = projections[projnum].get_attr('xform.projection')
		
		#Loop over all syms
		for sym in Symmetries.get(options.sym).get_syms():
			eulerangles = (sym*euler).get_rotation('eman')
		
			# Use has to count unique eulers
			hashkey = "%3.2f %3.2f"%(eulerangles['az'],eulerangles['alt'])
			if eulerhash.has_key(hashkey):
				eulerhash[hashkey] += 1
			else:
				eulerhash[hashkey] = 1
	
	# Now plot these eulers
	theta = []
	r = []
	size = []
	for euler, count in eulerhash.items():
		eulers = euler.split()
		theta.append(float(eulers[0]))
		r.append(float(eulers[1]))
		size.append(count*options.pointwidth)
		
	# Make the QT app and plot
	app = EMApp()
	plot = EMPolarPlot2DWidget()
	plot.set_yticklabels(not options.norticklabels)
	plot.set_xticklabels(not options.nothetaticklabels)
	plot.setAxisParms(False,False)
	plot.set_data((theta,r),linewidth=50,radcut=180)
	plot.setPointSizes(size)
	plot.show()
	
	app.exec_()
	
if __name__ == "__main__":
	main()