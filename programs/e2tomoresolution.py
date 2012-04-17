#!/usr/bin/env python

#
# Author: John Flanagan Oct 25th 2011 (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine
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


from EMAN2 import *
from emplot2d import EMPlot2DWidget
from emapplication import EMApp

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] 
	This is a program to compute the resolution of a n averaged subtomogram. Right now it is very simple simple divide the aligned
	subtomos into even/odd classes, average and then compute the FSC. In the future this program will be extended to compute 
	resolution of an averged subtomo vs a reference and hopefuly of a single sub/tomogram.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_pos_argument(name="tomodir",help="The refinement directory to use for tomoresolution.", default="", guitype='dirbox', dirbasename='spt_',  row=0, col=0,rowspan=1, colspan=2)
	parser.add_header(name="tomoresoheader", help='Options below this label are specific to e2tomoresolution', title="### e2tomoresolution options ###", row=1, col=0, rowspan=1, colspan=2)
	parser.add_argument("--averager",type=str,help="The averager used to generate the averages. Default is \'mean\'.",default="mean", guitype='combobox', choicelist='dump_averagers_list()', row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--sym",  type=str,help="The recon symmetry", default="c1", guitype='symbox', row=3, col=0, rowspan=1, colspan=2)
	parser.add_argument("--mask", type=str,help="The mask to apply before FSC calculation", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=4, col=0, rowspan=1, colspan=2)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	if options.mask: options.mask = parsemodopt(options.mask)
	
	logid=E2init(sys.argv,options.ppid)
	
	fscstrategy = EvenOddReso(args[0], options)
	fscstrategy.execute()
	
	results_db = db_open_dict("bdb:%s#convergence.results"%args[0])
	results_db["tomo_fsc"] = [fscstrategy.getFreq(),fscstrategy.getFSC(),fscstrategy.getError()]
	results_db.close()
	
	E2end(logid)
	
	# Plot FSC
	app = EMApp()
	plot = EMPlot2DWidget()
	plot.set_data((fscstrategy.getFreq(),fscstrategy.getFSC()))
	plot.show()
	app.exec_()

# Use strategy pattern here. Any new stategy needs to inherit this
class Strategy:
	def __init__(self, jobdir, options):
		self.jobdir = jobdir
		self.options = options
		self.fsc = None
		self.freq = None
		self.error = None
		
	def execute(self):
		raise NotImplementedError("Subclass must implement abstract method")
	
	def getFSC(self):
		return self.fsc
		
	def getFreq(self):
		return self.freq
		
	def getError(self):
		return self.error
		
class EvenOddReso(Strategy):
	def __init__(self, jobdir, options):
		Strategy.__init__(self, jobdir, options)
		
	def execute(self):
		tomo_db = db_open_dict("bdb:%s#class_ptcl"%self.jobdir)
		evenavgr = Averagers.get(self.options.averager)
		oddavgr = Averagers.get(self.options.averager)
		
		# Compute the even odd maps
		for tomo in xrange(0, len(tomo_db)/2, 2):
			evenavgr.add_image(tomo_db[tomo*2])
			oddavgr.add_image(tomo_db[tomo*2 + 1])
			
		# Now actually compute the FSC
		evenavg = evenavgr.finish()
		oddavg = oddavgr.finish()
		if self.options.sym != "c1":
			evenavg = evenavg.process('xform.applysym',{'sym':self.options.sym})
			oddavg = oddavg.process('xform.applysym',{'sym':self.options.sym})
			
		if self.options.mask:
			evenavg = evenavg.process(self.options.mask[0],self.options.mask[1])
			oddavg = oddavg.process(self.options.mask[0],self.options.mask[1])
		
		evenavg.write_image('bdb:%s#evenvol'%self.jobdir)
		oddavg.write_image('bdb:%s#oddvol'%self.jobdir)
		
		fscdata = evenavg.calc_fourier_shell_correlation(oddavg)
		
		size = len(fscdata)/3
		self.freq = fscdata[0:size]
		self.fsc = fscdata[size:size*2]
		self.error = fscdata[size*2:size*3]
		
		tomo_db.close()
		
if __name__ == "__main__":
    main()