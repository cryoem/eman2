#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya 03/2011, 
# (based on Steven Ludtke's initial implementation [02/15/2011] of Jesus's older scripts, from M.F.Schmid's methods).
# Last modification: July/08/2015
#
# Copyright (c) 2011 Baylor College of Medicine
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

from copy import deepcopy
import os
import sys

import random
from EMAN2jsondb import JSTask,jsonclasses
import datetime

from e2spt_hac import textwriter

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	This program runs different spt programs quickly, in testing mode, such that crashes
	can be identified more easily.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--testn",type=int,default=6,help="""default=6. size of dataset to run tests with; cannot be < 6, since initial model generation with HAC for gold-standard refinement requires at least 3 particles for the even set and 3 for the odd set.""")
	
	parser.add_argument("--path",type=str,default='spttests',help="""Default=spttests. Directory to store results in. The default is a numbered series of directories containing the prefix 'spttests'; for example, spttests_02 will be the directory by default if 'spttests_01' already exists.""")
	
	parser.add_argument("--parallel",type=str,default='',help="""the program will detect the number of cores available and use threaded parallelism by default. To use only one core, supply --parallel=thread:1. For MPI on clusters, see parallelism at http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")
	
	#parser.add_argument("--testsim",action='store_true',default=False,help="""default=False. If supplied, this option will test e2spt_simulation.py as well and use the generated simulated particles for subsequent tests, opposed to random volumes that do not have a missing wedge, noise or CTF.""")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
		
	(options, args) = parser.parse_args()

	logger = E2init(sys.argv, options.ppid)
	
	if not options.parallel:
		import multiprocessing
		nparallel = multiprocessing.cpu_count()
		options.parallel = 'thread:' + str(nparallel)
		print "\nfound %d cores" %(nparallel)
		print "setting --parallel to", options.parallel
	
	if options.testn < 6:
		print "\nERROR: --testn must be > 5."
		sys.exit()
	
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'spt_bt')
	
	from e2spt_classaverage import writeParameters
	writeParameters(options,'e2spt_test.py', 'spttests')
	
	for i in range( options.testn ):
		a = test_image_3d()
		t=Transform()
		if i > 0:
			az = random.randint(0,360)
			alt = random.randint(0,180)
			phi = random.randint(0,360)
			tx = random.randint(-5,5)
			ty = random.randint(-5,5)
			tz = random.randint(-5,5)
			t = Transform({'type':'eman','tx':tx,'ty':ty,'tz':tz,'alt':alt,'az':az,'phi':phi})
			a.transform(t)
		
		a.process_inplace('math.meanshrink',{'n':4})
		a['spt_randT'] = t
		a.write_image(options.path + '/testimgs.hdf',-1)
		
		if i==0:
			a.write_image(options.path + '/testimg_ref.hdf',0)
		
	cmds = []
	
	rootpath = os.getcwd()
	
	
	os.system('touch ' + options.path + '/output.txt')
	
	#input = rootpath + '/' + options.path + '/testimgs.hdf'
	#if options.testsim:
	
	simcmd = 'e2spt_simulation.py --input=' + options.path + '/testimg_ref.hdf --nptcls ' + str(options.testn) + ' --tiltrange 60 --nslices 25 --saveprjs --applyctf --snr 2' + ' --parallel=' +options.parallel + ' --path testsim'
	if options.verbose:
		simcmd += ' --verbose ' + str(options.verbose)
	
	cmds.append( simcmd )
	simcmd2 = 'mv testsim ' + options.path
	cmds.append( simcmd2 )

	input = rootpath + '/' + options.path +'/testsim/simptcls.hdf'
	
	
	btcmd = 'e2spt_binarytree.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --parallel=' +options.parallel + ' --path testbt'
	cmds.append( btcmd )
	btcmd2 = 'mv testbt ' + rootpath + '/' + options.path + '/'
	cmds.append( btcmd2 )
	
	haccmd = 'e2spt_hac.py  --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --parallel=' +options.parallel + ' --path testhac'
	cmds.append( haccmd )
	haccmd2 = 'mv testhac ' + rootpath + '/' + options.path + '/'
	cmds.append( haccmd2 )
	
	ssacmd = 'e2symsearch3d.py  --input ' + input + ' --sym icos --steps 2 --parallel=' +options.parallel + ' --path testssa'
	cmds.append( ssacmd )
	ssacmd2 = 'mv testssa ' + rootpath + '/' + options.path + '/'
	cmds.append( ssacmd2 )
	
	
	sptdefaultcmdgoldoff = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --goldstandardoff --parallel=' +options.parallel + ' --path testsptdefaultgoldoff'
	cmds.append( sptdefaultcmdgoldoff )
	sptdefaultcmdgoldoff2 = 'mv testsptdefaultgoldoff ' + rootpath + '/' + options.path + '/'
	cmds.append( sptdefaultcmdgoldoff2 )
	
	sptrefbtcmdgoldoff = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --goldstandardoff --btref 2 --parallel=' +options.parallel + ' --path testsptrefbtgoldoff'
	cmds.append( sptrefbtcmdgoldoff )
	sptrefbtcmdgoldoff2 = 'mv testsptrefbtgoldoff ' + rootpath + '/' + options.path + '/'
	cmds.append( sptrefbtcmdgoldoff2 )
	
	sptrefssacmdgoldoff = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --goldstandardoff --ssaref 2 --parallel=' +options.parallel + ' --path testsptrefssagoldoff'
	cmds.append( sptrefssacmdgoldoff )
	sptrefssacmdgoldoff2 = 'mv testsptrefssagoldoff ' + rootpath + '/' + options.path + '/'
	cmds.append( sptrefssacmdgoldoff2 )
	
	sptrefhaccmdgoldoff = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --goldstandardoff --hacref 3 --parallel=' +options.parallel + ' --path testsptrefhacgoldoff'
	cmds.append( sptrefhaccmdgoldoff )
	sptrefhaccmdgoldoff2 = 'mv testsptrefhacgoldoff ' + rootpath + '/' + options.path + '/'
	cmds.append( sptrefhaccmdgoldoff2 )
	
	
	sptdefaultcmdgoldon = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --parallel=' +options.parallel + ' --path testsptdefaultgoldon'
	cmds.append( sptdefaultcmdgoldon )
	sptdefaultcmdgoldon2 = 'mv testsptdefaultgoldon ' + rootpath + '/' + options.path + '/'
	cmds.append( sptdefaultcmdgoldon2 )
	
	sptrefbtcmdgoldon = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --btref 4 --parallel=' +options.parallel + ' --path testsptrefbtgoldon'
	cmds.append( sptrefbtcmdgoldon )
	sptrefbtcmdgoldon2 = 'mv testsptrefbtgoldon ' + rootpath + '/' + options.path + '/'
	cmds.append( sptrefbtcmdgoldon2 )
	
	sptrefssacmdgoldon = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --ssaref 4 --parallel=' +options.parallel + ' --path testsptrefssagoldon'
	cmds.append( sptrefssacmdgoldon )
	sptrefssacmdgoldon2 = 'mv testsptrefssagoldon ' + rootpath + '/' + options.path + '/'
	cmds.append( sptrefssacmdgoldon2 )
	
	sptrefhaccmdgoldon = 'e2spt_classaverage.py --input ' + input + ' --align rotate_symmetry_3d:sym=c1 --falign refine_3d_grid:range=3:delta=3 --hacref 6 --parallel=' +options.parallel + ' --path testsptrefhacgoldon'
	cmds.append( sptrefhaccmdgoldon )
	sptrefhaccmdgoldon2 = 'mv testsptrefhacgoldon ' + rootpath + '/' + options.path + '/'
	cmds.append( sptrefhaccmdgoldon2 )
	
	for cmd in cmds:
		runcmd( options, cmd )
	

	E2end(logger)
	sys.stdout.flush()
	
	
	
	return


def runcmd(options,cmd):
	if options.verbose > 9:
		print "(e2spt_classaverage)(runcmd) running command", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	lines = ['\n']
	for t in text:
		lines.append(t)
	lines.append('\n')
	
	outputfile = options.path +'/output.txt'
	f=open(outputfile,'a')
	f.writelines( lines )
	f.close()
	
	if options.verbose > 9:
		print "(e2spt_classaverage)(runcmd) done"
	return


if __name__ == '__main__':
	main()
