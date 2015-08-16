#!/usr/bin/env python
#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

#  VERSION 2  09/25/2014

from   EMAN2 import *
from   sparx import *
import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys
from   random import shuffle
import subprocess



"""
0   1   2   3   4   5
A   A   A   =   =   =
B   =   =   B   B   =
=   C   =   C   =   C
=   =   D   =   D   D

Rotations pair-wise
0-1  A
0-2  A
0-3  B
0-4  B
1-2  A
1-3  C
1-5  C
2-4  D
2-5  D
3-4  B
3-5  C
4-5  D
"""



#  data = [[alpha1,sx1,sy1], [alpha2,sx2,sy2], ...]
def average2dtransform(data, return_avg_pixel_error=False):

	from math import pi, sin, cos, radians, degrees

	L = len(data)

	sum_cosa = 0.0
	sum_sina = 0.0
	sx = 0.0
	sy = 0.0
	for j in xrange(L):
		sum_cosa += cos(radians(data[j][0]))
		sum_sina += sin(radians(data[j][0]))
		sx +=  data[j][1]
		sy +=  data[j][2]
	sx /= L
	sy /= L
	sqrtP = sqrt(sum_cosa**2+sum_sina**2)

	# Get ave transform params
	H = Transform({"type":"2D"})
	H.set_matrix([sum_cosa/sqrtP, sum_sina/sqrtP, 0.0, sx, -sum_sina/sqrtP, sum_cosa/sqrtP, 0.0, sy, 0.0, 0.0, 1.0, 0.0])
	dd = H.get_params("2D")

	H = Transform({"type":"2D","alpha":dd[ "alpha" ],"tx":dd[ "tx" ],"ty": dd[ "ty" ],"mirror":0,"scale":1.0})
	dd = H.get_params("2D")
	if return_avg_pixel_error:
		return sum(sqr_pixel_error)/N
	else:
		return [dd[ "alpha" ], dd[ "tx" ], dd[ "ty" ]]

'''
def rotate_angles(angleset1, angleset2, rot, indexes=None):
	"""
	  It returns list of angles describing differences between the given anglesets (rotations of anglesets don't matter).
	  The function works as follow:
	  1. use the rotation_between_anglesets function to find the transformation between the anglesets
	  2. apply the transformation to the first angleset (to fit it into the second angleset, which is the template)
	  3. calculate angles between corresponding projections directions from the anglesets 
	  INPUT: each of the parameters should be a list of pairs [phi, theta] (additional parameters (pdi, shifts) don't influence on the results)
	  OUTPUT: list of floats - angles in degrees (the n-th element of the list equals the angle between n-th projections directions from the anglesets)
	  The third parameter (indexes) is optional and may be set to list of indexes. In that case only elements from given list are taken into account.
	"""
	from EMAN2 import Transform, Vec2f
	
	"""
	if indexes != None:
		new_ang1 = []
		new_ang2 = []
		for i in indexes:
			new_ang1.append(angleset1[i])
			new_ang2.append(angleset2[i])
		angleset1 = new_ang1
		angleset2 = new_ang2
	"""
	rot = rotation_between_anglesets(angleset1, angleset2)
	T2 = Transform({"type":"spider","phi":rot[0],"theta":rot[1],"psi":rot[2]})
	angles = []
	angle_errors = []
	for i in xrange(len(angleset1)):
		T1 = Transform({"type":"spider","phi":angleset1[i][0],"theta":angleset1[i][1],"psi":angleset1[i][2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
		T1.set_trans(Vec2f(-angleset1[i][3], -angleset1[i][4]))
		T = T1*T2
		d = T.get_params("spider")
		d = T.get_params("spider")
		angleset1[i] = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]]
		#angle_errors.append( angle_between_projections_directions(angles[i], angleset2[i]) )
		#return angle_errors
'''

def apply_rotation(angleset1, qrot):
	from EMAN2 import Transform, Vec2f

	"""
	if indexes != None:
		new_ang1 = []
		new_ang2 = []
		for i in indexes:
			new_ang1.append(angleset1[i])
			new_ang2.append(angleset2[i])
		angleset1 = new_ang1
		angleset2 = new_ang2
	"""
	import types
	if(type(qrot) == types.ListType): rot = Transform({"type":"spider","phi":qrot[0],"theta":qrot[1],"psi":qrot[2]})
	else:                             rot = qrot
	for i in xrange(len(angleset1)):
		T1 = Transform({"type":"spider","phi":angleset1[i][0],"theta":angleset1[i][1],"psi":angleset1[i][2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
		T1.set_trans(Vec2f(-angleset1[i][3], -angleset1[i][4]))
		T = T1*rot
		d = T.get_params("spider")
		angleset1[i] = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]]

def apply_rotation_gaps(angleset1, qrot):
	from EMAN2 import Transform, Vec2f
	import types
	if(type(qrot) == types.ListType): rot = Transform({"type":"spider","phi":qrot[0],"theta":qrot[1],"psi":qrot[2]})
	else:                             rot = qrot
	for i in xrange(len(angleset1)):
		if(len(angleset1[i]) > 1):
			T1 = Transform({"type":"spider","phi":angleset1[i][0],"theta":angleset1[i][1],"psi":angleset1[i][2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
			T1.set_trans(Vec2f(-angleset1[i][3], -angleset1[i][4]))
			T = T1*rot
			d = T.get_params("spider")
			angleset1[i] = [d["phi"], d["theta"], d["psi"], -d["tx"], -d["ty"]]


def errors_per_image_original(nn, qt, asi, avgtrans, thresherr=1.0, radius = 1.0):
	#print "errors per image"
	#  We have average parameters and transformed parameters
	#  Find errors per image
	#  Find outliers, remove them, repeat calculations
	perr = [None]*nn
	for k in xrange(nn):
		lin = k//(nn//4)
		#  
		fifi = []
		for j in xrange(6):
			if(asi[lin][j] != -10):
				#print  " fff  ",k,lin,j,asi[lin][j],qt[j][k]
				fifi.append(qt[j][k])
		assert(len(fifi) == 3)
		nas = [0.0,0.0,0.0]
		r1,r2,r3 = getfvec(avgtrans[k][0],avgtrans[k][1])
		sacos = 0.0
		sder = 0.0
		ser3 = 0.0
		for i in xrange(3):
			#print i,fifi[i]
			d = max_3D_pixel_error(fifi[i], avgtrans[k], r=radius)
			#if(d>10.):  print  "  LARGE ERROR",k,i,d,fifi[i], avgtrans[k]
			ser3 += d
			n1,n2,n3 = getfvec(fifi[i][0],fifi[i][1])
			sacos += acos(min(1.0,n1*r1+n2*r2+n3*r3))
			sder += pixel_error_2D(fifi[i][2:], avgtrans[k][2:], r = 1.0)
		# average deviation in radians
		sacos /= 3.0
		sare   = tan(sacos)
		sder  /= 3.0
		ser3  /= 3
		#print  k, ser3, sare, sder
		perr[k] = [k, ser3, sacos, sare, sder]

	#write_text_row(perr, 'per.txt')
	perr = [perr[k][1] <= thresherr for k in xrange(nn)]
	return perr

def errors_per_image(params, avgtrans, thresherr=1.0, radius = 1.0):
	#print "errors per image"
	#  We have average parameters and transformed parameters
	#  Find errors per image
	ll = len(params)
	nn = len(params[0])
	perr = [None]*nn
	for k in xrange(nn):
		#  
		nas = [0.0,0.0,0.0]
		r1,r2,r3 = getfvec(avgtrans[k][0],avgtrans[k][1])
		sacos = 0.0
		sder = 0.0
		ser3 = 0.0
		for i in xrange(ll):
			#print i,fifi[i]
			d = max_3D_pixel_error(params[i][k], avgtrans[k], r=radius)
			#if(d>10.):  print  "  LARGE ERROR",k,i,d,fifi[i], avgtrans[k]
			ser3 += d
			n1,n2,n3 = getfvec(params[i][k][0],params[i][k][1])
			sacos += acos( min(1.0, n1*r1+n2*r2+n3*r3) )
			sder += pixel_error_2D(params[i][k][2:], avgtrans[k][2:], r = 1.0)
		# average deviation in radians
		sacos /= ll
		sare   = tan(sacos)
		sder  /= ll
		ser3  /= ll
		#print  k, ser3, sare, sder
		perr[k] = [ser3, sacos, sare, sder]

	#write_text_row(perr, 'per.txt')
	#perr = [perr[k][1] <= thresherr for k in xrange(nn)]
	return perr


def average_trans(params):
	#  Compute average projection params and pixel errors
	from utilities import getfvec
	from pixel_error import max_3D_pixel_error
	from math import sqrt, degrees, radians, acos
	nn = lem(params[0])
	avgtrans = [None]*nn
	pixer   = [0.0]*nn
	for j in xrange(nn):
		nas = [0.0,0.0,0.0]
		if( sym == 1):
			#pixer[q[0][0]][j] = max_3D_pixel_error(fifi[0], fifi[1], r=radius)
			for i in xrange(ll):
				n1,n2,n3 = getfvec(params[i][j][0],params[i][j][1])
				nas[0] += n1
				nas[1] += n2
				nas[2] += n3
		else:
			m1,m2,m3 = getfvec(params[0][j][0],params[0][j][1])
			nas[0] = m1
			nas[1] = m2
			nas[2] = m3
			for i in xrange(1,ll):
				qnom = -1.e10
				for j in xrange(-1,2,1):
					t1,t2,t3 = getfvec(params[i][j][0]+j*qsym,params[i][j][1])
					nom = t1*m1 + t2*m2 + t3*m3
					if(nom > qnom):
						qnom = nom
						n1=t1
						n2=t2
						n3=t3
					#if(k == 2):
					#	print '  t1,t2,t3 ',fifi[i][0]+j*qsym,fifi[i][1],t1,t2,t3,nom,qnom
				nas[0] += n1
				nas[1] += n2
				nas[2] += n3
				print qnom, n1,n2,n3,nas
			#  To get the correct pixer phi angle has to be taken from the above!!

		nom = sqrt(nas[0]**2 + nas[1]**2 + nas[2]**2)

		if(nom < 1.e-6):
			nphi   = 0.0
			ntheta = 0.0
		else:
			ntheta = degrees(acos(nas[2]/nom))%360.0
			if(sym>1 and ntheta>90.0):  nphi   = (degrees(atan2( nas[1], nas[0] ))-180.0)%qsym + 180.0
			else:                       nphi   =  degrees(atan2( nas[1], nas[0] ))%qsym

		#print   "FIFI     %4d     %7.2f     %7.2f    %7.2f    %7.2f     %7.2f     %7.2f    %7.2f    %7.2f"%(k,fifi[0][0],fifi[0][1],fifi[1][0],fifi[1][1],fifi[2][0],fifi[2][1],nphi,ntheta)
		twod = average2dtransform([params[ii][j][2:] for ii in xrange(ll)])
		avgtrans[j] = [nphi, ntheta, twod[0], twod[1], twod[2]]
	return avgtrans

def main():
	arglist = []
	for arg in sys.argv:	arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " stack outdir --phase=1 --ou=outer_radius|sxconsistency.py --phase=3 newlocal/main000 --ou=133 --thresherr=3.0 --params=paramsa outgrouparms"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--phase",     type= "int",         default= 1,       help="Phase =1 prepares resampled stacks, =2 analyzes consistency of orientation parameters")
	parser.add_option("--ou",        type= "int",         default= -1,      help="outer radius for calculation of pixel error")
	parser.add_option("--sym",       type="string",       default= "c1",    help="symmetry of the refined structure")
	parser.add_option("--thresherr", type="float",        default=1.0,      help="Threshold for accpetable orientation errors (in pixels)")
	parser.add_option("--ndigits",   type="int",          default=1,        help="Accuracy for checking whether parameters are identical")
	parser.add_option("--chunk",     type="string",       default="",       help="Root of of four chunk files with indeces")
	parser.add_option("--params",    type="string",       default="",       help="Root of of six parameter file names with refinement results")
	(options, args) = parser.parse_args(arglist[1:])
	global_def.BATCH = True
	if options.phase == 1 and len(args) == 2:
		inputbdb = args[0]
		outdir   = args[1]
		nn = EMUtil.get_image_count(inputbdb)
		t = range(nn)
		shuffle(t)
		chunks = []
		for i in xrange(4):
			#  I use the MPI function here just to easily get the balanced load
			j,k = MPI_start_end(nn, 4, i)
			chunks.append(t[j:k])
			chunks[i].sort()
			write_text_file(chunks[i],os.path.join(outdir,'chunk%01d.txt'%i))

		del t
		write_text_file([ len(chunks[i]) for i in xrange(4) ],os.path.join(outdir,'chunklengths.txt'))

		"""
		pt = [[None]]*6
		ll=0
		for i in xrange(3):
			for j in xrange(i+1,4):
				pt[ll] = chunks[i]+chunks[j]
				ll+=1

		for i in xrange(6):
			listfile = os.path.join(outdir,'lili%01d.txt'%i)
			write_text_file(pt[i],listfile)
			outbdb = "bdb:"+ os.path.join(outdir,"X%01d"%i)
			cmd = '{} {} {} {}'.format('e2bdb.py', inputbdb, '--makevstack='+outbdb, '--list='+listfile)
			subprocess.call(cmd, shell=True)
		"""
		#  Run 6 programs

	elif options.phase == 2 and len(args) == 2:
		outdir         = args[0]
		howmanythesame = args[1]
		ndigits = options.ndigits  #for chc5 1, for ribo 4#4.0
		prms = []
		for i in xrange(6):
			prms.append( read_text_row(os.path.join(outdir,options.params+"%01d.txt"%i)) )
			for j in xrange(len(prms[-1])):
				for k in xrange(5):
					prms[-1][j][k] = round(prms[-1][j][k], ndigits)

		nn = 2*len(prms[0])
		n4 = nn//4
		qt=[[[-1.0]]*nn for i in xrange(6)]
		ll=0
		for i in xrange(3):
			for j in xrange(i+1,4):
				qt[ll][i*n4:(i+1)*n4]=prms[ll][:n4]
				qt[ll][j*n4:(j+1)*n4]=prms[ll][n4:]
				ll+=1

		thesame = 0

		for ll in xrange(len(qt[0])):
			rw = []
			for j in xrange(6):
				if(len(qt[j][ll]) > 1):  rw.append(qt[j][ll])
			isame = True
			for j in xrange(3):
				if(rw[0][j] != rw[1][j]):
					isame = False
					#print  ll,rw[0][j], rw[1][j]
					break
				if(rw[0][j] != rw[2][j]):
					isame = False
					#print  ll,rw[0][j], rw[2][j]
					break
			if isame: thesame += 1
		qt = float(thesame)/nn
		print "Proportion of the same orientations ",qt
		write_text_file([qt], os.path.join(outdir,howmanythesame) )

	#########     PHASE 3
	elif options.phase == 3 and len(args) == 2:
		outdir      = args[0]
		outgrouparms= args[1]
		radius = options.ou
		thresherr = options.thresherr  #for chc5 1, for ribo 4#4.0
		sym = int(options.sym[1:])
		qsym = 360.0/sym

		blocks = ['A','B','C','D']
		params = {}
		ll = 0
		for i in xrange(3):
			for j in xrange(i+1,4):
				params[chr(65+i)+chr(48+ll)]=[]
				params[chr(65+j)+chr(48+ll)]=[]
				ll+=1

		chunks = {}
		chunklengths = {}
		for i in xrange(4):
			chunks[chr(65+i)] = map(int,read_text_file(os.path.join(outdir,options.chunk+"%01d.txt"%i)))
			chunklengths[chr(65+i)] = len(chunks[chr(65+i)])


		for i in xrange(6):
			prms = read_text_row(os.path.join(outdir,options.params+"%01d.txt"%i))
			for q in blocks:
				if params.has_key(q+chr(48+i)):
					params[q+chr(48+i)] = prms[:chunklengths[q]]
					del prms[:chunklengths[q]]


		pairs = [["A0","A1"],["B3","B4"],["C1","C5"],["D2","D4"]]
		lefts = ["A2","B0","C3","D5"]

		#  Compute average projection params and pixel errors
		avgtrans = {}
		pixer    = {}
		for i,q in enumerate(pairs):
			avgtrans[q[0][0]] = [0.0]*chunklengths[q[0][0]]
			pixer[q[0][0]]    = [0.0]*chunklengths[q[0][0]]
			for j in xrange(chunklengths[q[0][0]]):
				fifi = [ params[q[0]][j], params[q[1]][j] ]

				nas = [0.0,0.0,0.0]
				if( sym == 1):
					pixer[q[0][0]][j] = max_3D_pixel_error(fifi[0], fifi[1], r=radius)
					for i in xrange(2):
						n1,n2,n3 = getfvec(fifi[i][0],fifi[i][1])
						nas[0] += n1
						nas[1] += n2
						nas[2] += n3
				else:
					m1,m2,m3 = getfvec(fifi[0][0],fifi[0][1])
					nas[0] = m1
					nas[1] = m2
					nas[2] = m3
					#if(k == 2):
					#	print  "XXXX"
					#	print fifi[0],nas
					for i in xrange(1,2):
						qnom = -1.e10
						for j in xrange(-1,2,1):
							t1,t2,t3 = getfvec(fifi[i][0]+j*qsym,fifi[i][1])
							nom = t1*m1 + t2*m2 + t3*m3
							if(nom > qnom):
								qnom = nom
								n1=t1
								n2=t2
								n3=t3
							#if(k == 2):
							#	print '  t1,t2,t3 ',fifi[i][0]+j*qsym,fifi[i][1],t1,t2,t3,nom,qnom
						nas[0] += n1
						nas[1] += n2
						nas[2] += n3
						print qnom, n1,n2,n3,nas
					#  To get the correct pixer phi angle has to be taken from the above!!
					pixer[q[0][0]][j] = max_3D_pixel_error(fifi[0], fifi[1], r=radius)

				nom = sqrt(nas[0]**2 + nas[1]**2 + nas[2]**2)

				if(nom < 1.e-6):
					nphi   = 0.0
					ntheta = 0.0
				else:
					ntheta = degrees(acos(nas[2]/nom))%360.0
					if(sym>1 and ntheta>90.0):  nphi   = (degrees(atan2( nas[1], nas[0] ))-180.0)%qsym + 180.0
					else:                       nphi   = degrees(atan2( nas[1], nas[0] ))%qsym

				#print   "FIFI     %4d     %7.2f     %7.2f    %7.2f    %7.2f     %7.2f     %7.2f    %7.2f    %7.2f"%(k,fifi[0][0],fifi[0][1],fifi[1][0],fifi[1][1],fifi[2][0],fifi[2][1],nphi,ntheta)
				twod = average2dtransform([fifi[0][2:],fifi[1][2:]])
				avgtrans[q[0][0]][j] = [nphi, ntheta, twod[0], twod[1], twod[2]]


		perr = {}
		tgood = 0
		for q in blocks:
			perr[q] = [True]*chunklengths[q]
			for k in xrange(chunklengths[q]):
				if(pixer[q][k] > thresherr):  perr[q][k] = False
				if  perr[q][k]: tgood += 1
			if(tgood < 4):
				print "  No good images within the pixel error threshold specified"
				exit()
		print  " tgood ", tgood
		hi = hist_list([pixer[q][k] for q in blocks for k in xrange(chunklengths[q])  ],16)
		for i in xrange(len(hi[0])):
			print  "%4d   %12.3f    %12.0f "%(i,hi[0][i],hi[1][i])
		#  Finished, store average orientation params and table of good images



		#  store lists of good images for each group
		#  blocks is indexed by first letter
		good = [[] for i in xrange(4)]
		bad  = [[] for i in xrange(4)]
		for i,q in enumerate(blocks):
			for k in xrange(chunklengths[q]):
				#deprt = max_3D_pixel_error(params[lefts[i]][k],avgtrans[q][k],r=radius)
				if  perr[q][k]:
					good[i].append(chunks[q][k])
					#[chunks[q][k],pixer[q][k],deprt, params[pairs[i][0]][k],params[pairs[i][1]][k],params[lefts[i]][k],avgtrans[q][k] ])
				else:
					bad[i].append(chunks[q][k])
					#[chunks[q][k],pixer[q][k],deprt, params[pairs[i][0]][k],params[pairs[i][1]][k],params[lefts[i]][k],avgtrans[q][k]])
			write_text_file( good[i], os.path.join(outdir,"newgood%01d.txt"%i) )
			write_text_file( bad[i],  os.path.join(outdir,"newbad%01d.txt"%i) )

			#  write out parameters, for those in pairs write out average, for leftouts leave them as they were
			#  These parameters refer to the original X files.
			ll = 0
			for m in xrange(6):
				if params.has_key(q+chr(48+m)):
					prmsgood = []
					prmsbad  = []
					try:
						j = lefts.index(q+chr(48+m))
						for k in xrange(chunklengths[q]):
							if  perr[q][k]:
								prmsgood.append(params[q+chr(48+m)][k])
							else:
								prmsbad.append(params[q+chr(48+m)][k])
					except:
						for k in xrange(chunklengths[q]):
							if  perr[q][k]:
								prmsgood.append(avgtrans[q][k])
							else:
								prmsbad.append(avgtrans[q][k])
					write_text_row(prmsgood, os.path.join(outdir,"params-newgood%01d%01d.txt"%(i,ll)))
					write_text_row(prmsbad, os.path.join(outdir, "params-newbad%01d%01d.txt"%(i,ll)))
					ll += 1




		#  Generate newlili files from newgood, these contain original numbering of the total single file
		ll = 0
		for i in xrange(3):
			for j in xrange(i+1,4):
				write_text_file( good[i] + good[j], os.path.join(outdir,"newlili%01d.txt"%ll) )
				ll += 1

		#  write out parameters, for those in pairs write out average, for leftouts leave them as they were
		#  These parameters refer to the original X files.
		for i in xrange(6):
			prms = []
			for q in blocks:
				if params.has_key(q+chr(48+i)):
					try:
						j = lefts.index(q+chr(48+i))
						prms += params[q+chr(48+i)]
					except:
						prms += avgtrans[q]
			write_text_row(prms, os.path.join(outdir,outgrouparms+"%01d.txt"%i))


		#  Write chunklengths
		chunklengths = [len(good[i]) for i in xrange(4)]
		write_text_file(chunklengths, os.path.join(outdir,"chunklengths.txt") )
		"""
		#  We do not use consecutive numbering anymore
		#  Generate newx files from newgood, these contain consecutive (with gaps) numbering that allows to generate truncated X files from the previous X files
		ll = 0
		for i in xrange(3):
			firstblock = []
			l = len(perr[blocks[i]])
			for k in xrange(l):
				if  perr[blocks[i]][k]:
					firstblock.append(k)
			for j in xrange(i+1,4):
				secondblock = []
				for k in xrange(len(perr[blocks[j]])):
					if  perr[blocks[j]][k]:
						secondblock.append(k+l)

				write_text_file( firstblock + secondblock, os.path.join(outdir,"goodX%01d.txt"%ll) )
				ll += 1
		del good,bad,firstblock,secondblock,perr
		"""

		"""
		0   1   2   3   4   5

		A   A   A   =   =   =
		B   =   =   B   B   =
		=   C   =   C   =   C
		=   =   D   =   D   D

		"""

	elif options.phase == 4 and len(args) == 1:
		outdir      = args[0]
		bp = 'badparams'
		#outgrouparms= args[1]
		radius = options.ou
		thresherr = options.thresherr
		sym = int(options.sym[1:])
		qsym = 360.0/sym
		#params = [[None for i in xrange(3)] for j in xrange(4)]
		ll = 3 # this is hardwired as we have three groups.  however, I would like to keep the code general.
		for jj in xrange(4):
			params = [None for ii in xrange(ll)]
			newbad = map(int, read_text_file(options.params+"%01d.txt"%jj) )
			nn = len(newbad)
			for ii in xrange(ll):
				params[ii] = read_text_row(os.path.join(outdir,bp+"%01d%01d.txt"%(jj,ii)))
				assert(nn == len(params[ii]) )

			#  Compute average projection params and pixel errors
			avgtrans = [None]*nn
			pixer   = [0.0]*nn
			for j in xrange(nn):
				nas = [0.0,0.0,0.0]
				if( sym == 1):
					#pixer[q[0][0]][j] = max_3D_pixel_error(fifi[0], fifi[1], r=radius)
					for i in xrange(ll):
						n1,n2,n3 = getfvec(params[i][j][0],params[i][j][1])
						nas[0] += n1
						nas[1] += n2
						nas[2] += n3
				else:
					m1,m2,m3 = getfvec(params[0][j][0],params[0][j][1])
					nas[0] = m1
					nas[1] = m2
					nas[2] = m3
					for i in xrange(1,ll):
						qnom = -1.e10
						for j in xrange(-1,2,1):
							t1,t2,t3 = getfvec(params[i][j][0]+j*qsym,params[i][j][1])
							nom = t1*m1 + t2*m2 + t3*m3
							if(nom > qnom):
								qnom = nom
								n1=t1
								n2=t2
								n3=t3
							#if(k == 2):
							#	print '  t1,t2,t3 ',fifi[i][0]+j*qsym,fifi[i][1],t1,t2,t3,nom,qnom
						nas[0] += n1
						nas[1] += n2
						nas[2] += n3
						print qnom, n1,n2,n3,nas
					#  To get the correct pixer phi angle has to be taken from the above!!

				nom = sqrt(nas[0]**2 + nas[1]**2 + nas[2]**2)

				if(nom < 1.e-6):
					nphi   = 0.0
					ntheta = 0.0
				else:
					ntheta = degrees(acos(nas[2]/nom))%360.0
					if(sym>1 and ntheta>90.0):  nphi   = (degrees(atan2( nas[1], nas[0] ))-180.0)%qsym + 180.0
					else:                       nphi   =  degrees(atan2( nas[1], nas[0] ))%qsym

				#print   "FIFI     %4d     %7.2f     %7.2f    %7.2f    %7.2f     %7.2f     %7.2f    %7.2f    %7.2f"%(k,fifi[0][0],fifi[0][1],fifi[1][0],fifi[1][1],fifi[2][0],fifi[2][1],nphi,ntheta)
				twod = average2dtransform([params[ii][j][2:] for ii in xrange(ll)])
				avgtrans[j] = [nphi, ntheta, twod[0], twod[1], twod[2]]

			perr =  errors_per_image(params, avgtrans, thresherr, radius )
			rescued = []
			rejects = []
			for j in xrange(nn):
				#print  chr(65+jj),j,perr[j][0],[[params[m][j][i] for i in xrange(2)] for m in xrange(ll)]
				if( perr[j][0] <= thresherr ):  rescued.append([newbad[j],j])
				else:							rejects.append([newbad[j],j])
			if( len(rescued) == 0 ):    write_text_row([-1,-1],  os.path.join(outdir,"rescued%01d.txt"%jj))
			else:                  		write_text_row(rescued,  os.path.join(outdir,"rescued%01d.txt"%jj))
			#  We also have to write params.
			if( len(rescued) != 0 ):
				for ii in xrange(ll):  write_text_row([params[ii][rescued[k][1]] for k in xrange(len(rescued))]  ,  os.path.join(outdir,"params-rescued%01d%01d.txt"%(jj,ii)))
			if( len(rejects) == 0 ):    write_text_row([-1,-1],  os.path.join(outdir,'rejects'+"%01d.txt"%jj))
			else:                  		write_text_row(rejects,  os.path.join(outdir,'rejects'+"%01d.txt"%jj))
			if( len(rejects) != 0 ):
				for ii in xrange(ll):  write_text_row([params[ii][rejects[k][1]] for k in xrange(len(rejects))]  ,  os.path.join(outdir,"params-rejects%01d%01d.txt"%(jj,ii)))

			hi = hist_list([perr[j][0]  for j in xrange(nn)  ], 16)
			print  "Pixel errors for BAD GROUP  ",chr(65+jj)
			for ii in xrange(len(hi[0])):
				print  "%4d   %12.3f    %12.0f "%(ii,hi[0][ii],hi[1][ii])

	elif options.phase == 5 and len(args) == 1:
		#  This version is for meridien refinement.  There are simply three full sets of params.
		outdir      = args[0]
		bp = 'badparams'
		#outgrouparms= args[1]
		radius = options.ou
		thresherr = options.thresherr
		sym = int(options.sym[1:])
		qsym = 360.0/sym
		#params = [[None for i in xrange(3)] for j in xrange(4)]
		ll = 3 # this is hardwired as we have three groups.  however, I would like to keep the code general.
		for jj in xrange(1):
			params = [None for ii in xrange(ll)]
			#newbad = map(int, read_text_file(options.params+"%01d.txt"%jj) )
			#nn = len(newbad)
			for ii in xrange(ll):
				params[ii] = read_text_row(os.path.join(outdir,"params%01d.txt"%(ii)))
				#assert(nn == len(params[ii]) )
			nn = len(params[0])
			newbad = range(nn)
			#  Compute average projection params and pixel errors
			avgtrans = [None]*nn
			pixer   = [0.0]*nn
			for j in xrange(nn):
				nas = [0.0,0.0,0.0]
				if( sym == 1):
					#pixer[q[0][0]][j] = max_3D_pixel_error(fifi[0], fifi[1], r=radius)
					for i in xrange(ll):
						n1,n2,n3 = getfvec(params[i][j][0],params[i][j][1])
						nas[0] += n1
						nas[1] += n2
						nas[2] += n3
				else:
					m1,m2,m3 = getfvec(params[0][j][0],params[0][j][1])
					nas[0] = m1
					nas[1] = m2
					nas[2] = m3
					for i in xrange(1,ll):
						qnom = -1.e10
						for j in xrange(-1,2,1):
							t1,t2,t3 = getfvec(params[i][j][0]+j*qsym,params[i][j][1])
							nom = t1*m1 + t2*m2 + t3*m3
							if(nom > qnom):
								qnom = nom
								n1=t1
								n2=t2
								n3=t3
							#if(k == 2):
							#	print '  t1,t2,t3 ',fifi[i][0]+j*qsym,fifi[i][1],t1,t2,t3,nom,qnom
						nas[0] += n1
						nas[1] += n2
						nas[2] += n3
						print qnom, n1,n2,n3,nas
					#  To get the correct pixer phi angle has to be taken from the above!!

				nom = sqrt(nas[0]**2 + nas[1]**2 + nas[2]**2)

				if(nom < 1.e-6):
					nphi   = 0.0
					ntheta = 0.0
				else:
					ntheta = degrees(acos(nas[2]/nom))%360.0
					if(sym>1 and ntheta>90.0):  nphi   = (degrees(atan2( nas[1], nas[0] ))-180.0)%qsym + 180.0
					else:                       nphi   =  degrees(atan2( nas[1], nas[0] ))%qsym

				#print   "FIFI     %4d     %7.2f     %7.2f    %7.2f    %7.2f     %7.2f     %7.2f    %7.2f    %7.2f"%(k,fifi[0][0],fifi[0][1],fifi[1][0],fifi[1][1],fifi[2][0],fifi[2][1],nphi,ntheta)
				twod = average2dtransform([params[ii][j][2:] for ii in xrange(ll)])
				avgtrans[j] = [nphi, ntheta, twod[0], twod[1], twod[2]]

			write_text_row(avgtrans,  os.path.join(outdir,"avgtrans.txt"))
			perr =  errors_per_image(params, avgtrans, thresherr, radius )
			rescued = []
			rejects = []
			for j in xrange(nn):
				#print  chr(65+jj),j,perr[j][0],[[params[m][j][i] for i in xrange(2)] for m in xrange(ll)]
				if( perr[j][0] <= thresherr ):  rescued.append([newbad[j],j])
				else:							rejects.append([newbad[j],j])
			if( len(rescued) == 0 ):    write_text_row([-1,-1],  os.path.join(outdir,"rescued%01d.txt"%jj))
			else:                  		write_text_row(rescued,  os.path.join(outdir,"rescued%01d.txt"%jj))
			#  We also have to write params.
			if( len(rescued) != 0 ):
				for ii in xrange(ll):  write_text_row([params[ii][rescued[k][1]] for k in xrange(len(rescued))]  ,  os.path.join(outdir,"params-rescued%01d%01d.txt"%(jj,ii)))
			if( len(rejects) == 0 ):    write_text_row([-1,-1],  os.path.join(outdir,'rejects'+"%01d.txt"%jj))
			else:                  		write_text_row(rejects,  os.path.join(outdir,'rejects'+"%01d.txt"%jj))
			if( len(rejects) != 0 ):
				for ii in xrange(ll):  write_text_row([params[ii][rejects[k][1]] for k in xrange(len(rejects))]  ,  os.path.join(outdir,"params-rejects%01d%01d.txt"%(jj,ii)))

			hi = hist_list([perr[j][0]  for j in xrange(nn)  ], 16)
			print  "Pixel errors for BAD GROUP  ",chr(65+jj)
			for ii in xrange(len(hi[0])):
				print  "%4d   %12.3f    %12.0f "%(ii,hi[0][ii],hi[1][ii])


	else:
		print "Usage: "
		print """
		Phase 1:   sxconsistency.py  --phase=1  bdb:data  outdir
			output files are:
			  in directory outdir: lili0.txt to lili5.txt contain indices of images in resampled six groups
			  bdb:dataX0 to bdb:dataX5 are metafiles containing six resampled groups of images derived from bdb:data


		Phase 2    sxconsistency.py  --phase=2  outdir  --ndigits=1  --params=paramsb howmanythesame.txt
			  outdir - directory containing files lili0.txt to lili5.txt produced in phase 1
			  --params=master/main/params - Root of of six parameter file names with refinement results, the actual names should be
											  master/main/params0.txt  to master/main/params5.txt
			output files:
				howmanythesame.txt - contains one number, a ratio of number of images  that did not change orientations
				                                            to the total number of images


		Phase 3:   sxconsistency.py --phase=3 outdir --ou=133 --thresherr=3.0 --params=paramsa outgrouparms

			input files:
			  outdir - directory containing files lili0.txt to lili5.txt produced in phase 1
			  --params=master/main/params - Root of of six parameter file names with refinement results, the actual names should be
											  master/main/params0.txt  to master/main/params5.txt
			output files:
			  outgrouparms*.txt - Root of of six parameters files with average 3D orientation parameters computed from six runs, can be imported into bdb:data
                   The next two files contain the original image numbers refering to the top bdb.
			  outdir/newgood.txt    - Text file with indices of images whose orientation parameters arrors are below specified thresherr (to be kept)
			  outdir/bad.txt        - Text file with indices of images whose orientation parameters arrors are above specified thresherr (to be rejected)

		Phase 4:   sxconsistency.py --phase=4 outdir --ou=133 --thresherr=3.0 --params=master/main001/newbad

			input files:
			  outdir - directory containing files badparamsij.txt i= 0,3, j=0,2, which resulted from three alignments of four badchunk images
			  --params= - Root of four files with original indices of bad images, they will be read and a subset corresponding to rescued ones will be outputed
			output files:
			  rescued*.txt - four files with indices of accepted images from badchunk.
			  rejects*.txt - four files with indices of rejected images from badchunk.
			  					There are two columns: [number in original bdb:chunk, number in bdb:newbad] 
		"""
		print "Please run '" + progname + " -h' for detailed options"

	global_def.BATCH = False


main()
