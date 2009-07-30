#!/usr/bin/env python

# Author: Muthu Alagappan, 07/22/09
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#


import sys
from random import *
import random
from EMAN2 import *
from math import *
from numpy import *
import time
#from pylab import scatter, subplot, show, figure, colorbar
from emapplication import EMStandAloneApplication
from emplot3d import EMPlot3DModule


class E2FoldHunterStat:
	
	def __init__(self): pass

	def gen_data(self, mrc_file_name, pdb_file_name, trans_Num, iso_Thresh):

		target = EMData(mrc_file_name)

		apix_x = target["apix_x"]
		apix_y = target["apix_y"]
		apix_z = target["apix_z"]

		xMax = target.get_xsize()
		yMax = target.get_ysize()
		zMax = target.get_zsize()

		transNum = trans_Num

		s2iso = iso_Thresh
		dVol = .05

		################################################################


		######################### read in pdb file, store coordinates
		try:
			pdb_file = open(pdb_file_name, 'r')
			pdb_lines = pdb_file.readlines()
			pdb_file.close()
		except IOError:
			print "The file \"" + str(pdb_file_name) + "\" does not exist"
			sys.exit()

		atomNumber=[]
		atomCount = 0
		pixelValues =[]

		for item in pdb_lines:
			itemS = str(item[0:6].strip())
			if (itemS == "ATOM"):
				tempX=int((((float(item[30:38].strip()))/apix_x)+(xMax*.5))+.5)
				tempY=int((((float(item[38:46].strip()))/apix_y)+(yMax*.5))+.5)
				tempZ=int((((float(item[46:54].strip()))/apix_z)+(zMax*.5))+.5)	
				tempValue = target.get(tempX, tempY, tempZ)
				pixelValues.append(tempValue)
				atomCount=atomCount+1

		##############################################################

		######################## write into pdb file

		f_pdb = open (pdb_file_name, "r")
		entries = f_pdb.readlines()
		lastEntry = list(entries)
		f_pdb.close()

		for i in range (0, len(entries)):
			tempEntry = entries[i]
			entries[i] = tempEntry[0:61]
			tempLast = lastEntry[i]
			lastEntry[i] = tempLast[66:80]

		f = open (pdb_file_name, "w")
		for i in range (0, len(entries)):
			f.write(entries[i])
			if (i < len(pixelValues)):
				tempPixel = str(pixelValues[i])
				tempPixel = tempPixel[0:5]
				f.write(tempPixel + lastEntry[i]+ '\n')
			else: f.write('\n')

		f.close()

		################################################################

		xTrans = (.25*xMax)
		yTrans = (.25*yMax)
		zTrans = (.25*zMax)

		all_rots=[]

		s = Symmetries.get("c1")
		all_rots = s.gen_orientations("rand", {"n":transNum, "phitoo":1})
		for i in all_rots:
			i.set_trans(random.randint(-xTrans, xTrans),random.randint(-yTrans, yTrans),random.randint(-zTrans, zTrans))
	
		all_rots[0].set_params({"type":"eman", "alt":0, "az":0, "phi":0, 'ty':0.0, "tz":0.0, "tx":0.0})
		all_rots[0].set_trans(0,0,0)

		rotList = all_rots

		score1 =[]
		s1_sdv,s2_sdv,s3_sdv=0.0,0.0,0.0
		s1_mean,s2_mean,s3_mean=0.0,0.0,0.0
		score2 =[]
		score3 =[]

		a = PDBReader()
		a.read_from_pdb(pdb_file_name)


		for i in range(0,len(rotList)):

			t = Transform(rotList[i].get_params("eman"))
			b = a.copy()
			b.right_transform(t)  
			points = b.get_points()

			p = 0
			t_pixelValues = []

			for m in range(0, atomCount):
				tempX_t=int(((float(points[(p)])/apix_x)+(xMax*.5))+.5)
				tempY_t=int(((float(points[(p+1)])/apix_y)+(yMax*.5))+.5)
				tempZ_t=int(((float(points[(p+2)])/apix_z)+(zMax*.5))+.5)
				if (tempX_t>xMax): tempX_t = xMax	
				if (tempY_t>yMax): tempY_t = yMax
				if (tempZ_t>zMax): tempZ_t = zMax
				if (tempX_t<0): tempX_t = 0
				if (tempY_t<0): tempY_t = 0
				if (tempZ_t<0): tempZ_t = 0
				tempValue2 = target.get(tempX_t, tempY_t, tempZ_t)
				t_pixelValues.append(tempValue2)	
				p = p + 3

	
			###### score 1   #.001 seconds per transform
			sumValues = 0.0
			for x in t_pixelValues:
				sumValues = sumValues + x
			sumValues = (sumValues/atomCount) 
			####################

			###### score 2  #.001 seconds per transform
			includeValue = 0.0
			cutOff = (1/(s2iso*10.))
			isoValue = s2iso
			for x in t_pixelValues:
				if (x>=isoValue): includeValue = includeValue + 1.0
			percentAbove = (includeValue/atomCount)
			#######


			###### score 3 	
			pA = b.makePointArray(b)	
			pA.save_to_pdb("thisFile100.txt")
			probe_MRC = pA.pdb2mrc_by_summation(xMax,apix_x,4.0)
			probe_MRC.process_inplace("normalize.toimage",{"noisy":target,"keepzero":1}) 

			probe_MRC.process_inplace("threshold.binary",{"value":1.0})
			target_binary = target.process("threshold.binary",{"value":1.0})
			c = (target_binary - probe_MRC)
			c.process_inplace("threshold.binary",{"value":0.1}) 

			remainder_volume = c["mean"]*(xMax*yMax*zMax)
			MRC_volume = target_binary["mean"]*(xMax*yMax*zMax)

			excludePercent = (float(remainder_volume)/MRC_volume)
			volumePercent = 1-excludePercent

			#score3.append(volumePercent) #raw scores
			flag = 0
			if flag:
				if (volumePercent >= dVol and percentAbove>=cutOff): 
						score1.append(sumValues)
						score2.append(percentAbove)
						score3.append(volumePercent)
			
				if(percentAbove<cutOff or volumePercent <dVol):
					score1.append("none")
					score2.append("none")
					score3.append("none")		
			else: 
				score1.append(sumValues)
				score2.append(percentAbove)
				score3.append(volumePercent)
			
				####################
		'''
		print '\n'
		print score1
		print '\n'
		print score2
		print '\n'
		print score3
		print '\n'
		'''
		#########  Calculate standard deviations 
		calc1 = []  #calc holds only the good transform raw scores, only these are used to calculate mean and std
		calc2 = []
		calc3 = []
		for l in range(0,len(score1)):
			if (score1[l]!="none"): 
				calc1.append(score1[l])
				calc2.append(score2[l])
				calc3.append(score3[l])

		s1_mean = mean(calc1)
		s1_std = std(calc1)
		s2_mean = mean(calc2)
		s2_std = std(calc2)
		s3_mean = mean(calc3)
		s3_std = std(calc3)
		s1 =[]             #score1,2,3 holds raw scores of all transforms ("none" for bad transforms)-  s1,2,3 holds std scores of all transforms (0.0 for bad transforms)
		s2 =[]
		s3 =[]

		cCount = 0
		for i in range(0,len(score1)):
			if (score1[i]=="none"): 
				s1.append(0.0)
				s2.append(0.0)
				s3.append(0.0)
			else:
				s1.append(((calc1[cCount]-s1_mean)/s1_std))
				s2.append(((calc2[cCount]-s2_mean)/s2_std))
				s3.append(((calc3[cCount]-s3_mean)/s3_std))
				cCount=cCount+1
				#s3.append(calc3[cCount])


		print "number of transformations preformed: " + str(len(calc1))
		print "score 1: " + str(s1[0])
		print "score 2: " + str(s2[0])
		print "score 3: " + str(s3[0])

		##### Create Dictionary
		vals = {}
		tx = []
		ty = []
		tz = []
		alt = []
		az = []
		phi = []
		tNum = []
		q_num = 0
		for q in rotList:
			tempo = q.get_params("eman")
			tx.append(tempo['tx'])
			ty.append(tempo['ty'])
			tz.append(tempo['tz'])
			alt.append(tempo['alt'])
			az.append(tempo['az'])
			phi.append(tempo['phi'])
			tNum.append(q_num)
			q_num = q_num+1
		vals["tx"] = tx
		vals["ty"] = ty
		vals["tz"] = tz
		vals["alt"] = alt
		vals["az"] = az
		vals["phi"] = phi
		vals["score_1"] = s1
		vals["score_2"] = s2
		vals["score_3"] = s3
		vals["tNum"] = tNum

		s4 =[]
		for i in range(0,len(s1)):
			te = s1[i] + s2[i] + s3[i]
			s4.append(te)

		dist =[]
		for i in rotList:
			tempor = i.get_params("eman")
			d = (tempor['tx'])**2 + (tempor['ty'])**2 + (tempor['tz'])**2
			d = (d)**.5
			dist.append(d)

		blank = []
		for i in range(0,len(s1)):
			blank.append(0)


		mpl1 =[] #holds standard deviation scores for only good transforms
		mpl2 =[]
		mpl3 =[]
		mpl4 = []
		mplDist=[]
		for t in range (0, len(s1)):
			if (score1[t] != "none"):
				mpl1.append(s1[t])
				mpl2.append(s2[t])
				mpl3.append(s3[t])
				mpl4.append(s1[t]+s2[t]+s3[t])
				d = (vals['tx'][t])**2 + (vals['ty'][t])**2 + (vals['tz'][t])**2
				d = (d)**.5
				mplDist.append(d)

		
		############ emplot 3d
		data = []
		data.append(s1)
		data.append(s2)
		data.append(s3)

		initPoint = []
		initPoint.append(s1[:1])
		initPoint.append(s2[:1])
		initPoint.append(s3[:1])

		rotData = []
		rotData.append(alt)
		rotData.append(az)
		rotData.append(phi)
		rotData.append(s4)
		rotData.append(dist)

		rotData2 = []
		rotData2.append(alt)
		rotData2.append(az)
		rotData2.append(phi)
		rotData2.append(dist)
		b = a.copy()
		return vals, rotList, a, data, initPoint


if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window2 = E2FoldHunterStat()
	em_app.show()
	em_app.execute()

		
		
