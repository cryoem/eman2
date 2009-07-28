#!/usr/bin/env python
###	test.py	Muthu Alagappan June 2009


######################### 

import sys
from random import *
import random
from EMAN2 import *
from math import *
from sys import argv
from numpy import *
import time
#from pylab import scatter, subplot, show, figure, colorbar
from emapplication import EMStandAloneApplication
from emplot3d import EMPlot3DModule

def getRots():
	all_rots=[]
	s = Symmetries.get("c1")
	all_rots = s.gen_orientations("rand", {"n":transNum, "phitoo":1})
	for i in all_rots:
		i.set_trans(random.randint(-xTrans, xTrans),random.randint(-yTrans, yTrans),random.randint(-zTrans, zTrans))
	
	all_rots[0].set_params({"type":"eman", "alt":0, "az":0, "phi":0, 'ty':0.0, "tz":0.0, "tx":0.0})
	all_rots[0].set_trans(0,0,0)
	return all_rots


if (len(argv) < 3) :
	print "Usage: \n e2fhstat.py <mrc-file 1 > <probe 2 > <# of transformations 3 > <isosurface threshold 4 > <discard volume 5 > <binary threshold 6 > <apix 7 > "
	sys.exit()

target = EMData(argv[1])

try: 
	apix_x = float(argv[7])
	apix_y, apix_z  = apix_x, apix_x
except:
	apix_x = target["apix_x"]
	apix_y = target["apix_y"]
	apix_z = target["apix_z"]
xMax = target.get_xsize()
yMax = target.get_ysize()
zMax = target.get_zsize()
try:
	transNum = int(argv[3])
except:
	transNum = 25
try:
	s2iso = float(argv[4])
except: 
	s2iso = 5.0
try:
	dVol = float(argv[5])
except: 
	dVol = .05
################################################################


######################### read in pdb file, store coordinates
try:
	pdb_file = open(argv[2], 'r')
	pdb_lines = pdb_file.readlines()
	pdb_file.close()
except IOError:
	print "The file \"" + str(argv[2]) + "\" does not exist"
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

		
#printValues()
#getPixelForAtom()

##############################################################

######################## write into pdb file

f_pdb = open ("%s"%argv[2], "r")
entries = f_pdb.readlines()
lastEntry = list(entries)
f_pdb.close()

for i in range (0, len(entries)):
	tempEntry = entries[i]
	entries[i] = tempEntry[0:61]
	tempLast = lastEntry[i]
	lastEntry[i] = tempLast[66:80]

f = open ("%s"%argv[2], "w")
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


rotList = getRots()

score1 =[]  #.005
s1_sdv,s2_sdv,s3_sdv=0.0,0.0,0.0
s1_mean,s2_mean,s3_mean=0.0,0.0,0.0
score2 =[]
score3 =[]

a = PDBReader()
a.read_from_pdb("%s"%argv[2])


for i in range(0,len(rotList)):

	t = Transform(rotList[i].get_params("eman"))
	b = a.copy()
	b.right_transform(t)  
	points = b.get_points()

	p = 0
	t_pixelValues = []

	for m in range(0, atomCount):  #.025 seconds per transform
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
	#score1.append(sumValues) 
	####################

	###### score 2  #.001 seconds per transform
	includeValue = 0.0
	cutOff = (1/(s2iso*10.))
	isoValue = s2iso
	for x in t_pixelValues:
		if (x>=isoValue): includeValue = includeValue + 1.0
	percentAbove = (includeValue/atomCount)
	#score2.append(perecentAbove)
	#######


	###### score 3   #.1 second per transform	
	pA = b.makePointArray(b)	
	pA.save_to_pdb("thisFile100.txt")
	probe_MRC = pA.pdb2mrc_by_summation(xMax,apix_x,4.0)
	probe_MRC.process_inplace("normalize.toimage",{"noisy":target,"keepzero":1}) #.01 seconds

	probe_MRC.process_inplace("threshold.binary",{"value":1.0})#.02 seconds
	target_binary = target.process("threshold.binary",{"value":1.0})
	c = (target_binary - probe_MRC)
	c.process_inplace("threshold.binary",{"value":0.1}) #.015 seconds

	remainder_volume = c["mean"]*(xMax*yMax*zMax) #.018 seconds
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

print len(calc1)

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
	#print str(s1[i]) + '   ' + str(s2[i]) + '   ' + str(s3[i])	
	#print int(s1[i]+.5), int(s2[i]+.5), int(s3[i]+.5)

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

#for i in rotList:
	#print i.get_params("eman")["alt"]
	#print i.get_params("eman")["az"]
	#print '\n'

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


########################################################
'''
w=0
for i in range(0,(len(score1)/100)):
	me1 = mean(score1[0:(w+100)])
	st1 = std(score1[0:(w+100)])
	me2 = mean(score2[0:(w+100)])
	st2 = std(score2[0:(w+100)])
	me3 = mean(score3[0:(w+100)])
	st3 = std(score3[0:(w+100)])
	#print str(w+100) + ":   " + str(me1) + '   ' + str(me2) + '   ' + str(me3) + '     ' + str(st1) + '   ' + str(st2) + '   ' + str(st3) 
	#print '\n'
	w = w+100
'''

#####################################################

########  Graphing 1
'''
g1 = figure(1)
g1 = subplot(311)
g1.scatter(mpl1,mpl2, marker ="o", color ="green")
g1.scatter(mpl1[:1],mpl2[:1],color="red")
g1.set_xlabel("Real space correlation")
g1.set_ylabel("Atom Inclusion Percentage")

g1 = subplot(312)
g1.scatter(mpl3,mpl2, marker = "o", color="green")
g1.scatter(mpl3[:1],mpl2[:1],color="red")
g1.set_xlabel("Volume Overlap")
g1.set_ylabel("Atom Inclusion Percentage")
g1=subplot(313)

g1.scatter(mpl3,mpl1, marker = "o", color="green")
g1.scatter(mpl3[:1],s1[:1],color="red")
g1.set_xlabel("Volume Overlap")
g1.set_ylabel("Real space correlation")
#ga.axis([-3.5, 3.5, -3.5, 3.5])
#gb = subplot(111)
#ga.set_zlabel("Score 3")
#print s1[0], s2[0], rot_extra[0][4], rot_extra[0][5]
#show()

####################################################

######### Graphing 2 - individual lines

g2 = figure(2)
g2 = subplot(311)

g2.scatter(mpl1,blank, marker ="o", color ="green")
g2.scatter(mpl1[:1],blank[:1],color="red")
g2.set_xlabel("Real space correlation")
g2.set_ylabel(" ")

g2 = subplot(312)
g2.scatter(mpl2,blank, marker = "o", color="green")
g2.scatter(mpl2[:1],blank[:1],color="red")
g2.set_xlabel("Atom Inclusion Percentage")
g2.set_ylabel(" ")

g2=subplot(313)
g2.scatter(mpl3,blank, marker = "o", color="green")
g2.scatter(mpl3[:1],blank[:1],color="red")
g2.set_xlabel("Volume Overlap")
g2.set_ylabel(" ")
#show()

###################################################

########## Graphing 3

g3 = figure(3)
g3 = subplot(111)

# finds the best transform
#s4_s=0
#val = 0
#for i in range(0, len(s4)):
#	if (s4[i]>s4_s): 
#		s4_s = s4[i]
#		val = i
#valT = Transform(rotList[val].get_params("eman"))
#tArray = a.copy()
#tArray.right_transform(valT)  
#tArray.save_to_pdb("best.pdb")



g3.scatter(mpl4,blank, marker ="o", color ="green")
g3.scatter(mpl4[:1],blank[:1],color="red")
#g1.scatter(s4[val:(val+1)],blank[:1],color="blue")
g3.set_xlabel("Total z score")
g3.set_ylabel(" ")
#show()

##################################################

######## Graphing 4 - color

g4 = figure(4)
g4 = subplot(111)
g4.scatter(mpl1,mpl2, marker ="o", c=mpl3, faceted="false")
g4.scatter(mpl1[:1],mpl2[:1], marker=(3,2,0))
g4.set_xlabel("Real space correlation")
g4.set_ylabel("Atom Inclusion Percentage")
#show()

##################################################

######## Graphing 5 -  2d by distance

g5 = figure(5)
g5 = subplot(111)

g5.scatter(mpl1,mpl2, marker ="o", c=mplDist, faceted="false")
g5.scatter(mpl1[:1],mpl2[:1], marker=(3,2,0))
g5.set_xlabel("Real space correlation")
g5.set_ylabel("Atom Inclusion Percentage")
#show()

##################################################

######## Graphimg 6 - histogram

g6 = figure(6)
g6 = subplot(311)
g6.hist(mpl1,50, histtype='bar')
g6.scatter(mpl1[:1],blank[1:],color="red")
g6.set_xlabel("Real space correlation")

g6 = subplot(312)
g6.hist(mpl2,50, histtype='bar')
g6.scatter(mpl2[:1], blank[1:],color="red")
g6.set_xlabel("Atom Inclusion Percentage")

g6 = subplot(313)
g6.hist(mpl3,50, histtype='bar')
g6.scatter(mpl3[:1], blank[1:],color="red")
g6.set_xlabel("Volume inclusion")
#show()
#################################################

######## Graphing 7 -  2d distance by z score

g5 = figure(7)
g5 = subplot(111)

g5.scatter(mpl4,mplDist, marker ="o", c="green")
g5.scatter(mpl4[:1],mplDist[:1], c="red")
g5.set_xlabel("Cumulative Z Score")
g5.set_ylabel("Distance")
show()

##################################################
'''

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

em_app = EMStandAloneApplication()
window = EMPlot3DModule(application=em_app)
window.set_Vals(vals)
window.set_Rotations(rotList)
b = a.copy()
window.set_Probe(b)
#window.set_data(data,"Z-Scores Data")
#window.set_data(initPoint, "Original Probe", shape = "Cube")
window.set_data(rotData,"Rotation Angles with Z Score")         # origin is white
#window.set_data(rotData2,"Rotation Angles with Distance", shape = "Octahedron")       #origin is black
em_app.show()
em_app.execute()

'''
em2 = EMStandAloneApplication()
win2 = EMPDBViewer()
em2.show()
em2.execute()
'''
