#!/usr/bin/env python

#
# Author: Matthew Baker, 10/2005, modified 02/2006 by MFS  
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

import os
import sys
import string
import commands
import math
import EMAN2
#import Numeric
from math import *
from sys import argv


def tomoccf(targetMRC,probeMRC):
    ccf=targetMRC.calc_ccf(probeMRC)
    ccf.toCorner()
    return (ccf)

def updateCCF(bestCCF,bestALT,bestAZ,bestPHI,altrot,azrot,phirot,currentCCF,box,scalar):
    bestValue=-1.e10
    xbest=-10
    ybest=-10
    zbest=-10
    z=(box/2)-3
    while z < (box/2)+4:
        y=(box/2)-3
        while y < (box/2)+4:
            x=(box/2)-3
            while x < (box/2)+4:
                currentValue=currentCCF.get_value_at(x,y,z)/scalar
                if currentValue > bestValue:
                    bestValue = currentValue
                    xbest=x-(box/2)
                    ybest=y-(box/2)
                    zbest=z-(box/2)
                #print x,y,z, currentValue
                x=x+1
            y=y+1
        z=z+1
    inlist=0
    while inlist < 10:
        if  bestValue > bestCCF.get_value_at(inlist,1,1):
            swlist=9
            while swlist >inlist:
                #print swlist
                bestCCF.set_value_at(swlist,1,1,bestCCF.get_value_at(swlist-1,1,1))
                bestALT.set_value_at(swlist,1,1,bestALT.get_value_at(swlist-1,1,1))
                bestAZ.set_value_at(swlist,1,1,bestAZ.get_value_at(swlist-1,1,1))
                bestPHI.set_value_at(swlist,1,1,bestPHI.get_value_at(swlist-1,1,1))
                bestX.set_value_at(swlist,1,1,bestX.get_value_at(swlist-1,1,1))
                bestY.set_value_at(swlist,1,1,bestY.get_value_at(swlist-1,1,1))
                bestZ.set_value_at(swlist,1,1,bestZ.get_value_at(swlist-1,1,1))
                swlist=swlist-1
            bestCCF.set_value_at(inlist,1,1,bestValue)
            bestALT.set_value_at(inlist,1,1,altrot*180/3.14159)
            bestAZ.set_value_at(inlist,1,1,azrot*180/3.14159)
            bestPHI.set_value_at(inlist,1,1,phirot*180/3.14159)
            bestX.set_value_at(inlist,1,1,xbest)
            bestY.set_value_at(inlist,1,1,ybest)
            bestZ.set_value_at(inlist,1,1,zbest)
            break
        inlist=inlist+1
    #print "one"
    bestCCF.update()
    #print "two"
    bestALT.update()
    bestAZ.update()
    bestPHI.update()
    return(bestCCF)
    
def ccfFFT(currentCCF, thresh, box):
    tempCCF = currentCCF.copy()
    tempCOMPLEX=tempCCF.doFFT()
    tempCOMPLEX.gimmeFFT() # supposed to give me "ownership". whatever...
    tempCOMPLEX.ri2ap()
    tempCOMPLEX.realFilter(2,thresh)
    mapSum=tempCOMPLEX.Mean()*box*box*box/2.
    #print mapSum
    return(mapSum)

def peakSearch(bestCCF,Max_location,width,box):
    xmin=Max_location[0]-width
    xmax=Max_location[0]+width
    ymin=Max_location[1]-width
    ymax=Max_location[1]+width
    zmin=Max_location[2]-width
    zmax=Max_location[2]+width
    for x in (range(xmin,xmax)):
        for y in (range(ymin,ymax)):
            for z in (range(zmin,zmax)):
                print x,y,z
                bestCCF.set_value_at(x,y,z,-10000)
    bestCCF.update()
    return(bestCCF)        

if (len(argv)<3) :
    print "Usage:\ntomohunter.py <target mrc file> <probe mrc file> \n"
    print "Options: da=<search step (default=30)>, thresh=<threshold (default=0), maxpeaks=<number of results returned in log file (default=20)>, width=<peak width (default=2)>"
    sys.exit()
    
target=argv[1]
probe=argv[2]
thresh=0
da=60
maxPeaks=10
width=2

for a in argv[3:] :
    s=a.split('=')
    if (s[0]=='dal'):
        dal=float(s[1])
    elif (s[0]=='daz'):
        daz=float(s[1])
    elif (s[0]=='dap'):
        dap=float(s[1])
    elif (s[0]=='thresh'):
        thresh=float(s[1])
    elif (s[0]=='ral'):
        ral=float(s[1])
#    elif (s[0]=='raz'):
#        raz=float(s[1])
    elif (s[0]=='rap'):
        rap=float(s[1])
    else:
        print("Unknown argument "+a)
        exit(1)

print target, probe

targetMRC=EMAN2.EMData()
targetMRC.read_image(argv[1],-1)
targetMean=targetMRC.get_attr('mean')
targetSigma=targetMRC.get_attr('sigma')
print "Target Information"
print "   mean:       %f"%(targetMean)
print "   sigma:      %f"%(targetSigma)

target_xsize=targetMRC.get_xsize()
target_ysize=targetMRC.get_ysize()
target_zsize=targetMRC.get_zsize()
if (target_xsize!=target_ysize!=target_zsize) or (target_xsize%2==1):
    print "The density map must be even and cubic. Terminating."
    sys.exit()
box=target_xsize

probeMRC=EMAN2.EMData()
probeMRC.read_image(argv[2],-1)
probeMean=probeMRC.get_attr('mean')
probeSigma=probeMRC.get_attr('sigma')
print "Probe Information"
print "   mean:       %f"%(probeMean)
print "   sigma:      %f"%(probeSigma)


bestCCF=EMAN2.EMData()
bestCCF.set_size(box,box,box)
bestCCF.to_zero()

bestAZ=EMAN2.EMData()
bestAZ.set_size(box,box,box)
bestAZ.to_zero()

bestALT=EMAN2.EMData()
bestALT.set_size(box,box,box)
bestALT.to_zero()

bestPHI=EMAN2.EMData()
bestPHI.set_size(box,box,box)
bestPHI.to_zero()

bestX=EMAN2.EMData()
bestX.set_size(box,box,box)
bestX.to_zero()

bestY=EMAN2.EMData()
bestY.set_size(box,box,box)
bestY.to_zero()

bestZ=EMAN2.EMData()
bestZ.set_size(box,box,box)
bestZ.to_zero()


altarray=[]
alt=0
while alt <= ral:
    altarray.append(alt*3.14159/180)
    alt=alt+dal

azarray=[]
az=-180
while az < 180:
    azarray.append(az*3.14159/180)
    az=az+daz

#phiarray=[]
#phi=-180
#while phi < 180:
#    phiarray.append(phi*3.14159/180)
#    phi=phi+da
# I had to change this because the phi range varies depending on the azimuth (it is nearly the negative of it, within
# a range, when alt is near zero, which it will be at this point)
rarad=rap*3.14159/180.
darad=dap*3.14159/180.
maxfrac=0.
for altrot in altarray:
    print "Trying rotation %f"%(altrot*180./3.14159)
    if (altrot==0.):
        azrot=0.
        phirot=-azrot-rarad
        minnum=10000000
        maxnum=0
        #print "Trying rotation %f %f"%(altrot, azrot)
        while phirot <= -azrot+rarad:
            dMRC=EMAN2.EMData()
            dMRC = probeMRC.copy()
            dMRC.rotate(azrot, altrot, phirot)
            #print "Trying rotation %f %f %f"%(altrot, azrot, phirot)
            currentCCF=tomoccf(targetMRC,dMRC)
            scalar=ccfFFT(currentCCF,thresh,box)
            if scalar>maxnum:
                maxnum=int(scalar)
            if scalar<minnum:
                minnum=int(scalar)
            scalar=scalar/(box*box*box/2.)
            #scalar=1
            #scaledCCF=currentCCF/scalar
            #print "three"
            bestCCF=updateCCF(bestCCF,bestALT,bestAZ,bestPHI,altrot,azrot,phirot,currentCCF,box,scalar)
            phirot=phirot+darad
        #print minnum,maxnum, float(maxnum)/float(minnum)
    else:
        for azrot in azarray:
            phirot=-azrot-rarad
            #print "Trying rotation %f %f"%(altrot, azrot)
            while phirot <= -azrot+rarad:
                dMRC=EMAN2.EMData()
                dMRC = probeMRC.copy()
                dMRC.rotate(azrot, altrot, phirot)
                #print "Trying rotation %f %f %f"%(altrot, azrot, phirot)
                currentCCF=tomoccf(targetMRC,dMRC)
                scalar=ccfFFT(currentCCF,thresh,box)
                if scalar>maxnum:
                    maxnum=int(scalar)
                if scalar<minnum:
                    minnum=int(scalar)
                scalar=scalar/(box*box*box/2.)
                #scalar=1
                #scaledCCF=currentCCF/scalar
                #print "three"
                bestCCF=updateCCF(bestCCF,bestALT,bestAZ,bestPHI,altrot,azrot,phirot,currentCCF,box,scalar)
                phirot=phirot+darad
            #print minnum,maxnum, float(maxnum)/float(minnum)
print minnum,maxnum, float(maxnum)/float(minnum)
#outCCF="rl-%s"%(argv[1])
#outalt="alt-%s"%(argv[1])
#outaz="az-%s"%(argv[1])
#outphi="phi-%s"%(argv[1])

#bestCCF.write_image(outCCF)
#bestALT.write_image(outalt)
#bestAZ.write_image(outaz)
#bestPHI.write_image(outphi)

out=open("log-s3-%s%s.txt"%(argv[1],argv[2]),"w")
peak=0
while peak < 10:
    #Max_location=bestCCF.MinLoc()
    ALT=str(bestALT.get_value_at(peak,1,1))
    AZ=str(bestAZ.get_value_at(peak,1,1))
    PHI=str(bestPHI.get_value_at(peak,1,1))
    COEFF=str(bestCCF.get_value_at(peak,1,1))
    LOC=str( ( (bestX.get_value_at(peak,1,1)),(bestY.get_value_at(peak,1,1)),(bestZ.get_value_at(peak,1,1) ) ) )
    line="Peak %d rot=( %s, %s, %s ) trans= %s coeff= %s\n"%(peak,ALT,AZ,PHI,LOC,COEFF)
    out.write(line)
    #bestCCF=peakSearch(bestCCF,Max_location, width, box)
    peak=peak+1
out.close()
