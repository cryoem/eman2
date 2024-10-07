#!/usr/bin/env python

import os
import sys
import numpy as np
#from EMAN2 import EMData, EMArgumentParser, EMANVERSION, Transform
from EMAN2 import *


def project_map_and_compute_bispectrum(map_FN, NumProj, output_FN):
    """
    Projects a 3D map into a set of 2D projections, computes the bispectrum for each projection,
    and writes the average bispectrum to a file.

    Parameters:
    - map_file: path to the 3D map (MRC or HDF file)
    - NumProj: number of projections to generate
    - output_FN: file where the average bispectrum will be saved (MRC or HDF format)
    """
    # Load the 3D map
    volumeData = EMData(map_FN)

    # Array to store bispectra
    if 0:
        bispectra = []

    # Generate projections evenly distributed in angular space
    sym = Symmetries.get("c1"); # Create Projections.  Can skip this if they have already been created
    directionsEven = sym.gen_orientations("eman",{"n":NumProj,"inc_mirror":True});#  Change


    #nimg=EMUtil.get_image_count(args[i])

    rfp=4
    size=16;# was 16
    #kVal=11

    sig =1
    #BispecAvgr=Averagers.get("mean",{"sigma":sig})

    BispecAvgr=Averagers.get("mean")

    #prjArray = []
    print(len(directionsEven))
    for idirection,t in enumerate(directionsEven):

        if (idirection%100)==0: print(idirection)
        prj = volumeData.project("standard",t);
        #prj.process_inplace("normalize.edgemean");
        bispecOut =  prj.process("math.bispectrum.slice",{"rfp":rfp,"size":size});

        BispecAvgr.add_image(bispecOut)
        #if idirection==0:
        prj.write_image('prj.hdf',idirection)
        bispecOut.write_image('bispecOut.hdf',idirection)

    BispecAvg = BispecAvgr.finish()
    BispecAvg.write_image(output_FN)

    print(np.shape(bispecOut.numpy()))
    print(f"Average bispectrum (will have a symmetry) written to {output_FN}")

# Usage python project_bispectrum.py --input map.hdf --output avg_bispectrum.hdf --numproj 20





if __name__ == "__main__":
#def main():
    # Setup the command-line argument parser
    progname = os.path.basename(sys.argv[0])
    usage = f"""
    {progname} takes a 3D map and makes projections of it based on NumProj (a parameter passed into the program),
    then calculates the bispectrum on each projection. The average bispectrum is written out to a file.
    """

    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    parser.add_argument("--input", type=str, required=True, help="Path to the 3D map file (MRC or HDF format).")
    parser.add_argument("--output", type=str, required=True, help="Path to the output file for the average bispectrum (MRC or HDF format).")
    parser.add_argument("--numproj", type=int, default=1000, help="Number of projections to generate. Default is 1000.")

    # Parse the command-line arguments
    options, args = parser.parse_args()

    # Run the projection and bispectrum calculation
    project_map_and_compute_bispectrum(options.input, options.numproj, options.output)






























""""


def e3MakeBispectrumFromMap.py

mapSh = EMData('TwoCylindersAt60Degrees.hdf')

print(os.getcwd())

if 1:
    sym = Symmetries.get("c1"); # Create Projections.  Can skip this if they have already been created
    directionsEven600 = sym.gen_orientations("eman",{"n":600,"inc_mirror":False});#  Change

    directionsEven600min = []
    for idirection,t in enumerate(directionsEven600):
        trot = t.get_rotation()
        #taz  = trot['az']* np.pi/180.0
        #talt = trot['alt']* np.pi/180.0
        #tphi = - taz;
        if 0:
            trot['phi'] = - trot['az']
        t.set_rotation(trot);
        #tOut.set_rotation(
        directionsEven600min.append(t);



    XList = []
    YList = []
    ZList = []

    #seed(0)
    for idirection,t in enumerate(directionsEven600min):

        if (idirection%100)==0: print(idirection)
        if 1:
            prj = mapSh.project("standard",t);
            prj.process_inplace("normalize.edgemean");
            prj.write_image('EvenDistributionOf600_TwoCylindersAt60Degrees.hdf',idirection);
        #prj.translate( (4*(random.random()-0.5)) ,(4*(random.random()-0.5)),0);
        #prj.rotate(random.random()*360,0,0);    # Now apply noise
        #prj.process_inplace("normalize.edgemean");
        #prj.write_image(JustProjRT,idirection);
        trot = t.get_rotation()
        taz  = trot['az']* np.pi/180.0
        talt = trot['alt']* np.pi/180.0

        XNow = np.sin(talt) * np.sin(taz);      XList.append(XNow);
        YNow = np.sin(talt) * np.cos(taz);      YList.append(YNow);
        ZNow = np.cos(talt);                    ZList.append(ZNow);
    print(len(directionsEven600))


    Xarray = (np.array(XList)*1.0e6) ; Xarray = Xarray.astype(int)/1.0e6
Yarray = (np.array(YList)*1.0e6) ; Yarray = Yarray.astype(int)/1.0e6
Zarray = (np.array(ZList)*1.0e6) ; Zarray = Zarray.astype(int)/1.0e6
Zunique = np.unique(Zarray)
Rarray = np.sqrt(Xarray*Xarray+Yarray*Yarray);
PhiArray = np.arctan2(Yarray,Xarray);
ROutArray = np.sqrt(1 - np.sqrt(1-Rarray*Rarray));

XOutArray = ROutArray* np.cos(PhiArray);
YOutArray = ROutArray* np.sin(PhiArray);

Zunique,ZuniqueInds, ZuniqueCounts = np.unique(Zarray,  return_index=True, return_counts=True)

print(Zunique)
print(ZuniqueCounts)
print(ZuniqueInds[0])


plt.figure(figsize=(3,3))
plt.scatter(Xarray,Yarray,s=2)

plt.figure(figsize=(3,3))
plt.scatter(XOutArray,YOutArray,s=2)

plt.figure(figsize=(3,3))
plt.scatter(Xarray,Zarray,s=2)

#
if 0:
    XX = np.asarray(AllProjPoints[:,0])
    YY = np.asarray(AllProjPoints[:,1])

    XXB= XX.copy()
    YYB= YY.copy()

    XXB[IndsZZneg] *=-1
    YYB[IndsZZneg] *=-1
    CosThetaB[IndsZZneg] *=-1

    CosTheta2 = CosThetaB*CosThetaB
    SinThetaB = np.sqrt(1-CosTheta2)
    ThetaB = np.arcsin(SinThetaB)
    FuncThetaB = ThetaB/(np.pi/2)
    FuncThetaB = np.sin(ThetaB/2)/(np.sin(np.pi/4))

    XXFinal=XXB.copy()
    YYFinal=YYB.copy()

    GoodInds= np.where(np.abs(SinThetaB)>.0001)[0]
    XXFinal[GoodInds] = XXB[GoodInds]*FuncThetaB[GoodInds]/SinThetaB[GoodInds]
    YYFinal[GoodInds] = YYB[GoodInds]*FuncThetaB[GoodInds]/SinThetaB[GoodInds]


    PlotTheseProjections = 'EvenDistributionOf600_TwoCylindersAt90Degrees.hdf'
PlotTheseProjections = 'EvenDistributionOf600_TwoCylindersAt60Degrees.hdf'
prj.read_image(PlotTheseProjections,0);
XSize = prj.get_xsize();

XOutArrayby3 = XOutArray[::3];  YOutArrayby3 = YOutArray[::3]
XOutArrayby4 = XOutArray[::4];  YOutArrayby4 = YOutArray[::4]

print(len(XOutArray) , len(XOutArrayby3)  , len(XOutArrayby4) )

if 0:
    prj.read_image(PlotTheseProjections,0);
    prjNP = prj.numpy();
    print(np.mean(prjNP))

    prjNPs = np.shape(prjNP)[0]
    prjNP_N = prjNP + 0.1*np.random.randn(prjNPs,prjNPs)

    plt.imshow(prjNP,cmap='gray')


TwoDList    = np.vstack([XOutArray,YOutArray]).T
TwoDListBy3 = np.vstack([XOutArrayby3,YOutArrayby3]).T
TwoDListBy4 = np.vstack([XOutArrayby4,YOutArrayby4]).T
print(np.shape(TwoDList), TwoDList[100,1], YOutArrayby4[100])

ShiftXY = 1.05
TwoDList += [ShiftXY,ShiftXY]; TwoDList *= 300;
TwoDListBy3 += [ShiftXY,ShiftXY]; TwoDListBy3 *= 270;
TwoDListBy4 += [ShiftXY,ShiftXY]; TwoDListBy4 *= 300;

sTwoDList = np.shape(TwoDList)[0]
sTwoDListBy3 = np.shape(TwoDListBy3)[0]
sTwoDListBy4 = np.shape(TwoDListBy4)[0]

# Define the size of the canvas (change as needed)
canvas_size = 600;

# Create a larger canvas to place the images on
canvas = np.zeros((canvas_size, canvas_size))



# Loop through the list of coordinates and images to place them on the canvas
for jCoords  in range(sTwoDListBy3):
    coordsNow = TwoDListBy3[jCoords,:]
    #print(i)
    #if iImage>5: continue
    #print(iImage,coords)
    prj.read_image(PlotTheseProjections,3*jCoords);
    prjNP = prj.numpy();
    #print(np.mean(prjNP))

    sprjNP = np.shape(prjNP)[0]
    image = prjNP + 0.1*np.random.randn(sprjNP,sprjNP)


    try:
        place_image_on_canvas(canvas, image, int(coordsNow[0]+0.5), int(coordsNow[1]+0.5))
    except:
        #print(jCoords)
        continue

if 0:
    plt.figure(figsize=(28, 28))
    plt.imshow(canvas, cmap='gray')
    plt.title('Images Centered at Specified Coordinates')
    plt.grid(False)
    plt.show()


    #e2proc2d.py EvenDistributionOf600_12639.hdf  EvenDistributionOf600_12639_Np1.hdf --process math.addnoise=0.
rfp=4
size=22
size=16
kVal=11


print(sTwoDList,PlotTheseProjections)
EvenDistributionFileNamewNoise='EvenDistributionOf600_TwoCylindersAt90DegreesN6.hdf'
BSOutNoiseFileName='OutBispecNoise6.hdf'

EvenDistributionFileNamewNoise='EvenDistributionOf600_TwoCylindersAt60DegreesN6.hdf'
BSOutNoiseFileName='OutBispec60DegreesNoise6.hdf'

# Loop through the list of coordinates and images to place them on the canvas
np.random.seed(0)
prjNow = prj.copy()

for jProj  in range(sTwoDList):
    #coordsNow = TwoDList[jProj,:]
    #if jProj>1: continue
    prjNow.read_image(PlotTheseProjections,jProj);

    prjNoise =   prjNow.process("math.addnoise",{"noise":6});
    #a1.process_inplace("math.addnoise",{"noise":3})
    #prjNP1 = prj.numpy();
    #print(np.mean(prjNP))

    prjNoise.write_image(EvenDistributionFileNamewNoise,jProj)
    #sprjNP = np.shape(prjNP)[0]
    #bispecOutRFP=  prj.process("math.bispectrum.slice",{"rfp":rfp,"size":size});
    #image = prjNP + 0.1*np.random.randn(spr)
    #try:
    #    place_image_on_canvas(canvas, image, int(coordsNow[0]+0.5), int(coordsNow[1]+0.5))
    #except:
    #    print(jCoords)
    #    continue
    bispecOutNoise=  prjNoise.process("math.bispectrum.slice",{"rfp":rfp,"size":size});
    bispecOutNoise.write_image(BSOutNoiseFileName,jProj)
    #try:
    #    bispecOutNoise.write_image(BSOutNoiseFileName,jProj)
    #except RuntimeError as e:
    #    print(f"RuntimeError: {e}")

    if jProj==0:
        bispecAve = bispecOutNoise.copy()
    else:
        bispecAve += bispecOutNoise
    #bispecOutNoise.write_image(BSOutNoiseFileName, jProj)

"""
