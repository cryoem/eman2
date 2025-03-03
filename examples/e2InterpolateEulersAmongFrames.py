#!/usr/bin/env python
# coding: utf-8

# Usage python e2InterpolateEulersAmongFrames.py --subaverages StackOfSubaverages.hdf --fullframes FullFramesFN.hdf --numSA 30

import os
import sys
import numpy as np
from scipy.linalg import fractional_matrix_power
#from EMAN2 import EMData, EMArgumentParser, EMANVERSION, Transform
from EMAN2 import *



# js_open_dict(info_name(ptcl_name))
#
import random

# import glob
# from tqdm.auto import tqdm
# import scipy.ndimage as ndimage
from matplotlib import pylab, cm
import matplotlib.pyplot as plt 
import glob

#import datetime
#from time import time
#from datetime import date

import json





def e2InterpolateEulersEngine(SubAvgSeriesFN, FullFramesFN, AvgNumberFramesPerTilt):
    """
    Inspects an Image Stack consisting of subaveraged tilts for their transform objects.
           Interpolates these transforms to write to a set of much larger numbers of frames.

    Parameters:
    - SubAvgSeriesFN:         path to the stack of SubAverages
    - FullFramesFN:            path to the stack of Full Frames
    - AvgNumberFramesPerTilt: integer that was the number of frames used to average into each subaverage
    """


    # 0   Section 0  Set up   variables and functions

    #  ProjNow,tr = MatToHeader(FullFramesFN, PositionInFileName, Matrix4by4)
    def MatToHeader(FullFramesFN, PositionInFileName, Matrix4by4):
        FullFramesListFN = 'sets/'+FullFramesFN[:-4]+'.lst'
        lsx=LSXFile(FullFramesListFN)
        ProjNow =  EMData(FullFramesFN,PositionInFileName);
        tr= Transform(Matrix4by4.flatten()[:12].tolist())
        ProjNow.set_attr('xform.projection',tr)
        PN_attr_dict= ProjNow.get_attr_dict();
        lsx[PositionInFileName]=[PositionInFileName,FullFramesListFN,PN_attr_dict]
        #ProjNow.write_image(FullFramesListFN,PositionInFileName)
        #EMData.write_images(FullFramesListFN,[ProjNow],idxs=PositionInFileName,header_only=True)
        #return ProjNow,tr

    """
    >>> lsx=LSXFile("r3d_00/ptcls_01.lst")
    >>> print(len(lsx))
    399000
    >>> print(lsx[0])
    [0, 'particles/50ca-ND.hdf', {'class': 0, 'score': -0.128, 'xform.projection': Transform({'az':57.909,'alt':73.218,'phi':-85.855,'tx':-3.50,'ty':10.50,'tz':0.00,'mirror':0,'scale':1.0000,'type':'eman'})}]

    # this will update the header for the first particle (class=99)
    >>> lsx[0]=[0, 'particles/50ca-ND.hdf', {'class': 99, 'score': -0.128, 'xform.projection': Transform({'az':57.909,'alt':73.218,'phi':-85.855,'tx':-3.50,'ty':10.50,'tz':0.00,'mirror':0,'scale':1.0000,'type':'eman'})}]

    # note that assigning to [-1] will append to the end of the file:
    >>> lsx[-1]=[12,'particles/xyz.hdf',{"class":3}]
    >>> print(len(lsx))

    """

    if 0:
        NewContinuousTiltDataAt ='/home3/pbaldwin/Apollo/apollo_06222023/ContinuousTilt/ct_neurons_06282024/data/'
        CT_Directory = 'ct_7'
        #gainFileName = NewContinuousTiltDataAt+'/4Kby4K/GainFileMid500.hdf'
        os.chdir(NewContinuousTiltDataAt)
        os.chdir(CT_Directory)
        #
        SubAvgSeriesFN          =  'tiltseries/gainCorrectedFilewoShiftSA30_ct_7_.hdf'
        FullFramesFN            =  '20240630_77469_Neuron_continuousct_7__4Kby4K_Gain_Shifted.hdf'
        FullFramesFN            =  'Copy_CT7__4Kby4K_Gain_Shifted.hdf'
        AvgNumberFramesPerTilt  =    30


    nimgFrames              =EMUtil.get_image_count(FullFramesFN);
    print("Total number of Frames to be interpolated to="+str(nimgFrames))


    # Section  1 create a list of transforms from the SA frames

    JSFileName = info_name(SubAvgSeriesFN)
    JSFileData = js_open_dict(JSFileName)
    # print(JSFileData['tlt_params'])
    # output is   dict_keys(['ali_loss', 'apix_unbin', 'tlt_file', 'tlt_params'])
    # print(JSFileName,js_list_dicts(JSFileName))

    tltParams = np.array(JSFileData['tlt_params'])
    #  0 tx x shift 1 ty y shift  2  Z rotation       3 Y rotation    4   X rotation

    ShiftXArray, ShiftYArray, ZEulerArray, YEulerArray, XEulerArray = tltParams[:,0] , tltParams[:,1] , tltParams[:,2], tltParams[:,3], tltParams[:,4]

    NumTilts = len(ShiftXArray)

    print("Number of Subaverages (equivalently tilts)="+str(NumTilts))

    # create a list of transforms

    #return

    TransformList = [];

    for jTilt in range(NumTilts):
        trNow=Transform()
        trNow.set_rotation({"type":"xyz","xtilt":XEulerArray[jTilt], "ytilt":YEulerArray[jTilt],"ztilt":ZEulerArray[jTilt] })
        trNow.set_trans(ShiftXArray[jTilt],ShiftYArray[jTilt],0)
        TransformList.append(trNow)



    #     Section 2   create a List of the steps that one takes in heading from one tilt to another: MatStepArray (also MatSAList)
    Mat0 = np.array(TransformList[0].get_matrix_4x4()).reshape(4,4)
    MatPrior = Mat0.copy()
    MatSAList =[]
    MatStepArray = []

    MatSAList.append(Mat0)

    for jTilt in range(NumTilts-1):
        jTiltNext = jTilt+1;
        MatNext = np.array(TransformList[jTiltNext].get_matrix_4x4()).reshape(4,4)
        MatSAList.append(MatNext)
        MatBridge = np.dot(MatNext,np.linalg.inv(MatPrior))
        MatStep= np.real(fractional_matrix_power(MatBridge,1.0/AvgNumberFramesPerTilt))
        MatStepArray.append(MatStep)
        MatPrior =MatNext.copy()


    #       Section 3    Step between subaverages: write to full frames file
    # In each place where there is a print command, change to a write to the complete Stack

    PrintDiag=0; # If 1, print for debugging
    WriteToFile=1; # If 1,  Write to file

    # The First block all gets assigned to MatSAList[0]
    CountPrinted =0


    MatPrior     = MatSAList[0]
    for jStep in range(AvgNumberFramesPerTilt//2):
            vv = MatPrior.reshape(-1)
            if PrintDiag: print(f"{CountPrinted}, " +", ".join(f"{x:.4f}" for x in vv[:12]))
            PositionInFileName = CountPrinted
            if WriteToFile: MatToHeader(FullFramesFN, PositionInFileName, MatNext)
            CountPrinted +=1



    for jTilt in range(NumTilts-1):
        #if jTilt>1: continue
        #print(jTilt)
        MatPrior     = MatSAList[jTilt]
        MatStepNow   = MatStepArray[jTilt]
        if PrintDiag: print(CountPrinted,MatPrior.reshape(1,16));
        PositionInFileName = CountPrinted
        if WriteToFile: MatToHeader(FullFramesFN, PositionInFileName, MatPrior)
        #print(f"{CountPrinted}, " +", ".join(f"{x:.4f}" for x in MatPrior.reshape(1,16)))
        CountPrinted += 1
        for jStep in range(AvgNumberFramesPerTilt-1):
            MatNext = np.dot(MatStepNow,MatPrior)

            vv = MatNext.reshape(-1)
            if PrintDiag: print(f"{CountPrinted}, " +", ".join(f"{x:.4f}" for x in vv[:12]))
            PositionInFileName = CountPrinted
            if WriteToFile: MatToHeader(FullFramesFN, PositionInFileName, MatNext)
            CountPrinted +=1
            MatPrior= MatNext.copy()

    jTilt+=1;
    MatPrior     = MatSAList[jTilt]
    if PrintDiag: print(CountPrinted,MatSAList[jTilt].reshape(1,16))
    PositionInFileName = CountPrinted
    if WriteToFile: MatToHeader(FullFramesFN, PositionInFileName, MatPrior)
    CountPrinted += 1

    for jStep in range(AvgNumberFramesPerTilt-1):
        if CountPrinted >= nimgFrames: continue
        MatNext = MatPrior
        vv = MatNext.reshape(-1)
        if PrintDiag: print(f"{CountPrinted}, " +", ".join(f"{x:.4f}" for x in vv[:12]))
        PositionInFileName = CountPrinted
        if WriteToFile: MatToHeader(FullFramesFN, PositionInFileName, MatNext)
        CountPrinted +=1
        MatPrior= MatNext.copy()

    if PrintDiag: print(jTilt, (jTilt)*30,CountPrinted)



if __name__ == "__main__":
#def main():
    # Setup the command-line argument parser
    progname = os.path.basename(sys.argv[0])
    usage = f"""
    {progname} Inspects an Image Stack consisting of subaveraged tilts for their transform objects.
        Interpolates these transforms to write to a set of much larger numbers of frames.
        A list file (.lst) is written to the sets directory, based on the fullrames file name.
        """

    #print('Hi')

    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    parser.add_argument("--subaverages", type=str, required=True, help="path to the stack of SubAverages")
    parser.add_argument("--fullframes", type=str, required=True, help="path to the stack of Full Frames.")
    parser.add_argument("--numSA", type=int, default=30, help=" integer that was the number of frames used to average into each subaverage")

    # Parse the command-line arguments
    options, args = parser.parse_args()

    # Run the main loop
    e2InterpolateEulersEngine(options.subaverages, options.fullframes, options.numSA)

    #python ./e2InterpolateEulers.py --subaverages gainCorrectedFilewoShiftSA30_ct_7_.hdf  --fullframes Copy_CT7__4Kby4K_Gain_Shifted.hdf --numSA 30



