#!/usr/bin/env python
# coding: utf-8

# Usage python e2InterpolateJsonEulersAmongFrames.py --JsonIn StackOfSubaverages.json --JsonOut FullFramesFN.json --numFrames 1495

# Usage python e2InterpolateJsonEulersAmongFrames.py --JsonIn stage2_full_options.json --JsonOut FullFramesFN.json --numFrames 1495

#      These are the fields for the json file
#		# Save selected options required for Stage 2
#		required_options = {
#			"tlt_params":tlt_params,
#			"pks":pks,
#			"ali_loss": loss0.tolist(),
#			"apix_init": options.apix_init,
#			"autoclipxy": options.autoclipxy,
#			"basename": options.basename,
#			"clipz": options.clipz,
#			"compressbits": options.compressbits,
#			"ctf": options.ctf,
#			"extrapad": options.extrapad,
#			"filterto": options.filterto,
#			"inputname": options.inputname,
#			"moretile": options.moretile,
#			"normslice": options.normslice,
#			"reconmode": options.reconmode,
#			"threads": options.threads,
#			"tltkeep": options.tltkeep
#		}# filterto and reconmde are new

# This will take in the JsonIn ,
#  keep all the parameters the same
#   except for tlt_params
#  it will reinterpolate the data from the first  small number which is the len(tlt_params)  to the final number numFrames



import os
import sys
import numpy as np
from scipy.linalg import fractional_matrix_power
#from EMAN2 import EMData, EMArgumentParser, EMANVERSION, Transform
from EMAN2 import *

from matplotlib import pylab, cm
import matplotlib.pyplot as plt



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



# Implementing the correct version that guarantees:
# - First element = n//2
# - Last element = N - n//2
# - Total of n elements
# - Gaps are only base_step or base_step+1
# - Evenly distributed +1 steps to avoid large last gap

def generate_final_correct_array(start, end,n, N):
    #start = n // 2
    #end = N - start
    total_range = end - start
    steps_needed = n - 1

    base_step = total_range // steps_needed
    extra_steps = total_range % steps_needed

    # Spread the +1s as evenly as possible
    step_pattern = []
    for i in range(steps_needed):
        if (i * extra_steps) % steps_needed < extra_steps:
            step_pattern.append(base_step + 1)
        else:
            step_pattern.append(base_step)

    # Build the array
    values = [start]
    for step in step_pattern:
        values.append(values[-1] + step)

    return values, step_pattern


def e2InterpolateEulersEngine(tlt_params_In, nImgFrames):
    """
    Inspects an Image Stack consisting of subaveraged tilts for their transform objects.
           Interpolates these transforms to write to a set of much larger numbers of frames.

    Parameters:
    - SubAvgSeriesFN:         path to the stack of SubAverages
    - FullFramesFN:            path to the stack of Full Frames
    - AvgNumberFramesPerTilt: integer that was the number of frames used to average into each subaverage
    """

    #nImgFrames              =EMUtil.get_image_count(FullFramesFN);
    print("Total number of Frames to be interpolated to="+str(nImgFrames))


    # Section  1 create a list of transforms from the SA frames
    #print(JsonIn)

    #JSFileData = js_open_dict(JsonIn)
    #print(JSFileData['tlt_params'])
    #print(JSFileData['tlt_params'])
    # output is   dict_keys(['ali_loss', 'apix_unbin', 'tlt_file', 'tlt_params'])
    #print(JSFileData)
    #print(JsonIn,js_list_dicts(JsonIn))
    #print("Keys in JSON:", JSFileData.keys())

    tltParams = np.array(tlt_params_In)
    #  0 tx x shift 1 ty y shift  2  Z rotation       3 Y rotation    4   X rotation

    ShiftXArray, ShiftYArray, ZEulerArray, YEulerArray, XEulerArray = tltParams[:,0] , tltParams[:,1] , tltParams[:,2], tltParams[:,3], tltParams[:,4]
    #print(tltParams)
    NumTilts = len(ShiftXArray)

    #tltParams_Out = tltParams


    AvgNumberFramesPerTilt = nImgFrames//NumTilts

    print("Number of Subaverages (equivalently tilts)="+str(NumTilts), '.  AvgNumberFramesPerTilt = ', AvgNumberFramesPerTilt )

    # Section 1.  Create a list of transforms


    #return

    TransformList = [];

    for jTilt in range(NumTilts):
        trNow=Transform();#  PRB Initialize  with next 2 lines
        trNow.set_rotation({"type":"xyz","xtilt":XEulerArray[jTilt], "ytilt":YEulerArray[jTilt],"ztilt":ZEulerArray[jTilt] })
        trNow.set_trans(ShiftXArray[jTilt],ShiftYArray[jTilt],0)
        TransformList.append(trNow)

    print('Length Transform List, ' ,len(TransformList))


    PrintDiag=1; # If 1, print for debugging
    WriteToFile=1; # If 1,  Write to file


    for array_Ind in range(NumTilts):
        MatPrior = TransformList[array_Ind].get_matrix_4x4(); MatPrior4by4 = np.array(MatPrior).reshape(4,4)
        values_str = ' '.join(f"{x:.4f}" for x in MatPrior)
        if PrintDiag: print(f"{array_Ind} [{values_str}]")




    # Section 2.  Create a list of integers.  The first is at offset, and the last is at -offset

    offset=  AvgNumberFramesPerTilt//2;
    MatrixFlatOut = np.zeros([nImgFrames,12])

    final_array, step_sizes = generate_final_correct_array(offset,nImgFrames-offset, NumTilts, nImgFrames)

    print(np.shape(final_array),final_array)
    print(step_sizes)
    print('np.shape(final_array),np.shape(step_sizes')
    print(np.shape(final_array), np.shape(step_sizes))


    MatrixFlatFixed =  TransformList[0].get_matrix_4x4()
    for jBegin in range(offset):
        MatrixFlatOut[jBegin,:] =MatrixFlatFixed[:12]
    #MatrixFlatOut[0,:]  = MatrixFlatFixed[:12]

    if PrintDiag:
        print('np.shape(MatrixFlatFixed)',np.shape(MatrixFlatFixed))
        print(MatrixFlatFixed)


    MatrixFlatFixed = TransformList[-1].get_matrix_4x4()
    for jEnd in range(offset):
        MatrixFlatOut[-jEnd-1,:] =MatrixFlatFixed[:12]

    for final_array_Ind in range(NumTilts-1):
        final_Ind_Now = final_array[final_array_Ind]
        MatPrior = TransformList[final_array_Ind].get_matrix_4x4(); MatPrior4by4 = np.array(MatPrior).reshape(4,4)
        MatNext   = TransformList[final_array_Ind+1].get_matrix_4x4(); MatNext4by4 =  np.array(MatNext).reshape(4,4)
        MatBridge4by4 = np.dot(MatNext4by4,np.linalg.inv(MatPrior4by4))# apply this to MatPrior to get MatNext
        step_sizeNow = step_sizes[final_array_Ind]
        MatStep4by4= np.real(fractional_matrix_power(MatBridge4by4,1.0/step_sizeNow))
        for jStepNow in range(step_sizeNow):
            MatStepTotalNow = np.linalg.matrix_power(MatStep4by4,jStepNow)
            MatFinalNow = np.dot(MatStepTotalNow,MatPrior4by4)
            MatrixFlatOut[final_Ind_Now+jStepNow,:] = MatFinalNow.flatten().tolist()[:12]

    if PrintDiag:
        print(MatFinalNow)
        print(MatFinalNow.flatten().tolist())

    for jMat in range(nImgFrames):
        MatFlatNow = MatrixFlatOut[jMat, :]
        values_str = ' '.join(f"{x:.4f}" for x in MatFlatNow)
        if PrintDiag: print(f"{jMat} [{values_str}]")

    #return

    print('Conversion to tlt_params_Out  ', nImgFrames)

    tltParams_Out=[]

    for jMatNow in range(nImgFrames):
        tltParams_Out_Now =  MatrixFlatOut[jMatNow].flatten()
        #print(jMatNow,tltParams_Out_Now)
        trNow = Transform(tltParams_Out_Now.tolist()[:12])

        XYZtilts = trNow.get_rotation("xyz")
        XYZshifts = trNow.get_trans()

        Xtilt = XYZtilts['xtilt'];  Ytilt = XYZtilts['ytilt'] ;  Ztilt = XYZtilts['ztilt']
        xshift = XYZshifts[0];  yshift = XYZshifts[1] ;

        #ShiftXArray, ShiftYArray, ZEulerArray, YEulerArray, XEulerArray = tltParams[:,0] , tltParams[:,1] , tltParams[:,2], tltParams[:,3], tltParams[:,4]

        tlt_paramsNow = [   xshift, yshift, Ztilt, Ytilt, Xtilt  ]
        tlt_paramsNow_rounded = [round(v, 3) for v in tlt_paramsNow]
        tltParams_Out.append(tlt_paramsNow_rounded)

    tltparamsNP = np.array(tltParams_Out)
    plt.figure()
    plt.plot(tltparamsNP[:,3])
    plt.savefig('tltParams.jpg')

    return tltParams_Out


if __name__ == "__main__":
#def main():
    # Setup the command-line argument parser
    progname = os.path.basename(sys.argv[0])
    usage = f"""
    {progname} Takes a json file describing a stack of subaveraged tilts.
        Interpolates thes transform objects for a small stack  to a set of much larger numbers of frames.
        A list file (.lst) is written to the sets directory, based on the fullframes file name.
        """

    #print('Hi')

    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    parser.add_argument("--JsonIn", type=str, required=True, help="path to the json with initial eulers")
    parser.add_argument("--JsonOut", type=str, required=True, help="path to the json with final eulers")
    parser.add_argument("--finalstackName", type=str, required=True, help=" name of the final image stack")

    # Parse the command-line arguments
    options, args = parser.parse_args()

    # Initialize logging for this command
    logid=E2init(sys.argv)


    # --- your EMAN2 code here ---
    print("Doing some EMAN2 stuff")

    # tltParams = np.array(JSFileData['tlt_params'])
    # #  0 tx x shift 1 ty y shift  2  Z rotation       3 Y rotation    4   X rotation

    # ShiftXArray, ShiftYArray, ZEulerArray, YEulerArray, XEulerArray = tltParams[:,0] , tltParams[:,1] , tltParams[:,2], tltParams[:,3], tltParams[:,4]
    # print(tltParams)


    # Load original JSON from file
    with open(options.JsonIn, "r") as f:
        data = json.load(f)

    ali_loss_In = data['ali_loss']
    tlt_params_In = data['tlt_params']

    #basename = data['basename']
    inputname = data['inputname']
    numFramesIn = EMUtil.get_image_count(inputname);


    fileNameFrames = options.finalstackName
    numFramesOut  = EMUtil.get_image_count(fileNameFrames);

    basename = os.path.basename(fileNameFrames)
    data['basename']  =  basename
    data['inputname'] =  fileNameFrames


    # Create evenly spaced original and target indices
    x_old = np.linspace(0, 1, numFramesIn)
    x_new = np.linspace(0, 1, numFramesOut)

    # Interpolate
    ali_loss_Out = np.interp(x_new, x_old, ali_loss_In)
    ali_loss_Out = ali_loss_Out.tolist()


    data['ali_loss']  = ali_loss_Out

    if 0:
        tlt_params_Out = tlt_params_In
    else:
        # Run the main loop
        tlt_params_Out = e2InterpolateEulersEngine(tlt_params_In, numFramesOut)
        print("This is tlt_params_Out \n")
        print(tlt_params_Out)

    print('type(tlt_params_Out)')
    print(type(tlt_params_Out))

    # Replace the 'tlt_params' key with your new value
    data['tlt_params'] = tlt_params_Out
    #data['tlt_params'] = tlt_params_Out.tolist() if hasattr(tlt_params_Out, "tolist") else tlt_params_Out

    # Write updated data to new JSON file
    with open(options.JsonOut, "w") as f:
        json.dump(data, f, indent=2)

    # End logging
    E2end(logid)
    #python ./e2InterpolateEulers.py --subaverages gainCorrectedFilewoShiftSA30_ct_7_.hdf  --fullframes Copy_CT7__4Kby4K_Gain_Shifted.hdf --numSA 30
    #python ~/src/eman2/examples/e2InterpolateEulersAmongFrames.py --subaverages gainCorrectedFilewoShiftSA30_ct_7_.hdf  --fullframes Copy_CT7__4Kby4K_Gain_Shifted.hdf --numSA 30 >Out.txt
    #python ./e2InterpolateEulers.py --subaverages gainCorrectedFilewoShiftSA30_ct_7_.hdf  --fullframes Copy_CT7__4Kby4K_Gain_Shifted.hdf --numSA 30


#      e2InterpolateJsonEulersAmongFrames.py --JsonIn gainCorrectedFilewShiftSA60_ct_7_wTransformSh4_EN_recon_tomo_final.json --JsonOut FullFrames1495FN.json --finalstackName 20240630_77469_Neuron_continuousct_7__4Kby4K_Gain_ShiftedSh4.hdf  > TempNewAngs.txt
# Usage python e2InterpolateJsonEulersAmongFrames.py --JsonIn StackOfSubaverages.json --JsonOut FullFramesFN.json --finalstackName 20240630_77469_Neuron_continuousct_7__4Kby4K_Gain_ShiftedSh4.hdf

#      These are the fields for the json file
#		# Save selected options required for Stage 2
#		required_options = {
#			"tlt_params":tlt_params,
#			"pks":pks,
#			"ali_loss": loss0.tolist(),
#			"apix_init": options.apix_init,
#			"autoclipxy": options.autoclipxy,
#			"basename": options.basename,
#			"clipz": options.clipz,
#			"compressbits": options.compressbits,
#			"ctf": options.ctf,
#			"extrapad": options.extrapad,
#			"filterto": options.filterto,
#			"inputname": options.inputname,
#			"moretile": options.moretile,
#			"normslice": options.normslice,
#			"reconmode": options.reconmode,
#			"threads": options.threads,
#			"tltkeep": options.tltkeep
#		}# filterto and reconmde are new

# This will take in the JsonIn ,
#  keep all the parameters the same
#   except for tlt_params
#  it will reinterpolate the data from the first  small number which is the len(tlt_params)  to the final number numFrames


