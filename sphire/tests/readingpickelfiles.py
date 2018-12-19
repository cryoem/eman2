import numpy
import pickle
import os
import EMAN2_cppwrap as e2cpp

ABSOLUTE_PATH =  os.path.dirname(os.path.realpath(__file__))
filepath = os.path.join(ABSOLUTE_PATH, "files/PICKLE.ornq")

with open(filepath, 'rb') as rb:
      (image,crefim,xrng,yrng,step,mode,numr,cnx,cny,deltapsi) = pickle.load(rb)
      print(numpy.shape(image.get_3dview()))
      print(numpy.shape(crefim.get_3dview()))
      print(xrng)
      print(yrng)
      print(step)
      print(mode)
      print(numpy.shape(numr))
      print(cnx)
      print(cny)
      print(deltapsi)
