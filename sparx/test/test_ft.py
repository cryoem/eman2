#!/bin/env python
from EMAN2  import *
from sparx  import *


myft = EMData()
myft.set_size(6,4)
myft.set_complex(True) # "True" in python, "true" in C++
myft.set_ri(True)
myft.set_fftpad(True)
myft.set_fftodd(False)
myft.print_image()
for ix in range(3):
      for iy in range(4):
          myft[ix,iy] = 2*ix+iy

myft.center_origin_fft()
myft.print_image()
myft.do_ift_inplace()
myft.print_image()
myft.do_fft_inplace()
myft.print_image()
