#! /usr/bin/env python
# testImage.py  05/08/2005  Grant Tang
#this script is used for create all kinds of test images,
#it can take a command line option "-i inputimage" or run without any option

from EMAN2 import *

def create_param_gau():
    params = {}
    params["sigma"] = 200
    params["axis"] = 'y'
    params["c"] = 120
    return params

def create_param_sin():
    params = {}
    params["wave_length"] = 30
    params["axis"] = 'y'
    params["c"] = 50
    params["phase"] = 20
    return params
    
def create_param_squ_cube():
    params = {}
    params["edge_length"] = 300
    params["axis"] = 'y'
    params["odd_edge"] = 500
    params["fill"] = 'yes'
    return params
    
def create_param_cir_sph():
    params = {}
    params["radius"] = 200
    params["axis"] = 'y'
    params["c"] = 100
    params["fill"] = 'yes'
    return params

def create_param_noi_gau():
    params = {}
    params["noise_level"] = 0.4
    return params

def create_test_image( img ):
    #create a gaussian blob test image
    params = create_param_gau()
    img.process("testimage.gaussian", params)
    img.write_image("img_gau.mrc")
    print "Write test image img_gau.mrc..."
    
    #create a sine wave test image
    params = create_param_sin()
    img.process("testimage.sinewave", params)
    img.write_image("img_sin.mrc")
    print "Write test image img_sin.mrc..."
    
    #create a square/cube test image
    params = create_param_squ_cube()
    img.process("testimage.squarecube", params)
    img.write_image("img_squ.mrc")
    print "Write test image img_squ.mrc..."
    
    #create a circle/sphere test image
    params = create_param_cir_sph()
    img.process("testimage.circlesphere", params)
    img.write_image("img_cir.mrc")
    print "Write test image img_cir.mrc..."
    
    #create a uniform noise image
    img.process("testimage.noise.uniform.rand")
    img.write_image("img_noi_uni.mrc")
    print "Write test image img_noi_uni.mrc..."

    #create a Gaussian distributed noise
    params = create_param_noi_gau()
    img.process("testimage.noise.gauss", params)
    img.write_image("img_noi_gau.mrc")
    print "Write test image img_noi_gau.mrc..."

def testImage_main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--inputfile",action="store", type="string", dest="input_filename")
    
    from sys import argv
    (options, argv) = parser.parse_args(argv)
    if options.input_filename != None:
        img = EMData()
        img.read_image( options.input_filename )
        create_test_image( img )
    else:
        img = EMData()
        img.set_size(1024,1024,1)
        img.to_zero()
        create_test_image( img )

if __name__ == '__main__':
    testImage_main()
