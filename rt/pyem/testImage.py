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
    
def get_img( file_name ):
    img = EMData()
    img.read_image( file_name )
    return img
    
def apply_gaussian( input_file_name, file_name, file_ext ):
    #apply a gaussian blob to input image
    print "Apply gaussian blob to input image", input_file_name, "..."
    img = get_img( input_file_name )
    params = create_param_gau()
    img.process("testimage.gaussian", params)
    img2 = get_img( input_file_name )
    img3 = img * img2
    output_filename = file_name+'_gau.'+file_ext
    img3.write_image( output_filename )
    print "Write processed image to--", output_filename
    
def apply_sinewave( input_file_name, file_name, file_ext ):
    #apply a sine wave to input image
    print "Apply sine wave to input image", input_file_name, "..."
    img = get_img( input_file_name )
    params = create_param_sin()
    img.process("testimage.sinewave", params)
    img2 = get_img( input_file_name )
    img3 = img * img2
    output_filename = file_name+'_sin.'+file_ext
    img3.write_image( output_filename )
    print "Write processed image to--", output_filename
    
def apply_square( input_file_name, file_name, file_ext ):
    #apply a square/cube to input image
    print "Apply square/cube to input image", input_file_name, "..."
    img = get_img( input_file_name )
    params = create_param_squ_cube()
    img.process("testimage.squarecube", params)
    img2 = get_img( input_file_name )
    img3 = img * img2
    output_filename = file_name+'_squ.'+file_ext
    img3.write_image( output_filename )
    print "Write processed image to--", output_filename
    
def apply_circle( input_file_name, file_name, file_ext ):
    print "Apply circle/sphere to input image", input_file_name, "..."
    img = get_img( input_file_name )
    params = create_param_cir_sph()
    img.process("testimage.circlesphere", params)
    img2 = get_img( input_file_name )
    img3 = img * img2
    output_filename = file_name+'_cir.'+file_ext
    img3.write_image( output_filename )
    print "Write processed image to--", output_filename
    
def apply_uniform_noise( input_file_name, file_name, file_ext ):
    print "Apply uniform noise to input image", input_file_name, "..."
    img = get_img( input_file_name )
    img.process("testimage.noise.uniform.rand")
    img2 = get_img( input_file_name )
    img3 = 0.1 * img2.get_attr("mean") * img + img2
    output_filename = file_name+'_noiuni.'+file_ext
    img3.write_image( output_filename )
    print "Write processed image to--", output_filename
    
def apply_gaussian_noise( input_file_name, file_name, file_ext ):
    print "Apply gaussian noise to input image", input_file_name, "..."
    img = get_img( input_file_name )
    params = create_param_noi_gau()
    img.process("testimage.noise.gauss", params)
    img2 = get_img( input_file_name )
    img3 = 0.1 * img2.get_attr("mean") * img + img2
    output_filename = file_name+'_noigau.'+file_ext
    img3.write_image( output_filename )
    print "Write processed image to--", output_filename
    
def apply_test_image( input_file_name ):
    file_ext = input_file_name[-3:]
    file_name = input_file_name[0:-4]
    
    apply_gaussian( input_file_name, file_name, file_ext )
    apply_sinewave( input_file_name, file_name, file_ext )
    apply_square( input_file_name, file_name, file_ext )
    apply_circle( input_file_name, file_name, file_ext )
    apply_uniform_noise( input_file_name, file_name, file_ext )
    apply_gaussian_noise( input_file_name, file_name, file_ext )
    
def testImage_main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--inputfile",action="store", type="string", dest="input_filename")
    
    from sys import argv
    (options, argv) = parser.parse_args(argv)
    fname = options.input_filename
    if fname != None:
        apply_test_image( fname )
    else:
        img = EMData()
        img.set_size(1024,1024,1)
        img.to_zero()
        create_test_image( img )

if __name__ == '__main__':
    testImage_main()