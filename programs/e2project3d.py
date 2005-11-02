#!/bin/env python

# initial version of project3d; more features/documentation to come
# try "e2project3d --help"

import sys, math, os, random
from EMAN2 import *
from optparse import OptionParser
deg2rad = math.pi / 180.0
rad2deg = 180.0 / math.pi

def main():

    (options, files) = load_args()
    
    logger=E2init(sys.argv)

    data = EMData()

    try:
        data.read_image(files[0])
    except:
        sys.stderr.write("ERROR: Input file: %s does not exist\n" % files[0])
        exit(1)

    if (options.dump_angles):
        try:
            os.unlink("proj_angles.txt")
        except OSError :
            pass
        try:
            options.proj_angle_file = open("proj_angles.txt", "w")
            if (options.angletype == "EMAN"):
                options.proj_angle_file.write("x\talt\taz\tphi\n")
            if (options.angletype == "SPIDER"):
                options.proj_angle_file.write("x\ttheta\tphi\tpsi\n")
        except:
            sys.stderr.write("ERROR: Output list file: %s could not be opened\n"
                             % options.list_file)
            options.proj_angle_file = ""
    else:
        options.proj_angle_file = ""
    

    try:
        os.unlink(options.outfile[:len(options.outfile) - 3] + "img")
        os.unlink(options.outfile[:len(options.outfile) - 3] + "hed")
    except OSError :
        pass
    
    
    # setup alt1, az1, phi1 based on sym
    az1 = 2.0 * math.pi
    alt1 = math.pi / 2.0
    alt0, az0 = 0.0, 0.0

    if options.sym:
        if (options.sym[0] in ["c","d"]):
            val = int(options.sym[1:])
            az1 = az1 / (val * 1.0)
            if (options.sym[0]=="d"):
                az1 = az1 / 2.0
        elif (options.sym[0] == "h"):
            val = int(options.sym[1:])
            az1 = val * math.pi / 180.
            alt0 = 80 * math.pi / 180.
        elif (options.sym == "icos"):
            val = 60
            az1 = 10.
            alt1 = 37.3773681406497 * math.pi / 180.
        elif (options.sym == "oct"):
            val = 24
            az1 = az1 / 4.0
            alt1 = 54.736 * math.pi / 180.
        else:
            sys.stderr.write("ERROR: Invalid symmetry type\n")
        
    if options.prop:
        #        prop(data, options, val, az0, alt0, az1, alt1)
        prop(data, options, val, az0, az1)
    elif options.grid:
        grid(data, options, az0, alt0, az1, alt1)
    elif options.user_angles:
        user_angles(data, options)
    else:  #create single projection
        if (options.pad):
            data = data.clip((data.get_xsize() - options.pad) / 2,
                             (data.get_ysize() - options.pad) / 2,
                             pad, pad)
            
        if options.euler == []:
            p=Projectors.get(options.projector, {"alt":0.0,
                                                 "az":0.0,
                                                 "phi":0.0,
                                                 "angletype":options.angletype,
                                                 "mode":options.mode
                                                 })
            angs = [0, 0, 0]

        else:
            p=Projectors.get(options.projector, {"az":options.euler[0],
                                                 "alt":options.euler[1],
                                                 "phi":options.euler[2],
                                                 "angletype":options.angletype,
                                                 "mode":options.mode})

            angs = options.euler
            

        sys.stdout.write("Angle: %4.2f\t%4.2f\t%4.2f\n" % (angs[1], angs[0], angs[2]))
        if (options.proj_angle_file):
            options.proj_angle_file.write("%d\t%4.2f\t%4.2f\t%4.2f\n" %
                                          (1 ,angs[1], angs[0], angs[2]))


            
        output = p.project3d(data)
        output.set_rotation(angs[0] * deg2rad, angs[1] * deg2rad, angs[2] * deg2rad)
        output.process("eman1.mask.sharp", {"outer_radius":options.mask})

        try:
            output.write_image(options.outfile, -1)
        except:
            sys.stderr.write("ERROR: Cannot write to file %s" % options.outfile)
            exit(1)

    print "Output File: %s" % options.outfile

    if (options.proj_angle_file):
        try:
            options.proj_angle_file.close()
            print "Angle Projection File: proj_angles.txt"
        except:
            sys.stderr.write("ERROR: Could not close file proj_angles.txt")

    E2end(logger)          # enable this once things are more stable

def user_angles(data, options):
    count = 0
    for x in read_angles(options.user_angles, options.angletype):
        sys.stdout.write("%d\t%4.2f\t%4.2f\t%4.2f\n" % (count, alt, az, phi))

        if (options.proj_angle_file):
            options.proj_angle_file.write("%d\t%4.2f\t%4.2f\t%4.2f\n" % (count, alt, az, phi))
                
        p = Projectors.get(options.projector, {"alt" : alt,
                                               "az" : az, 
                                               "phi" : phi,
                                               "angletype":options.angletype,
                                               "mode": options.mode})
        q = p.project3d(data)
        q.set_rotation(az, alt, phi)
        q.process("eman1.mask.sharp", {"outer_radius":options.mask})

        try:
            q.write_image(options.outfile, -1)
        except:
            sys.stderr.write("Error: Cannot write to file %s" % options.outfile)
            exit(1)

        count = count + 1


def grid(data, options, az0, alt0, az1, alt1):
    print "Projection grid %d x %d" % (options.grid, options.grid)

    phi = az = 0
    i, k, alt = 0, 0, alt0

    while (i < options.grid):

        j, az = 0, az0
        while (j < options.grid):

            if (options.randomphi):
                phi = random.random() * 2.0 * math.pi
            if (options.phic):
                phi = -az * math.cos(alt)

            if (options.proj_angle_file):
                options.proj_angle_file.write("%d\t%4.2f\t%4.2f\t%4.2f\n" %
                                              (j + i * options.grid, alt * rad2deg,
                                               az * rad2deg, phi * rad2deg))
            sys.stdout.write("%d\t%f\t%f\t%f\n" % (j + i * options.grid, alt * rad2deg,
                                                   az  * rad2deg, phi * rad2deg))
            
            p = Projectors.get(options.projector, {"alt" : alt * rad2deg,
                                                   "az" : az * rad2deg, 
                                                   "phi" : phi * rad2deg,
                                                   "angletype":options.angletype,
                                                   "mode": options.mode})
            q = p.project3d(data)
            q.set_rotation(az, alt, phi)
            q.process("eman1.mask.sharp", {"outer_radius":options.mask})

            if (options.pad > 0.0):
                q = q.get_clip( Region(q.get_xsize() - options.pad / 2.0,
                                       q.get_ysize() - options.pad / 2.0,
                                       options.pad, options.pad))
                    
            try:
                q.write_image(options.outfile, k)
            except:
                sys.stderr.write("Error: Cannot write to file %s" % options.outfile)
                exit(1)

            j = j + 1
            k = k + 1
            az = az + az1 / float(options.grid)        

        i = i + 1
        alt = alt + alt1 / (options.grid - 1.0)
        
    

    

#def prop(data, options, val,  az0, alt0, az1, alt1):
def prop(data, options, val, az0, az1):
    alt0 = options.prop_start * deg2rad
    alt1 = options.prop_end * deg2rad
    
    phi = az = 0
    if (val == 60):
        val = -5
    if (val == 24):
        val = -4

    i, j, alt2 = 1, 0, alt0
    while (alt2 <= alt1 + options.prop * math.pi / 360.):
        alt = alt2
        h = math.floor(360. / (options.prop * 1.1547))
        h = int(math.floor(h * math.sin(alt) + .5))
        if (h == 0):
            h = 1
        h = abs(val) * math.floor(h / float(abs(val)) + .5)

        if h==0:
            h=2.0**30                       # ok to do this?
        else:
            h = math.pi * 2.0 / h
        if (alt > 0) and ((az1 - az0) / h < 2.8):
            h = (az1 - az0) / 2.1
        if (alt == 0):
            h = az1
            
        az2 = az0
        if (j % 1):                        # when is this ever true??
            az2 = az2 + h / 2.0

        while (az2 < az1 - h / 4):
            alt = alt2
            az = az2
     
            phi = 0
                
            if (options.phitoo):
                limit = math.pi * 2. - options.phitoo * math.pi / 360.
                step = options.phitoo * math.pi / 180.
            else:
                limit = .0001
                step = math.pi
            while (phi < limit):

                if (az > math.pi and alt > math.pi / 2 -.001 and alt < math.pi / 2.0 + .001):
                    continue
                if (az > math.pi / 4.0):
                    tmp_check = math.pi / 2.0 - az
                else:
                    tmp_check = az
                if (val == -4 and math.tan(alt) * math.cos(tmp_check) > 1.0):
                    i = i -1 #423
                    continue

                if (options.randomphi):
                    phi = random.random() * 2.0 * math.pi
                if (options.phic):
                    phi = -az * math.cos(alt)
              
                if (val == -5):
                    az = az + 3. * math.pi / 2.

                sys.stdout.write("%d\t%4.2f\t%4.2f\t%4.2f\n" % (i, alt * rad2deg,
                                                                az * rad2deg,
                                                                phi * rad2deg))

                q = EMData()
                q2 = EMData()
                if (options.phitoo and options.smear):   #haven't checked these options yet
                    p = Projectors.get(options.projector, {"alt" : alt * rad2deg,
                                                           "az" : az * rad2deg, 
                                                           "phi" : phi * rad2deg,
                                                           "angletype":options.angletype,
                                                           "mode": options.mode})
                    q = p.project3d(data)

                    smr = phi + options.phitoo * math.pi / 1800.

                    while (smr < phi + options.phitoo * math.pi / 180):
                        p = Projectors.get(options.projector, {"alt" : alt * rad2deg,
                                                               "az" : az * rad2deg, 
                                                               "phi" : smr * rad2deg,
                                                               "angletype":options.angletype,
                                                               "mode": options.mode})
                        q2 = p.project3d(data)
                        q.add(q2)
                            
                        smr = smr + options.phitoo * math.pi / 1800.
                else:
                    p = Projectors.get(options.projector, {"alt" : alt * rad2deg,
                                                           "az" : az * rad2deg, 
                                                           "phi" : phi * rad2deg,
                                                           "angletype":options.angletype,
                                                           "mode": options.mode
                                                           })
                    q = p.project3d(data)

                    if (options.proj_angle_file):
                        options.proj_angle_file.write("%d\t%4.2f\t%4.2f\t%4.2f\n" %
                                                      (i, alt * rad2deg,
                                                       az * rad2deg,
                                                       phi * rad2deg))
                        
                # q.process() #figure this out later
                q.set_rotation(az, alt, phi)
                q.process("eman1.mask.sharp", {"outer_radius":options.mask})

                if (options.pad > 0.0):
                    q = q.get_clip( Region(q.get_xsize() - options.pad / 2.0,
                                           q.get_ysize() - options.pad / 2.0,
                                           options.pad, options.pad))
                    
                try:
                    q.write_image(options.outfile, i-1)
                except:
                    sys.stderr.write("Error: Cannot write to file %s" % options.outfile)
                    exit(1)

                phi = phi + step
                i = i + 1

                    
            az2 = az2 + h

        alt2 = alt2 + options.prop * math.pi / 180.
        j = j + 1


def load_args():
    parser=OptionParser(usage = "%prog <input file> [options]", version = "%prog 2.0"+chr(223))
    parser.add_option("--out", dest = "outfile", default = "e2proj.img",
                      help = "Output file. Default is 'e2proj.img'")
    parser.add_option("--projector", dest = "projector", default = "standard",
                      help = "Projector to use")
    parser.add_option("--prop", dest = "prop", type = "float", default = 0.0,
                      help = "Generates projections with a relatively uniform projection density in the unit triangle.")
    parser.add_option("--prop_start", dest = "prop_start", type = "float", default = 0.0,
                      help = "Start projections at the specified alt (in degrees). Default=0.0")
    parser.add_option("--prop_end", dest = "prop_end", type = "float", default = 90.0,
                      help = "Create projections up to the specified alt (in degrees). Default=90.0")
    parser.add_option("--euler", dest = "euler", type = "float", default = [],
                      action = "store", nargs = 3, metavar = "<az> <alt> <phi>",
                      help = "Generate a single projection with the given orientation")
    parser.add_option("--sym", dest = "sym",
                      help = "Set the symmetry; choices are: c<n>, d<n>, h<n>, i, t, icos, or oct")
    parser.add_option("--mode", dest = "mode", type = "int", default = 2,
                      help = "Default is real-space projection, this specifies various Fourier modes")
    parser.add_option("--angletype", dest = "angletype", default = "EMAN",
                      help = "Angle convention to use: [EMAN, SPIDER].  EMAN is the default")
    parser.add_option("--dump_angles", dest = "dump_angles", default = False,
                      action = "store_true",
                      help = "Dumps Euler angles to a text file")
    parser.add_option("--mask", dest = "mask", type = "float", default = 0.0,
                      help = "Specify a circular mask radius for the projections")
    parser.add_option("--randomphi", dest = "randomphi", default = False, action = "store_true",
                      help = "Randomize phi")
    parser.add_option("--phicomp", dest = "phic", default = False, action = "store_true",
                      help = "Roughly compensate for in-plane rotation from az with phi")
    parser.add_option("--phitoo", dest="phitoo", type = "float", default = 0,
                      help = "This will also vary phi in the final file. Warning: This works only with '--prop=' and generates a LOT of projections.")    
    parser.add_option("--smear", dest = "smear", default = False, action = "store_true",
                      help="Used in conjunction with '--phitoo=', this will rotationally smear between phi steps.")
    parser.add_option("--grid", dest = "grid", type = "int", default = 0.0,
                      help="Generate projections on an rectangular alt/az mesh, nonuniform projection density")
    parser.add_option("--angle_input_file", dest = "user_angles", metavar = "FILE_NAME",
                      help="Us an input file to specify what angles to project")
    parser.add_option("--pad", dest = "pad", type = "int", default = 0,
                      help="Pad image")


    (opt,args)=parser.parse_args()

    if len(args) < 1:
        sys.stderr.write("ERROR: No input file given\n")
        exit(1)

    opt.angletype = opt.angletype.upper()

    
    try:                             # check for valid projector type
        Projectors.get(opt.projector)
    except:
        sys.stderr.write("ERROR: %s is not a valid projector - use one of the following\n"%
                         opt.projector)
        dump_projectors()
        exit(1)
        

    if not(opt.sym):             # remove this once unknown symmetry is handled
        sys.stderr.write("ERROR: No symmetry type given\n");
        exit(1)


    opt.sym=opt.sym.lower()           # check for valid symmetry type
    if (opt.sym[0] in ["c","d"]):
        if not(opt.sym[1:].isdigit()):
            sys.stderr.write("ERROR: %s is an invalid symmetry type\n"%opt.sym)


    # test for mutual exclusion   
    if (opt.randomphi and opt.phic):
        sys.stderr.write("WARNING: --phicomp and --randomphi are mutually exclusive\n")
    if (opt.grid and opt.prop):  #doesn't work right anymore
        sys.stderr.write("ERROR: --grid, --prop, --angle_input_file are mutually exclusive\n")
        exit(1)

    return (opt, args)


def read_angles(filename, angletype):
    
    if (angletype == "EMAN"):
        words = ["alt", "az", "phi"]
    elif (angletype == "SPIDER"):
        words = ["phi", "psi", "theta"]
    else:
        sys.stderr.write("ERROR: unsupported angle type - use EMAN or SPIDER\n")
        exit(1)

    try:
        infile = open(filename)
        informat = infile.readline()
    except:
        sys.stderr.write("ERROR: could not read file: %s\n" % filename)
        exit(1)

    informat = informat.lower()
    informat = informat.split()
    for x in words:
        if (informat.count(x) > 1):
            sys.stderr.write("ERROR: %s is listed %d times in the angle-input file format string\n"
                             % (x, informat.count(x)))
            exit(1)
        if (informat.count(x) < 1):
            sys.stderr.write("ERROR: %s is no listed in the angle-input file format string\n"
                             % x)
            exit(1)

    spacing = {}
    for x in range(len(informat)):
        if (informat[x] in words):
            spacing[informat[x]] = x

    for x in infile.readlines():
        if x == "\n":
            continue
        values = x.split()
        for y in spacing:
            exec(y+"="+values[spacing[y]]) in globals()
        yield 1

    
    try:
        infile.close()
    except:
        sys.stderr.write("ERROR: could not properly close file: %s\n" % filename)

        

if __name__=="__main__":
    main()
