#!/usr/bin/env python

def myexec( cmd ):
    import os
    print  "         ", cmd
    r = os.system( cmd )
    if r != 0:
        print "Command: "
        print "         ", cmd 
        print "Failed!"

        print "if it is due to no internet connection, try download the file from other machine,"
        print "copy it to the current directory and restart install_mpi.py."
        print "otherwise, check the file log and try to resolve it"
        os.exit(-1)
 


def get_pythonroot( options ) :
    import commands
    import os
    from string import split

    if options.pythonroot is None:
        result = commands.getoutput( "which python" )
        if result.find( "no python" ) != -1:
            print "Error: pythonroot was not given, and python exectuables cannot be found"
            return None,None
        
        pythonroot = os.path.dirname( os.path.dirname(result) )
    else:
        pythonroot = options.pythonroot


    python = pythonroot + "/bin/python"
    if not os.path.exists(python):
        print "Error: python executable ", python, " does not exist."
        return None,None

    pythonver = commands.getoutput( python + " -V" )
    pythonver = split( pythonver )[1]

    header = pythonroot + "/include/python" + pythonver[0:3] + "/Python.h"
    if not os.path.exists(header):
        print "Error: we cannot find file ", header, "."
        print "       This usually is because you did not install the python development package."
        print "       the python itself, you have to install the python development package"
        print "       such as python-devel."
        return None,pythonver

    return pythonroot,pythonver



def get_mpiroot(options):
    import os
    import commands
    if options.mpiroot is None:
        result = commands.getoutput( "which mpirun" )
        if result.find( " " ) != -1 :
            print "Error: mpiroot was not given, and mpirun exectuable cannot be found"
            print "       install mpi first. or use the --force the program will install"
            print "       openmpi for you"
            return None

        mpiroot = os.path.dirname( os.path.dirname(result) )
    else:
        mpiroot = options.mpiroot

    header = mpiroot + "/include/mpi.h"
    if not os.path.exists(header):
        print "Error: header ", mpiinc + "/mpi.h does not exist"
        print "This is usually because you did not install mpi development package"
    
    return mpiroot

def get_Numeric(pythonroot, pythonver):
 
    try:
        import Numeric
	return pythonroot 

    except Exception, inst:
        
        print "problem import Numeric:", inst
        return None

def install_python( pythonver ):
    import os
    root = os.environ["HOME"] + "/EMAN2"

    pwd = os.getcwd()

    if pythonver is None:
        pythonver = "2.5.2"
  
    print "Installing python ", pythonver
    file = "Python-" + pythonver + ".tgz"

    if not( os.path.exists(file) ):
        myexec( "wget http://www.python.org/ftp/python/" + file )
    
    myexec( "tar -zxf " + file )
    os.chdir( "Python-" + pythonver )
    myexec( "./configure --enable-shared --preifx=" + root + ">& log" )
    myexec( "make >& log" )
    myexec( "make install >& log" )

    os.chdir( pwd )
    return root,pythonver

def install_openmpi( ):
    import os
    root = os.environ["HOME"] + "/EMAN2"

    pwd = os.getcwd()

    file = "openmpi-1.2.8.tar.gz"

    print "Installing openmpi-1.2.8"
    if not( os.path.exists(file) ):
        myexec( "wget http://www.open-mpi.org/software/ompi/v1.2/downloads/" + file )
    
    myexec( "tar -zxf " + file )
    os.chdir( "openmpi-1.2.8" )
    myexec( "./configure --enable-shared --preifx=" + root + ">& log" )
    myexec( "make >& log" )
    myexec( "make install >& log" )
    os.chdir( pwd )
    return root

       
def install_numeric(pythonroot, pythonver):
    import os
    root = os.environ["HOME"] + "/EMAN2"
    python = pythonroot + "/bin/python"
   
    file = "Numeric-24.2.tar.gz"

    print "Installing Numeric-24.2"
    if not( os.path.exists(file) ):
        myexec( "wget http://internap.dl.sourceforge.net/sourceforge/numpy/" + file )

    myexec( "tar -zxf " + file )
    os.chdir( "Numeric-24.2" )

    if root==pythonroot:
        myexec( python + " setup.py install --prefix=" )
    else:
        myexec( python + " setup.py install --prefix=" + root )
        print "Info: please add " + root + "/lib/python" + pythonver[0:3] + "/site-packages/Numeric to you PYTHONPATH"
    os.chdir( pwd )
    return root




from optparse import OptionParser

usage = "install_mpi.py --cc=c_compiler --pythonroot --mpiroot --Numeric"
parser = OptionParser(usage)
parser.add_option( "--cc", type="string", default="gcc", help="the compiler, default is gcc" )
parser.add_option( "--pythonroot", type="string", help="the directory containing python header files" )
parser.add_option( "--mpiroot", type="string", help="the directory containing mpi header and mpi libraries" )
parser.add_option( "--Numeric", type="string", help="the directory containing Numeric header" )
parser.add_option( "--force", action="store_true", default=False, help="if can't find some packages, install it forcefully" )
options,args = parser.parse_args()


pythonroot,pythonver = get_pythonroot( options )
mpiroot = get_mpiroot( options )
numeric = get_Numeric( pythonroot, pythonver )
numeric = None
if pythonroot is None:
     if options.force:
         pythonroot,pythonver = install_python( pythonver )
     else:
         print "Error: cannot find python executables or headers, cannot proceed"
         exit(-1)


if mpiroot is None:
    if options.force:
        mpiroot = install_openmpi()
    else:
        print "Error: cannot find mpi header file, cannot proceed"
        exit(-2)

if numeric is None:
     if options.force:
         numeric = install_numeric(pythonroot, pythonver )
     else:
         print "Error: cannot find Numeric header file, cannot proceed"
         exit(-3)

import os
import commands

eman2 = os.environ["HOME"] + "/EMAN2"
if not os.path.exists( eman2 ):
    myexec( "mkdir " + eman2 )

if not os.path.exists( eman2 + "/src" ):
    myexec( "mkdir " + eman2 + "/src" )

if not os.path.exists( eman2 + "/lib" ):
    myexec( "mkdir " + eman2 + "/lib" )


print "Builing the mpi python binding"
pythoninc = pythonroot + "/include/python" + pythonver[0:3] 
pythonlib = "-lpython" + pythonver[0:3] + " -L" + pythonroot + "/lib/" + " -L" + pythonroot + "/lib64/"

mpiinc = mpiroot + "/include/"
mpilib = "-lmpi -L" + mpiroot + "/lib"

numericinc = numeric + "/include/python" + pythonver[0:3]

cmd = "gcc -c -fPIC mympimodule.c -I%s -I%s -I%s >& log" % ( pythoninc,  mpiinc, numericinc )
myexec(cmd)

uname = commands.getoutput( "uname" ) 
if uname == "Darwin":
    cmd = options.cc + " -dynamiclib -single_module -o mpi.so mympimodule.o >& log " + mpilib + " " + pythonlib
else:
    cmd = options.cc + " -shared -o mpi.so mympimodule.o >& log " + mpilib + " " + pythonlib

myexec(cmd)

cmd = "cp mpi.so " + eman2 + "/lib"
myexec(cmd)

print "Installation complete successfully!"
print "make sure you have %s in your PYTHONPATH" % ( eman2 + "/lib" )



