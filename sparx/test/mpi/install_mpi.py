#!/usr/bin/env python

def myexec( cmd ):
    import os
    import sys
    print  "         ", cmd
    r = os.system( cmd )
    if r != 0:
        print "Failed!"

        print "if it failed at wget or curl due to no internet connection, try download the file from other machine,"
        print "copy it to the current directory and restart install_mpi.py."
        print "otherwise, check the file log and try to resolve it"
        sys.exit(-1)

def chdir( dir ):
    import os
    print  "          cd ", dir
    os.chdir(dir)
 
def macos():
    import commands
    r = commands.getoutput( "uname" )
    return r=="Darwin"

def geturl( url, file ):
    import os
    import commands

    if macos():
	myexec( "curl " + url + " -o " + file + " >& log" )
    else:
        myexec( "wget " + url + " >& log" )
    


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
            print "       openmpi-1.2.8 for you"
            return None

        mpiroot = os.path.dirname( os.path.dirname(result) )
    else:
        mpiroot = options.mpiroot

    header = mpiroot + "/include/mpi.h"
    if not os.path.exists(header):
        print "Error: header ", mpiinc + "/mpi.h does not exist"
        print "       This is usually because you did not install mpi development package"
        print "       You can either install the mpi-devel package by yourself or use "
        print "       the --force option so that the program install openmpi-1.2.8 for you"

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

    print ""
    print "Installing python ", pythonver
    file = "Python-" + pythonver + ".tgz"
    url= "http://www.python.org/ftp/python/" + pythonver + "/" + file

    if not( os.path.exists(file) ):
        geturl( url, file )
    
    myexec( "tar -zxf " + file )
    chdir( "Python-" + pythonver )
    myexec( "./configure --enable-shared --prefix=" + root + ">& log" )
    myexec( "make >& log" )
    myexec( "make install >& log" )

    os.chdir( pwd )
    return root,pythonver

def install_openmpi( ):
    import os
    root = os.environ["HOME"] + "/EMAN2"

    pwd = os.getcwd()

    file = "openmpi-1.2.8.tar.gz"
    url = "http://www.open-mpi.org/software/ompi/v1.2/downloads/" + file

    print ""
    print "Installing openmpi-1.2.8"
    if not( os.path.exists(file) ):
        geturl( url, file )
    
    myexec( "tar -zxf " + file )
    chdir( "openmpi-1.2.8" )
    myexec( "./configure --enable-shared --prefix=" + root + ">& log" )
    myexec( "make >& log" )
    myexec( "make install >& log" )
    os.chdir( pwd )
    return root

       
def install_numeric(pythonroot, pythonver):
    import os
    root = os.environ["HOME"] + "/EMAN2"
    python = pythonroot + "/bin/python"
   
    pwd = os.getcwd()
    file = "Numeric-24.2.tar.gz"
    url = "http://internap.dl.sourceforge.net/sourceforge/numpy/" + file

    print ""
    print "Installing Numeric-24.2"
    if not( os.path.exists(file) ):
        geturl( url, file )

    myexec( "tar -zxf " + file )
    chdir( "Numeric-24.2" )

    if macos():
        myexec( "patch -p1 < ../Numeric_macos.patch >& log" )


    if root==pythonroot:
        myexec( python + " setup.py install >& log" )
    else:
        myexec( python + " setup.py install --prefix=" + root  + ">& log")

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

import os
import commands

eman2 = os.environ["HOME"] + "/EMAN2"
if not os.path.exists( eman2 ):
    myexec( "mkdir " + eman2 )

if not os.path.exists( eman2 + "/src" ):
    myexec( "mkdir " + eman2 + "/src" )

if not os.path.exists( eman2 + "/lib" ):
    myexec( "mkdir " + eman2 + "/lib" )


pythonroot,pythonver = get_pythonroot( options )
mpiroot = get_mpiroot( options )
numeric = get_Numeric( pythonroot, pythonver )
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
        numericpath = numeric + "/lib/python" + pythonver[0:3] + "/site-packages/Numeric"
    else:
        print "Error: cannot find Numeric header file, cannot proceed"
        exit(-3)
else:
    numericpath = None

print ""
print "Builing the mpi python binding"
pythoninc = pythonroot + "/include/python" + pythonver[0:3] 
pythonlib = "-lpython" + pythonver[0:3] + " -L" + pythonroot + "/lib/" + " -L" + pythonroot + "/lib64/"

mpiinc = mpiroot + "/include/"
mpilib = "-lmpi -L" + mpiroot + "/lib"

numericinc = numeric + "/include/python" + pythonver[0:3]

cmd = "gcc -c -fPIC mympimodule.c -I%s -I%s -I%s >& log" % ( pythoninc,  mpiinc, numericinc )
myexec(cmd)

uname = commands.getoutput( "uname" ) 
if macos():
    cmd = options.cc + " -dynamiclib -single_module -o mpi.so mympimodule.o >& log " + mpilib + " " + pythonlib
else:
    cmd = options.cc + " -shared -o mpi.so mympimodule.o >& log " + mpilib + " " + pythonlib

myexec(cmd)

e2lib = eman2 + "/lib"
cmd = "cp mpi.so " + e2lib
myexec(cmd)

print "Installation complete successfully!"

if macos():
    ldkey = "DYLD_LIBRARY_PATH"
else:
    ldkey = "LD_LIBRARY_PATH"

ldlibpath = None
if os.environ[ldkey].find(e2lib) == -1:
    ldlibpath = e2dir

pythonpath = None
if os.environ["PYTHONPATH"].find(e2lib) == -1:
    pythonpath = e2dir

if not(numericpath is None):
    if pythonpath is None:
        pythonpath = numericpath
    else:
        pythonpath = pythonpath + ":" + numericpath


if not(ldlibpath is None) or not(pythonpath is None):
    import os
    print ""
    print "Setting up your environment"

    print "  1. Bash user: the following line(s) has been added to ~/.bashrc"
    if not( ldlibpath is None) : 
        line = "        export %s=%s:$%s" % ( ldkey, ldlibpath, ldkey )
        print line
        os.system( "echo '%s' >> ~/.bashrc" % line )

    if not(pythonpath is None) : 
        line = "        export PYTHONPATH=%s:$PYTHONPATH" % (pythonpath)
        print line
        os.system( "echo '%s' >> ~/.bashrc" % line )

    print "     Please run"
    print "           source ~/.bashshrc"
    print ""

    print "  2. Csh  user: the following line(s) has been added to ~/.cshrc"
    if not( ldlibpath is None) : 
        line = "        setenv LD_LIBRARY_PATH %s:$LD_LIBRARY_PATH" % (ldkey, eman2+"/lib", ldkey)
        print line
        os.system( "echo '%s' >> ~/.cshrc" % line )
       
    if not(pythonpath is None) : 
        line = "        setenv PYTHONPATH %s:$PYTHONPATH" % (pythonpath)
        print line
        os.system( "echo '%s' >> ~/.cshrc" % line )

    print "     Please run"
    print "           source ~/.cshrc"
    print ""

