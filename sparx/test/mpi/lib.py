
def myexec( cmd ):
    import os
    import sys
    print  "         ", cmd
    r = os.system( cmd + ">& install.log" )
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
    


