#!/usr/bin/env python

from optparse import OptionParser
import os
from sys import exit


def myexec(cmd):
	print  "	 ", cmd
	r = os.system(cmd)
	if r != 0:
		print "Command execution failed!"
		print "If it failed at wget or curl due to no Internet connection, try download the file from other machine,",
		print "copy it to the current directory and restart install_mpi.py.",
		print "Otherwise, check the log file and try to resolve it."
		sys.exit(-1)

def chdir(dir):
	print  "	  cd ", dir
	os.chdir(dir)

def get_mpiroot(options):
	print "Checking mpicc"
		
	r = os.system("mpicc --version")
	if r != 0:
		print "Cannot find mpicc"
		return False
		
	return True

default_version_of_open_mpi_to_istall = "1.10.2"

parser = OptionParser()
parser.add_option("--openmpi_ver", type="string",  default=default_version_of_open_mpi_to_istall, help="version of openmpi to forcefully install, default = %s"%default_version_of_open_mpi_to_istall)
options,args = parser.parse_args()

if not get_mpiroot(options):
	print "You need MPI environment (both runtime and developer packages) and gcc compiler to continue. "
	print "If you work on professional HPC cluster, in all likelihood both are already installed. "
	print "In this case read the user guide - you have to probably load appriopriate module by \"module load\" command."
	print "You can also run this script again with the --force option - it will download and install MPI (openmpi-%s) for you."%default_version_of_open_mpi_to_istall
	exit(-1)

print ""
print "=====> Configuring the mpi python binding"

myexec("./configure --prefix=%s" % os.environ['SP_DIR'])

print ""
print "=====> Building the mpi python binding"
myexec("make clean >> log.txt")	
myexec("sed -i.bak 's/\(^LDFLAGS.*$\)/\\1 -lfftw3_mpi -lfftw3/' src/Makefile")
myexec("make all >> log.txt")


print ""
print "=====> Install the mpi python binding"
myexec("make install >> log.txt")	
