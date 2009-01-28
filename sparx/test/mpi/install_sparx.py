#!/usr/bin/env python

from os import environ
from lib import *

def get_pyver():
	import commands
	from string import split
	r = commands.getoutput( "python -V" )
	r = split( r )[1]
	return r[0:3], r[3:]



ROOT = environ["HOME"] + "/EMAN2"
PYVER,PYSUB = get_pyver()
environ["PYTHONPATH"] = ROOT + "/lib/python" + PYVER + "/site-packages"

info = {}

info["fftw"] = {"key":"/include/fftw3.h", 
                "url":"http://www.fftw.org/fftw-3.2.tar.gz",
                "src":"fftw-3.2.tar.gz",
                "dir":"fftw-3.2" }

info["gsl"]  = {"key":"/include/gsl",
		"url":"ftp://ftp.gnu.org/gnu/gsl/gsl-1.10.tar.gz",
		"src":"gsl-1.10.tar.gz",
		"dir":"gsl-1.10"}

info["boost"]= {"key":"/include/boost-1_36",
		"url":"http://superb-west.dl.sourceforge.net/sourceforge/boost/boost_1_36_0.tar.gz",
		"src":"boost_1_36_0.tar.gz",
		"dir":"boost_1_36_0"}

info["hdf5"] = {"key":"/include/hdf5.h",
		"url":"ftp://ftp.hdfgroup.org/HDF5/current16/src/hdf5-1.6.8.tar.gz",
		"src":"hdf5-1.6.8.tar.gz",
		"dir":"hdf5-1.6.8"}

info["db4"]  = {"key":"/include/db4.h",
		"url":"http://download.oracle.com/berkeley-db/db-4.7.25.NC.tar.gz",
		"src":"db-4.7.25.NC.tar.gz",
		"dir":"db-4.7.25.NC/build_unix"}

info["cmake"]= {"key":"/bin/cmake",
		"url":"http://www.cmake.org/files/v2.6/cmake-2.6.2.tar.gz",
		"src":"cmake-2.6.2.tar.gz",
		"dir":"cmake-2.6.2"}

info["python"] = {"key":"/include/python%s/Python.h"%PYVER,
		"url":"http://www.python.org/ftp/python/%s%s/Python-%s%s.tgz" %(PYVER,PYSUB,PYVER,PYSUB),
		"src":"Python-%s%s.tgz" %(PYVER,PYSUB),
		"dir":"Python-%s%s" %(PYVER,PYSUB)}

info["numpy"]= {"key":"numpy",
		"url":"http://superb-west.dl.sourceforge.net/sourceforge/numpy/numpy-1.0.4.tar.gz",
		"src":"numpy-1.0.4.tar.gz",
		"dir":"numpy-1.0.4"}

info["setuptools"] =     {"key":"setuptools.pth",
		"url":"http://pypi.python.org/packages/source/s/setuptools/setuptools-0.6c8.tar.gz",
		"src":"setuptools-0.6c8.tar.gz",
		"dir":"setuptools-0.6c8"}

info["bdb3"] = {"key":"bsddb3-4.7.2-py2.4-linux-i686.egg",
		"url":"http://pypi.python.org/packages/source/b/bsddb3/bsddb3-4.7.2.tar.gz",
		"src":"bsddb3-4.7.2.tar.gz",
		"dir":"bsddb3-4.7.2"}

info["ipython"] = {"key":"IPython",
		"url":"http://ipython.scipy.org/dist/ipython-0.9.1.tar.gz",
		"src":"ipython-0.9.1.tar.gz",
		"dir":"ipython-0.9.1"}


def install_clib( name, opts=None, confdir="." ):
	import os

	dict = info[name]

	if os.path.exists( ROOT+dict["key"] ):
		return

	print "Installing ", name
	geturl( dict["url"], dict["src"] )
	myexec( "tar -zxvf " + dict["src"] )
	chdir( dict["dir"] )

	config = confdir + "/configure --prefix=" + ROOT 
	if not(opts is None):
		config += opts;

	myexec( config )
	myexec( "make" )
	myexec( "make install" )
	chdir( ".." )


def install_plib( name ):
	import os

	dict = info[name]
	if os.path.exists( ROOT+"/lib/python"+PYVER+"/site-packages/"+dict["key"] ):
		return

	print "Installing ", name
	geturl( dict["url"], dict["src"] )
	myexec( "tar -zxvf " + dict["src"] )
	chdir( dict["dir"] )

	install = "python setup.py install --prefix=" + ROOT 
	myexec( install )
	chdir( ".." )


install_clib( "python", "--enable-shared" )
install_clib( "fftw", "--enable-float --enable-shared" )
install_clib( "gsl" )
install_clib( "boost" )
install_clib( "db4" )
install_clib( "hdf5" )
install_clib( "cmake" )
install_plib( "numpy" )
install_plib( "setuptools" )
install_plib( "bdb3" )
install_plib( "ipython" )
