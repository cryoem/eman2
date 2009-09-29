#!/usr/bin/env python
# This was designed to install most of the EMAN2 dependencies from source in a specified
# target directory. If you don't have root permissions, you can specify your home directory,
# but you will need to install python in your home directory before doing this. The URL links are up to date as of 5/6/2008
# In theory if one of the URL's is out of date, you should only need to replace
# the URL. The rest of the script should still function properly. Note that this
# script installs dependencies, but not EMAN2 itself. Also, it is safe to rerun
# this script, as it will do basic checks to avoid unnecessary recompilations.
#
# You will need to install Qt-Mac:
#
# ftp://ftp.trolltech.com/qt/source/qt-mac-opensource-4.4.0.dmg
# 
# Manually before running this script

from urllib import urlopen
from os import *

ftg={ "setup":"http://peak.telecommunity.com/dist/ez_setup.py",
#"png":"http://superb-east.dl.sourceforge.net/sourceforge/libpng/libpng-1.2.28.tar.bz2",
"jpeg":"http://www.ijg.org/files/jpegsrc.v6b.tar.gz",
"gpg":"ftp://ftp.gnupg.org/gcrypt/gnupg/gnupg-1.4.9.tar.bz2",
"sip":"http://www.riverbankcomputing.co.uk/static/Downloads/sip4/sip-4.9.tar.gz",
#"pyqt":"http://www.riverbankcomputing.co.uk/static/Downloads/PyQt4/PyQt-mac-gpl-4.6.tar.gz",
"pyqt":"http://www.riverbankcomputing.co.uk/static/Downloads/PyQt4/PyQt-x11-gpl-4.6.tar.gz",
"fftw":"http://www.fftw.org/fftw-3.2.1.tar.gz",
"gsl":"http://mirror.anl.gov/pub/gnu/gsl/gsl-1.12.tar.gz",
"jam":"http://internap.dl.sourceforge.net/sourceforge/boost/boost-jam-3.1.17.tgz",
"boost":"http://internap.dl.sourceforge.net/sourceforge/boost/boost_1_39_0.tar.bz2",
"cmake":"http://www.cmake.org/files/v2.6/cmake-2.6.4.tar.gz",
"hdf5":"ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-1.8.3.tar.gz",
"tiff":"ftp://ftp.remotesensing.org/pub/libtiff/tiff-3.8.2.tar.gz",
"pyopengl":"http://downloads.sourceforge.net/project/pyopengl/PyOpenGL/3.0.0/PyOpenGL-3.0.0.tar.gz"}

from sys import argv
if len(argv)<2 : 
	print """Please specify an installation prefix location. On a mac, /usr/local is good (assuming
you have administrative permissions. If you lack root access, use your home directory."""
prefix=argv[1]

fsp={}
for i in ftg: fsp[i]=ftg[i].split("/")[-1]

path=getenv("HOME")+"/EMAN2/src"
try: makedirs(path)
except: pass
system("chown -R %s ~/EMAN2"%getenv("SUDO_USER"))

chdir(path)
print "Running in ",path

for i in ftg:
	if not access(fsp[i],R_OK) : 
		print "Retrieving ",i
		file(fsp[i],"w").write(urlopen(ftg[i]).read())
	else: print "Already have ",i

# easy setup
if system("which ipython") :
	system("python ez_setup.py")
	system("easy_install ipython")
	
# fftw
#if not access("%s/lib/libfftw3f.3.dylib"%prefix,R_OK) :
if not access("%s/lib/libfftw3f.3.so"%prefix,R_OK) :
	system("tar xvzf "+fsp["fftw"])
	system("cd %s; ./configure --enable-float --enable-shared --prefix=%s; make; make install"%(fsp["fftw"][:-7],prefix))

# GSL
#if not access("%s/lib/libgsl.dylib"%prefix,R_OK):
if not access("%s/lib/libgsl.so"%prefix,R_OK):
	system("tar xvzf "+fsp["gsl"])
	system("cd %s; ./configure --prefix=%s; make; make install"%(fsp["gsl"][:-7],prefix))

# jpg
if not access("%s/lib/libjpeg.a"%prefix,R_OK):
	system("tar xvzf "+fsp["jpeg"])
	system("cd jpeg-6b; cp /usr/share/libtool/config.sub .; cp /usr/share/libtool/config.guess .; ./configure --enable-shared --enable-static --prefix=%s; make; make install;ranlib %s/lib/libjpeg.a"%(prefix,prefix))

# tiff
#if not access("%s/lib/libtiff.dylib"%prefix,R_OK):
if not access("%s/lib/libtiff.so"%prefix,R_OK):
	system("tar xvzf "+fsp["tiff"])
	system("cd %s; ./configure --prefix=%s; make; make install"%(fsp["tiff"][:-7],prefix))

# boost
# this version is for OSX
#if system("which bjam") :
	#system("tar xvzf "+fsp["jam"])
	#system("cd %s; ./build.sh; cp bin.macosxx86/bjam %s/bin/"%(fsp["jam"][:-4],prefix))
#	
	#system("tar xvjf "+fsp["boost"])
	#system("cd %s; bjam --toolset=darwin install"%fsp["boost"][:-8])

# This version is for Linux x86_64
if system("which bjam") :
	system("tar xvzf "+fsp["jam"])
	system("cd %s; ./build.sh; cp bin.linuxx86_64//bjam %s/bin/"%(fsp["jam"][:-4],prefix))
	
	system("tar xvjf "+fsp["boost"])
	system("cd %s; bjam --toolset=gcc install"%fsp["boost"][:-8])

# HDF5
#if not access("%s/lib/libhdf5.dylib"%prefix,R_OK):
if not access("%s/lib/libhdf5.so"%prefix,R_OK):
	system("tar xvzf "+fsp["hdf5"])
	system("cd %s; ./configure --prefix=%s --with-default-api-version=v16; make; make install"%(fsp["hdf5"][:-7],prefix))

# cmake
if not access("%s/bin/cmake"%prefix,R_OK):
	system("tar xvzf "+fsp["cmake"])
	system("cd %s; ./configure --prefix=%s; make; make install"%(fsp["cmake"][:-7],prefix))

# SIP
try: import sip
except:
	system("tar xvzf "+fsp["sip"])
	system("cd %s; python configure.py; make; make install"%fsp["sip"][:-7])

# PyQt
try: import PyQt4
except:
	system("tar xvzf "+fsp["pyqt"])
	system("cd %s; echo 'yes' | python configure.py; make; make install"%fsp["pyqt"][:-7])

# PyOpenGL
try: import OpenGL
except:
	system("tar xvzf "+fsp["pyopengl"])
	system("cd %s; python setup.py install"%fsp["pyopengl"][:-7])

try: import matplotlib
except:
	system("easy_install matplotlib")

# GPG
if not access("%s/bin/gpg"%prefix,R_OK) :
	system("tar xvjf "+fsp["gpg"])
	system("cd %s; ./configure --prefix=%s; make; make install"%(fsp["gpg"][:-8],prefix))

