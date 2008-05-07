#!/usr/bin/env python
# This is designed for OSX Leopard. The URL links are up to date as of 5/6/2008
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
"sip":"http://www.riverbankcomputing.co.uk/static/Downloads/sip4/sip-4.7.4.tar.gz",
"pyqt":"http://www.riverbankcomputing.co.uk/static/Downloads/PyQt4/PyQt-mac-gpl-4.3.3.tar.gz",
"fftw":"http://www.fftw.org/fftw-3.1.2.tar.gz",
"gsl":"http://mirror.anl.gov/pub/gnu/gsl/gsl-1.11.tar.gz",
"jam":"http://internap.dl.sourceforge.net/sourceforge/boost/boost-jam-3.1.16.tgz",
"boost":"http://internap.dl.sourceforge.net/sourceforge/boost/boost_1_34_1.tar.bz2",
"cmake":"http://www.cmake.org/files/v2.4/cmake-2.4.8.tar.gz",
"hdf5":"ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-1.8.0.tar.gz",
"tiff":"ftp://ftp.remotesensing.org/pub/libtiff/tiff-3.8.2.tar.gz",
"pyopengl":"http://superb-east.dl.sourceforge.net/sourceforge/pyopengl/PyOpenGL-3.0.0b2.tar.gz"}

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
if not access("/usr/local/lib/libfftw3f.3.dylib",R_OK) :
	system("tar xvzf "+fsp["fftw"])
	system("cd %s; ./configure --enable-float --enable-shared --prefix=/usr/local; make; make install"%fsp["fftw"][:-7])

# GSL
if not access("/usr/local/lib/libgsl.dylib",R_OK):
	system("tar xvzf "+fsp["gsl"])
	system("cd %s; ./configure --prefix=/usr/local; make; make install"%fsp["gsl"][:-7])

# jpg
if not access("/usr/local/lib/libjpeg.a",R_OK):
	system("tar xvzf "+fsp["jpeg"])
	system("cd jpeg-6b; cp /usr/share/libtool/config.sub .; cp /usr/share/libtool/config.guess .; ./configure --enable-shared --enable-static --prefix=/usr/local; make; make install;ranlib /usr/local/lib/libjpeg.a")

# tiff
if not access("/usr/local/lib/libtiff.dylib",R_OK):
	system("tar xvzf "+fsp["tiff"])
	system("cd %s; ./configure --prefix=/usr/local; make; make install"%fsp["tiff"][:-7])

# boost
if system("which bjam") :
	system("tar xvzf "+fsp["jam"])
	system("cd %s; ./build.sh; cp bin.macosxx86/bjam /usr/local/bin/"%fsp["jam"][:-4])
	
	system("tar xvjf "+fsp["boost"])
	system("cd %s; bjam --toolset=darwin install"%fsp["boost"][:-8])

# HDF5
if not access("/usr/local/lib/libhdf5.dylib",R_OK):
	system("tar xvzf "+fsp["hdf5"])
	system("cd %s; ./configure --prefix=/usr/local --with-default-api-version=v16; make; make install"%fsp["hdf5"][:-7])

# cmake
if not access("/usr/local/bin/cmake",R_OK):
	system("tar xvzf "+fsp["cmake"])
	system("cd %s; ./configure --prefix=/usr/local; make; make install"%fsp["cmake"][:-7])

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
if not access("/usr/local/bin/gpg",R_OK) :
	system("tar xvjf "+fsp["gpg"])
	system("cd %s; ./configure --prefix=/usr/local; make; make install"%fsp["gpg"][:-8])

