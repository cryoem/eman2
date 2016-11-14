#!/usr/bin/env python
# Ian Rees, 2013

# DEPENDENCY ORDERS:
#   1) python -> boost, and other python packages
#   2) db -> bsddb3
#   3) freetype -> ftgl
#   4) qt4 -> pyqt4, matplotlib
#   5) tiff, libpng -> matplotlib
#   6) everything else

import sys
import os
import re
import shutil
import subprocess
import glob
import datetime
import argparse

from build import *

def cmd(*popenargs, **kwargs):
    print "Running:", 
    print " ".join(*popenargs)
    process = subprocess.Popen(*popenargs, **kwargs)
    process.wait()  
    
def join(*path):
    # Fix issues with Windows path escaping.
    p = os.path.join(*path)
    return p.replace("\\","/")

##### Base Dependencies #####

class Dependency(object):
    priority = None
    registered = {}
    
    def __init__(self, args, src):
        self.args = args
        self.src = src
        self.package = self._get_package(src)
        self.package_dir = join(os.getcwd(), self.package)
        
    def _get_package(self, src):
        packagename = os.path.basename(self.src)
        removeexts = [".tar", ".tgz", ".gz", ".bz2", ".zip"]
        for i in removeexts:
            packagename = packagename.replace(i, "")
        return packagename
        
    @classmethod
    def register(cls, name):
        def f(o):
            cls.registered[name] = o
            return o
        return f
        
    def untar(self):
        print "Retree: '%s' "%self.package_dir
        retree(self.package_dir)
        args = ['tar', '-x', '-v']
        if self.src.endswith("gz"):
            args.append("-z")
        elif self.src.endswith("bz2"):
            args.append("-j")
        args.extend(["-f", self.src])
        args.extend(["-C", self.package_dir])
        args.extend(['--strip-components', '1'])
        cmd(args)
        
    def configure(self):
        args = []
        args += self.configure_script()
        args += self.configure_prefix()
        args += self.configure_args()
        cmd(args, cwd=self.package_dir)
    
    def configure_script(self):
        # Windows is not bright -- needs 'bash'
        args = ['bash', join(self.package_dir, 'configure')]
        return args
        
    def configure_prefix(self):
        if self.args.prefix:
            return ['--prefix=%s'%self.args.prefix]
        return []

    def configure_args(self):
        return []
        
    def build(self):
        pass
    
    def install(self):
        pass

    def postinstall(self):
        pass

class AutoConfDependency(Dependency):
    def configure_args(self):
        args = []
        args += ['--enable-shared']
        # if 'apple' in self.args.target:
        args += ['--disable-dependency-tracking']
        #args += ['--disable-libtool-lock']
        return args

    def build(self):
        args = ['make']
        cmd(args, cwd=self.package_dir)
    
    def install(self):
        args = ['make']
        args += ['install']
        cmd(args, cwd=self.package_dir)
        
class DistUtilsDependency(Dependency):
    def configure(self):
        args = [sys.executable]
        args += [join(self.package_dir, 'setup.py')]
        args += ['build', 'install']
        args += self.configure_prefix()
        args += self.configure_args()
        cmd(args, cwd=self.package_dir)
    
    def configure_prefix(self):
        args = []
        if self.args.prefix:
            args += ['--root', self.args.prefix]
            args += ['--install-lib', 'site-packages']
            args += ['--install-headers', 'include']
            args += ['--install-scripts', 'bin']
            args += ['--install-data', 'share']
        return args

##### Unique Build System Dependencies #####

@Dependency.register('^boost')
class dep_boost(Dependency):
    # On Windows, you need to run boostrap.bat manually, and
    # patch: 
    def configure_script(self):
        # Windows is not bright -- needs 'bash'
        args = ['bash', join(self.package_dir, 'bootstrap.sh')]
        return args
        
    def configure_args(self):
        args = ['--with-libraries=python'] # ,system,filesystem,thread
        if 'linux' in self.args.target:
            args += ['--with-python=%s'%join(self.args.prefix, 'python', 'bin', 'python')]
        # elif 'darwin' in self.args.target:
        #    args += ['toolset=darwin', 'architecture=x86', 'address-model=32_64'] # Mac
        return args
        
    def build(self):
        args = [join(self.package_dir, 'b2')]
        args += self.configure_prefix()
        if 'darwin' in self.args.target:
            args += ['toolset=darwin', 'architecture=x86', 'address-model=32_64'] # Mac        
        cmd(args, cwd=self.package_dir)
        
    def install(self):
        args = [join(self.package_dir, 'b2')]
        args += self.configure_prefix()
        args += ['install']
        if 'darwin' in self.args.target:
            args += ['toolset=darwin', 'architecture=x86', 'address-model=32_64'] # Mac        
        cmd(args, cwd=self.package_dir)

@Dependency.register('^db')
class dep_db(Dependency):
    def configure(self):
        # We have to run in a different directory...
        cwd = join(self.package_dir, 'build_unix')
        args = ['bash', join(self.package_dir, 'dist', 'configure')]
        args += self.configure_prefix()
        args += self.configure_args()
        cmd(args, cwd=cwd)

    def configure_args(self):
        if 'mingw' in self.args.target:
            return ['--enable-mingw']
        return []
        
    def build(self):
        cwd = join(self.package_dir, 'build_unix')
        args = ['make']
        cmd(args, cwd=cwd)
    
    def install(self):
        cwd = join(self.package_dir, 'build_unix')    
        args = ['make']
        args += ['install']
        cmd(args, cwd=cwd)

@Dependency.register('^qt-everywhere-opensource')
class dep_qt4(AutoConfDependency):
    def configure_script(self):
        if 'mingw' in self.args.target:
            return [join(self.package_dir, 'configure.exe')]
        else:
            return [join(self.package_dir, 'configure')]            
            
    def configure_args(self):
        # Don't build debug libraries
        args = []
        args += ['-release']
        args += ['-opensource']
        args += ['--confirm-license']
        
        # Exclude various items we don't need to save time
        args += ['-no-libtiff', '-no-libmng', '-no-webkit', '-no-phonon', '-no-multimedia', '-no-script', '-no-xmlpatterns']
        args += ['-nomake','examples', '-nomake','demos', '-nomake','docs', '-nomake', 'translations']
        
        # Use built in plugins to save complexity
        args += ['-qt-libpng', '-qt-libjpeg']

        if 'linux' in self.args.target:
            pass
        elif 'darwin' in self.args.target:
            args += ['-cocoa', '-framework', '-arch','i386', '-arch', 'x86_64']
        return args    
    
    def configure_prefix(self):
        if self.args.prefix:
            return ['-prefix', join(self.args.prefix, 'qt4')]
        return []

class SIPDependency(Dependency):
    def configure_script(self):
        return [sys.executable, join(self.package_dir, 'configure.py')]
    
    def build(self):
        args = ['make']
        cmd(args, cwd=self.package_dir)
    
    def install(self):
        args = ['make']
        args += ['install']
        cmd(args, cwd=self.package_dir)            

@Dependency.register('^sip')    
class dep_sip(SIPDependency):
    def configure_prefix(self):
        args = []
        if self.args.prefix:
            args += ['--bindir', join(self.args.prefix, 'bin')]
            args += ['--destdir', join(self.args.prefix, 'site-packages')]
            args += ['--incdir', join(self.args.prefix, 'include')]
            args += ['--sipdir', join(self.args.prefix, 'share', 'sip')]
        return args
    
    def configure_args(self):
        args = []
        if 'apple' in self.args.target:
            args += ['--arch=i386', '--arch=x86_64']
        if 'darwin10' in self.args.target:
            args += ['--sdk=MacOSX10.6.sdk']
        elif 'darwin11' in self.args.target:
            args += ['--sdk=MacOSX10.7.sdk']            
        if 'mingw' in self.args.target:
            args += ['-p', 'win32-g++']
        return args

@Dependency.register('^PyQt')    
class dep_pyqt4(SIPDependency):
    def configure_prefix(self):
        # NOTE: sip and pyqt4 have slightly different prefix-related options,
        # so this method is mostly-the-same-but-different.
        args = []
        if self.args.prefix:
            args += ['-q', join(self.args.prefix, 'qt4','bin', 'qmake')]
            args += ['--bindir', join(self.args.prefix, 'bin')]
            args += ['--destdir', join(self.args.prefix, 'site-packages')]
            args += ['--sipdir', join(self.args.prefix, 'share', 'sip')]
        return args

    def configure_args(self):
        args = []
        args += ['--confirm-license']
        args += ['--no-designer-plugin']
        args += ['-e', 'QtGui', '-e', 'QtCore', '-e', 'QtOpenGL']
        return args

##### AutoConf Dependencies #####

@Dependency.register('^cmake')
class dep_cmake(AutoConfDependency):
    def configure_args(self):
        # CMake bootstrap doesnt' like options.
        return []

@Dependency.register('^tiff')
class dep_tiff(AutoConfDependency):
    pass

@Dependency.register('^jpegsrc')
class dep_jpeg(AutoConfDependency):
    pass    

@Dependency.register('^libpng')
class dep_png(AutoConfDependency):
    pass

@Dependency.register('^gsl')
class dep_gsl(AutoConfDependency):
    pass

@Dependency.register('^hdf5')
class dep_hdf5(AutoConfDependency):
    pass

@Dependency.register('^freetype')
class dep_freetype(AutoConfDependency):
    pass

@Dependency.register('^szip')
class dep_szip(AutoConfDependency):
    pass

@Dependency.register('^zlib')
class dep_zlib(AutoConfDependency):
    # Needed by Windows. Needs manual build:
    # make -f win32/Makefile.gcc
    # make -f win32/Makefile.gcc install SHARED_MODE=1 BINARY_PATH=$PREFIX/bin LIBRARY_PATH=$PREFIX/lib INCLUDE_PATH=$PREFIX/include 
    pass

@Dependency.register('^fftw-3')
class dep_fftw3(AutoConfDependency):
    def configure_args(self):
        args = super(dep_fftw3, self).configure_args()
        args += ['--enable-float', '--enable-sse2']
        if 'mingw' in self.args.target:
            args += ['--with-our-malloc16', '--with-windows-f77-mangling']
        return args        

@Dependency.register('^ftgl')
class dep_ftgl(AutoConfDependency):
    def configure_args(self):
        args = super(dep_ftgl, self).configure_args()
        if self.args.prefix:
            args += ['--with-ft-prefix=%s'%self.args.prefix]
        return args

@Dependency.register('^Python')        
class dep_python(AutoConfDependency):
	# Argh, note, on Linux I need to set CPPFLAGS=CFLAGS and LD_LIBRARY_PATH=$PREFIX/lib
	# to find OpenSSL.
    def configure_prefix(self):
        if self.args.prefix:
            return ['--prefix', join(self.args.prefix, 'python')]
        return []
        
    def configure_args(self):
        return ['--enable-unicode=ucs4']

##### DistUtils Packages #####

@Dependency.register('^ipython')
class dep_ipython(DistUtilsDependency):
    pass
    
@Dependency.register('^readline')
class dep_readline(DistUtilsDependency):
    pass
    
@Dependency.register('^PyOpenGL')        
class dep_pyopengl(DistUtilsDependency):
    pass

@Dependency.register('^numpy')
class dep_numpy(DistUtilsDependency):
    pass
    
@Dependency.register('^bsddb3')
class dep_bsddb3(DistUtilsDependency):
    def configure_args(self):
        if self.args.prefix:
            return ['--berkeley-db=%s'%self.args.prefix]
        return []
    
@Dependency.register('^matplotlib')    
class dep_matplotlib(DistUtilsDependency):
    pass
    
##### For EMEN2 #####

@Dependency.register('^setuptools')
class dep_setuptools(DistUtilsDependency):
    pass

@Dependency.register('^Twisted')
class dep_twisted(DistUtilsDependency):
    pass    

@Dependency.register('^Markdown')
class dep_markdown(DistUtilsDependency):
    pass
    
@Dependency.register('^Mako')    
class dep_mako(DistUtilsDependency):
    pass
    
@Dependency.register('^jsonrpc')
class dep_jsonrpc(DistUtilsDependency):
    pass

@Dependency.register('^python-dateutil')
class dep_dateutil(DistUtilsDependency):
    pass

@Dependency.register('^argparse')
class dep_dateutil(DistUtilsDependency):
    pass

@Dependency.register('^py-bcrypt')
class dep_bcrypt(DistUtilsDependency):
	pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('packages',    help='Build packages', nargs='+')
    parser.add_argument('--target',    help='Build target. Checks $TARGET: %s'%os.getenv('TARGET'), default=os.getenv('TARGET') or 'linux')
    parser.add_argument('--prefix',    help="Installation prefix. Checks $PREFIX: %s"%os.getenv('PREFIX'), default=os.getenv('PREFIX'))
    parser.add_argument('--skip',  action="append", help="Skip build step: untar, configure, build, install, postinstall")
    args = parser.parse_args()
    args.skip = args.skip or []

    build = []
    for package in args.packages:
        found = False
        for k,v in Dependency.registered.items():
            if re.match(k, os.path.basename(package)):
                # print "Found builder for %s: %s"%(package, v)
                found = True
                build.append(v(args, package))
        if not found:
            print "Error: Could not find builder for %s: %s"%(package, v) 

    for b in build:
        print "===== Installing: %s ====="%b.src        
        if 'untar' not in args.skip:
            b.untar()
        if 'configure' not in args.skip:
            b.configure()
        if 'build' not in args.skip:
            b.build()
        if 'install' not in args.skip:
            b.install()
        if 'postinstall' not in args.skip:
            b.postinstall()



