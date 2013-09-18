===== Overview =====

Using "/build", as the root, it will contain the following:
	/build/var.sh					Shell script to set various environment variables (TARGET, CFLAGS, PREFIX, LD_LIBRARY_PATH, etc.)
    /build/build.py                 The build script
	/build/src						Build area for dependencies
	/build/local					Install prefix for dependencies (e.g. --prefix=/build/local)

And the following directories with a subdirectory for each release, e.g. "eman2.pre2-1":
	/build/extlib					"Slimmed down" copy of /build/local with just necessities to run EMAN2. This gets copied to the distribution.
	/build/co						EMAN2 source code
	/build/build					EMAN2 CMake build directories
	/build/stage					EMAN2 install directories
	/build/images					Final .tar.gz and .dmg distributions that are uploaded to website.

===== Installing dependencies =====

The following dependencies are installed, in the following order.

1. Python (Windows, Linux)
    Python-2.7.3.tgz

2. Base
    cmake-2.8.10.2.tar.gz
    cmake-2.8.10.2-win32-x86.exe (Windows)
    boost_1_53_0.tar.gz
    fftw-3.3.3.tar.gz
    freetype-2.4.11.tar.gz
    ftgl-2.1.3-rc5.tar.gz
    gsl-1.15.tar.gz
    db-5.3.21.NC.tar

3. File formats
    zlib-1.2.7.tar.gz (Windows, Linux)
    szip-2.1.tar.gz
    hdf5-1.8.10-patch1.tar.gz
    libpng-1.6.1.tar.gz
    jpegsrc.v9.tar.gz
    tiff-4.0.3.tar.gz

4. Qt
    qt-everywhere-opensource-src-4.8.4.tar.gz
    sip-4.14.5.tar.gz
    PyQt-mac-gpl-4.10.tar.gz (Mac)
    PyQt-win-gpl-4.10.zip (Windows)
    PyQt-x11-gpl-4.10.tar.gz (Linux)

5. Python packages
    setuptools-0.6c11.tar.gz
    ipython-0.10.2.tar.gz
    numpy-1.7.0.tar.gz
    bsddb3-5.3.0.tar.gz
    matplotlib-1.2.1.tar.gz
    PyOpenGL-3.0.2.tar.gz
    readline-6.2.4.1.tar.gz

6. Additional Python (EMEN2)
    Markdown-2.3.1.tar.gz
    Mako-0.7.3.tar.gz
    Twisted-12.3.0.tar.bz2
    jsonrpc-1.2.tar.gz
    pyOpenSSL-0.13.tar.gz
    python-dateutil-1.5.tar.gz
    bcrypt

The script /build/shared/scripts/dep.py *MAY* help automate the process of compiling all of EMAN2's dependencies. 

    /build/shared/scripts/dep.py <dependency.tar.gz>

It tries to be smart enough to know what package you're trying to compile, find a class I've defined to help manage the dependency build, compile it, and install it. This sets all of the various configure flags correctly, detect platform specific weirdness, runs the correct build command, etc.

It works reliably on Mac and Linux, but no guarantees are provided. You may need to fall back and manually compile a dependency, making sure to set flags like --prefix=/build/local, --enable-shared, etc.

===== Build script =====

The main script for managing the build process is:
	build.py <command>
	
The build.py script defines a number of target platforms and build commands.

A target defines some platform specific variables, such as compile flags, python version, bashrc, etc. The target can be specific using the --target argument, or using the $TARGET environment variable. Available targets are i686-apple-darwin11, i686-redhat-linux, and x86_64-redhat-linux.

A build command executes one step of the process. Available build commands are checkout, build, and upload.

The complete process is as follows:

	# Load platform specific environment variables.
	source vars.sh

	# Run cvs checkout command to get current trunk source.
	build.py checkout

	# Run build command to run CMake and install into /build/stage/...
	build.py build

	# Create distribution, package, and upload
	build.py install package upload
	
===== Command: checkout =====

The "checkout" command performs a clean checkout of the EMAN2 source code into /build/co. It uses the following arguments:

	--cvsroot			eman@blake.grid.bcm.edu:/usr/local/CVS/CVS
	--cvsmodule			eman2
	--cvstag			daily
	
The new checkout will be in /build/co/<cvsmodule>.<cvstag>. With the defaults, this would be /build/co/eman2.daily.

===== Command: build =====

The "build" command runs cmake, make clean, make, make install, and other post-processing steps necessary to prepare the EMAN2 binary distribution.

You can skip the "make clean" step by using --clean=0 if you are just debugging.

===== Command: install =====

The "install" command performs post-processing steps.

The post-processing steps are platform specific. For instance, the "build" command on Mac OS X runs:
	 CopyExtlib			Copy /build/extlib to EMAN2/extlib
	 CopyShrc			Write EMAN2/eman2.bashrc and EMAN2/eman2.cshrc
	 FixInterpreter		Change the script interpreter line to #!/usr/bin/python2.7
	 FixLinks			Create symlinks from libEM.dylib to libEM.so, etc.
	 FixInstallNames	Use the Mac OS X "install_name_tool" to fix library locations; this is used instead of LD_LIBRARY_PATH.

The Linux target omits FixInterpreter, FixLinks, and FixInstallNames. Other platforms may have similar variations.

===== Command: package =====

Prepares a .tar.gz distribution on Linux, or a .dmg distribution on Mac, of the /build/stage/eman2.daily directory.

===== Upload =====

The "upload" command uploads the release to the web server using the following arguments:

	--scphost	eman@10.10.9.104
	--scpdest	/home/zope-extdata/reposit/ncmi/software/counter_222/software_86

===== Notes on Virtual Machine Setup =====

The build is currently done using a number of Virtual Machines running in VMWare.

	Mac OS X 10.6 Server
	Mac OS X 10.7
	CentOS 5.9 32-bit
	CentOS 5.9 64-bit
	Partially working: Windows 7

User accounts on the host machine and guest VMs are "eman" with the standard password.

In each VM, a nightly cron job executes to prepare and upload the builds:

	0 12 * * * source /build/vars.sh; /build/build.py checkout
	5 12 * * * source /build/vars.sh; /build/build.py build install package
	0 17 * * * source /build/vars.sh; /build/build.py upload

The cron scripts assume that public key SSH authentication is working to both the CVS server and the web server.

===== Notes on /build/local and /build/extlib ======

/build/local is a normal "PREFIX" directory, similar to /usr/local. Dependencies are passed --prefix=/build/local. This places shared libraries in /build/local/lib, binaries in /build/local/bin, etc. 

Python and Qt are also installed with the /build/local prefix. Previously they used local/python and local/qt4.

Python packages are installed in /build/local/site-packages. On Linux, extlib/lib/python2.7/site-packages becomes symlinked to this directory.

/build/extlib is a copy of /build/local with extraneous items removed. For instance, it omits /build/share, /build/man, /build/include, debug versions of libraries, etc. This is done to save space in the distribution.

