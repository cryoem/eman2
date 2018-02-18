### 1. download/install `eman2.21a` tar from the NCMI website

1.1) go to the [eman2.21a application form](http://ncmi.bcm.edu/ncmi/software/counter_222/software_138/manage_addProduct/NCMI/attendee_factory?myname=eman2.21.MacOS.sh) on the NCMI website

1.2) execute the following commands in your terminal to extract and move the contents into a `~/.eman2` directory
```
$ cd ~/Downloads
$ tar xopf eman2-2.2.tar.gz  # extract the contents of the tar into a folder
$ mv eman2-2.2 ~/.eman2  # move extracted contents into an eman2 directory in your HOME directory
$ cd ~/.eman2  # cd into the new directory
```

1.3) open the `~/.eman2/INSTALL` file for further instructions on installing dependencies

### 2. Install the dependencies

2.1) fttw2

2.1.1) download the following tar: http://www.fftw.org/fftw-2.1.5.tar.gz

2.1.2) execute the following commands in your terminal to extract and move the contents into a `~/.fftw2` directory
```
$ cd ~/Downloads
$ tar xopf fftw-2.1.5.tar.gz  # extract the contents of the tar into a folder
$ mv fftw-2.1.5 ~/.fftw2  # move extracted contents into an fftw2 directory in your HOME directory
$ cd ~/.fftw2  # cd into the new directory
```

2.1.3) run the following commands to install fftw2
```
$ ./configure --enable-static=no --enable-shared=yes --enable-float --enable-type-prefix
$ make
```

2.2) install GNU Scientific Library (`gsl`) via `brew install gsl`

2.3) install python via `brew install python`

2.4) install numpy via `brew install numpy`

2.5) install boost via `brew install boost`

2.6) install cmake via `brew install cmake`

2.7) install PyQT4 via `brew install cartr/qt4/pyqt`

2.8) install ipython via `brew install ipython`

### 3) run the following commands from `~/.eman2/INSTALL` and hope for the best:

3.1) create a `build` directory and `cmake` into it
```
$ cd ~/.eman2/
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install
```

3.2) setup login shell for csh/tcsh
```
$ cat >> ~/.cshrc << CSHRC
setenv EMAN2DIR ${HOME}/EMAN2
setenv PATH ${EMAN2DIR}/bin:${PATH}
setenv LD_LIBRARY_PATH  ${LD_LIBRARY_PATH}:${EMAN2DIR}/lib
setenv PYTHONPATH .:${HOME}/EMAN2/lib:${PYTHONPATH}
CSHRC
```

3.3) setup login shell for bash
```
$ cat >> ~/.bashrc << BASHRC
export EMAN2DIR=${HOME}/EMAN2
export PATH=${PATH}:${EMAN2DIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${EMAN2DIR}/lib
export PYTHONPATH=${PYTHONPATH}:${HOME}/EMAN2/lib
BASHRC
```

3.4) close and reopen your terminal to ensure that your environmental variables are set correctly


### 4) you're done maybe??

4.1) return to your `~/.eman2` directory and run some test `programs/`
```
$ cd ~/.eman2/programs/
$ python e2version.py
$ python e2speedtest.py
$ python e2display.py
$ python e2proc2d.py :64:64:1 test.hdf --process mask.sharp:outer_radius=24
```

4.2) if any of those don't work, something has gone wrong and i don't know how to help

