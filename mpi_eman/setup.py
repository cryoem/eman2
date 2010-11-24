from distutils.core import setup, Extension

module1 = Extension("mpi_eman", sources = ["mpi_eman.c"])

setup (name = "mpi_eman", version = "0.1", description = "A simple MPI wrapper for use from python in EMAN2",ext_modules = [module1])

