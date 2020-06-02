#!/usr/bin/env python
import numpy as np
from EMAN2 import *
import EMAN2
import mpi
import sp_utilities
import ctypes
import numpy
import EMAN2

from future import standard_library
standard_library.install_aliases()
mpi.mpi_init(0, [])

Blockdata = {
    "main_node" : 0,
    "nproc" : mpi.mpi_comm_size(mpi.MPI_COMM_WORLD),
    "shared_comm": mpi.mpi_comm_split_type(
    mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED, 0, mpi.MPI_INFO_NULL),
    "myid_on_node": mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD),
    }

nx = 4
ny = 4
nz = 4
data = []
if( Blockdata["myid_on_node"] == 0 ):
    for i in range(nz):
        data.append(EMData(nx, ny) + i)

size = nx * ny * nz
disp_unit = 4
nimastack = nz
target_nx = nx
size_of_one_image = nx * ny


win_sm, base_ptr = mpi.mpi_win_allocate_shared( size*disp_unit , disp_unit, mpi.MPI_INFO_NULL, Blockdata["shared_comm"])

if( Blockdata["myid_on_node"] != 0 ):
    base_ptr, = mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)

# buffer = np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype = 'f4')
# buffer = buffer.reshape(nimastack, target_nx, target_nx)

ptr = ctypes.cast(base_ptr, ctypes.POINTER(ctypes.c_int * size))
buffer = numpy.frombuffer(ptr.contents, dtype="f4")
buffer = buffer.reshape(nimastack, target_nx, target_nx)

emnumpy2 = EMNumPy()
bigbuffer = emnumpy2.register_numpy_to_emdata(buffer)

# bigbuffer = EMNumPy.numpy2em(buffer)


if( Blockdata["myid_on_node"] == 0 ):
    #  read data on process 0 of each node
    #print "  READING DATA FIRST :",Blockdata["myid"],Blockdata["stack_ali2d"],len(plist)
    for i in range(nimastack):
        bigbuffer.insert_clip( data[i],(0,0,i) )

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
alldata =  [None] * nimastack
emnumpy3 = [None] * nimastack

msk = sp_utilities.model_blank(target_nx, target_nx,1,1)
# for i in range(nimastack):
#     pointer_location = base_ptr + i*size_of_one_image*disp_unit
#     img_buffer  = np.frombuffer(np.core.multiarray.int_asbuffer(pointer_location, size_of_one_image*disp_unit), dtype = 'f4')
#     img_buffer  = img_buffer.reshape(target_nx, target_nx)
#     emnumpy3[i] = EMNumPy()
#     alldata[i]  = emnumpy3[i].register_numpy_to_emdata(img_buffer)


for i in range(nimastack):
    newpoint = base_ptr + i * size_of_one_image * disp_unit
    pointer_location = ctypes.cast(newpoint, ctypes.POINTER(ctypes.c_int * size_of_one_image))
    img_buffer = numpy.frombuffer(pointer_location.contents, dtype="f4")
    img_buffer  = img_buffer.reshape(target_nx, target_nx)
    emnumpy3[i] = EMNumPy()
    alldata[i]  = emnumpy3[i].register_numpy_to_emdata(img_buffer)

mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

print(Blockdata['myid_on_node'], [i.get_2dview() for i in alldata])
mpi.mpi_finalize()
