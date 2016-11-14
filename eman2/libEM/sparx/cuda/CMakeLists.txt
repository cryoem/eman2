
# Add current directory to the nvcc include line.
#CUDA_INCLUDE_DIRECTORIES(
#  ${CMAKE_CURRENT_SOURCE_DIR}
#)

##############################################################################
# Use one executable only.
#CUDA_ADD_EXECUTABLE(sparx_cuda_test
# sparx_cuda_test.cu
# cuda_ccf.cu
# )


##############################################################################
# Build a shared library
CUDA_ADD_LIBRARY(EM2SparxCuda
  cuda_ccf.cu
  cuda_mpi_kmeans.cu
)

  
# Specify the dependency.
TARGET_LINK_LIBRARIES(EM2SparxCuda
  ${CUDA_LIBRARIES}
)

##############################################################################

INSTALL_TARGETS(/lib EM2SparxCuda)
#INSTALL_TARGETS(/bin sparx_cuda_test)

