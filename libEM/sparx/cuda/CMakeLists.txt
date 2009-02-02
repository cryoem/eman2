
# Add current directory to the nvcc include line.
#CUDA_INCLUDE_DIRECTORIES(
#  ${CMAKE_CURRENT_SOURCE_DIR}
#)

##############################################################################
# Use one executable only.
#CUDA_ADD_EXECUTABLE(sparx_cuda_test
# sparx_cuda_test.cu
# other_code.c
# )

##############################################################################
# Build a shared library
#CUDA_ADD_LIBRARY(EM2SparxCuda SHARED
#  sparx_blas.cu
#  )
  

# Specify the dependency.
#TARGET_LINK_LIBRARIES(EM2SparxCuda
#  cufft
#  cudart
#  )

##############################################################################

#INSTALL_TARGETS(/lib EM2SparxCuda)
#INSTALL_TARGETS(/bin sparx_cuda_test)