
# Add current directory to the nvcc include line.
#CUDA_INCLUDE_DIRECTORIES(
#  ${CMAKE_CURRENT_SOURCE_DIR}
#)

##############################################################################
# And an executable
#CUDA_ADD_EXECUTABLE(proj3d
#	proj3d.cu
# )

##############################################################################
# Build a shared library
CUDA_ADD_LIBRARY(EM2Cuda
  cuda_projector.cu
  cuda_util.cu
  cuda_emfft.cu
  cuda_processor.cu
  cuda_cmp.cu
 )


# Specify the dependency.
TARGET_LINK_LIBRARIES(EM2Cuda
  ${CUDA_LIBRARIES}
  )

##############################################################################
# Use the CUDA_COMPILE macro.
#CUDA_COMPILE(CUDA_FILES test_bin.cu)
# 

INSTALL_TARGETS(/lib EM2Cuda)
#INSTALL_TARGETS(/bin proj3d)

