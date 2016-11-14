
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
  cuda_util.cu
  cuda_emfft.cu
 )

# Specify the dependency.
TARGET_LINK_LIBRARIES(EM2Cuda
  ${CUDA_LIBRARIES}
  )

##############################################################################
# Use the CUDA_COMPILE macro.
#CUDA_COMPILE(CUDA_FILES test_bin.cu)
# 

INSTALL(TARGETS EM2Cuda 
        DESTINATION lib)
#INSTALL_TARGETS(/bin proj3d)

