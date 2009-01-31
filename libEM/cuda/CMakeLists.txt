
# Add current directory to the nvcc include line.
#CUDA_INCLUDE_DIRECTORIES(
#  ${CMAKE_CURRENT_SOURCE_DIR}
#)

##############################################################################
# Use one executable only.
#CUDA_ADD_EXECUTABLE(test
#  test_bin.cu
#  main.cc
#  )

##############################################################################
# Link cuda code in a library to something else.
CUDA_ADD_LIBRARY(EM2Cuda SHARED
  cuda_projector.cu
  cuda_util.cu
  cuda_emfft.cu
  )
  
# Something else to link to.
#ADD_EXECUTABLE(lib_example
#  main.cc
#  )

# Specify the dependency.
TARGET_LINK_LIBRARIES(EM2Cuda
  cufft
  )

##############################################################################
# Use the CUDA_COMPILE macro.
#CUDA_COMPILE(CUDA_FILES test_bin.cu)

#ADD_EXECUTABLE(cuda_compile_example 
#  ${CUDA_FILES} 
#  main.cc
#  )
#TARGET_LINK_LIBRARIES(cuda_compile_example
#  ${CUDA_TARGET_LINK}
#  )

INSTALL_TARGETS(/lib EM2Cuda)