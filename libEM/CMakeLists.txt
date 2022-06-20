add_subdirectory(pmconfig)

if(WIN32)
	add_library(EM2 STATIC "")
else()
	add_library(EM2 SHARED "")
endif()

target_sources(EM2 PRIVATE
			   byteorder.cpp
			   boxingtools.cpp
			   emobject.cpp
			   emfft.cpp
			   log.cpp
			   io/imageio.cpp
			   util.cpp
			   emutil.cpp
			   interp.cpp
			   quaternion.cpp
			   transform.cpp
			   symmetry.cpp
			   geometry.cpp
			   pointarray.cpp
			   pdbreader.cpp
			   io/mrcio.cpp
			   io/eerio.cpp
			   io/dm3io.cpp
			   io/dm4io.cpp
			   io/tifio.cpp
			   io/hdfio.cpp
			   io/hdfio2.cpp
			   io/jpegio.cpp
			   io/renderer.cpp
			   emdata.cpp
			   emdata_io.cpp
			   emdata_core.cpp
			   emdata_cuda.cpp
			   emdata_modular.cpp
			   emdata_metadata.cpp
			   emdata_transform.cpp
			   io/pifio.cpp
			   io/v4l2io.cpp
			   io/vtkio.cpp
			   io/spiderio.cpp
			   io/sspiderio.cpp
			   io/pgmio.cpp
			   io/fitsio.cpp
			   io/icosio.cpp
			   io/lstio.cpp
			   io/lstfastio.cpp
			   io/pngio.cpp
			   io/salio.cpp
			   io/amiraio.cpp
			   io/xplorio.cpp
			   io/gatan2io.cpp
			   io/emio.cpp
			   io/imagicio.cpp
			   io/imagicio2.cpp
			   io/df3io.cpp
			   io/omapio.cpp
			   io/situsio.cpp
			   io/serio.cpp
			   emcache.cpp
			   ctf.cpp
			   xydata.cpp
			   processor.cpp
			   aligner.cpp
			   projector.cpp
			   cmp.cpp
			   averager.cpp
			   reconstructor.cpp
			   reconstructor_tools.cpp
			   exception.cpp
			   testutil.cpp
			   analyzer.cpp
			   randnum.cpp
			   io/hdf_filecache.cpp
			   polardata.cpp
			   tomoseg.cpp
			   )

add_subdirectory(gorgon)
add_subdirectory(sparx)
add_subdirectory(sphire)

if(WIN32)
	target_compile_definitions(EM2 PUBLIC _CRT_SECURE_NO_WARNINGS _SCL_SECURE_NO_WARNINGS)
endif()

target_link_libraries(EM2 HDF5::HDF5 GSL::gsl GSL::gslcblas)

install(TARGETS EM2
		DESTINATION ${SP_DIR}
		)

if(ENABLE_FFTW3)
	find_package(FFTW3 REQUIRED)
	
	if(ENABLE_NFFT2)
		find_package(NFFT2 REQUIRED)
		target_link_libraries(NFFT2 FFTW3D::FFTW3D)
	endif()
endif()

if(WIN32 OR ENABLE_FFTW3)
    target_link_libraries(EM2 FFTW3::FFTW3)
endif()

option(ENABLE_TIFF "enable TIFF support" ON)
if(ENABLE_TIFF)
	include(${CMAKE_SOURCE_DIR}/cmake/TIFF.cmake)
    target_link_libraries(EM2 TIFF::TIFF)
endif()

option(ENABLE_JPEG "enable JPEG support" ON)
if(ENABLE_JPEG)
	include(${CMAKE_SOURCE_DIR}/cmake/JPEG.cmake)
    target_link_libraries(EM2 JPEG::JPEG)
endif()

option(ENABLE_PNG "enable PNG support" ON)
if(ENABLE_PNG)
	include(${CMAKE_SOURCE_DIR}/cmake/PNG.cmake)
    target_link_libraries(EM2 PNG::PNG)
endif()

if(ENABLE_SPARX_CUDA)
    target_link_libraries(EM2 EM2SparxCuda FFTW3::FFTW3)
endif()

if(ENABLE_EMAN_CUDA)
    target_link_libraries(EM2 EM2Cuda FFTW3::FFTW3)
endif()

if(ENABLE_NFFT2)
	target_sources(EM2 PRIVATE pointarray.cpp)
	target_link_libraries(EM2 NFFT2::NFFT2 FFTW3D::FFTW3D)
endif()

option(ENABLE_DJBFFT "enable djbfft support" OFF)
mark_as_advanced(ENABLE_DJBFFT)
if(ENABLE_DJBFFT)
	option(ENABLE_NFFT2 "enable nfft2 support" OFF)
	
	find_package(DJBFFT REQUIRED)
	target_link_libraries(EM2 DJBFFT::DJBFFT)
endif()

option(ENABLE_OPTPP "enable Opt++ support" OFF)
if(ENABLE_OPTPP)
	find_package(OPTPP REQUIRED)
	target_link_libraries(EM2 OPTPP::OPTPP)
endif()

if(ENABLE_ACML_FFT)
	find_package(ACML REQUIRED)
	target_link_libraries(EM2 ACML::ACML)
endif()

option(ENABLE_OPENGL "enable EMAN2 OpenGL usage" ON)
if(ENABLE_OPENGL)
	include(${CMAKE_SOURCE_DIR}/cmake/OpenGL.cmake)

	option(ENABLE_STATIC_FTGL "If build EMAN with libftgl.a" OFF)
	mark_as_advanced(ENABLE_STATIC_FTGL)
	
	find_package(FTGL REQUIRED)
	
	add_subdirectory(fonts)
endif()
