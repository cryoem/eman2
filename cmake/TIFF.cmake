CHECK_REQUIRED_LIB(TIFF tiff tiffio.h "" "")
add_definitions(-DUSE_TIFF)
CHECK_LIB_ONLY(JPEG jpeg)
