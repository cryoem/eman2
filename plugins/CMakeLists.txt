SET(LIB_EXT_SOURCE
    io_template.cpp
    aligner_template.cpp
    averager_template.cpp
    cmp_template.cpp
    filter_template.cpp
    projector_template.cpp
    reconstructor_template.cpp
   )

IF(IS_WINDOWS)
    ADD_LIBRARY(EMext STATIC ${LIB_EXT_SOURCE})
ELSE(IS_WINDOWS)
    ADD_LIBRARY(EMext SHARED ${LIB_EXT_SOURCE})
ENDIF(IS_WINDOWS)

TARGET_LINK_LIBRARIES(EMext EM2)
INSTALL_TARGETS(/lib EMext)
