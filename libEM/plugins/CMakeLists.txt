INSTALL_FILES(/include/plugins "^.+\\.h$")

#SET(LIB_EXT_SOURCE
#    io_template.cpp
#    aligner_template.cpp
#    averager_template.cpp
#    cmp_template.cpp
#    processor_template.cpp
#    projector_template.cpp
#    reconstructor_template.cpp
#   )

#IF(WIN32)
#    ADD_LIBRARY(EMext STATIC ${LIB_EXT_SOURCE})
#ELSE(WIN32)
#    ADD_LIBRARY(EMext SHARED ${LIB_EXT_SOURCE})
#ENDIF(WIN32)

#TARGET_LINK_LIBRARIES(EMext EM2)
#INSTALL_TARGETS(/lib EMext)
