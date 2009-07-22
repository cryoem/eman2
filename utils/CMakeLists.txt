LINK_LIBRARIES(GLEM2)
# ADD_EXECUTABLE(e2iminfo iminfo.cpp)
# ADD_EXECUTABLE(shared_ptr baldwinweights.cpp)
ADD_EXECUTABLE(e2speedtest speedtest.cpp)
ADD_EXECUTABLE(pdbreader pdbreader.cpp)
# ADD_EXECUTABLE(e2ptcltrace e2ptcltrace.cpp)


INSTALL_TARGETS(/bin e2speedtest)
