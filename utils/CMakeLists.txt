#ADD_EXECUTABLE(e2iminfo iminfo.cpp)
#ADD_EXECUTABLE(shared_ptr baldwinweights.cpp)
#ADD_EXECUTABLE(pdbreader pdbreader.cpp)
#ADD_EXECUTABLE(e2ptcltrace e2ptcltrace.cpp)
ADD_EXECUTABLE(e2speedtest speedtest.cpp)
TARGET_LINK_LIBRARIES(e2speedtest EM2)
INSTALL(TARGETS e2speedtest 
    DESTINATION    bin
)
