LINK_LIBRARIES(EM2)
ADD_EXECUTABLE(e2iminfo iminfo.C)
ADD_EXECUTABLE(e2proc2d proc2d.C)
ADD_EXECUTABLE(e2speedtest speedtest.C)
INSTALL_TARGETS(/bin e2iminfo e2proc2d e2speedtest)

