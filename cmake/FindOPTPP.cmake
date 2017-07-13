FIND_PATH(OPTPP_INCLUDE_PATH Opt.h /usr/include /usr/include/opt++
		  /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
		  $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
FIND_PATH(NEWMAT_INCLUDE_PATH newmat.h /usr/include /usr/include/opt++
		  /usr/include/newmat /usr/local/include/newmat
		  /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
		  $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
FIND_LIBRARY(OPTPP_LIBRARY NAMES opt-linux PATHS /usr/lib /usr/local/lib)
FIND_LIBRARY(NEWMAT_LIBRARY NAMES newmat-linux PATHS /usr/lib /usr/local/lib)

ADD_DEFINITIONS(-DUSE_OPTPP)
INCLUDE_DIRECTORIES(${OPTPP_INCLUDE_PATH} ${NEWMAT_INCLUDE_PATH})
