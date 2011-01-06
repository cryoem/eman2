INSTALL(DIRECTORY . DESTINATION mpi_eman 
        PATTERN "CVS" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN "*" 
        )