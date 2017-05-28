if(NOT WIN32)
    INSTALL(PROGRAMS uninstall_openmpi.sh build_and_install_openmpi.sh build_pydusa_numpy.sh install_pydusa_numpy.sh build_eman_numpy.sh install_eman_numpy.sh install_pydusa.sh
        DESTINATION    utils
    )
endif()
