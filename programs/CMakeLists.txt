FILE(GLOB e2programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

FILE(GLOB exclusion_list
			e2plotFSC.py
			e2refine_evenodd.py
			e2refine_evenodd.py
			e2refmerge.py
			e2resolution.py
			e2montecarlorecon.py
			e2markbadparticles.py
			e22.py
			e2boxer21.py
			e2ali2img.py
			e2avg_ffts.py
			e2fileinfo.py
			e2flick.py
			e2findsubunit.py
			e2helical_recons.py
			e2modeleval.py
			e2preferences.py
			e2scp.py
			e2tomoallvall.py
			e2tomoaverage.py
			e2tomohunter.py
			e2remoted.py
			e2tomosim.py
		)

add_custom_target(Programs
				  COMMAND ${CMAKE_COMMAND} -DCOMPONENT=Programs -P cmake_install.cmake
				  )

LIST(REMOVE_ITEM e2programs ${exclusion_list})
INSTALL(PROGRAMS ${e2programs}
		DESTINATION    bin 
		COMPONENT Programs
)

IF(WIN32)
	INSTALL(FILES e2boxer.py
		DESTINATION bin
		RENAME	e2boxer.pyw
		)
		
	INSTALL(FILES e2display.py
		DESTINATION bin
		RENAME	e2display.pyw
		)
ENDIF(WIN32)

find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
	execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
					WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
					OUTPUT_VARIABLE EMAN_GITHASH
					OUTPUT_STRIP_TRAILING_WHITESPACE
					)
	# git-describe output: <tag>-<num-of-commits-since-tag>-g<hash>
	string(REGEX REPLACE
		   "^.*-.*-g" ""
		   EMAN_GITHASH ${EMAN_GITHASH}
		   )
endif()

string(TIMESTAMP EMAN_TIMESTAMP "%Y-%m-%d %H:%M")
configure_file(${CMAKE_SOURCE_DIR}/libpyEM/EMAN2_meta.py.in
				${SP_DIR}/EMAN2_meta.py
				)
