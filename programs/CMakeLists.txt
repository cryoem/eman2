FILE(GLOB e2programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

INSTALL(PROGRAMS ${e2programs}
		DESTINATION    bin 
		COMPONENT PythonFiles
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
		   EMAN_GITHASH "${EMAN_GITHASH}"
		   )
endif()

string(TIMESTAMP EMAN_TIMESTAMP "%Y-%m-%d %H:%M")
configure_file(${CMAKE_SOURCE_DIR}/libpyEM/EMAN2_meta.py.in
				${SP_DIR}/EMAN2_meta.py
				)
