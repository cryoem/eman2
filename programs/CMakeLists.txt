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
	# Get git repository's directory
	execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --git-dir
					WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
					OUTPUT_VARIABLE EMAN_GITDIR
					OUTPUT_STRIP_TRAILING_WHITESPACE
					)
	get_filename_component(EMAN_GITDIR "${EMAN_GITDIR}/.." ABSOLUTE BASE_DIR ${CMAKE_SOURCE_DIR})
	
	# If the source directory is missing .git, but EMAN_GITDIR is still non-empty,
	# the git repository found is, probably, a parent repository.
	# So, we test if the directories found by git and cmake are the same. 
	cmake_print_variables(EMAN_GITDIR)
	cmake_print_variables(CMAKE_SOURCE_DIR)
	
	if(NOT "${EMAN_GITDIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
		set(EMAN_GITHASH "NOT-INSTALLED-FROM-GIT-REPO")
	else()
		if("$ENV{CONDA_BUILD_STATE}" STREQUAL "BUILD" )
			set(EMAN_GITHASH $ENV{GIT_FULL_HASH})
		else()
			execute_process(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
							WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
							OUTPUT_VARIABLE EMAN_GITHASH
							OUTPUT_STRIP_TRAILING_WHITESPACE
							)
		endif()
		# git-describe output: <tag>-<num-of-commits-since-tag>-g<hash>
		# We are insterested in extracting the hash part of the git-describe output
		string(REGEX REPLACE
			   "^.*-.*-g" ""
			   EMAN_GITHASH "${EMAN_GITHASH}"
			   )
	endif()
endif()

string(TIMESTAMP EMAN_TIMESTAMP "%Y-%m-%d %H:%M")
configure_file(${CMAKE_SOURCE_DIR}/libpyEM/EMAN2_meta.py.in
				${SP_DIR}/EMAN2_meta.py
				)
