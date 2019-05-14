enable_testing()

add_custom_target(test-verbose
		COMMAND ${CMAKE_COMMAND} -P cmake_install.cmake
		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release
		)

add_test(NAME imports
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_imports.py
		 )

add_test(NAME test-EMAN2DIR
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_EMAN2DIR.py
		 )

if(NOT WIN32)
	add_test(NAME nose-tests
			COMMAND ${NOSETESTS_EXECUTABLE} -vv --exe -m "^test_*" -e "^test_image_" -e "test_main" -e "test_result" -e "test_boxing" -a \!broken
					${CMAKE_SOURCE_DIR}/rt/pyem/ 
			)
else()
	set(test_methods-win
			test_imageio.py:TestHdfIO.test_read_image
			test_imageio.py:TestHdfIO.test_write_image
			test_imageio.py:TestHdfIO.test_read_write_hdf
			test_imageio.py:TestPNGIO.test_write_png
		)

	foreach(t ${test_methods-win})
		add_test(NAME ${t}
				COMMAND ${NOSETESTS_EXECUTABLE} -v -m "^test_*" -a \!broken
						${CMAKE_SOURCE_DIR}/rt/pyem/${t}
				)
	endforeach()
endif()

add_custom_target(test-rt
		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R nose-tests
		DEPENDS PythonFiles
		)

add_custom_target(test-verbose-broken
		COMMAND ${NOSETESTS_EXECUTABLE} -v -m "^test_*" -a broken ${CMAKE_SOURCE_DIR}/rt/pyem/*
		)

if(NOT WIN32)
	add_custom_target(test-progs
					  COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R progs
					  DEPENDS PythonFiles
					  )
	
	add_test(NAME progs
			 COMMAND bash ${CMAKE_SOURCE_DIR}/tests/run_prog_tests.sh
			 )
endif()

add_test(NAME py-compile
		COMMAND ${PYTHON_EXECUTABLE} -m compileall -q ${CMAKE_SOURCE_DIR}
		)

add_custom_target(test-py-compile
		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R py-compile
		DEPENDS PythonFiles
		)
