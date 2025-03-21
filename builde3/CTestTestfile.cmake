# CMake generated Testfile for 
# Source directory: /home/stevel/pro/eman2
# Build directory: /home/stevel/pro/eman2/builde3
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(py-compile "/home2/stevel/miniforge3/envs/e3/bin/python3.12" "-m" "compileall" "-q" "-x" ".git" "-x" "sparx" "-x" "sphire" "/home/stevel/pro/eman2")
set_tests_properties(py-compile PROPERTIES  _BACKTRACE_TRIPLES "/home/stevel/pro/eman2/cmake/testing.cmake;8;add_test;/home/stevel/pro/eman2/cmake/testing.cmake;0;;/home/stevel/pro/eman2/CMakeLists.txt;298;include;/home/stevel/pro/eman2/CMakeLists.txt;0;")
add_test(imports "/home2/stevel/miniforge3/envs/e3/bin/python3.12" "/home/stevel/pro/eman2/tests/test_imports.py")
set_tests_properties(imports PROPERTIES  _BACKTRACE_TRIPLES "/home/stevel/pro/eman2/cmake/testing.cmake;17;add_test;/home/stevel/pro/eman2/cmake/testing.cmake;0;;/home/stevel/pro/eman2/CMakeLists.txt;298;include;/home/stevel/pro/eman2/CMakeLists.txt;0;")
add_test(test-EMAN2DIR "/home2/stevel/miniforge3/envs/e3/bin/python3.12" "/home/stevel/pro/eman2/tests/test_EMAN2DIR.py")
set_tests_properties(test-EMAN2DIR PROPERTIES  _BACKTRACE_TRIPLES "/home/stevel/pro/eman2/cmake/testing.cmake;21;add_test;/home/stevel/pro/eman2/cmake/testing.cmake;0;;/home/stevel/pro/eman2/CMakeLists.txt;298;include;/home/stevel/pro/eman2/CMakeLists.txt;0;")
add_test(nose-tests "NOSETESTS_EXECUTABLE-NOTFOUND" "-vv" "--exe" "-m" "^test_*" "-e" "^test_image_" "-e" "test_main" "-e" "test_result" "-e" "test_boxing" "-a" "!broken" "/home/stevel/pro/eman2/rt/pyem/")
set_tests_properties(nose-tests PROPERTIES  _BACKTRACE_TRIPLES "/home/stevel/pro/eman2/cmake/testing.cmake;26;add_test;/home/stevel/pro/eman2/cmake/testing.cmake;0;;/home/stevel/pro/eman2/CMakeLists.txt;298;include;/home/stevel/pro/eman2/CMakeLists.txt;0;")
add_test(progs "/home2/stevel/miniforge3/envs/e3/bin/python3.12" "/home/stevel/pro/eman2/tests/run_prog_tests.py")
set_tests_properties(progs PROPERTIES  _BACKTRACE_TRIPLES "/home/stevel/pro/eman2/cmake/testing.cmake;61;add_test;/home/stevel/pro/eman2/cmake/testing.cmake;0;;/home/stevel/pro/eman2/CMakeLists.txt;298;include;/home/stevel/pro/eman2/CMakeLists.txt;0;")
subdirs("libEM")
subdirs("libpyEM")
subdirs("sparx")
subdirs("sphire")
subdirs("utils")
subdirs("examples")
subdirs("programs")
subdirs("doc")
subdirs("images")
subdirs("rt")
