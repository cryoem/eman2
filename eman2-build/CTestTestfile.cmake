# CMake generated Testfile for 
# Source directory: /Users/landang/Lan_EMAN/eman2
# Build directory: /Users/landang/Lan_EMAN/eman2/eman2-build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(py-compile "/Users/landang/miniconda3/envs/clone_folk/bin/python3.9" "-m" "compileall" "-q" "-x" ".git" "-x" "sparx" "-x" "sphire" "/Users/landang/Lan_EMAN/eman2")
set_tests_properties(py-compile PROPERTIES  _BACKTRACE_TRIPLES "/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;8;add_test;/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;0;;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;298;include;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;0;")
add_test(imports "/Users/landang/miniconda3/envs/clone_folk/bin/python3.9" "/Users/landang/Lan_EMAN/eman2/tests/test_imports.py")
set_tests_properties(imports PROPERTIES  _BACKTRACE_TRIPLES "/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;17;add_test;/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;0;;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;298;include;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;0;")
add_test(test-EMAN2DIR "/Users/landang/miniconda3/envs/clone_folk/bin/python3.9" "/Users/landang/Lan_EMAN/eman2/tests/test_EMAN2DIR.py")
set_tests_properties(test-EMAN2DIR PROPERTIES  _BACKTRACE_TRIPLES "/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;21;add_test;/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;0;;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;298;include;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;0;")
add_test(nose-tests "/Users/landang/miniconda3/envs/clone_folk/bin/nosetests" "-vv" "--exe" "-m" "^test_*" "-e" "^test_image_" "-e" "test_main" "-e" "test_result" "-e" "test_boxing" "-a" "!broken" "/Users/landang/Lan_EMAN/eman2/rt/pyem/")
set_tests_properties(nose-tests PROPERTIES  _BACKTRACE_TRIPLES "/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;26;add_test;/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;0;;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;298;include;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;0;")
add_test(progs "/Users/landang/miniconda3/envs/clone_folk/bin/python3.9" "/Users/landang/Lan_EMAN/eman2/tests/run_prog_tests.py")
set_tests_properties(progs PROPERTIES  _BACKTRACE_TRIPLES "/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;61;add_test;/Users/landang/Lan_EMAN/eman2/cmake/testing.cmake;0;;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;298;include;/Users/landang/Lan_EMAN/eman2/CMakeLists.txt;0;")
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
