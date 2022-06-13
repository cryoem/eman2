add_executable(test_eer test_eer.cpp)
target_link_libraries(test_eer EM2)
add_test(test-eer test_eer)
