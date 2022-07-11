add_executable(test_eer test_eer.cpp)
target_link_libraries(test_eer EM2)
add_test(test-eer test_eer)

add_custom_target(test-eer
        COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R test-eer
        DEPENDS test_eer
        )

if (ENABLE_CATCH2)
    Include(FetchContent)

    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG        v3.0.1 # or a later release
    )

    FetchContent_MakeAvailable(Catch2)


    add_executable(test_renderer test_renderer.cpp)
    target_link_libraries(test_renderer PRIVATE EM2 Catch2::Catch2WithMain)
    add_test(test-renderer test_renderer)

    add_custom_target(test-renderer
            COMMAND ${CMAKE_CTEST_COMMAND} -V --output-on-failure -C Release -R test-renderer
            DEPENDS test_renderer
            )
endif ()
