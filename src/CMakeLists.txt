add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE smoke.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE vendor)
target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Release>:
        $<$<CXX_COMPILER_ID:MSVC>:/O2 /arch:AVX2 /fp:fast>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-O3 -march=native -ffast-math>
    >
)
