add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE smoke.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE vendor)
