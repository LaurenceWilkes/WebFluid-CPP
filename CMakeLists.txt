cmake_minimum_required(VERSION 3.24)
project(Smoke)

set(CMAKE_CXX_STANDARD 20)

add_executable(smoke src/smoke.cpp)
target_sources(smoke PRIVATE src/smoke.cpp)
