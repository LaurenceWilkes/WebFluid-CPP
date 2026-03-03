# WebFluid-CPP

WebFluid-CPP is a C++/SDL3 implementation of the smoke simulation from [WebFluid](https://github.com/LaurenceWilkes/WebFluid).

## Goal and status

The goal of this project is to explore the performance of a lower-level version of the same stable fluid simulation outside of the browser.
Currently, the implementation is a fairly direct translation of the smoke simulation from WebFluid but with double precision instead of single.
The performance is slightly better than the browser version when compiled in release mode. 
Hopefully, there will be ways of optimising/parallelising so that this version runs even faster.

#### Build and Run (Windows)
`cmake --build build --config Release && .\build\src\Release\Smoke.exe`

## Features
- Implemented in C++20
- Uses SDL3 for display and rendering instead of HTML canvas
- Uses double precision for simulation state
- Built with CMake


The implementation lives in `src/smoke.cpp` with SDL3 as a submodule in `vendor/SDL `.
