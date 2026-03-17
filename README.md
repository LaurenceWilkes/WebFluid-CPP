# WebFluid-CPP

WebFluid-CPP is a C++/SDL3 implementation of the smoke simulation from [WebFluid](https://github.com/LaurenceWilkes/WebFluid).

## Goal and status

The goal of this project is to explore the performance of a lower-level version of the same stable fluid simulation outside of the browser.
Currently, the implementation is a fairly direct translation of the smoke simulation from WebFluid but including the possibility of displaying the movement of the fluid with particles in a similar way to the fluid visualisation in [Animations](https://github.com/LaurenceWilkes/Animations).
The performance is slightly better than the browser version when compiled in release mode. 
Hopefully, there will be ways of optimising/parallelising so that this version runs even faster.

#### Build and Run (Windows)
`cmake --build build --config Release`

`.\build\src\Release\Smoke.exe`

Run with the flag `--particles` for the particle display of the fluid motion.

## Features
- Uses SDL3 for display and rendering instead of HTML canvas
- Uses single precision for simulation state
- Built with CMake
- Fluid motion is displayed with either colourful smoke or particle trails.


The implementation lives in `src/smoke.cpp` with SDL3 as a submodule in `vendor/SDL `.
