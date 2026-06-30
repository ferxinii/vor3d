Single include in your code:
    #include "indirect_predicates.h"
That one include pulls the entire chain. Just keep all 8 files in the same directory.

One runtime call needed before using any predicates (once per thread):
    initFPU();
This sets up the floating-point environment (defined in numerics.hpp).

Compiler flags required:
    - GCC/Clang: -frounding-math (needed for interval arithmetic's FPU rounding switches to work correctly)
    - MSVC: /fp:strict
