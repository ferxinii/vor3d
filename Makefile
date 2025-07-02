

make: 
	clang test_voronoi.c voronoi.c vd_3d.c bpoly.c algebra.c geometry.c dt_3d_incremental.c simplical_complex.c convhull_3d.c ./predicates/build/Bin/libpredicates.dylib -o test_voronoi -Wpedantic -Wextra -Wall -Werror -Wl,-rpath,./predicates/build/Bin/ -O3 -g
