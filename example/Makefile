CC      = clang
CFLAGS  = -Wall -Wextra -Wpedantic -I../include 

# default target
all: test_voronoi

test_voronoi: test_voronoi.o ../voronoi3d.a
	$(CC) -o $@ $^ $(CFLAGS)
	rm test_voronoi.o

# compile source
test_voronoi.o: test_voronoi.c
	$(CC) -c $< $(CFLAGS)

