CC      = gcc
CFLAGS  = -Wall -Wextra -Wpedantic -I../include \
		  -Wl,--gc-sections,--print-gc-sections

# default target
all: test_voronoi

test_voronoi: test_voronoi.o ../voronoi3d.a
	$(CC) -o $@ $^ $(CFLAGS)
	rm test_voronoi.o

# compile source
test_voronoi.o: test_voronoi.c
	$(CC) -c $< $(CFLAGS)

