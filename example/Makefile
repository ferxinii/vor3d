CC      = gcc
CFLAGS  = -Wall -Wextra -Wpedantic -I../include

# default target
all: example

example: example.o ../vor3d.a
	$(CC) -o $@ $^ $(CFLAGS)
	rm example.o

# compile source
example.o: example.c
	$(CC) -c $< $(CFLAGS)

