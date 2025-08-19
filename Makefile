CC      = clang
CFLAGS  = -Wall -Wextra -Wpedantic -Iinclude -Igeometric_predicates/include
AR      = ar

SRC     = $(wildcard src/*.c)
OBJ     = $(SRC:.c=.o)

GEOM_PRED_LIB = geometric_predicates/build/Bin/libpredicates.a

LIB     = voronoi3d.a

all: $(LIB)

# link program into a single archive
$(LIB): $(OBJ) $(GEOM_PRED_LIB)
	@mkdir -p tmp_objs
	@cd tmp_objs && ar -x ../$(GEOM_PRED_LIB)
	@$(AR) rcs $@ $(OBJ) tmp_objs/*.o
	@rm -rf tmp_objs

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

