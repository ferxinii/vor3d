CC      = clang
CFLAGS  = -Wall -Wextra -Wpedantic -O3
INCLUDE = -Iinclude -Iexternal/geometric_predicates/include -Iexternal/convhull_3d
AR      = ar

SRC     = $(wildcard src/*.c)
OBJ_DIR = obj
OBJ     = $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(SRC))

GEOM_PRED_LIB = external/geometric_predicates/build/Bin/libpredicates.a
CONVHULL_OBJ = external/convhull_3d/convhull_3d.o

LIB     = vor3d.a

all: $(LIB)

# compile source files into obj/
$(OBJ_DIR)/%.o: src/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# link program into a single archive
$(LIB): $(OBJ)
	@mkdir -p tmp_objs
	@cd tmp_objs && ar -x ../$(GEOM_PRED_LIB)
	@$(AR) rcs $@ $(OBJ) $(CONVHULL_OBJ) tmp_objs/*.o
	@rm -rf tmp_objs

