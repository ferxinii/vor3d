CC      = clang
CFLAGS  = -Wall -Wextra -Wpedantic -Iinclude -Igeometric_predicates/include
AR      = ar

SRC     = $(wildcard src/*.c)
OBJ_DIR = obj
OBJ     = $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(SRC))

GEOM_PRED_LIB = geometric_predicates/build/Bin/libpredicates.a
LIB     = vor3d.a

all: $(LIB)

# compile source files into obj/
$(OBJ_DIR)/%.o: src/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# link program into a single archive
$(LIB): $(OBJ) $(GEOM_PRED_LIB)
	@mkdir -p tmp_objs
	@cd tmp_objs && ar -x ../$(GEOM_PRED_LIB)
	@$(AR) rcs $@ $(OBJ) tmp_objs/*.o
	@rm -rf tmp_objs

