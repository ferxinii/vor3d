CC      = clang
AR      = ar
CFLAGS  = -Iinclude -I. -Wall -Wextra -Wpedantic -O3


# EXTERNAL DEPENDENCIES:
EXT_INC = -Iexternal/gnuplotC/include \
          -Iexternal/convhull_3d \
          -Iexternal/geometric_predicates/include
EXT_SRC = external/gnuplotC/gnuplotc.c external/convhull_3d/convhull_3d.c
EXT_OBJ = 
EXT_LIB = external/geometric_predicates/build/Bin/libpredicates.a


SRC     = $(wildcard src/*.c)
OBJ_DIR = obj
LIB     = vor3d.a

OBJ           = $(patsubst src/%.c,$(OBJ_DIR)/%.o,$(SRC))
EXT_OBJ_BUILT = $(addprefix $(OBJ_DIR)/,$(notdir $(EXT_SRC:.c=.o)))
ALL_OBJS = $(OBJ) $(EXT_OBJ_BUILT) $(EXT_OBJ)


all: $(LIB)

# compile source files into obj/
$(OBJ_DIR)/%.o: src/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(EXT_INC) -c $< -o $@


# Compile external sources into obj/
$(OBJ_DIR)/%.o: external/%/*.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(EXT_INC) -c $< -o $@

# Link everything into one archive
$(LIB): $(ALL_OBJS) $(EXT_LIB)
	@mkdir -p tmp_objs
	@cd tmp_objs && ar -x ../$(EXT_LIB)
	@$(AR) rcs $@ $(ALL_OBJS) tmp_objs/*.o
	@rm -rf tmp_objs

clean:
	rm -rf $(OBJ_DIR) $(LIB)


