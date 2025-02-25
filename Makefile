

make: 
	clang test.c -o test -Wpedantic -Wextra -Wall -Werror -g -O0 -fsanitize=address
