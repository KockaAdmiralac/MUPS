CC = gcc

CC_FLAGS = -fopenmp -O3 -lm
CC_FLAGS += -Wall -Wextra

ifeq ($(DEBUG), 1)
CC_FLAGS += -DDEBUG
endif

all: dz1z1/prime dz1z3/feynman

dz1z1/prime: dz1z1/prime.c util.c
	$(CC) $(CC_FLAGS) $(^) -o $(@)

dz1z3/feynman: dz1z3/feynman.c util.c
	$(CC) $(CC_FLAGS) $(^) -o $(@)

clean:
	rm -f dz1z1/prime
	rm -f dz1z3/feynman
