CC = gcc

CC_FLAGS = -fopenmp -O3
CC_FLAGS += -Wall -Wextra
LIBS = -lm

ifeq ($(DEBUG), 1)
CC_FLAGS += -DDEBUG
endif

all: dz1z1/prime dz1z3/feynman dz1z5/moldyn

dz1z1/prime: dz1z1/prime.c util.c
	$(CC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

dz1z3/feynman: dz1z3/feynman.c util.c
	$(CC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

dz1z5/moldyn: dz1z5/moldyn.c util.c
	$(CC) $(CC_FLAGS) $(^) -o $(@) $(LIBS)

clean:
	rm -f dz1z1/prime
	rm -f dz1z3/feynman
