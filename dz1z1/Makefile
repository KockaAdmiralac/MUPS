CC = gcc

CC_FLAGS = -fopenmp -O3
CC_FLAGS += -Wall -Wextra

prime: prime.c ../util.c
	$(CC) $(CC_FLAGS) prime.c ../util.c -o prime

clean:
	rm -f prime
