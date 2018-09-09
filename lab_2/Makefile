CC = mpicc
CFLAGS =
LIBS =
TARGETS = grid_4_4 grid_512_512

default: clean $(TARGETS)
all: default

.PHONY: clean

grid_4_4:
	$(CC) -o grid_4_4 grid_4_4.c $(CFLAGS) -std=c99 $(LIBS)

grid_512_512:
	$(CC) -o grid_512_512 grid_512_512.c $(CFLAGS) -std=c99 $(LIBS)

clean:
	-rm -f *.o
	-rm -f $(TARGETS)