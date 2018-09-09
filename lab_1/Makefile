CC = gcc
CFLAGS =
LIBS = -fopenmp
TARGETS = rectify pool convolve

default: clean $(TARGETS)
all: default

.PHONY: clean

rectify:
	$(CC) -o rectify rectification.c lodepng.c $(CFLAGS) -std=c99 -Wall $(LIBS)

pool:
	$(CC) -o pool pooling.c lodepng.c $(CFLAGS) -std=c99 -Wall $(LIBS)

convolve:
	$(CC) -o convolve convolution.c lodepng.c $(CFLAGS) -std=c99 -Wall $(LIBS)

clean:
	-rm -f *.o
	-rm -f $(TARGETS)