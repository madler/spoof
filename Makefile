CFLAGS=-O3 -Wall -Wextra -Wcast-qual -std=c99 -pedantic
CXXFLAGS=-O3 -Wall -Wextra -std=c++11 -pedantic
all: spoof flip ruse
fline.o: fline.c fline.h
spoof: spoof.c fline.o
flip: flip.c
ruse: ruse.cc
clean:
	@rm -f *.o spoof flip ruse
