# Make file for Mosaic

CC = gcc
CFLAGS =  -O2 -pedantic
#CFLAGS = -g
LIB = -lm

OBJ = tools.o seqtools.o

all: mosaic

mosaic: $(OBJ) mosaic.o tools.h seqtools.h mosaic_fb.h
	$(CC) $(CFLAGS) -o mosaic $(OBJ) mosaic.o $(LIB)

mosaic.o: $(OBJ) mosaic_fb.c tools.h seqtools.h mosaic_fb.h
	$(CC) $(CFLAGS) -c -o mosaic.o mosaic_fb.c

tools.o: tools.c tools.h
	$(CC) $(CFLAGS) -c -o tools.o tools.c

seqtools.o: seqtools.c seqtools.h tools.h
	$(CC) $(CFLAGS) -c -o seqtools.o seqtools.c 

clean:
	rm -rf *.o
	rm -rf *~
