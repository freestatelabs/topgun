
CC = gcc 
CFLAGS = -march=native 
CFILES = main.c vector.c matrix.c blas2.c 
OFILE = main.o

neon: 
	$(CC) $(CFLAGS) $(CFILES) -o $(OFILE)
	./$(OFILE)