CC=gcc
CFLAGS=-I ~/src/OpenBLAS/ -L ~/src/OpenBLAS/ -lopenblas -lpthread -lgfortran -lm -lf77blas

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
	
	
