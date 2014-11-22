CC=gcc
# CFLAGS=-I ~/src/OpenBLAS/ -L ~/src/OpenBLAS/ -lopenblas -lpthread -lgfortran -lm -lf77blas
CFLAGS = -g ../../OpenBLAS/libopenblas.a -lm -lpthread
%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
	
	
