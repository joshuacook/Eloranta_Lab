CC=gcc
CFLAGS=-I ~/src/OpenBLAS/ -L ~/src/OpenBLAS/ -lopenblas -lpthread -lgfortran

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
	
	
