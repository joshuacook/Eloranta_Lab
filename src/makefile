CC=gcc
CFLAGS=-I ~/bin/OpenBLAS/ -L ~/bin/OpenBLAS/ -lopenblas -lpthread -lgfortran

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
	
	