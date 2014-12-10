CC=gcc
CFLAGS=-I /usr/include/openblas -L /usr/include/openblas -lopenblas -lpthread -lgfortran -llapack -lm 
# CFLAGS = -g ../../OpenBLAS/libopenblas.a -lm -lpthread
%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
	
	
