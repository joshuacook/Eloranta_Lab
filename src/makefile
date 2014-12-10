UNAME := $(shell uname)

CC=gcc

ifeq ($(UNAME), Darwin)
  CFLAGS = -g ../../OpenBLAS/libopenblas.a -lm -lpthread
else
  CFLAGS=-I /usr/include/openblas -L /usr/include/openblas -lopenblas -lpthread -lgfortran -llapack -lm 
endif

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)
	
	
