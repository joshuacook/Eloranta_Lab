UNAME := $(shell uname)

CC=gcc

ifeq ($(UNAME), Darwin)
  CFLAGS = -g /usr/local/Cellar/openblas/0.2.12/lib/libopenblas.a -lm -lpthread
else
  CFLAGS=-I /usr/lib/openblas-base -L /usr/lib/openblas-base -lopenblas -lpthread -lgfortran -llapack -lm 
endif

%.o: %.c $(CC) -c -o $@ $< $(CFLAGS)
  
  