include ../make.conf

all: cgrd2veloc grd2cdf vector2cdf view

vector2cdf: vector2cdf.o
	$(CC) $(CFLAGS) -o vector2cdf vector2cdf.o $(LDFLAGS) -lnetcdf

vector2cdf.o: vector2cdf.c
	$(CC) $(CFLAGS) -c vector2cdf.c

cgrd2veloc: cgrd2veloc.o
	$(CC) $(CFLAGS) -o cgrd2veloc cgrd2veloc.o $(LDFLAGS) -lnetcdf

grd2cdf: grd2cdf.o
	$(CC) $(CFLAGS) -o grd2cdf grd2cdf.o $(LDFLAGS) -lnetcdf

grd2cdf.o: grd2cdf.c
	$(CC) $(CFLAGS) -c grd2cdf.c

cgrd2veloc.o: cgrd2veloc.c
	$(CC) $(CFLAGS) -c cgrd2veloc.c

view: view.o
	$(CC) $(CFLAGS) -o view view.o -lgrace_np

view.o: view.c
	$(CC) $(CFLAGS) -c view.c

clean:
	-rm grd2cdf cgrd2veloc *~ *.o vector2cdf view

