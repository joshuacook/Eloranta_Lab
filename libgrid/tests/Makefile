
include ../make.conf

all: test_wf1d test_wf1d_itp_linear test_wf1d_itp_nonlinear \
     test_wf2d test_wf2d_itp_linear test_wf2d_itp_nonlinear \
     test_wf3d test_wf3d_itp_linear test_wf3d_itp_nonlinear cfft-2d cfft-3d

cfft-2d: cfft-2d.o
	$(CC) $(CFLAGS) -o cfft-2d cfft-2d.o $(LDFLAGS)

cfft-2d.o: cfft-2d.c
	$(CC) $(CFLAGS) -c cfft-2d.c

cfft-3d: cfft-3d.o
	$(CC) $(CFLAGS) -o cfft-3d cfft-3d.o $(LDFLAGS)

cfft-3d.o: cfft-3d.c
	$(CC) $(CFLAGS) -c cfft-3d.c

test_wf1d: test_wf1d.o
	$(CC) $(CFLAGS) -o test_wf1d test_wf1d.o $(LDFLAGS)

test_wf1d.o: test_wf1d.c

test_wf1d_itp_linear: test_wf1d_itp_linear.o
	$(CC) $(CFLAGS) -o test_wf1d_itp_linear test_wf1d_itp_linear.o $(LDFLAGS)

test_wf1d_itp_linear.o: test_wf1d_itp_linear.c

test_wf1d_itp_nonlinear: test_wf1d_itp_nonlinear.o
	$(CC) $(CFLAGS) -o test_wf1d_itp_nonlinear test_wf1d_itp_nonlinear.o $(LDFLAGS)

test_wf1d_itp_nonlinear.o: test_wf1d_itp_nonlinear.c

test_wf2d: test_wf2d.o
	$(CC) $(CFLAGS) -o test_wf2d test_wf2d.o $(LDFLAGS)

test_wf2d.o: test_wf2d.c

test_wf2d_itp_linear: test_wf2d_itp_linear.o
	$(CC) $(CFLAGS) -o test_wf2d_itp_linear test_wf2d_itp_linear.o $(LDFLAGS)

test_wf2d_itp_linear.o: test_wf2d_itp_linear.c

test_wf2d_itp_nonlinear: test_wf2d_itp_nonlinear.o
	$(CC) $(CFLAGS) -o test_wf2d_itp_nonlinear test_wf2d_itp_nonlinear.o $(LDFLAGS)

test_wf2d_itp_nonlinear.o: test_wf2d_itp_nonlinear.c

test_wf3d.o: test_wf3d.c

test_wf3d: test_wf3d.o
	$(CC) $(CFLAGS) -o test_wf3d test_wf3d.o $(LDFLAGS)

test_wf3d_itp_linear: test_wf3d_itp_linear.o
	$(CC) $(CFLAGS) -o test_wf3d_itp_linear test_wf3d_itp_linear.o $(LDFLAGS)

test_wf3d_itp_linear.o: test_wf3d_itp_linear.c

test_wf3d_itp_nonlinear: test_wf3d_itp_nonlinear.o
	$(CC) $(CFLAGS) -o test_wf3d_itp_nonlinear test_wf3d_itp_nonlinear.o $(LDFLAGS)

test_wf3d_itp_nonlinear.o: test_wf3d_itp_nonlinear.c

clean:
	-rm *.o *~ 
	-rm test_wf1d test_wf1d_itp_linear test_wf1d_itp_nonlinear
	-rm test_wf2d test_wf2d_itp_linear test_wf2d_itp_nonlinear 
	-rm test_wf3d test_wf3d_itp_linear test_wf3d_itp_nonlinear 
	-rm cfft-3d cfft-2d
