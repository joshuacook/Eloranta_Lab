
include ../make.conf

OBJS = cgrid1d.o cgrid1d_fft.o rgrid1d.o rgrid1d_fft.o grid1d_wf.o grid1d_itp.o cgrid2d.o cgrid2d_cyl.o cgrid2d_fft.o rgrid2d.o rgrid2d_cyl.o rgrid2d_fft.o grid2d_wf.o grid2d_wf_cyl.o grid2d_itp.o cgrid3d.o cgrid3d_cyl.o cgrid3d_fft.o rgrid3d.o rgrid3d_cyl.o rgrid3d_fft.o grid3d_wf.o grid3d_wf_cyl.o grid3d_itp.o grid2d_cfft.o threads.o linalg.o ctimer.o other.o

libgrid.a: prototypes $(OBJS)
	$(AR) cr libgrid.a $(OBJS)
	$(RANLIB) libgrid.a

prototypes:
	egrep -h \^EXPORT *.c | tr \{ \; | sed -e "s/EXPORT //g" > proto.h

install:
	cp libgrid.a $(ROOT)/lib
	$(RANLIB) $(ROOT)/lib/libgrid.a
	-mkdir $(ROOT)/include/grid
	cp *.h $(ROOT)/include/grid
	cp ../make.conf $(ROOT)/include/grid

clean:
	-rm libgrid.a *.o *~
