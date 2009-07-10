CC=gfortran
CFLAGS=-g -I/usr/include/ -O2 -Wuninitialized
MUMPS=$(HOME)/tmp/MUMPS_4.8.4/
LDFLAGS= -L$(HOME)/local/lib \
         $(MUMPS)/lib/libmumps_common.a \
         $(MUMPS)/lib/libzmumps.a \
         -lm  -lsuperlu
INCLUDEPATH=-I$(HOME)/tmp/eigen2/ -I/usr/include/superlu \
            -I$(HOME)/local/include -I$(HOME)/tmp/MUMPS_4.8.4/include

OPTFLAGS=-DNDEBUG -msse2 -O3 -fweb -fwhole-program -ffast-math -fassociative-math -freciprocal-math  -ffinite-math-only
DEBUGFLAGS=-Wextra -g -O -fbounds-check -fstack-protector-all
# GPPFLAGS=$(INCLUDEPATH) $(LDFLAGS) -ltaucs -lsuperlu -DEIGEN_SUPERLU_SUPPORT -DEIGEN_TAUCS_SUPPORT
GPPFLAGS=$(INCLUDEPATH) $(LDFLAGS) -DEIGEN_SUPERLU_SUPPORT 

all: cppspin

inv_general_complex_mat.o: Makefile inv_general_complex_mat.f
	$(CC) -c $(CFLAGS) inv_general_complex_mat.f

#spin_optics1.o:
#	$(CC) -c $(CFLAGS) spin_optics1.f.f

spin: nano0903c.f inv_general_complex_mat.o Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) -fbounds-check -o spin  inv_general_complex_mat.o nano0903c.f DBL_COMBO_generatorjs1.f
#	$(CC) $(CFLAGS) $(LDFLAGS) -o spin  inv_general_complex_mat.o nano0903c.f

cppspin: spin.cpp Makefile math-utils.h 
	g++ $(GPPFLAGS) -Wall $(DEBUGFLAGS) -o cppspin spin.cpp
#	g++ $(GPPFLAGS) -Wall $(OPTFLAGS) -o cppspin spin.cpp

vspin: spin.cpp Makefile math-utils.h visualize.h
	g++ $(GPPFLAGS) -Wall $(DEBUGFLAGS) -DVISUALIZE `imlib2-config --libs` -o vspin spin.cpp

t: test.cpp Makefile math-utils.h visualize.h
	g++ $(GPPFLAGS) `imlib2-config --libs` -Wall $(DEBUGFLAGS) -o t test.cpp
clean:
	rm -f *.o *~ spin cppspin

.PHONY: clean

# vim: tw=0
