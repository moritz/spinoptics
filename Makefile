CC=gfortran
CFLAGS=-g -I/usr/include/ -O2 -Wuninitialized
LDFLAGS= -L$(HOME)/local/lib
INCLUDEPATH=-I$(HOME)/tmp/eigen2-git/ -I/usr/include/superlu -I$(HOME)/local/include
OPTFLAGS=-DNDEBUG -msse2 -O3 -fweb -fwhole-program -ffast-math -fassociative-math -freciprocal-math  -ffinite-math-only
DEBUGFLAGS=-Wextra -g -O -fbounds-check -fstack-protector-all
# GPPFLAGS=$(INCLUDEPATH) $(LDFLAGS) -ltaucs -lsuperlu -DEIGEN_SUPERLU_SUPPORT -DEIGEN_TAUCS_SUPPORT
GPPFLAGS=$(INCLUDEPATH) $(LDFLAGS) -lsuperlu -DEIGEN_SUPERLU_SUPPORT

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

t: test.cpp Makefile math-utils.h
	g++ $(GPPFLAGS) -Wall $(DEBUGFLAGS) -o t test.cpp
clean:
	rm -f *.o *~ spin cppspin

.PHONY: clean

# vim: tw=0
