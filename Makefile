CC=gfortran
CFLAGS=-g -I/usr/include/ -O2 -Wuninitialized

MUMPS=$(HOME)/src/lib/MUMPS_4.8.4/
LDFLAGS= -L$(HOME)/local/lib \
	 -L$(MUMPS)libseq \
	 -L$(MUMPS)PORD/lib \
         $(MUMPS)/lib/libzmumps.a \
         $(MUMPS)/lib/libmumps_common.a \
         -lm  -lsuperlu -lmpiseq -lpthread -lblas -lpord
INCLUDEPATH=-I$(HOME)/tmp/eigen2/ -I/usr/include/superlu \
            -I$(HOME)/local/include -I$(MUMPS)/include
OPTFLAGS=-DNDEBUG -msse2 -O3 -fweb -fwhole-program -ffast-math -fassociative-math -freciprocal-math  -ffinite-math-only
DEBUGFLAGS=-Wextra -g -O -fbounds-check -fstack-protector-all
# GPPFLAGS=$(INCLUDEPATH) $(LDFLAGS) -ltaucs -lsuperlu -DEIGEN_SUPERLU_SUPPORT -DEIGEN_TAUCS_SUPPORT
GPPFLAGS=$(INCLUDEPATH) -Wall -DEIGEN_SUPERLU_SUPPORT 

all: cppspin

inv_general_complex_mat.o: Makefile inv_general_complex_mat.f
	$(CC) -c $(CFLAGS) inv_general_complex_mat.f

#spin_optics1.o:
#	$(CC) -c $(CFLAGS) spin_optics1.f.f

spin: nano0903c.f inv_general_complex_mat.o Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) -fbounds-check -o spin  inv_general_complex_mat.o nano0903c.f DBL_COMBO_generatorjs1.f
#	$(CC) $(CFLAGS) $(LDFLAGS) -o spin  inv_general_complex_mat.o nano0903c.f

cppspin.o: spin.cpp Makefile math-utils.h 
	g++ $(GPPFLAGS) $(DEBUGFLAGS) -c -o cppspin.o spin.cpp

cppspin: cppspin.o Makefile
	g++ $(GPPFLAGS) $(DEBUGFLAGS) -o cppspin cppspin.o $(LDFLAGS)
#	g++ $(GPPFLAGS) $(DEBUGFLAGS) -o cppspin cppspin.o $(LDFLAGS)
#	g++ $(GPPFLAGS) -Wall $(OPTFLAGS) -o cppspin spin.cpp

cppspin2.o: spin2.cpp Makefile math-utils2.h 
	g++ $(GPPFLAGS) -c -o cppspin2.o spin2.cpp

cppspin2: cppspin2.o Makefile
	g++ $(GPPFLAGS) $(DEBUGFLAGS) -o cppspin2 cppspin2.o -lm -lsuperlu
#	g++ $(GPPFLAGS) $(DEBUGFLAGS) -o cppspin2 cppspin2.o $(LDFLAGS)
#	g++ $(GPPFLAGS) -Wall $(OPTFLAGS) -o cppspin2 spin2.cpp

vspin: spin.cpp Makefile math-utils.h visualize.h
	g++ $(GPPFLAGS) -Wall $(DEBUGFLAGS) -DVISUALIZE `imlib2-config --libs` -o vspin spin.cpp

t: test.cpp Makefile math-utils.h visualize.h
	g++ $(GPPFLAGS) `imlib2-config --libs` -Wall $(DEBUGFLAGS) -o t test.cpp $(LDFLAGS)
clean:
	rm -f *.o *~ spin cppspin

.PHONY: clean

# vim: tw=0
