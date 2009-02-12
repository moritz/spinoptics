CC=gfortran
CFLAGS=-g -I/usr/include/ -O2 -Wuninitialized
LDFLAGS=
INCLUDEPATH=-I/home/mlenz/tmp/eigen2/
GPPFLAGS=$(INCLUDEPATH) -lsuperlu -DEIGEN_SUPERLU_SUPPORT

all: spin cppspin

inv_general_complex_mat.o: Makefile inv_general_complex_mat.f
	$(CC) -c $(CFLAGS) inv_general_complex_mat.f

#spin_optics1.o:
#	$(CC) -c $(CFLAGS) spin_optics1.f.f

spin: nano0903c.f inv_general_complex_mat.o Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) -fbounds-check -o spin  inv_general_complex_mat.o nano0903c.f DBL_COMBO_generatorjs1.f
#	$(CC) $(CFLAGS) $(LDFLAGS) -o spin  inv_general_complex_mat.o nano0903c.f

cppspin: spin.cpp Makefile
	g++ $(GPPFLAGS) -Wall -DNDEBUG -msse2 -O3 -o cppspin spin.cpp
#	g++ $(GPPFLAGS) -Wall -g -O -fbounds-check -o cppspin spin.cpp

clean:
	rm -f *.o *~ spin cppspin

.PHONY: clean
