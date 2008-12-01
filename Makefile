CC=g77
CFLAGS=-g -I/usr/include/ -O2 -Wuninitialized
LDFLAGS=

all: spin cppspin

inv_general_complex_mat.o: Makefile inv_general_complex_mat.f
	$(CC) -c $(CFLAGS) inv_general_complex_mat.f

#spin_optics1.o:
#	$(CC) -c $(CFLAGS) spin_optics1.f.f

spin: nano0903c.f inv_general_complex_mat.o Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) -o spin  inv_general_complex_mat.o nano0903c.f DBL_COMBO_generatorjs1.f
#	$(CC) $(CFLAGS) $(LDFLAGS) -o spin  inv_general_complex_mat.o nano0903c.f

cppspin: spin.cpp
	g++ -O2 -g -Wall -o cppspin spin.cpp
