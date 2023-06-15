.SUFFIXES: .out .o .e .r .f .for .y .l .s .p .e

.f.o: ; f77 -O -c $*.f

all:
	make scalefree.i

scalefree.i: scalefree.o
	f77 -O -o scalefree.e scalefree.o
