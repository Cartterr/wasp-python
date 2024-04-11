## If you want to change the size of GF bank, please just modify the size.h
FFLAGS=-O
MY_BIN=.
F77=gfortran -ffixed-line-length-132 -mcmodel=medium -O2 -Wall -g
# linux or Mac
modules = constants.o retrieve_gf.o bessel2.o vel_model_data.o wave_travel.o fk_source.o haskell.o prop.o fk_kernel.o fk.o bessel.o layer.o fk_st_kernel.o st_haskell.o fft.o Complex.o 

SUBS = fft.o Complex.o 

all: green_bank_fk

green_bank_fk: $(modules) green_bank_fk.o
	$(F77) -o $(MY_BIN)/green_bank_f95 $(modules) green_bank_fk.o

%.o: %.f95
	$(F77) -c $<

%.o: %.f
	$(F77) -c $<


%.o: %.c
	cc -c $<

# This generic rule handles the creation of .mod files implicitly.
# Explicit module compilation rules are not necessary unless special compilation flags are needed.

#clean:
#	rm -f *.mod
#	rm -f *.o
#	rm -f $(MY_BIN)/green_bank_f95

.PHONY: all clean
