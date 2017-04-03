### EDIT AS NEEDED ###

# fortran compiler
FC = gfortran

# compiler flags
FCFLAGS = -O3

###


# programs to compile
PROGRAMS = mol2_conf_ensemble

# make (all)
all: $(PROGRAMS)

mol2_conf_ensemble.o: boxio.o mol2io_min.o
mol2_conf_ensemble: boxio.o mol2io_min.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

install:
	mv $(PROGRAMS) ../bin

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD

