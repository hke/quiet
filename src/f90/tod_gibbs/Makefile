TARGETS := tod_gibbs

all : $(TARGETS)

% : %.o
	$(MPF90) -o $@ $< $(LINK) 

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(TARGETS)
