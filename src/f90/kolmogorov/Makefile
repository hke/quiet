TARGETS := kolmogorov

all : $(TARGETS)

%: %.o ../include/libquiet.a
	$(MPF90) -o $@ $< $(LINK)

%.o : %.f90
	$(MPF90) -c $(F90COMP) $(OPENMP) $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(TARGETS)
