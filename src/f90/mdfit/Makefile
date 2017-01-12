TARGETS := mdfit

all : $(TARGETS)

mdfit  : mdfit.o
#	$(MPF90) -o $@ $^  $(LINK)
	ifort -o $@ $^  $(LINK)

mdfit.o : mdfit.f90 ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(TARGETS) 
