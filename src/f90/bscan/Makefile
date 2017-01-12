TARGETS := bscan_skip bscan_viking

all : $(TARGETS)

bscan_skip  : bscan_skip.o
	$(MPF90) -o $@ $^  $(LINK)

bscan_viking  : bscan_viking.o
	$(MPF90) -o $@ $^  $(LINK)

bscan_skip.o : bscan_skip.f90 ../include/libquiet.a

bscan_viking.o : bscan_viking.f90 ../include/libquiet.a

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
