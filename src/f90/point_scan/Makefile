TARGETS := point_scan point_sim

all: $(TARGETS)

point_scan : point_scan.o ../include/libquiet.a
	$(MPF90) -o $@ $< $(LINK)

point_sim : point_sim.o ../include/libquiet.a
	$(MPF90) -o $@ $< $(LINK)

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ $(TARGETS)
