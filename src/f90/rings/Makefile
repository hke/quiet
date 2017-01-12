all : rings statring

rings  : rings.o 
	$(MPF90) -o rings rings.o $(LINK)
statring  : statring.o 
	$(MPF90) -o statring statring.o $(LINK)

rings.o : rings.f90 ../include/libquiet.a

statring.o : statring.f90 ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ rings

