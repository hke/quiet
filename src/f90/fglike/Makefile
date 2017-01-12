all : fglike_multi 

fglike  : fglike.o 
	$(MPF90) -o fglike fglike.o $(LINK)

fglike.o : fglike.f90 ../include/libquiet.a


fglike_coreplus  : fglike_coreplus.o 
	$(MPF90) -o fglike_coreplus fglike_coreplus.o $(LINK)

fglike_coreplus.o : fglike_coreplus.f90 ../include/libquiet.a

fglike_multi  : fglike_multi.o 
	$(MPF90) -o fglike_multi fglike_multi.o $(LINK)

fglike_multi.o : fglike_multi.f90 ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ fglike_multi

