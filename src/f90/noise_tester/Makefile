all : noise_tester gen gen2

noise_tester : noise_tester.o ../include/libquiet.a 
	$(MPF90) -o noise_tester noise_tester.o $(LINK) 

gen : gen.o ../include/libquiet.a
	$(MPF90) -o gen gen.o $(LINK) 

gen2 : gen2.o ../include/libquiet.a
	$(MPF90) -o $@ $< $(LINK) 

noise_tester.o   : ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	icc $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ noise_tester

