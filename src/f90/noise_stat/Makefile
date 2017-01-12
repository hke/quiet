targets: noise_stat cross_stat

noise_stat : noise_stat.o ../include/libquiet.a
	$(MPF90) -o noise_stat noise_stat.o $(LINK)
cross_stat : cross_stat.o ../include/libquiet.a
	$(MPF90) -o cross_stat cross_stat.o $(LINK)

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	icc $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ noise_stat
