OBJS = newton_minimizing_mod.o sindex_maxlike_mod.o sindex.o

all : sindex multisindex sindex_utils convert_WMAP9_noise convert_Planck_noise

sindex  : $(OBJS)
	$(MPF90) -o sindex $(OBJS) $(LINK)

multisindex  : multisindex.o 
	$(MPF90) -o multisindex multisindex.o $(LINK)

sindex_utils  : sindex_utils.o sindex_maxlike_mod.o 
	$(MPF90) -o sindex_utils sindex_utils.o sindex_maxlike_mod.o $(LINK)

convert_WMAP9_noise  : convert_WMAP9_noise.o 
	$(MPF90) -o convert_WMAP9_noise convert_WMAP9_noise.o $(LINK)

convert_Planck_noise  : convert_Planck_noise.o 
	$(MPF90) -o convert_Planck_noise convert_Planck_noise.o $(LINK)

sindex.o : sindex.f90 ../include/libquiet.a

multisindex.o : multisindex.f90 ../include/libquiet.a

sindex_utils.o : sindex_utils.f90 ../include/libquiet.a

convert_WMAP9_noise.o : convert_WMAP9_noise.f90 ../include/libquiet.a
convert_Planck_noise.o : convert_Planck_noise.f90 ../include/libquiet.a

%.o : %.F90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ 

