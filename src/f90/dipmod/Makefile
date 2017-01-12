all : dipmod dipmod_postproc dipmod_gencov

dipmod : libdipmod.a dipmod.o ../include/libquiet.a
	$(MPF90) -o dipmod dipmod.o libdipmod.a $(LINK) $(COMMON_FLAGS)
libdipmod.a : dipmod_utils.o dipmod_powspec.o dipmod_chisq.o dipmod_mcmc.o dipmod_evidence.o 
	ar crs $@ $^
dipmod_postproc : dipmod_postproc.o ../include/libquiet.a
	$(MPF90) -o dipmod_postproc dipmod_postproc.o $(LINK) $(COMMON_FLAGS)
dipmod_gencov : dipmod_gencov.o ../include/libquiet.a
	$(MPF90) -o dipmod_gencov dipmod_gencov.o $(LINK) $(COMMON_FLAGS)

dipmod.o   : libdipmod.a
dipmod_chisq.o : dipmod_evidence.o
dipmod_evidence.o : dipmod_utils.o

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $< -Os

clean :
	@rm -f *.o *.mod *.MOD *.a *~ dipmod dipmod_gencov dipmod_postproc

