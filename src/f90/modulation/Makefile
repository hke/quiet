TARGETS := modulation gen_covar_mat

all : $(TARGETS)

modulation  : modulation_utils.o modulation_mcmc.o modulation_evidence.o modulation_chisq.o modulation.o
	$(MPF90) -o $@ $^ $(LINK)
#	$(MPF90) -o $@ $^ $(SCALAPACK_LIBS) $(LINK)

gen_covar_mat  : gen_covar_mat.o 
	$(MPF90) -o $@ $^ $(LINK)

scalapost.o : scalapost.f90 ../include/libquietscala.a ../include/libquiet.a

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
