OBJS = gibcl_utils.o gibcl_noise_mod.o gibcl.o

all : gibcl

gibcl : $(OBJS) ../include/libquiet.a
	$(MPF90) -o gibcl $(OBJS) $(LINK) 
#gibcl : gibcl.o ../include/libquiet.a
#	$(MPF90) -pg -o gibcl gibcl.o $(LINK) 

gibcl.o : gibcl.f90 ../include/libquiet.a

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<
#	$(MPF90) $(F90COMP) $(OPENMP) -I../map2cl/ -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ gibcl

