all : comap

comap : libcomap.a comap.o ../include/libquiet.a
	$(MPF90) -o comap comap.o libcomap.a $(LINK) $(COMMON_FLAGS)

libcomap.a : comap_utils.o 
	ar crs $@ $^

comap.o : libcomap.a comap_utils.o

%.o : %.f90
	$(MPF90) $(F90COMP) $(OPENMP) -c $< -Os

clean :
	@rm -f *.o *.mod *.MOD *.a *~ comap

