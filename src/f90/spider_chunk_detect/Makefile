all: spider_chunk_detect

spider_chunk_detect : spider_chunk_detect.o ../include/libquiet.a
	$(MPF90) -o spider_chunk_detect spider_chunk_detect.o $(LINK)

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ spider_chunk_detect
