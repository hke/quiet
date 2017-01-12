TARGETS := minisindex
OBJS := minisindex.o 

all: $(TARGETS)

minisindex : $(OBJS) ../include/libquiet.a
	$(F90) -o $@ $(OBJS)  $(LINK)

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ minisindex

