TARGETS := point_fit
OBJS := point_fit_mod.o point_fit_params.o point_fit.o

all: $(TARGETS)

point_fit : $(OBJS) ../include/libquiet.a
	$(MPF90) -o $@ $(OBJS)  $(LINK)

point_fit_mod.o: point_fit_params.o

%.o : %.f90
	$(MPF90) $(F90COMP) -c $<

%.o : %.f
	$(MPF77) $(FCOMP) -c $<

%.o : %.c
	$(MPCC) $(CCOMP) -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ point_fit

