OBJS := rms_mod.o find_pts_mod.o plot_patch_mod.o patch_sel.o

all: patch_sel

patch_sel: $(OBJS)
	$(MPF90) -o $@ $(OBJS) $(LINK) ../include/libquiet.a

%.o: %.f90 ../include/libquiet.a
	$(MPF90) $(F90COMP) -c $<

clean :
	@rm -f *.o *.mod *.MOD *.a *~ 
