TOOLS    = sla.o fits_helper.o ds_precision.o ds_utils.o ds_types.o ds_maptools.o ds_mapping.o ds_fitstools.o \
           random.o inifile.o

OBJ      = ds_multidetector.o ds_simple_prior.o ds_solver.o \
           ds_solver.o fl_lists.o ds_oslo_options.o ds_oslotools.o 

all : descart_oslo

descart_oslo : libdescart.a ds_oslo.o
	$(MPF90) -o descart_oslo ds_oslo.o libdescart.a  $(LINK) 

quick:
	$(MPF90) -o descart_oslo $(F90COMP) sla.f fits_helper.f90 ds_precision.f90 ds_utils.F90 ds_types.f90 ds_maptools.f90 ds_mapping.f90 ds_fitstools.f90 \
	           random.f90 inifile.f90 ds_multidetector.f90 ds_simple_prior.F90 ds_solver.f90  fl_lists.f90 ds_oslo_options.f90 ds_oslotools.F90 $(LINK)

test : ../include/libquiet.a ../tod2map/libtod2map.a libdescart.a test.o
	 $(MPF90) -o test test.o libdescart.a ../tod2map/libtod2map.a $(LINK)

libdescart.a : $(TOOLS) $(OBJ)
	$(AR) $(ARFLAGS) libdescart.a $(TOOLS) $(OBJ)
	$(RANLIB) libdescart.a

fits_helper.f90: make_helper.py
	python make_helper.py

ds_oslo_options.f90: build_oslo_options.py
	python build_oslo_options.py

fl_lists.f90: build_list_type.py
	python build_list_type.py > fl_lists.f90

%.o : %.F90
	$(MPF90) $(F90COMP)  -c  $<

%.o : %.f90
	$(MPF90) $(F90COMP)  -c $<

%.o : %.f
	$(MPF77) $(FCOMP)  -c $<

%.o : %.c
	$(MPCC) $(CCOMP)  -c $< 

clean :
	@rm -f *.o *.mod *.MOD *.a *~ descart_oslo
	@rm -f fits_helper.f90 ds_oslo_options.f90 fl_lists.f90 default_params.ini
