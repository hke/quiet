#
# Adapted from SVN setup written by Ted Kisner
#
# QUIET Makefile for use with EITHER the planck
# module build system OR a set of stand-alone
# config files.  Default make target prints help.
#
# This makefile requires the use of GNU make or
# compatible tools that allow the substitution
# operator ":=".
#
#

NOW := $(shell date +%Y%m%d-%H%M)
DIR := quiet_$(NOW)

TOPDIR := $(shell pwd)
export TOPDIR
export PYTHONPATH := "$(TOPDIR)/src/python:$(PYTHONPATH)"

ifdef TARGET
	ifndef PREFIX
		PREFIX := ./$(TARGET)
	endif
	INSTALL := $(PREFIX)/quiet
else
	ifdef QUIET
		include $(TOPDIR)/config/config.$(QUIET)
		ifndef INSTALL
	INSTALL := $(TOPDIR)/install_$(QUIET)
		endif
	else
		$(error QUIET undefined)UNDEFINED
	endif
	ifndef MAKE
		export MAKE := make
	endif
	ifndef AR
		export AR := ar
	endif
	ifndef ARFLAGS
		export ARFLAGS := crs
	endif
	ifndef RANLIB
		export RANLIB := ranlib
	endif
	ifndef CTAR
		export CTAR := tar czvf
	endif
	ifndef F90
		export F90 := f90
	endif
	ifndef MPF90
		export MPF90 := mpif90
	endif
	ifndef F90FLAGS
		export F90FLAGS := -g -O2
	endif
	ifndef MPF77
		export MPF77 := mpif77
	endif
	ifndef FFLAGS
		export FFLAGS := -g -O2
	endif
	ifndef MPCC
		export MPCC := cc
	endif
	ifndef CFLAGS
		export CFLAGS := -O2
	endif
	ifndef LDFLAGS
		export LDFLAGS := -lm
	endif
	ifndef FORTRAN_UPPER
		export FORTRAN_UPPER := 0
	endif
	ifndef CFITSIO_LINK
		export CFITSIO_LINK := -L/usr/local/lib -lcfitsio
	endif
	ifndef LAPACK_LINK
		export LAPACK_LINK := -L/usr/local/lib -llapack -lblas
	endif
	ifndef F90FLAGS_DEBUG
		export F90FLAGS_DEBUG := -g -O0
	endif
	ifndef FFLAGS_DEBUG
		export FFLAGS_DEBUG := -g -O0
	endif
	ifndef CFLAGS_DEBUG
		export CFLAGS_DEBUG := -g -O0
	endif
	ifndef LDFLAGS_DEBUG
		export LDFLAGS_DEBUG := -lm
	endif
endif

ifdef DEBUG
	export F90FLAGS := $(F90FLAGS_DEBUG)
	export FFFLAGS  := $(FFFLAGS_DEBUG)
	export CFLAGS   := $(CFLAGS_DEBUG)
	export LDFLAGS  := $(LDFLAGS_DEBUG)
endif

export CCOMP := $(CFLAGS)  -I$(TOPDIR)/src/f90/include $(HEALPIX_INCLUDE) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE)

export F90COMP := $(F90FLAGS) -I$(TOPDIR)/src/f90/include $(HEALPIX_INCLUDE) $(LAPACK_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE) $(HDF_INCLUDE)
export FCOMP := $(FFLAGS) -I$(TOPDIR)/src/f90/include $(LAPACK_INCLUDE) $(HEALPIX_INCLUDE) $(CFITSIO_INCLUDE) $(FFTW_INCLUDE)
export LINK := -L$(TOPDIR)/src/f90/include -lquiet $(HEALPIX_LINK) $(CFITSIO_LINK) $(LAPACK_LINK) $(FFTW_LINK) $(NOVAS_LINK) $(LDFLAGS) $(HDF_LINK) $(LDFLAGS) $(OPENMP)
export TEMPITA := "$(TOPDIR)/src/python/tempita_proc.py"

all : libquiet ces_detect spider_chunk_detect l2gen l3gen postmap tod2map comap map_editor ephcom libutil utils_f90 point_fit ces_validate noise_tester point_scan dipmod maptool #modulation l3planck utils #glitch sindex descart print_powspec smooth_covar_mat patch_sel

full : all libquietscala scalapost map2cl

jz : full install

help :
	@echo ' '
	@echo '  This Makefile is used to build Quiet in a way that is'
	@echo '  customized to your system.  You must export the QUIET'
	@echo '  environment variable and set it to the name of your platform.'
	@echo '  Then you must create a config file named config/config.<$QUIET>.'
	@echo '  I suggest copying the example config file and modifying the'
	@echo '  parameters for your platform.'
	@echo ' '
	@echo '  The following make targets are supported:'
	@echo ' '
	@echo '    make         : build everything'
	@echo '    make help    : print this help screen'
	@echo '    make install : install everything'
	@echo '    make clean   : remove build files'
	@echo '    make dist    : construct a date-stamped source tarball'
	@echo ' '

install : all
	@mkdir -p $(INSTALL)/lib
	@mkdir -p $(INSTALL)/include
	@mkdir -p $(INSTALL)/bin
	@cp src/f90/include/libquiet.a $(INSTALL)/lib
	@cp src/f90/utils/convert_covar_mat $(INSTALL)/bin
	@cp src/f90/point_fit/minisindex $(INSTALL)/bin
	@cp src/f90/point_fit/point_fit $(INSTALL)/bin
	@cp src/f90/point_scan/point_scan $(INSTALL)/bin
	@cp src/f90/dipmod/dipmod $(INSTALL)/bin
	@cp src/f90/smooth_covar_mat/smooth_covar_mat $(INSTALL)/bin
	@cp src/f90/postmap/postmap $(INSTALL)/bin
	@cp src/f90/map2cl/map2cl $(INSTALL)/bin
	@cp src/f90/scalapost/scalapost $(INSTALL)/bin
	@cp src/f90/map_editor/map_editor $(INSTALL)/bin
	@cp src/f90/data_editor/data_editor $(INSTALL)/bin
	@cp src/f90/tod2map/tod2map $(INSTALL)/bin
	@cp src/f90/comap/comap $(INSTALL)/bin
	@cp src/f90/print_powspec/print_powspec $(INSTALL)/bin
	@cp src/tcl/quiet_scan_viewer/quiet_scan_viewer.tcl $(INSTALL)/bin/quiet_scan_viewer
	@cp src/cpp/utils/{l2txt,l2patch,mjd,eph,idate,map2png,map2tga_pol,mjd2unix,plot,scantime,transform,fft,powspec,tcat,gp2tab,patch2img,rmdiode,cut2accept} $(INSTALL)/bin
	@if test $(FORTRAN_UPPER) = 1; then \
	cp src/f90/include/*.MOD $(INSTALL)/include; \
	else \
	cp src/f90/include/*.mod $(INSTALL)/include; \
	fi
	@cp src/f90/ces_validate/ces_validate $(INSTALL)/bin
	@cp src/f90/noise_tester/noise_tester $(INSTALL)/bin
	@cp src/f90/l3gen/l3gen $(INSTALL)/bin
	@cp src/f90/l3planck/l3planck $(INSTALL)/bin
	@cp src/f90/ces_detect/ces_detect $(INSTALL)/bin
	@cp src/f90/spider_chunk_detect/spider_chunk_detect $(INSTALL)/bin	
	@cp src/f90/descart/descart_oslo $(INSTALL)/bin

docs :
	@cd docs; $(MAKE)

libquiet :
	@cd src/cpp/ephcom;  $(MAKE)
	@cd src/f90/include; $(MAKE)

libquietscala:
	@cd src/f90/include; $(MAKE) libquietscala.a

tester :
	@cd src/f90/tester; $(MAKE)

gibcl :
	@cd src/f90/gibcl; $(MAKE)

noise_stat :
	@cd src/f90/noise_stat; $(MAKE)

minisindex :
	@cd src/f90/minisindex; $(MAKE)

point_fit :
	@cd src/f90/point_fit; $(MAKE)

point_scan :
	@cd src/f90/point_scan; $(MAKE)

dipmod :
	@cd src/f90/dipmod; $(MAKE)

smooth_covar_mat :
	@cd src/f90/smooth_covar_mat; $(MAKE)

postmap :
	@cd src/f90/postmap; $(MAKE)

fglike :
	@cd src/f90/fglike; $(MAKE)

map2cl :
	@cd src/f90/map2cl; $(MAKE)

scalapost :
	@cd src/f90/scalapost; $(MAKE)

sharptest :
	@cd src/f90/sharptest; $(MAKE)

rmscomp :
	@cd src/f90/rmscomp; $(MAKE)

map_editor :
	@cd src/f90/map_editor; $(MAKE)

data_editor :
	@cd src/f90/data_editor; $(MAKE)

comap :
	@cd src/f90/comap; $(MAKE)

ces_validate :
	@cd src/f90/ces_validate; $(MAKE)

noise_tester :
	@cd src/f90/noise_tester; $(MAKE)

l3gen :
	@cd src/f90/l3gen; $(MAKE)

l3planck :
	@cd src/f90/l3planck; $(MAKE)

ces_detect :
	@cd src/f90/ces_detect; $(MAKE)

spider_chunk_detect :
	@cd src/f90/spider_chunk_detect; $(MAKE)

l2gen :
	@cd src/f90/l2gen; $(MAKE)

print_powspec :
	@cd src/f90/print_powspec; $(MAKE)

descart : 
	@cd src/f90/descart; $(MAKE)

utils: ephcom libutil
	@cd src/cpp/utils; $(MAKE)

libutil:
	@cd src/cpp/libutil; $(MAKE)

utils_f90: 
	@cd src/f90/utils; $(MAKE)

ephcom:
	@cd src/cpp/ephcom; $(MAKE)

l2info :
	@cd src/f90/l2info; $(MAKE)

sindex :
	@cd src/f90/sindex; $(MAKE)

rings :
	@cd src/f90/rings; $(MAKE)

l1condense :
	@cd src/f90/l1condense; $(MAKE)

tod_gibbs :
	@cd src/f90/tod_gibbs; $(MAKE)

kolmogorov: 
	@cd src/f90/kolmogorov; $(MAKE)

crosscorr: 
	@cd src/f90/crosscorr; $(MAKE)


patch_sel:
	@cd src/f90/patch_sel; $(MAKE)

maptool:
	@cd src/f90/maptool; $(MAKE)

modulation:
	@cd src/f90/modulation; $(MAKE)

mdfit:
	@cd src/f90/mdfit; $(MAKE)

bscan:
	@cd src/f90/bscan; $(MAKE)

glitch:
	@cd src/f90/glitch; $(MAKE)

clean : clean_libquiet clean_point_fit clean_minisindex clean_point_scan clean_dipmod clean_postmap clean_map2cl clean_scalapost clean_map_editor clean_tod2map clean_comap clean_ephcom clean_libutil clean_utils clean_descart clean_ces_validate clean_noise_tester clean_l3gen clean_l3planck clean_ces_detect clean_spider_chunk_detect clean_l2gen clean_utils_f90 clean_maptool clean_modulation clean_mdfit clean_rmscomp clean_bscan clean_fglike

clean_minisindex :
	@cd src/f90/minisindex; $(MAKE) clean

clean_point_fit :
	@cd src/f90/point_fit; $(MAKE) clean

clean_point_scan :
	@cd src/f90/point_scan; $(MAKE) clean

clean_dipmod :
	@cd src/f90/dipmod; $(MAKE) clean

clean_smooth_covar_mat :
	@cd src/f90/smooth_covar_mat; $(MAKE) clean

clean_postmap :
	@cd src/f90/postmap; $(MAKE) clean

clean_fglike :
	@cd src/f90/fglike; $(MAKE) clean

clean_map2cl :
	@cd src/f90/map2cl; $(MAKE) clean

clean_scalapost :
	@cd src/f90/scalapost; $(MAKE) clean

clean_rmscomp :
	@cd src/f90/rmscomp; $(MAKE) clean

clean_map_editor :
	@cd src/f90/map_editor; $(MAKE) clean

clean_data_editor :
	@cd src/f90/data_editor; $(MAKE) clean

clean_sindex :
	@cd src/f90/sindex; $(MAKE) clean

clean_tod2map :
	@cd src/f90/tod2map; $(MAKE) clean

clean_comap :
	@cd src/f90/comap; $(MAKE) clean

clean_utils_f90 :
	@cd src/f90/utils; $(MAKE) clean

clean_ces_validate :
	@cd src/f90/ces_validate; $(MAKE) clean 

clean_noise_tester :
	@cd src/f90/noise_tester; $(MAKE) clean 

clean_l3gen :
	@cd src/f90/l3gen; $(MAKE) clean 

clean_l3planck :
	@cd src/f90/l3planck; $(MAKE) clean 

clean_ces_detect :
	@cd src/f90/ces_detect; $(MAKE) clean 

clean_spider_chunk_detect :
	@cd src/f90/spider_chunk_detect; $(MAKE) clean 

clean_l2gen :
	@cd src/f90/l2gen; $(MAKE) clean 

clean_print_powspec :
	@cd src/f90/print_powspec; $(MAKE) clean

clean_descart :
	@cd src/f90/descart; $(MAKE) clean; 

clean_libquiet :
	@cd src/f90/include; $(MAKE) clean

clean_utils:
	@cd src/cpp/utils; $(MAKE) clean

clean_libutil :
	@cd src/cpp/libutil; $(MAKE) clean

clean_ephcom:
	@cd src/cpp/ephcom; $(MAKE) clean

clean_tod_gibbs :
	@cd src/f90/tod_gibbs; $(MAKE) clean

clean_kolmogorov :
	@cd src/f90/kolmogorov; $(MAKE) clean

clean_crosscorr :
	@cd src/f90/crosscorr; $(MAKE) clean

clean_patch_sel :
	@cd src/f90/patch_sel; $(MAKE) clean

clean_maptool :
	@cd src/f90/maptool; $(MAKE) clean

clean_modulation :
	@cd src/f90/modulation; $(MAKE) clean

clean_mdfit :
	@cd src/f90/mdfit; $(MAKE) clean

clean_bscan :
	@cd src/f90/bscan; $(MAKE) clean

clean_glitch :
	@cd src/f90/glitch; $(MAKE) clean

clean_tester :
	@cd src/f90/tester; $(MAKE) clean

dist : clean
	@mkdir $(DIR)
	@mkdir -p $(DIR)/src/f90/include
	@mkdir -p $(DIR)/src/f90/minisindex
	@mkdir -p $(DIR)/src/f90/point_fit
	@mkdir -p $(DIR)/src/f90/point_scan
	@mkdir -p $(DIR)/src/f90/dipmod
	@mkdir -p $(DIR)/src/f90/smooth_covar_mat
	@mkdir -p $(DIR)/src/f90/postmap
	@mkdir -p $(DIR)/src/f90/fglike
	@mkdir -p $(DIR)/src/f90/map2cl
	@mkdir -p $(DIR)/src/f90/scalapost
	@mkdir -p $(DIR)/src/f90/rmscomp
	@mkdir -p $(DIR)/src/f90/map_editor
	@mkdir -p $(DIR)/src/f90/data_editor
	@mkdir -p $(DIR)/src/f90/tod2map
	@mkdir -p $(DIR)/src/f90/comap
	@mkdir -p $(DIR)/src/f90/ces_validate
	@mkdir -p $(DIR)/src/f90/noise_tester
	@mkdir -p $(DIR)/src/f90/l3gen
	@mkdir -p $(DIR)/src/f90/l3planck
	@mkdir -p $(DIR)/src/f90/l2gen
	@mkdir -p $(DIR)/src/f90/ces_detect
	@mkdir -p $(DIR)/src/f90/spider_chunk_detect
	@mkdir -p $(DIR)/src/f90/print_powspec
	@mkdir -p $(DIR)/src/f90/descart
	@mkdir -p $(DIR)/src/f90/utils
	@mkdir -p $(DIR)/src/cpp/utils
	@mkdir -p $(DIR)/src/cpp/ephcom
	@cp -r config Makefile $(DIR)
	@cp src/f90/include/*.f90 src/f90/include/Makefile $(DIR)/src/f90/include
	@cp src/f90/minisindex/*.f90 src/f90/minisindex/Makefile $(DIR)/src/f90/minisindex
	@cp src/f90/point_fit/*.f90 src/f90/point_fit/Makefile $(DIR)/src/f90/point_fit
	@cp src/f90/point_scan/*.f90 src/f90/point_scan/Makefile $(DIR)/src/f90/point_scan
	@cp src/f90/dipmod/*.f90 src/f90/dipmod/Makefile $(DIR)/src/f90/dipmod
	@cp src/f90/smooth_covar_mat/*.f90 src/f90/smooth_covar_mat/Makefile $(DIR)/src/f90/smooth_covar_mat
	@cp src/f90/postmap/*.f90 src/f90/postmap/Makefile $(DIR)/src/f90/postmap
	@cp src/f90/fglike/*.f90 src/f90/fglike/Makefile $(DIR)/src/f90/fglike
	@cp src/f90/map2cl/*.f90 src/f90/map2cl/Makefile $(DIR)/src/f90/map2cl
	@cp src/f90/scalapost/*.f90 src/f90/scalapost/Makefile $(DIR)/src/f90/scalapost
	@cp src/f90/rmscomp/*.f90 src/f90/rmscomp/Makefile $(DIR)/src/f90/rmscomp
	@cp src/f90/map_editor/*.f90 src/f90/map_editor/Makefile $(DIR)/src/f90/map_editor
	@cp src/f90/data_editor/*.f90 src/f90/data_editor/Makefile $(DIR)/src/f90/data_editor
	@cp src/f90/tod2map/*.f90 src/f90/tod2map/Makefile $(DIR)/src/f90/tod2map
	@cp src/f90/comap/*.f90 src/f90/comap/Makefile $(DIR)/src/f90/comap
	@cp src/f90/ces_validate/*.f90 src/f90/ces_validate/Makefile $(DIR)/src/f90/ces_validate
	@cp src/f90/noise_tester/*.f90 src/f90/noise_tester/Makefile $(DIR)/src/f90/noise_tester
	@cp src/f90/l3gen/*.f90 src/f90/l3gen/Makefile $(DIR)/src/f90/l3gen
	@cp src/f90/l3planck/*.f90 src/f90/l3planck/Makefile $(DIR)/src/f90/l3planck
	@cp src/f90/ces_detect/*.f90 src/f90/ces_detect/Makefile $(DIR)/src/f90/ces_detect
	@cp src/f90/spider_chunk_detect/*.f90 src/f90/spider_chunk_detect/Makefile $(DIR)/src/f90/spider_chunk_detect	
	@cp src/f90/l2gen/*.f90 src/f90/l2gen/Makefile $(DIR)/src/f90/l2gen
	@cp src/f90/print_powspec/*.f90 src/f90/print_powspec/Makefile $(DIR)/src/f90/print_powspec
	@cp src/f90/descart/*.f90 src/f90/descart/Makefile $(DIR)/src/f90/descart
	@cp src/f90/utils/*.f90 src/f90/utils/Makefile $(DIR)/src/f90/utils
	@cp src/cpp/utils/{*.cpp,*.h,plot,scantime,mjd,transform,patchmont,patchmont,patch2img,rmdiode,Makefile} $(DIR)/src/cpp/utils
	@cp src/cpp/ephcom/{*.c,*.h,Makefile,README.txt,LICENSE.txt} $(DIR)/src/cpp/ephcom
	@rm -rf $(DIR)/config/.svn
	@$(CTAR) $(DIR).tar.gz $(DIR)
	@rm -rf $(DIR)
