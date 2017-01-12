program test_oslo

use ds_types
use ds_fitstools
use ds_multidetector
use ds_solver
use inifile
!use mpi
use ds_types
use ds_utils
use quiet_utils
use l2_fileutils
use quiet_db_mod
use quiet_fileutils
use quiet_calib_mod
use quiet_noise_mod
use ds_oslo_fitstools

implicit none

type(ds_detpointing), dimension(:), allocatable :: pointing
type(ds_timestream), dimension(:), allocatable :: timestream
integer rank,nproc
logical(lgt) ::  CES
character(len=128) :: level2_dir, target_name, target_type, accepted_scan_list
character(len=128) :: parfile_tod2map, parfile_tod_proc, database_file, parfile_common

character(len=256) ::  module_list_file
integer ierror, nside
real(dp), allocatable, dimension(:) :: alpha, sigma,fknee
integer, allocatable, dimension(:) :: ntod
integer(i4b) unit
type(l2_filelist),   pointer, dimension(:) :: target_list
integer(i4b), parameter :: num_subsets_per_target=1
call MPI_Init(ierror)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierror)
call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)


  call getarg(1,parfile_common)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD2MAP',   par_string=parfile_tod2map)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD_PROC',  par_string=parfile_tod_proc)
  call get_parameter(unit, parfile_common,  'SCAN_DATABASE',        par_string=database_file)
  call get_parameter(unit, parfile_common,  'LEVEL2_DIR',           par_string=level2_dir)
  call get_parameter(unit, parfile_tod2map, 'TARGET_TYPE',          par_string=target_type)
  call get_parameter(unit, parfile_tod2map, 'TARGET_NAME',          par_string=target_name)
  call get_parameter(unit, parfile_common,  'ACCEPTED_SCANS',         par_string=accepted_scan_list)

  call get_parameter(unit, parfile_common,  'MODULE_LIST',          par_string=module_list_file)
!  call get_parameter(unit, parfile_descart, 'NSIDE_OUT',            par_int=nside)

  ! Get parameters

  call initialize_quiet_noise_mod(unit, parfile_tod2map, parfile_tod_proc)
!  call initialize_temperature_map_maker(unit, parfile_tod2map)  
 ! call initialize_polarization_map_maker(unit, parfile_tod2map)  
  call initialize_quiet_calib_mod(unit, parfile_common)
  call initialize_quiet_db_mod(unit, database_file, level2_dir)
  call initialize_accepted_scan_list(unit, accepted_scan_list)

  call get_target_list(CES, num_subsets_per_target, target_type, target_name, level2_dir, target_list)
write(*,*) size(target_list), lbound(target_list)
write(*,*) target_list%numfiles, size(target_list(1)%filenames), lbound(target_list(1)%filenames)
write(*,*) trim(target_list(1)%filenames(1))

!call readDataAssignWork(pointing,timestream,rank,nproc,CES,level2_dir, target_name, target_type, module_list_file,nside,alpha,sigma,fknee,ntod)

!call prepareFFTQ(correlator,offsetLength,ntod,ndet,MPI_COMM_WORLD,sigma,fknee,alpha)


end program test_oslo
