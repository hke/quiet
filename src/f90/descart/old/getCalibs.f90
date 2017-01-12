program tod2map
  use healpix_types
  use pix_tools
  use fitstools
  use quiet_utils
  use l2_fileutils
  use quiet_fileutils
  use quiet_calib_mod
  use quiet_noise_mod
  use ds_types
  use ds_oslo_fitstools
  use ds_multidetector
  implicit none

!  include "mpif.h"

  integer(i4b)       :: i, j, k, l, m, n, q, p
  integer(i4b)       :: num_scans, num_targets, file_count
  integer(i4b)       :: unit, myid, numprocs, ierr, root
  integer(i4b)       :: nside, npix, nmaps, ordering
  integer(i4b)       :: nside_lowres, npix_lowres, frac, num_subsets_per_target
  logical(lgt)       :: reprocess, exist, segment_exist, produce_CES_maps
  logical(lgt)       :: include_temperature, include_polarization, temperature
  logical(lgt)       :: assume_white_noise, skip_scan, internal_filtering
  real(dp)           :: max_nhits_ratio, eigenvalue_threshold
  character(len=3)   :: file_count_text
  character(len=128) :: parfile_common, parfile_tod2map, parfile_tod_proc
  character(len=128) :: level2_dir, target_name, target_type, excluded_scan_list, accepted_scan_list
  character(len=256) :: infile, outfile, database_file, module_list_file

  type(ds_detpointing), dimension(:), allocatable :: pointing
  type(ds_timestream), dimension(:), allocatable :: timestream

  logical(lgt), allocatable, dimension(:,:)     :: my_mask, mask
  integer(i4b), allocatable, dimension(:)       :: numpix, map2mask
  integer(i4b), pointer,     dimension(:)       :: modules, module_list
  integer(i4b), pointer,     dimension(:,:,:)   :: pix
  real(dp),     allocatable, dimension(:,:)     :: map, nhits, my_nhits
  real(dp),     allocatable, dimension(:,:)     :: covar, my_covar, eigenvectors
  real(dp),     allocatable, dimension(:)       :: rhs, my_rhs, eigenvals
  real(dp),     pointer,     dimension(:)       :: time
  real(dp),     pointer,     dimension(:,:,:)   :: psi
  real(dp),     pointer,     dimension(:,:,:)   :: tod
  logical(lgt),              dimension(9)       :: output_object
  character(len=256), pointer, dimension(:)     :: excluded_scans

  type(l2_filelist),   dimension(:), pointer :: target_list
  type(module_struct), dimension(:), pointer :: data, data_red

  real(dp), allocatable, dimension(:) :: alpha, sigma,fknee
  real(dp), allocatable, dimension(:,:) :: rho
  integer, allocatable, dimension(:) :: ntod

  real(dp), parameter :: nyquist = 50.0d0 / 2
  integer ndet_total, offsetLength
  type(ds_correlator) :: correlator
  character(*), parameter :: pixelFile="pixels.dat"
  integer(i4b), allocatable, dimension(:) :: originalIndices
  integer :: totalNpix

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0
  unit        = 20+myid

  if (iargc() /= 1) then
     if (myid == root) then
        write(*,*) 'Usage: mpirun -n N tod2map [paramfile]'
     end if
     call mpi_finalize(ierr)
     stop
  else if (myid == root) then
     
     write(*,*) '************ DESCART mapmaker ***************'
     write(*,*)

  end if

!  write(*,*) myid, "hit barrier"
!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!  if (myid == root) write(*,*) ' All processors have reported in. Commencing work.'
!  write(*,*) "I crossed the barrier!", myid, ierr
!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!  write(*,*) "I crossed the barrier again!", myid
  ! Get name of parameter file
  call getarg(1,parfile_common)
!write(*,*) "Z:", myid
  call get_parameter(unit, parfile_common, 'PARFILE_TOD2MAP',   par_string=parfile_tod2map)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD_PROC',  par_string=parfile_tod_proc)
 !write(*,*) "1 one: A", myid
  ! Get parameters
  call get_parameter(unit, parfile_common,  'SCAN_DATABASE',          par_string=database_file)
  call get_parameter(unit, parfile_common,  'LEVEL2_DIR',             par_string=level2_dir)
  call get_parameter(unit, parfile_common,  'EXCLUDED_SCANS',         par_string=excluded_scan_list)
  call get_parameter(unit, parfile_common,  'ACCEPTED_SCANS',         par_string=accepted_scan_list)
  call get_parameter(unit, parfile_common,  'REPROCESS_ALL_FILES',    par_lgt=reprocess)
 !write(*,*) "1 one: B", myid
  call get_parameter(unit, parfile_common,  'INCLUDE_TEMPERATURE',    par_lgt=include_temperature)
  call get_parameter(unit, parfile_common,  'INCLUDE_POLARIZATION',   par_lgt=include_polarization)
  call get_parameter(unit, parfile_common,  'MODULE_LIST',            par_string=module_list_file)
!write(*,*) "1 one: C", myid
 
 call get_parameter(unit, parfile_tod2map, 'NSIDE_OUT',              par_int=nside)
  call get_parameter(unit, parfile_tod2map, 'ORDERING_OUT',           par_int=ordering)
  call get_parameter(unit, parfile_tod2map, 'PRODUCE_CES_MAPS',       par_lgt=produce_CES_maps)
  call get_parameter(unit, parfile_tod2map, 'TARGET_TYPE',            par_string=target_type)
  call get_parameter(unit, parfile_tod2map, 'TARGET_NAME',            par_string=target_name)
  call get_parameter(unit, parfile_tod2map, 'OFFSET_LENGTH',            par_int=offsetLength)
  
  
  


  npix  = 12*nside**2
  nmaps = 3

 ! write(*,*) myid, "read params"

  ! Read list of modules to process
  call read_module_list(unit, module_list_file, module_list)
  if (myid == root) write(*,*) '  Number of modules considered in analysis = ', size(module_list)
  !write(*,*) myid, "read modules"

  ! Read list of runs to exclude
 ! call read_excluded_scan_list(unit, excluded_scan_list, excluded_scans)

  ! Set up auxilliary modules
  call initialize_quiet_noise_mod(unit, parfile_tod2map, parfile_tod_proc)
  call initialize_quiet_calib_mod(unit, parfile_common)

 ! write(*,*) myid, "init1"
 ! write(*,*) trim(database_file)
  !write(*,*) trim(level2_dir)
  call initialize_quiet_db_mod(unit, database_file, level2_dir)
  call initialize_accepted_scan_list(unit, accepted_scan_list)
!write(*,*) trim(target_type), "  ", trim(target_name), "  ", trim(level2_dir)
num_subsets_per_target = 1
  !write(*,*) myid, "init2"



!call get_target_list(produce_CES_maps, num_subsets_per_target,target_type, target_name, level2_dir, target_list)
!write(*,*) trim(target_list(1)%filenames(1))
!write(*,*) myid, "Starting load"
correlator%comm = MPI_COMM_WORLD
call readDataAssignWork(pointing,timestream,correlator,produce_CES_maps,level2_dir, target_name, target_type, module_list,nside,alpha,sigma,fknee)
ndet_total = correlator%ndet
!pointings and TOD are now set up 
allocate(rho(0:ndet_total-1,0:ndet_total-1))
rho=0.0
do i=0,ndet_total-1
   rho(i,i)=1.0
enddo
write(*,*) "Loaded all data.  Renumbering Pointing", myid
totalNpix = 12*nside*nside
call makeSharedPointingContinuous(correlator,pointing,totalNpix,pixelFile,originalIndices)
npix = pointing(0)%npix
if (myid==0) call ds_log("Shared pointing complete.",ds_feedback_quiet)
if (myid==0) write(*,*) "Total pixels hit = ", npix
if (myid==0) call prepareFFTQ2(correlator,offsetLength,ndet_total,MPI_COMM_WORLD,sigma,fknee,alpha,rho,nyquist,.false.)
write(*,*) "prepareFFTQ complete"
call MPI_Finalize(ierr)


end program tod2map
