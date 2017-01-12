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
  use ds_azimuth
  use ds_oslo_fitstools
  use ds_simple_prior
  use ds_multidetector
  use ds_solver
  use ds_fitstools
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
  character(len=256) :: infile, outfile, database_file, module_list_file, message

  type(ds_detpointing), dimension(:), allocatable :: pointing
  type(ds_timestream), dimension(:), allocatable :: timestream
  type(ds_timestream),dimension(:),allocatable :: tempTimestream
  type(ds_offsets), dimension(:), allocatable :: offsetTarget, offsetSolution
  integer :: nd, d
  logical :: naive_mode
  character(*), parameter :: outputName = "map.fits"

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

  real(dp), allocatable, dimension(:) :: alpha, sigma,fknee,sigma2
  real(dp), allocatable, dimension(:,:) :: rho
  integer, allocatable, dimension(:) :: ntod

  real(dp), parameter :: sampling_rate = 50.0d0
  real(dp), parameter :: nyquist = sampling_rate / 2
  integer ndet_total, offsetLength
  type(ds_correlator) :: correlator
  character(*), parameter :: pixelFile="pixels.dat"
  integer(i4b), allocatable, dimension(:) :: originalIndices
  integer :: totalNpix
  type(ds_map) :: qmap, umap, qmap2, umap2
character(128) :: cwd
  logical traditionalMode
  integer azimuth_correlation_length
  integer azimuth_bins
real(dp) azimuth_correlation_time
logical feedback_set

type(ds_correlation_data) :: correlation_data

	logical should_cut_timestreams
	real(dp) cut_fraction

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
call getcwd(cwd)
  root = 0
  unit        = 20+myid

  ds_global_feedback=ds_feedback_debug

  if (iargc() /= 1) then
     if (myid == root) then
        call ds_log('Usage: mpirun -n N descart_oslo [paramfile]',ds_feedback_silent) !always prints
     end if
     call mpi_finalize(ierr)
     stop
  else if (myid == root) then
     
     call ds_log('************ DESCART mapmaker ***************',ds_feedback_silent) ! always prints

  end if

!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Get name of parameter file
  call getarg(1,parfile_common)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD2MAP',   par_string=parfile_tod2map)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD_PROC',  par_string=parfile_tod_proc)
  call get_parameter(unit, parfile_common,  'SCAN_DATABASE',          par_string=database_file)
  call get_parameter(unit, parfile_common,  'LEVEL2_DIR',             par_string=level2_dir)
  call get_parameter(unit, parfile_common,  'EXCLUDED_SCANS',         par_string=excluded_scan_list)
  call get_parameter(unit, parfile_common,  'ACCEPTED_SCANS',         par_string=accepted_scan_list)
  call get_parameter(unit, parfile_common,  'REPROCESS_ALL_FILES',    par_lgt=reprocess)
  call get_parameter(unit, parfile_common,  'INCLUDE_TEMPERATURE',    par_lgt=include_temperature)
  call get_parameter(unit, parfile_common,  'INCLUDE_POLARIZATION',   par_lgt=include_polarization)
  call get_parameter(unit, parfile_common,  'MODULE_LIST',            par_string=module_list_file)
 call get_parameter(unit, parfile_tod2map, 'NSIDE_OUT',              par_int=nside)
  call get_parameter(unit, parfile_tod2map, 'ORDERING_OUT',           par_int=ordering)
  call get_parameter(unit, parfile_tod2map, 'PRODUCE_CES_MAPS',       par_lgt=produce_CES_maps)
  call get_parameter(unit, parfile_tod2map, 'TARGET_TYPE',            par_string=target_type)
  call get_parameter(unit, parfile_tod2map, 'TARGET_NAME',            par_string=target_name)
  call get_parameter(unit, parfile_tod2map, 'OFFSET_LENGTH',            par_int=offsetLength)
  call get_parameter(unit, parfile_tod2map, 'NAIVE_MODE',            par_lgt=naive_mode)
  call get_parameter(unit, parfile_tod2map, 'TRADITIONAL_MODE',            par_lgt=traditionalMode)
  call get_parameter(unit, parfile_tod2map, 'AZ_CORR_TIME',            par_dp=azimuth_correlation_time)
  call get_parameter(unit, parfile_tod2map, 'AZ_NBINS',            	  par_int=azimuth_bins)


  call get_parameter(unit, parfile_tod2map, 'CUT_FRACTION',            par_dp=cut_fraction, par_present=should_cut_timestreams)
  call get_parameter(unit, parfile_tod2map, 'DESCART_FEEDBACK',            par_int=ds_global_feedback, par_present=feedback_set)

  if (.not. feedback_set) ds_global_feedback=ds_feedback_noisy
azimuth_correlation_length = azimuth_correlation_time * sampling_rate

  call read_module_list(unit, module_list_file, module_list)
if (myid == root) call ds_log('Loaded module list from ' //trim(module_list_file),ds_feedback_quiet)

	call initialize_quiet_noise_mod(unit, parfile_tod2map, parfile_tod_proc)
	if (myid == root) call ds_log('Loaded noise module ',ds_feedback_quiet)
    call initialize_quiet_calib_mod(unit, parfile_common)
	if (myid == root) call ds_log('Loaded calibration module ',ds_feedback_quiet)
	call initialize_quiet_l2_db_mod(unit, database_file, level2_dir, cwd,accepted_scan_list)

if (myid == root) then
	write(message,*) 'Loaded database from ' //trim(database_file)
	call ds_log(message,ds_feedback_quiet)
endif

  npix  = 12*nside**2
  nmaps = 3
 correlator%traditionalMode = traditionalMode
if (naive_mode) traditionalMode=.true.
  ! Read list of modules to process

num_subsets_per_target = 1

if (myid==0 .and. ds_global_feedback .ge. ds_feedback_quiet) call ds_logTime("Started loading at:")

if (should_cut_timestreams) then
	call readDataAssignWork(pointing,timestream,correlator,&
	produce_CES_maps,level2_dir, target_name, target_type, module_list,accepted_scan_list, &
	nside,alpha,sigma,fknee,offsetLength,offsetTarget,azimuth_correlation_length,azimuth_bins,cut_fraction)
else
	call readDataAssignWork(pointing,timestream,correlator,&
	produce_CES_maps,level2_dir, target_name, target_type, module_list,accepted_scan_list, &
	nside,alpha,sigma,fknee,offsetLength,offsetTarget,azimuth_correlation_length,azimuth_bins)
endif



!write(*,*) myid, " finished loading data"
write(message,*) myid, "Loaded all data.  ndet = ", correlator%ndet
call ds_log(message,ds_feedback_noisy) !more verbose output

ndet_total = correlator%ndet
!pointings and TOD are now set up 
!allocate(rho(0:ndet_total-1,0:ndet_total-1))
!rho=0.0
!do i=0,ndet_total-1
!   rho(i,i)=1.0
!enddo
if (myid==0) then 
	write(message,*) "Renumbering Pointing", myid
	call ds_log(message,ds_feedback_quiet)
endif

totalNpix = 12*nside*nside
call makeSharedPointingContinuous(correlator,pointing,totalNpix,pixelFile,originalIndices)
npix = pointing(0)%npix

if (myid==0) then 
	call ds_log("Shared pointing complete.",ds_feedback_quiet)
	write(message,*) "Total pixels hit = ", npix
	call ds_log(message,ds_feedback_quiet)
	if (ds_global_feedback .ge. ds_feedback_quiet) call ds_logTime("Started work at:")
endif

!#warning !setting noise info by hand
!sigma= 130.0_8
!fknee= 0.1_8
!alpha= -1.0_8

!call prepareFFTQ2(correlator,offsetLength,ndet_total,MPI_COMM_WORLD,sigma,fknee,alpha,correlation_data,nyquist,.true.)
call d_prepare_fftq(correlator,offsetLength,ndet_total,MPI_COMM_WORLD,sigma,fknee,alpha,nyquist,offsetTarget,.true.)



if (myid==0) call ds_log("completed covariance preparations",ds_feedback_quiet)
nd = correlator%my_ndet

!do d=0,nd-1
!	timestream(d)%timestream = timestream(d)%timestream*1e6
!enddo

!allocate(tempTimestream(0:nd-1))
allocate(offsetSolution(0:nd-1))
!jz offsettarget is now allocated and populated in the readDataAssignWork, above.

do d=0,nd - 1
   call ds_assert(correlator%ntod(correlator%my_det(d))==timestream(d)%nt, "Length mismatch timestream data")
!   call prepareOffsets(offsetTarget(d), correlator%ntod(correlator%my_det(d))/offsetLength, offsetLength)
!Not really	!JZ Here is where we add azimuth information
   call prepareOffsets(offsetSolution(d), offsetTarget(d))

enddo


call prepareMap(qmap,npix)
call prepareMap(umap,npix)
call prepareMap(qmap2,npix)
call prepareMap(umap2,npix)
if(myid==0) then
    call ds_assert(size(qmap%indices)==size(originalIndices),"map size not equal to hit pixels size")
    qmap%indices=originalIndices
    umap%indices=originalIndices
    qmap2%indices=originalIndices
    umap2%indices=originalIndices

    deallocate(originalIndices)
endif


call MPI_Barrier(correlator%comm,ierr)


write(message,*) myid, " Removing Naive Signal"
if (myid==0) then
	call ds_log(message,ds_feedback_noisy)
else
	call ds_log(message,ds_feedback_debug)
endif
call removeSignalNaive(timestream,pointing,qmap,umap,correlator)
!A nice side effect of this is that the qmap and umap variables now contain naive maps

write(message,*) myid, " done removal"
call ds_log(message,ds_feedback_debug)

do d=0,nd-1
    call deprojectTimestreamOntoOffset(timestream(d),offsetTarget(d))
    offsetTarget(d)%values= offsetTarget(d)%values / (sigma(correlator%my_det(d))**2)
!    call destroyTimestream(tempTimestream(d))
enddo

write(message,*) myid, " done projection"
call ds_log(message,ds_feedback_debug)

if (naive_mode) then
	if (myid==0) call ds_log("IN NAIVE MODE! NO PCG.",ds_feedback_quiet)
else  !Normal operation mode - PCG
	if (myid==0) call ds_log("Ready for PCG.", ds_feedback_quiet)

	allocate(sigma2(0:correlator%ndet-1))
	sigma2=sigma*sigma*regularisationParameter

	call PCG(correlator, pointing,sigma2, npix, offsetTarget, offsetSolution)

call mpi_barrier(correlator%comm,ierr)
	if (myid==0) call ds_log("Completed PCG",ds_feedback_quiet)


!At this point the timestreams have had the naive map removed.
!We now subtract the offsets so they have had both the offsets and the naive map removed.
	do d=0,nd-1
!       call prepareTimestream(tempTimestream(1),correlator%ntod(correlator%my_det(d)))
!       call projectOffsetOntoTimestream(offsetSolution(d),tempTimestream(1))
!       do i= 1,timestream(d)%nt
!          timestream(d)%timestream(i) = timestream(d)%timestream(i) - tempTimestream(1)%timestream(i)
!       enddo
!       call destroyTimestream(tempTimestream(1))

    !DWPS: we don't need the tempTimestream here because projectOffsetOntoTimestream() now either
    !adds or subtracts offsets to/from the input timestream, depending on the supplied string 
!    call projectOffsetOntoTimestream(offsetSolution(d),timestream(d),'subtract')
	call subtractOffsetFromTimestream(offsetSolution(d), timestream(d))
    enddo
    
	if (myid==0) call ds_log("Subtracted offsets", ds_feedback_quiet)
    
 endif
 
!call destroy_azimuth_cache()

!naive(timestream-offsets) = naive(timestream-project(naive(timestream)) - offsets) + naive(timestream)
call makeNaiveMap(timestream,qmap2,umap2,pointing,correlator)
!call d_test_naive(timestream,qmap2,umap2,pointing,correlator) !remember to comment out removeSignalNaive
call addMap(qmap2,qmap)  !add the naive maps together, into qmap
call addMap(umap2,umap)  !add the naive map together, into umap


if (correlator%proc == 0) then
	call savePolMaps(outputName, qmap,umap,nside) 
endif

if (myid==0 .and. ds_global_feedback .ge. ds_feedback_quiet) call ds_logTime("Ended work at:")



call MPI_Finalize(ierr)


end program tod2map
