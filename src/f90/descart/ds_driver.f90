program descart_oslo

!use mpi
use ds_types
use ds_fitstools
use ds_azimuth
use ds_multidetector
!use ds_covmat
use ds_solver
use inifile
use ds_types
use ds_utils
use ds_simple_prior
use quiet_utils
use l2_fileutils
use quiet_db_mod
use quiet_fileutils
use quiet_calib_mod
use quiet_noise_mod
use ds_oslo_fitstools
implicit none


  integer(i4b)       :: i!, j, k, l, m, n, q, p
!  integer(i4b)       :: num_scans, num_targets, file_count
  integer(i4b)       :: unit, myid, numprocs, ierr, root
  integer(i4b)       :: nside, npix, nmaps
!  integer(i4b)       :: nside_lowres, npix_lowres, frac, num_subsets_per_target
!  logical(lgt)       :: reprocess, exist, segment_exist, produce_CES_maps
!  logical(lgt)       :: include_temperature, include_polarization!, temperature
!  logical(lgt)       :: assume_white_noise, skip_scan, internal_filtering
!  real(dp)           :: max_nhits_ratio, eigenvalue_threshold
!  character(len=3)   :: file_count_text
  character(len=512) :: parfile
!  character(len=512) :: level2_dir, target_name, target_type
!character(len=512), dimension(1) :: accepted_scan_lists
  character(len=256) :: database_file, module_list_file, message, covFilename, pixelsFile, output_filename
  logical :: output_filename_present

	type(ds_oslo_options) :: options

character(*), parameter :: MC_TYPE_SIGNAL = "signal"
character(*), parameter :: MC_TYPE_NOISE = "noise"
character(*), parameter :: MC_TYPE_SIGNAL_NOISE = "signalnoise"


  !descart specific objects
  type(ds_modulescan), pointer, dimension(:) :: offsetTarget, offsetSolution
  type(ds_map) :: data_qmap, data_umap, offset_qmap, offset_umap
  type(ds_correlator) :: correlator
!  type(ds_covmatrix) :: matrix
  
  integer :: d
  logical :: naive_mode, partial_out
  
  !which of these do WE actually need?
!  logical(lgt), allocatable, dimension(:,:)     :: my_mask, mask
!  integer(i4b), allocatable, dimension(:)       :: numpix, map2mask
  integer(i4b), pointer,     dimension(:)       :: module_list!,modules
!  integer(i4b), pointer,     dimension(:,:,:)   :: pix
!  real(dp),     allocatable, dimension(:,:)     :: map, nhits, my_nhits
!  real(dp),     allocatable, dimension(:,:)     :: covar, my_covar, eigenvectors
 ! real(dp),     allocatable, dimension(:)       :: rhs, my_rhs, eigenvals
!  real(dp),     pointer,     dimension(:)       :: time
!  real(dp),     pointer,     dimension(:,:,:)   :: psi
!  real(dp),     pointer,     dimension(:,:,:)   :: tod
!  logical(lgt),              dimension(9)       :: output_object
!  character(len=256), pointer, dimension(:)     :: excluded_scans

!  type(l2_filelist),   dimension(:), pointer :: target_list
!  type(module_struct), dimension(:), pointer :: data, data_red


  real(dp), parameter :: sampling_rate = 50.0d0
  real(dp), parameter :: nyquist = sampling_rate / 2


	logical :: save_offsets, save_offsets_set, offset_dir_set
	character(256) :: offset_dir
	
	
	real(dp) :: sigma_cut
	integer sigma_when
	logical :: sigma_cut_set, sigma_when_set


  integer offsetLength
  integer(i8b), allocatable, dimension(:) :: originalIndices
!  integer :: totalNpix
character(128) :: cwd
  logical traditionalMode
  integer azimuth_correlation_length
  integer azimuth_bins
real(dp) azimuth_correlation_time
logical feedback_set, covFilePresent,pixelsFilePresent
logical ring_ordering
integer ordering


!logical should_cut_timestreams
type(ds_noiseinfo), pointer, dimension(:) :: noiseInfo
logical data_prior, data_prior_set



!#############################!
!##							##!
!##		Initial set-up		##!
!##							##!
!#############################!

! Initialize MPI environment
call ds_init_milestone()
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




!#####################################!
!##									##!
!##		Read QUIET parameters		##!
!##									##!
!#####################################!

  ! Get name of parameter file
  call getarg(1,parfile)
  call get_parameter(unit, parfile,  'L2_DATABASE',          par_string=database_file)
!  call get_parameter(unit, parfile,  'LEVEL2_DIR',             par_string=level2_dir)
!  call get_parameter(unit, parfile,  'ACCEPTED_SCANS1',         par_string=accepted_scan_list)
!  call get_parameter(unit, parfile,  'REPROCESS_ALL_FILES',    par_lgt=reprocess)
!  call get_parameter(unit, parfile,  'INCLUDE_TEMPERATURE',    par_lgt=include_temperature)
!  call get_parameter(unit, parfile,  'INCLUDE_POLARIZATION',   par_lgt=include_polarization)
  call get_parameter(unit, parfile,  'MODULE_LIST',            par_string=module_list_file)
  call get_parameter(unit, parfile, 'NSIDE_OUT',              par_int=nside)
  call get_parameter(unit, parfile, 'ORDERING_OUT',           par_int=ordering)
  if (ordering==1) then
	ring_ordering=.true.
  else if (ordering==2) then
	ring_ordering=.false.
  else
	stop 'ORDERING_OUT parameters must be 1 (RING) or 2 (NEST)'
	endif
	


!  call get_parameter(unit, parfile, 'TARGET_TYPE',            par_string=target_type)
!  call get_parameter(unit, parfile, 'TARGET_NAME',            par_string=target_name)
  call get_parameter(unit, parfile, 'OFFSET_LENGTH',          par_int=offsetLength)
  call get_parameter(unit, parfile, 'NAIVE_MODE',             par_lgt=naive_mode)
  call get_parameter(unit, parfile, 'TRADITIONAL_MODE',       par_lgt=traditionalMode)
  call get_parameter(unit, parfile, 'AZ_CORR_TIME',           par_dp=azimuth_correlation_time)
  call get_parameter(unit, parfile, 'AZ_NBINS',            	  par_int=azimuth_bins)
  call get_parameter(unit, parfile, 'PARTIAL_OUT',            par_lgt=partial_out)
  call get_parameter(unit, parfile, 'COV_FILENAME',           par_string=covFilename, &
  		par_present=covFilePresent)
  call get_parameter(unit, parfile, 'PIXELS_FILENAME',        par_string=pixelsFile, &
  		par_present=pixelsFilePresent)
  

!  call get_parameter(unit, parfile, 'CUT_FRACTION',            par_dp=cut_fraction, par_present=should_cut_timestreams)
!  call get_parameter(unit, parfile, 'RING_ORDERING',            par_lgt=ring_ordering, par_present=ring_ordering_set)
  call get_parameter(unit, parfile, 'DESCART_FEEDBACK',            par_int=ds_global_feedback, par_present=feedback_set)
  call get_parameter(unit, parfile, 'OUTPUT_MAP_FILENAME',  par_string=output_filename,par_present=output_filename_present)
  call get_parameter(unit, parfile, 'DESCART_SAVE_OFFSETS',  par_lgt = save_offsets, par_present=save_offsets_set)
  call get_parameter(unit, parfile, 'OFFSET_DIR',            par_string=offset_dir, par_present=offset_dir_set)

  call get_parameter(unit, parfile, 'PCG_CUT_SIGMA',            par_dp=sigma_cut, par_present=sigma_cut_set)
  call get_parameter(unit, parfile, 'PCG_CUT_ITERATION',            par_int=sigma_when, par_present=sigma_when_set)

	call get_parameter(unit, parfile, 'DESCART_DATA_PRIOR',            par_lgt=data_prior, par_present=data_prior_set)
	if (.not. data_prior_set) data_prior=.false.
	

  if (.not. offset_dir_set) offset_dir = "saved_offsets"
  if (.not. save_offsets_set) save_offsets=.false.

  
!  if (.not. ring_ordering_set) ring_ordering = .false.
  if (.not. output_filename_present) output_filename = 'map.fits'

  if (.not. feedback_set) ds_global_feedback=ds_feedback_noisy
azimuth_correlation_length = azimuth_correlation_time * sampling_rate

	call read_options(parfile,options)
	call initialize_modules

  call read_module_list(unit, module_list_file, module_list)
  if (myid == root) call ds_log('Loaded module list from ' //trim(module_list_file),ds_feedback_quiet)


	if (myid == root) call ds_log('Loaded calibration module ',ds_feedback_quiet)

if (myid == root) then
	write(message,*) 'Loaded database from ' //trim(database_file)
	call ds_log(message,ds_feedback_quiet)
endif

  npix  = 12*nside**2
  nmaps = 3
 correlator%traditionalMode = traditionalMode
if (naive_mode) traditionalMode=.true.
  ! Read list of modules to process


if (myid==root) call ds_log_milestone("MODULE_SETUP")


!#############################################!
!##											##!
!##		Process Data -> Destriped Map		##!
!##											##!
!#############################################!


call readDataAssignWork(offsetTarget,noiseInfo,correlator,parfile,originalIndices,data_qmap,data_umap)
if (myid==root) call ds_log_milestone("DATA_LOAD_COMPLETE")


!optionally output map pixel and noise information
if(correlator%proc==0) then
	if(covFilePresent .and. covFilename .ne. '') call saveCovariance(covFilename)
	if(pixelsFilePresent .and. pixelsFile .ne. '') call writePixelsFile(pixelsFile,nside,data_qmap)
endif

if(.not. naive_mode) then

	!this routine prepares prior. At the moment it also determines mpi context - why?
	if (.not. data_prior) call prepare_prior(correlator,offsetLength,nyquist,noiseinfo,offsetTarget)
	if (myid==root) call ds_log_milestone("PRIOR_PREPARED")
	
	
	

	if (myid==0) call ds_log("Ready for PCG.", ds_feedback_quiet)

	!this routine iterates to find the max-posterior offsets. Target is replaced with residuals
	if (sigma_cut_set) then
		if (sigma_when_set) then
			call PCG(correlator,offsetTarget,offsetSolution,sigma_cut,sigma_when)
		else 
			call PCG(correlator,offsetTarget,offsetSolution,sigma_cut)
		endif
	else
		call PCG(correlator,offsetTarget,offsetSolution)
	endif
	if (myid==root) call ds_log_milestone("PCG_COMPLETE")
	
	
	call mpi_barrier(correlator%comm,ierr)
	if (myid==0) call ds_log("Completed PCG",ds_feedback_quiet)

!	if (save_offsets) call saveAllOffsetsToFiles(offsetSolution,offset_dir)

	!this routine finds naive map of the projected best offsets.
	call makeNaiveMap(correlator,offsetSolution,offset_qmap,offset_umap)

	!The destriped maps is the difference between the naive map and the naive offsets map.
	call subtractMap(data_qmap,offset_qmap)
	call subtractMap(data_umap,offset_umap)

else
	if (myid==0) call ds_log("IN NAIVE MODE! NO PCG.",ds_feedback_quiet)
endif



!output the map
if(correlator%proc==0) then
	if (ring_ordering) then
		call ds_log("Saving map in RING ordering",ds_feedback_quiet)
	else
		call ds_log("Saving map in NESTED ordering",ds_feedback_quiet)		
	endif
   call savePolMaps(output_filename, data_qmap,data_umap,nside, partial_out,units='mK      ',isRing=ring_ordering)
endif

if (myid==0) call ds_log("Starting MCs", ds_feedback_quiet)
if (options%mc_signal) call mc_simulations(correlator,offsetTarget,MC_TYPE_SIGNAL,options,originalIndices=originalIndices)
if (options%mc_noise) call mc_simulations(correlator,offsetTarget,MC_TYPE_NOISE,options,noise=noiseInfo)
if (options%mc_signalnoise) call mc_simulations(correlator,offsetTarget,MC_TYPE_SIGNAL_NOISE,options,originalIndices,noiseInfo)



	!clean up offset residual vector. Add optional output?
call destroy_modulescan(offsetTarget,deallocateoffsets=.true.)

!!optional calculation of the covariance matrix
!!call pixel_pixel_covariance(corr,nside,nside_dg,modulescans,map,matrix)
!if(correlator%proc==0) print*,'calling pixel_pixel_covariance'
!call pixel_pixel_covariance(correlator,nside,nside,offsetTarget,data_qmap,matrix)
!if(correlator%proc==0) then
!	print*,'ouputting covariance'
!	open(unit=10234,file="covmat.bin",status="replace",form="binary")
!		write(10234) matrix%matrix
!	close(10234)
!endif

!if (myid==0 .and. ds_global_feedback .ge. ds_feedback_quiet) call ds_logTime("Ended work at:")
call MPI_Finalize(ierr)
if (myid==root) call ds_log_milestone("END")

contains 

subroutine saveAllOffsetsToFiles(moduleScans,dir)
        type(ds_modulescan), pointer, dimension(:) :: moduleScans
        character(*) dir
    type(ds_modulescan), pointer :: moduleScan
    integer i
    character(512) :: filename
    do i=0,size(moduleScans)-1
       moduleScan=>moduleScans(i)
                do d=1,ndiodes_max
                filename=filenameForSavedOffset(dir,moduleScan%id(1),moduleScan%id(2),moduleScan%id(3),d)
                        call saveOffsetsToFile(moduleScan%offsets(d),filename)
       enddo
    enddo
end subroutine saveAllOffsetsToFiles

function filenameForSavedOffset(directory,run,seg,module,diode) result(f)
	character(512) :: f
	character(*) :: directory
	character(*), parameter :: fmt = '( A, "/savedOffset_", I5.5, "_", I3.3, "_", I3.3, "_", I2.2, ".off"  )'
	integer :: run, seg, module, diode
	write(f,fmt) trim(directory),run,seg,module,diode
end function filenameForSavedOffset

subroutine initialize_modules()
	use quiet_array_mod
	use quiet_validation_mod
	
!  integer(i4b)     :: base_seed
!  type(planck_rng) :: rng_handle
!  real(dp)         :: mem
	character(512) :: focalplane_info
    call get_parameter(unit, parfile,  'FOCALPLANE_INFO',        par_string=focalplane_info)
if (myid == root) call ds_log_milestone("START_LIBQUIET")
  call initialize_noise_estimation_mod
  call initialize_quiet_array_mod(unit, focalplane_info)
  call initialize_quiet_assembly_mod(unit, parfile)
  call initialize_quiet_validation_mod(unit, parfile)
  call initialize_quiet_noise_mod(unit, parfile)
  call initialize_quiet_pointing_mod(unit, parfile, .false.)
!  call initialize_map_maker(parfile)
  call initialize_quiet_calib_mod(unit, parfile)
  call initialize_gain_mod(unit, parfile)
  call initialize_target_mod(parfile)
  call initialize_quiet_tod_filter_mod(unit, parfile)
  if (myid==root)   call ds_log("Libquiet Initialized",ds_feedback_quiet)
if (myid == root) call ds_log_milestone("END_LIBQUIET")


!  if (apply_i2qu) then
!     if (debug) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', myid,' time: ', t1-t2, ' mem: ', mem, ' -- i2qu module'
!     call initialize_i2qu(unit, parfile)
!     call cpu_time(t1)
!     mem = get_mem_use()/1024d0**3
!  end if
!  if (analyze_simulated_data) then 
!     if (debug) write(*,fmt="(a,i3,a,f8.4,a,f8.4,a)") 'myid = ', myid,' time: ', t1-t2, ' mem: ', mem, ' -- tod sim module'
!     call get_parameter(unit, parfile, 'BASE_SEED', par_int=base_seed)
!     call initialize_random_seeds(MPI_COMM_WORLD, base_seed, rng_handle)
!     base_seed = nint(rand_uni(rng_handle)*1.d7)
!     call initialize_todsim_mod(unit, base_seed, parfile)
!     call cpu_time(t1)
!     mem = get_mem_use()/1024d0**3
!  end if
! call get_targets(categories, objects, splits, targets)
end subroutine initialize_modules

subroutine mc_simulations(correlator,offsetTarget,mc_type,options,originalIndices,noise)
	use ds_montecarlo
	use quiet_fileutils, only : read_map
	
	character(*) :: mc_type
	type(ds_map) :: naive_Q, naive_U
	type(ds_map) :: signal_Q, signal_U
	type(ds_map) :: offset_Q, offset_U
	type(ds_correlator) :: correlator
	type(ds_modulescan), pointer, dimension(:) :: offsetTarget, offsetSolution
	type(ds_oslo_options) :: options
 	character(512) :: output_filename, input_filename
  	integer(i8b), optional, dimension(:) :: originalIndices
	real(dp), pointer, dimension(:,:) :: signal_data
	integer number_maps
	integer signal_nside, signal_ordering
	type(ds_noiseinfo), optional, dimension(:)  :: noise		
	integer npix_scanned, signal_npix
	integer nside
	logical signal_ring_ordering
	nside = options%nside
	npix_scanned = p_npix
	
	if (mc_type==MC_TYPE_NOISE) call ds_assert(present(noise), 'Must pass noise to simulate noise-only')
	if (mc_type==MC_TYPE_SIGNAL_NOISE) call ds_assert(present(noise), 'Must pass noise to simulate signal+noise')
	if (mc_type==MC_TYPE_SIGNAL) call ds_assert(present(originalIndices), 'Must pass original indices to simulate signal')
	
	do i=1,options%mc_niterations

		!If needed, load signal maps.
		if (mc_type==MC_TYPE_SIGNAL .or. mc_type==MC_TYPE_SIGNAL_NOISE) then
			input_filename=filename_for_input_signal(options%signal_input_format,i)
			if (correlator%proc==0) call ds_log("Reading signal map:"//trim(input_filename), ds_feedback_noisy)
			call read_map(signal_data, signal_ordering, input_filename,nside=signal_nside,nmap=number_maps)
			signal_ring_ordering=.false.
			if (signal_ordering==1) signal_ring_ordering=.true.
			signal_npix=nside2npix(signal_nside)
			call ds_assert(signal_nside==nside,"Must use the same resolution signal MC maps as output maps")
			call ds_assert(signal_ring_ordering .eqv. options%ring_ordering,"Must use the same ordering signal MC maps as output maps")

			call prepareMap(signal_q, signal_npix, zero_based=.true.)
			call prepareMap(signal_u, signal_npix, zero_based=.true.)
			
			if (number_maps==2) then
				signal_q%map = signal_data(:,1)
				signal_u%map = signal_data(:,2)
			else if (number_maps==2) then
				signal_q%map = signal_data(:,2)
				signal_u%map = signal_data(:,3)			
			else
				call ds_assert(.false.,'Do not know how to get signal from loaded maps of given size')
			endif
			deallocate(signal_data)
		endif
		
		!Not sure if you can use select case with strings.
		if (mc_type==MC_TYPE_SIGNAL) then
!			subroutine simulate_signal_scans(scans,naive_q,naive_u,comm,nside,originalIndices,signal_q,signal_u)
			write(message,*) "MC signal simulation # ",i
			if (correlator%proc==0) call ds_log(message, ds_feedback_noisy)
			call simulate_signal_scans(offsetTarget,naive_q,naive_u,correlator%comm,nside,originalIndices,signal_q,signal_u)
		else if (mc_type==MC_TYPE_NOISE) then
			write(message,*) "MC noise simulation # ",i
			if (correlator%proc==0) call ds_log(message, ds_feedback_noisy)			
			call simulate_noise_scans(offsetTarget,naive_q,naive_u,correlator%comm,noise)
		else if (mc_type==MC_TYPE_SIGNAL_NOISE) then
			write(message,*) "MC signal+noise simulation # ",i
			if (correlator%proc==0) call ds_log(message, ds_feedback_noisy)
			call simulate_signal_noise_scans(offsetTarget,naive_q,naive_u,correlator%comm,nside,originalIndices,signal_q,signal_u,noise)
		else
			call ds_assert(.false., "Unknown MC type in mc_simulations")
		endif	
		
				
		call PCG(correlator,offsetTarget,offsetSolution)
		call makeNaiveMap(correlator,offsetSolution,offset_Q,offset_U)
		call subtractMap(naive_Q,offset_Q)
		call subtractMap(naive_U,offset_U)
		output_filename=filename_for_saved_mc(options%mc_dir,i,mc_type)
		write(message,*) "Saving MC map to: ",trim(output_filename)
		if (myid==0) call ds_log(message, ds_feedback_noisy)
		if (myid==0) call savePolMaps(trim(output_filename), naive_Q,naive_U, nside, options%partial_output,units='mK      ',isRing=options%ring_ordering)
		call destroyMap(naive_Q)
		call destroyMap(naive_U)
	enddo
	
end subroutine mc_simulations

function filename_for_input_signal(fmt,i) result(f)
	character(512) :: f
	character(*) :: fmt
	integer i
	write(f,fmt) i
end function 

function filename_for_saved_mc(directory,i,mc_type) result(f)
	character(512) :: f
	character(*) :: directory
	character(*), parameter :: fmt = '( A, "/mc_",A,"_", I5.5, ".fits"  )'
	integer i
	character(*) :: mc_type
	write(f,fmt) trim(directory),trim(mc_type),i
end function filename_for_saved_mc


end program
