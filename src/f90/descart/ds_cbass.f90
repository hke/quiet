program descart_cbass

!use mpi
use ds_types
use ds_fitstools
use ds_multidetector
use ds_solver
use inifile
use ds_types
use ds_utils
use ds_simple_prior
use ds_cbass_fitstools
use ds_cbass_option_utils
implicit none

character(256) :: arg
integer :: arg_count
  integer(i4b)       :: i, j, k, l, m, n, q, p
  integer(i4b)       :: myid, numprocs, ierr, root
  integer(i4b)       :: nside, nmaps
  character(len=512) :: message
  type(ds_cbass_options) :: options

  !descart specific objects
  type(ds_modulescan), pointer, dimension(:) :: offsetTarget, offsetSolution
  type(ds_trimap) :: data_maps, offset_maps
  type(ds_correlator) :: correlator
  type(ds_covariance) :: map_covariance
  integer :: nd, d
  
  !which of these do WE actually need?


  real(dp), parameter :: sampling_rate = 50.0d0
  real(dp), parameter :: nyquist = sampling_rate / 2


	character(256) :: offset_dir
	
	


  integer ndet_total, offsetLength
  integer(i8b), allocatable, dimension(:) :: originalIndices
  integer :: totalNpix
character(128) :: cwd
  integer azimuth_correlation_length
  integer azimuth_bins
real(dp) azimuth_correlation_time
logical feedback_set

!type(ds_correlation_data) :: correlation_data

		character(256), dimension(:), allocatable :: file_list
type(fl_string80_list) :: fits_output_cards
character(len=80), dimension(:), pointer :: fits_output_cards_array

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

!Set the default feedback parameter
ds_global_feedback=ds_feedback_debug


!#############################################!
!##		Setup parameters					##!
!#############################################!


! Read the command line.
! If no args just print some help.
arg_count=iargc()
if (arg_count==0) then
	write(*,*) 'Syntax (version 1): descart_cbass scan1.fits scan2.fits ...'
	write(*,*) 'Syntax (version 2): descart_cbass params.ini'
	write(*,*) 'For version 1, all parameters will have their default value.'
	write(*,*) 'For version 2, the list of files will be given by the parameter file_list from the ini file'
	write(*,*) 'Version 2 should be used for everything except very basic testing.'
	write(*,*) 'See default_params.ini for an explanation of the parameters'
	stop
endif
call get_command_argument(1,arg)



!Load parameters
if (arg_count>1 .or. string_ends_with(arg,'.fits') ) then  
	!Use default arguments and assume the command line contains a list of FITS files to process
	!This is probably not a good idea in most cases, because you cannot change any of the parameters.
	call ds_log("Using default options",ds_feedback_noisy)
	call set_default_options(options)
	options%file_list=trim(arg)
	allocate(file_list(arg_count))
	do i=1,arg_count
		call get_command_argument(i,file_list(i))
	enddo
else
	!The usual case - read parameters from a parameter file.
	!One of these parameters is a file_list, which should contain a list of (presumably) FITS files to process.
	!Although if your system stores the data in a different way (e.g. a database) then that is fine too - 
	!just use that to record what data to look up.
	call ds_log("Reading options from: "//trim(arg),ds_feedback_quiet)
	call read_options(arg,options)
	call read_file_list(options%file_list, file_list)
endif

!Set the global feedback parameter.
!This parameter controls how many log messages are sent to the stdout.
ds_global_feedback = options%verbosity

!Setup the MPI communicator.
correlator%comm = MPI_COMM_WORLD
correlator%traditionalMode = options%traditional_mode

if (myid==root) call ds_log_milestone("MODULE_SETUP")


!#############################################!
!##		Load data							##!
!#############################################!


call readDataAssignWork(file_list,offsetTarget,correlator,originalIndices,data_maps, map_covariance,options, fits_output_cards)
if (myid==root) call ds_log_milestone("DATA_LOAD_COMPLETE")

!optionally output map pixel and noise information
if(correlator%proc==0) then
	if(options%save_covariance) call saveCovariance(map_covariance,options%covariance_filename)
endif


!#############################################!
!##		Destripe the data					##!
!#############################################!

if(options%naive_mode) then
	if (myid==0) call ds_log("IN NAIVE MODE! NO PCG.",ds_feedback_quiet)
else 
	! Run destriping
	
	! Set up the noise prior.  If traditional mode is on this does nothing.
	if (.not. options%data_prior) call prepare_prior(correlator,offsetLength,nyquist,offsetTarget)
	if (myid==root) call ds_log_milestone("PRIOR_PREPARED")
	if (myid==root) call ds_log("Ready for PCG.", ds_feedback_quiet)

	! Find the maximum-probability offsets.
	! We may have set up a cut removing any particularly divergent scans after a few iterations
	! but probably not.
	! This subroutine fills in the offsetSolution structure.
	if (options%do_sigma_cut) then
		call PCG(correlator,offsetTarget,offsetSolution,map_covariance,options%pcg_tolerance,options%sigma_cut_value, options%sigma_cut_time)
	else
		call PCG(correlator,offsetTarget,offsetSolution,map_covariance,options%pcg_tolerance)
	endif
	if (myid==root) call ds_log_milestone("PCG_COMPLETE")
	
	! Reclaim the memory used in the target data since we now have the solution.
	call destroy_modulescan(offsetTarget,deallocateoffsets=.true.)

	! We need to wait for all the processes to finish destriping to make maps.
	call mpi_barrier(correlator%comm,ierr)
	if (myid==0) call ds_log("Completed PCG",ds_feedback_quiet)

	! If desired, save the offsets to files.
	! This can be nice for diagnostics and visulaization
	if (options%save_offsets) call saveAllOffsetsToFiles(offsetSolution,options%offset_dir)

	! Make a naive map from all the calculated offsets.
	! This map represents the contribution to the naive maps of correlated noise.
	call makeNaiveMap(correlator,offsetSolution,offset_maps,map_covariance,options%do_temperature,options%do_polarization)

	!The destriped map is the difference between the naive map and the offset map.
	!Subtract to get the destriped maps.
	call subtractTriMap(data_maps,offset_maps)
endif



!#############################################!
!##		Save the result 					##!
!#############################################!
if(correlator%proc==0) then
	if (options%ring_ordering) then
		call ds_log("Saving map in RING ordering",ds_feedback_quiet)
	else
		call ds_log("Saving map in NESTED ordering",ds_feedback_quiet)		
	endif

	!Add the options used to the FITS cards to be put in the output map.
	call options_to_fits_cards(options,fits_output_cards)
	fits_output_cards_array => fl_string80_list_to_array(fits_output_cards) 
	
	
	if (options%plane_map) then
		call ds_write_planar_map(options%output_filename,data_maps,options%plane_size,options%plane_size,options%do_temperature, options%do_polarization)
	else
		!Write the map as a healpix map
		call ds_write_fits_map(options%output_filename, data_maps,options%nside, options%partial_out,units='unknown',isRing=options%ring_ordering, extra_headers=fits_output_cards_array)
	endif
	deallocate(fits_output_cards_array)
endif

!Finish up.
call MPI_Finalize(ierr)
if (myid==root) call ds_log_milestone("END")

contains 

subroutine read_file_list(filename, file_list)
	!Read a list of strings from a file into an array.
	character(*) :: filename
	character(256), dimension(:), allocatable :: file_list
	character(512) :: line
	integer unit,io
	
	
	!Count the number of lines in the file by reading each in turn until end of file
	!Count only lines the are 
	unit=ds_get_lun()
	open(unit=unit,file=filename)
	n=0
	do
		read(unit,'(A)',iostat=io) line
		if (io<0) exit
		if (trim(line) .ne. "") then
			line = adjustl(line)
			if(line(1:1) .ne. '#') n=n+1
		endif
	enddo
	rewind(unit)
	
	call ds_assert(n>0,"No uncommented written lines found in file: "//trim(filename))
	
	!
	allocate(file_list(n))
	n=1
	do
		read(unit,'(A)',iostat=io) line
		if (io<0) exit
		if (trim(line) .ne. "") then
			line = adjustl(line)
			if(line(1:1) .ne. '#') then
				file_list(n)=trim(line)
				n=n+1
			endif
		endif
	enddo

end subroutine read_file_list


end program
