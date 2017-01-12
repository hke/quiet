program descart_oslo

!use mpi
use ds_types
use ds_fitstools
use ds_multidetector
use ds_solver
use inifile
use ds_types
use ds_utils
use ds_simple_prior
use ds_oslotools
use ds_oslo_option_utils
use quiet_target_mod
use quiet_acceptlist_mod
implicit none

character(256) :: arg
integer :: arg_count
  integer(i4b)       :: i, j, k, l, m, n, q, p
  integer(i4b)       :: myid, numprocs, ierr, root
  integer(i4b)       :: nside, nmaps
  character(len=512) :: message
  type(ds_oslo_options) :: options

  !descart specific objects
  type(ds_modulescan), pointer, dimension(:) :: offsetTarget, offsetSolution, allScans
  type(quiet_target), dimension(:), allocatable :: splitTargets
  type(ds_trimap) :: data_maps, offset_maps
  type(ds_correlator) :: correlator
  type(ds_covariance) :: map_covariance
  integer :: nd, d
  real(dp), dimension(:,:), allocatable :: ces_statistics
  real(dp), dimension(:,:,:), allocatable :: diode_statistics


  !which of these do WE actually need?


  real(dp), parameter :: sampling_rate = 50.0d0
  real(dp), parameter :: nyquist = sampling_rate / 2


	character(256) :: offset_dir
	
	


  integer ndet_total, offsetLength
  integer(i8b), allocatable, dimension(:) :: originalIndices
  integer :: totalNpix, npix
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
write(*,*) "Processor", myid+1, " of ", numprocs, " running."
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
	write(*,*) 'Syntax (version 1): descart_oslo scan1.fits scan2.fits ...'
	write(*,*) 'Syntax (version 2): descart_oslo params.ini'
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
!	call read_file_list(options%file_list, file_list)
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


call read_data_assign_work(file_list, allScans, correlator, originalIndices, npix, options, ces_statistics, diode_statistics, fits_output_cards)
!subroutine read_data_assign_work(fileList,scans,correlator,originalIndices,npix,opt,output_fits_cards)

if (myid==root) call ds_log_milestone("DATA_LOAD_COMPLETE")



call build_jackknives(ces_statistics, diode_statistics, splitTargets, options)
!Loop through jack-knives
write(message, *) "Running ",size(splitTargets)," split targets plus one main run."
if (myid==root) call ds_log(message,ds_feedback_quiet)

!First make the unsplit maps

call copy_all_scans(allScans, offsetTarget)
if (myid==root) call ds_log("Copied scans",ds_feedback_quiet)
call generate_naive_and_subtract(correlator%comm, offsetTarget, data_maps, npix, map_covariance, options%covariance_cut, options%do_temperature, options%do_polarization)
if (myid==root) call ds_log("Naive map generated and subtracted",ds_feedback_quiet)
call add_indices_to_map(data_maps, originalIndices)
if(myid==root) call save_maps(data_maps,"naive_"//options%output_filename,options)	
if (myid==root) call ds_log("Naive map saved",ds_feedback_quiet)

if (myid==root) then
open(file='hitcount.txt', unit=51)
do i=1,npix
	write(51,*) map_covariance%QQ(i)
enddo
close(51)
endif

if(.not. options%naive_mode) then
	! Destripe the data.  Hurrah!
	call destripe_data(offsetTarget, data_maps, map_covariance, options)
	if (myid==root) call ds_log("Destriping complete.",ds_feedback_quiet)
	!save the results.
	if (myid==root) then
		call add_indices_to_map(data_maps, originalIndices)
		call save_maps(data_maps, options%output_filename, options)
	endif
endif
call destroy_modulescan(offsetTarget, deallocateoffsets = .true., deallocatepointing=.true.)
if (myid==root) call ds_log("Full map complete.",ds_feedback_quiet)

do i=1,size(splitTargets)
	write(message, *) "Running split-target number ",i
	if (myid==root) call ds_log(message,ds_feedback_quiet)
	!For each null test:
	!construct the array of copies of each scan (copy since we will be molesting the offsets) - set duplicateoffsets to True
	call filter_scans_for_target(splitTargets(i),allScans, offsetTarget)
	
	!call this function to get the new naive map for that half, the covariance, and subtract the naive map from the offsets	
	!so that no signal remains in offset space.  
	call generate_naive_and_subtract(correlator%comm, offsetTarget, data_maps, npix, map_covariance,  options%covariance_cut, options%do_temperature, options%do_polarization)


	if(options%save_covariance) call saveCovariance(map_covariance,options%covariance_filename)
	
	if(correlator%proc==0) call save_maps(data_maps,"naive_"//options%output_filename,options)	
	
	if(.not. options%naive_mode) then
		! Destripe the data.  Hurrah!
		call destripe_data(offsetTarget, data_maps, map_covariance, options)

		!save the results.
		call add_indices_to_map(data_maps, originalIndices)
		call save_maps(data_maps, options%output_filename, options)
		
	endif

	!Clean up.
	call destroy_modulescan(offsetTarget, deallocateoffsets = .true., deallocatepointing=.true.)
enddo

if (size(splitTargets)==0) then
	if (myid==root) call ds_log("No splits to consider.  Run complete.",ds_feedback_quiet)
endif


!










!Finish up.
call MPI_Finalize(ierr)
if (myid==root) call ds_log_milestone("END")

contains 

subroutine add_indices_to_map(maps,indices)
	type(ds_trimap) :: maps
	integer(i8b), dimension(0:) :: indices
	integer p
	
	p=1
	if (maps%has_t) then
		do i=0,size(indices)-1
			if (indices(i) .gt. 0) then
				maps%T%indices(p) = i
				p=p+1
			endif
		enddo
	endif
!	call ds_assert(p-1==maps%T%npix,"Indices not right size for assigned map in add_indices_to_map")

	p=1
	if (maps%has_p) then
		do i=0,size(indices)-1
			if (indices(i) .gt. 0) then
				maps%Q%indices(p) = i
				maps%U%indices(p) = i
				p=p+1
			endif
		enddo
!		call ds_assert(p-1==maps%Q%npix,"Indices not right size for assigned map in add_indices_to_map")
	endif

end subroutine

subroutine destripe_data(offset_target, data_maps, map_covariance, options)
	type(ds_modulescan), pointer, dimension(:) :: offset_target
	type(ds_trimap) :: data_maps
	type(ds_covariance) :: map_covariance
	type(ds_oslo_options) :: options

	type(ds_modulescan), pointer, dimension(:) :: offset_solution
	type(ds_trimap) :: offset_maps

	if (myid==0) call ds_log_milestone("STARTING_PCG")
!	write(*,*) "Size of offsetTarget = ", size(offsetTarget)
	!Destripe!
	if (options%do_sigma_cut) then
		call PCG(correlator,offsetTarget,offset_solution,map_covariance,options%pcg_tolerance,options%sigma_cut_value, options%sigma_cut_time)
	else
		call PCG(correlator,offsetTarget,offset_solution,map_covariance,options%pcg_tolerance)
	endif
	
	! Save the offsets to file, if desired.
	if (options%save_offsets) call saveAllOffsetsToFiles(offset_solution,options%offset_dir)
	
	!Make a naive map from the offsets.  We will subtract this from the naive map of the data to get the destriped naive map.
	call makeNaiveMap(correlator,offset_solution,offset_maps,map_covariance,options%do_temperature,options%do_polarization)
	
	!The destriped map is the difference between the naive map and the offset map.
	!Subtract to get the destriped maps.
	call subtractTriMap(data_maps,offset_maps)
	
	!Clean up.
	call destroyTriMap(offset_maps)
	call destroy_modulescan(offset_solution)
end subroutine destripe_data
	


subroutine save_maps(maps, filename, options)
	type(ds_trimap) :: maps
	character(*) :: filename
	type(ds_oslo_options) :: options
	
!	if (options%ring_ordering) then
!		call ds_log("Saving map in RING ordering",ds_feedback_quiet)
!	else
		call ds_log("Saving map in NESTED ordering",ds_feedback_quiet)		
!	endif

	!Add the options used to the FITS cards to be put in the output map.
	call options_to_fits_cards(options,fits_output_cards)
	fits_output_cards_array => fl_string80_list_to_array(fits_output_cards) 
	
	
!	if (options%plane_map) then
!		call ds_write_planar_map(filename,maps,options%plane_size,options%plane_size,options%do_temperature, options%do_polarization)
!	else !Write the map as a healpix map
		call ds_write_fits_map(filename, maps,options%nside, options%partial_out,units='unknown ',isRing=.false., extra_headers=fits_output_cards_array)
!	endif
	deallocate(fits_output_cards_array)
end subroutine


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


subroutine build_jackknives(ces_statistics, diode_statistics, targets, opt)
	real(dp), dimension(:,:) :: ces_statistics
	real(dp), dimension(:,:,:) :: diode_statistics
	type(quiet_target) :: target
	type(quiet_target),    dimension(:),     allocatable :: targets
	type(acceptlist) :: alist
	type(swiss_knife),     dimension(:),     allocatable :: knife_defs
	integer(i4b),          dimension(:,:,:), allocatable :: knife_res
  
	type(ds_oslo_options) :: opt
	character(512) :: acceptfile, jackknives, target_name
	
	call get_parameter(0, opt%quiet_parameter_file, 'TARGET_NAME',par_string=target_name)
    call get_parameter(0, opt%quiet_parameter_file, 'JACKKNIVES', par_string=jackknives)
	call get_parameter(0, opt%quiet_parameter_file, 'ACCEPTLIST', par_string=acceptfile)

	if (myid==0) call ds_log("Jackkinfe specification: "//trim(jackknives), ds_feedback_noisy)
	call initialize_accept_list(acceptfile, alist)
	call init_target(target, alist)
	call filter_object(target, target_name)
	!call get_accepted_ceses(target%alist, cnums)


  call jackknife(target, jackknives, ces_statistics, diode_statistics, targets, knife_defs, knife_res)
!  call print_knives(knife_defs, knife_res, cid_list)
  deallocate(knife_res)


!	call jackknife(target, jackknives, ces_statistics, targets)
	
	!We now have a set of targets with different accept lists in targets(:)%alist.
end subroutine

subroutine copy_all_scans(all_scans, scans)
	type(ds_modulescan), dimension(0:)  :: all_scans
	type(ds_modulescan), pointer, dimension(:)  :: scans
	integer i
	
	allocate(scans(0:size(all_scans)-1))
	do i=0,size(all_scans)-1
		call setup_moduleScan_ndiodes(scans(i), all_scans(i)%ndiodes_t, all_scans(i)%ndiodes_p)
		call copy_modulescan(all_scans(i),scans(i),duplicateoffsets=.true., duplicatemaps=.true., duplicatepointing=.true.)
	enddo
	
end subroutine copy_all_scans

subroutine filter_scans_for_target(target, all_scans, scans)
	type(quiet_target) :: target
	type(ds_modulescan), dimension(0:)  :: all_scans
	type(ds_modulescan), pointer, dimension(:)  :: scans
	
	integer i
	integer n
	
	!Count the number of scans accepted by this target.
	n=0
	do i=0,size(all_scans)-1
		if (is_accepted(target%alist,all_scans(i)%id(1),all_scans(i)%id(3))) n=n+1
	enddo
	allocate(scans(0:n-1))
	n=0
	do i=0,size(all_scans)-1
		if (is_accepted(target%alist,all_scans(i)%id(1),all_scans(i)%id(3))) then
			call setup_moduleScan_ndiodes(scans(n), all_scans(i)%ndiodes_t, all_scans(i)%ndiodes_p)
			call copy_modulescan(all_scans(i),scans(n),duplicateoffsets=.true., duplicatemaps=.true., duplicatepointing=.true.)
			n=n+1
		endif
	enddo
	
	
end subroutine filter_scans_for_target





end program
