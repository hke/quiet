! This file is the main one to change when creating a new executable for a new experiment.
!
! Its job is to read in the data and get it in the format used by the rest of the code.
! It also makes a naive map as it goes along.  Since there are a number of pieces of auxiliary information
! that need to be loaded in this code is a little complicated.

module ds_oslotools
	use ds_types
	use ds_utils
	use ds_simple_prior
	use healpix_types
	use pix_tools
	use ds_multidetector
	use ds_oslo_option_utils
	use quiet_lx_mod
#ifndef __GFORTRAN__
	use ieee_arithmetic
#endif
	use fl_lists
	use quiet_fileutils
	implicit none
	
	
	!JZ These are temporary, till we read what all the numbers are.
	integer, parameter :: QUIET_NDIODES_MAX = 4  !maximum number of diodes per module.
	integer, parameter :: QUIET_N_MODULES = 91  !maximum number of diodes per module.
	integer, parameter :: QUIET_NDIODES_P = 4  !maximum number of diodes per module.
	integer, parameter :: QUIET_NDIODES_T = 0  !maximum number of diodes per module.
	

	real(dp), parameter :: DEGRA = 0.0174532925
	integer, parameter :: FITS_READ_ONLY = 0
	character(8), parameter :: NAXIS2_KEYWORD = "NAXIS2  "
	
	integer rank

	
#if 0
#pragma mark -
#pragma mark DATA TYPES
#endif	
	type ds_moduleScanInfo
		integer id1,id2,id3
		integer ndiodes_T, ndiodes_P
		integer owner
		integer partner_module
		integer ndiodes
		integer, dimension(QUIET_NDIODES_MAX) :: diode_flag
		logical needs_azfilter
		character(len=256) filename
		integer modScan !The local index into moduleScans
		type(ds_moduleScanInfo),pointer :: next
		integer n
	end type ds_moduleScanInfo




	type ds_moduleScanList
		type(ds_moduleScanInfo), pointer :: first,last
		integer length
	end type ds_moduleScanList
	

	
	
	

	character(len=125), private :: message

	contains



#if 0	
#pragma mark -
#pragma mark PRIMARY MODULE FUNCTION
#pragma mark -
#endif
	subroutine read_data_assign_work(fileList,scans,correlator,originalIndices,npix,opt,ces_statistics, diode_statistics, output_fits_cards)
#if 0	
#pragma mark -
#endif
	! fileList - input, array of strings of length 256
	! scans - output, array of (module,scan) pairs (ds_modulescan)
	! correlator - input/output.  object containing data on which processor owns what, and noise correlations.  On input, should have %comm set.  Other stuff set on output
	! originalIndices - output.  contains map from our pixel numbering to healpix.
	! opt - input.  The options file.
	! output_fits_cards - output.  List of cards to propagate to the FITS header of the output map.
	

	!This is the only subroutine that you need to re-write to use 
		use quiet_ces_mod

		implicit none
		!Subroutine arguments
		type(ds_moduleScan), pointer, dimension(:) :: scans
		type(ds_correlator) :: correlator
		integer(i8b), allocatable, dimension(:) :: originalIndices
		type(ds_oslo_options) :: opt
		integer(i8b) :: full_map_pixels
		character(256), dimension(:) :: fileList
		type(fl_string80_list) :: output_fits_cards

		logical, parameter ::  do_correlations = .false.
		integer progress

		!Information about the data to be loaded.
		integer nmodules
		integer(i4b) npix

		!Work assignment information
		integer nproc
		integer nmin,nextra
		integer nScan,my_nScan

		!Things needed to load the data from file
		type(Lx_struct) :: full_file_data
		logical haveBadScans, thisScanBad, needsAz
		integer status
		integer ordering
		integer scan_number_in_file
		character(256) :: current_filename
		integer file_number
		integer my_nfiles, current_filecount
		character(256), dimension(:), allocatable :: my_files

		!Things needed to process the data and put it in the moduleScan structure.
		real(dp) :: gain(QUIET_NDIODES_MAX)
		real(dp), dimension(:,:), allocatable :: ces_statistics
		real(dp), dimension(:,:,:), allocatable :: diode_statistics
		
		integer ntod,na,n_az
		real(dp), pointer, dimension(:,:)  :: az
		type(ds_moduleScan), pointer :: scan
		type(ds_moduleScanInfo), pointer :: scanInfo
		type(ds_trimap) :: simulation

		!Iteration variables and misc
		integer(i4b) :: unit,one, ierr,dummy
		integer ms,i,j,t,diode,diode2
		real(dp) :: mjd
		
		type(ds_moduleScanList) :: scanList
!		type(ces_index_set) :: targets




		!Set up the size and rank information and save it in the correlator too.
		!This requires that correlator%comm has already been set in the driver.
		unit = ds_get_lun()
		call MPI_Comm_size(correlator%comm, nproc,ierr)
		call MPI_Comm_rank(correlator%comm, rank, ierr)
		correlator%nproc = nproc
		correlator%proc = rank
		call initialize_libquiet_modules(opt)

#if 0
#pragma mark - CHECK PARAMS
#endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!PARAMETER VALIDITY CHECKS
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (opt%should_cut_timestreams) then
			call ds_assert(opt%cutfraction>0 .and. opt%cutfraction<=1,"Invalid cutfraction specified (in parameter file)")
			write(*,*) "Cutting timestreams.  Fraction = ",opt%cutfraction
		endif

		call ds_assert(.not. associated(scans),"Module scans should not already be allocated in readDataAssignWork")
		if (opt%plane_map) then
			npix=opt%plane_size**2
		else
			npix=nside2npix(opt%nside)
		endif
		
		call ds_assert(npix>0,"Invalid nside in readDataAssignWork (in parameter file)")
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		!Work around an occasional titan bug where certain processes do not start until a collective MPI call
		!by doing a simple reduction operation
		!If this fails then you have big problems.
		dummy=1
		call MPI_Allreduce(MPI_IN_PLACE, dummy, 1, MPI_INTEGER,MPI_SUM,correlator%comm,ierr)


#if 0
#pragma mark - DIVIDE SCANS	
#endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!GET AND DIVIDE LIST OF SCANS
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!Turn the target_list data into a form slightly more useful for us.
		!The moduleScanList wraps a linked list of ds_moduleScanInfo types, each of which 
		!contains the module number, run, segment, etc.
		
!		if (rank==0 .and. opt%cut_extension_name/="") call ds_log("Reading scans from extension: "//trim(opt%cut_extension_name), ds_feedback_quiet)
		
		call buildModuleScanList(scanList, opt, output_fits_cards)
!		if (opt%save_name_mapping) call save_name_mapping(scanList,opt%name_mapping_file)

		!This is the total number of moduleScans for all the processors

		nscan = scanlist%length

		if(rank==0) then
			write(message,*) "Total number of (module,scans) over all = ", nScan
			call ds_log(message,ds_feedback_quiet)
		endif

		if (nScan==0) then
		    write(*,*) "NO SCANS!  CANNOT MAKE A MAP WITHOUT SCANS!"
		    call MPI_Finalize(ierr)
		    stop
		endif		


		call allocate_scans_compute_stats(correlator, opt%nside, npix, scanList, ces_statistics, diode_statistics, originalIndices, my_nscan)

#if 0
#pragma mark -  ALLOCATE SPACE
#endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Setup space for this proc
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!Now we have that info we check how many moduleScans this process has and allocate space
		!and initialize the pointers inside
		allocate(scans(0:my_nScan-1))



		if ((.not. do_correlations) .and. correlator%proc==0) write(*,*) "Cross-correlation deactivated"


		!JAZ - not any more - now we do these separately for different null tests outside this function
		!Allocate the maps
!		call prepareTriMap(maps, npix, opt%do_temperature, opt%do_polarization)

		if (opt%simulate) then
			call ds_log("REPLACING TIMESTREAMS WITH SIMULATION FROM "//trim(opt%simulation_filename),ds_feedback_quiet)
			call load_fits_maps(simulation,opt%do_temperature, opt%do_polarization, opt%simulation_filename)
		endif
		
		!Set up the numbering of the scans.
		call setupScans(scans,scanList,rank)

		if(rank==0 .and. opt%data_prior) call ds_log("Data-prior mode active: priors will be computed direct from the data. This may be slower.", ds_feedback_quiet)
		if(opt%subtractMeans .and. (rank==0) ) call ds_log("subtracting means",ds_feedback_noisy)

		!Finally, having built the information about the scans and modules, we can load the data.

		if (rank==0) call ds_log_milestone("DATA_SELECTION")



#if 0
#pragma mark - LOAD LOOP
#endif

		haveBadScans=.false.
		progress=0
		
		my_nfiles = count_my_filenames(scanList, rank)
		call collect_my_filenames(scanList, rank, my_files)
		
		do file_number=1,my_nfiles
			current_filename = my_files(file_number)
			write(message,*) "L3 FILE DOES NOT EXIST: ", trim(current_filename)
			call ds_assert(file_exists(trim(current_filename)), message)
			write(message,'(I0, " Starting loading file: ",A, " (", I0,"/",I0,")" )') ,rank, trim(current_filename), file_number, my_nfiles
			call ds_log(message,ds_feedback_noisy)
			call read_l3_file(trim(current_filename), full_file_data)
			write(message,'(I0, " Loaded file: ",A, " (", I0,"/",I0,")" )') ,rank, trim(current_filename), file_number, my_nfiles
			call ds_log(message,ds_feedback_noisy)
			if (rank==0) then
				write(message,'(A,F7.1,A)') "Approximate load progress: ", (file_number*100.0)/my_nfiles,"%"
				call ds_log(message,ds_feedback_noisy)
				call log_time_memory(trim(current_filename),.true.)
			endif
!			subroutine process_scans(full_file_data, current_filename, scanInfoList, scans, originalIndices, correlator, opt)
			call process_scans(full_file_data,current_filename,scanList,scans,originalIndices, correlator, npix, opt)
			
		call free_lx_struct(full_file_data)	
		enddo
		


#if 0
#pragma mark - FINALIZE loMAP
#endif
		write(message,*) "Rank ", rank," loaded all data."
		call ds_log(message,ds_feedback_debug)

		!If any of the scans were bad and did not load, quit here
		!The user should supply a better accepted list!
		!We may want to reconsider this behaviour.
		call ds_assert(.not. haveBadScans, "Bad scans reported")

		if (rank==0) call ds_log_milestone("DATA_LOADED")




		!JAZ we no longer do this here because we will be generating various combinations of data for different null tests
		!Repixelize both the pointings and the accumulated maps
!		if (opt%plane_map) then
!			full_map_pixels = opt%plane_size**2
!		else
!			full_map_pixels = nside2npix(opt%nside)
!		endif
!		write(*,*) "Starting repixelize", rank
!        if (opt%save_hits) then
!           call repixelizeData(correlator,scans,full_map_pixels,originalIndices,maps,hitMapFilename=opt%hits_filename)           
!        else
!           call repixelizeData(correlator,scans,full_map_pixels,originalIndices,maps)
!        endif

!write(*,*) "Done repixelize", rank

		
! 		if (opt%do_temperature) then
! 			npix = maps%T%npix
! 		else 
! 			npix = maps%Q%npix
! 		endif
		
		
! 		!Sum the accumulated maps accross all the processes, thereby summing for all modules for all scans
! 		if (opt%do_temperature)  call MPI_Allreduce(MPI_IN_PLACE, maps%T%map, maps%T%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
! 		if (opt%do_polarization) call MPI_Allreduce(MPI_IN_PLACE, maps%Q%map, maps%Q%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
! 		if (opt%do_polarization) call MPI_Allreduce(MPI_IN_PLACE, maps%U%map, maps%U%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)

		!We do not generate_cov_subtract_naive here, because our varied null tests mean that we do not do that until later.
		
		call destroyModuleScanList(scanList)



	end subroutine read_data_assign_work




!test

subroutine generate_naive_and_subtract(comm, scans, maps, npix, covariance, covariance_cut, do_T, do_P)
	integer :: comm
	type(ds_modulescan), dimension(0:) :: scans
	type(ds_trimap) :: maps
	type(ds_covariance) :: covariance
	logical :: do_T, do_P
	integer npix, ms, diode
	integer nscan, ierr, p
	real(dp) :: covariance_cut
	type(ds_modulescan), pointer :: scan	
	real(dp) :: covariance_limit

	nscan = size(scans)
!	write(*,*) "NSCAN = ", nscan
	
	call destroyTriMap(maps)
	call prepareTriMap(maps, npix, do_T, do_P)

	!Build the new naive map by summing the naive maps of the separate scans.
	!The separate scans should all have been pre-multiplied by C_w^{-1} so this builds the whole RHS.
	if (do_T) then
		maps%T%map=0
		do ms=0,nscan-1
			do p=1,npix
				maps%T%map(p) = maps%T%map(p) + scans(ms)%maps%T%map(p)
			enddo
		enddo
	endif
	if (do_P) then
		maps%Q%map=0
		maps%U%map=0
		do ms=0,nscan-1
			do p=1,npix
				maps%Q%map(p) = maps%Q%map(p) + scans(ms)%maps%Q%map(p)
				maps%U%map(p) = maps%U%map(p) + scans(ms)%maps%U%map(p)
				if (isnan(scans(ms)%maps%Q%map(p)) .or. isnan(scans(ms)%maps%U%map(p))) write(*,*) "NAN",ms,p
			enddo
		enddo
	endif
	
	!Sum the naive map across processors.
	if (do_T) call MPI_Allreduce(MPI_IN_PLACE, maps%T%map, maps%T%npix, MPI_REAL8, MPI_SUM, comm, ierr)
	if (do_P) call MPI_Allreduce(MPI_IN_PLACE, maps%Q%map, maps%Q%npix, MPI_REAL8, MPI_SUM, comm, ierr)
	if (do_P) call MPI_Allreduce(MPI_IN_PLACE, maps%U%map, maps%U%npix, MPI_REAL8, MPI_SUM, comm, ierr)

	!Get rid of the old covariance.
	call destroy_covariance(covariance)
	
	
	!Prepare and build the new covariance, and invert it.
	call make_naivecov(npix,covariance,comm,scans, do_T, do_P)
	
!	if (rank==0) then
!		write(*,*) "NPIX = ", npix
!		write(*,*) "MAX QQ = ", maxval(covariance%QQ), " at ", maxloc(covariance%QQ)
!		write(*,*) covariance%QQ(1:10)
!		write(*,*) "PP:"
!		write(*,*) scans(0)%inv_Cw_PP
!	endif
	
	
	!Multiply by the covariance
	call cov_mult(maps,covariance)
!	if (rank==0) write(*,*) "Map Q Max = ", maxval(maps%q%map)
		
!	if (rank==0) then
!		do p=1,npix
!			if (isnan(covariance%QQ(p))) write(*,*) "NAN", p
!		enddo
!	endif
	covariance_limit = covariance_cut * best_covariance(covariance)
	if (rank==0) call report_covariance_cut(covariance,covariance_limit)

	!Subtract the newly-made naive map from the offsets
	do ms=0,nscan-1
		call remove_missed_pixels(scans(ms), covariance, covariance_limit)
		
		do diode=1,QUIET_NDIODES_MAX
			if (scans(ms)%flags(diode)==1) call prepareTimestream(scans(ms)%timestreams(diode),scans(ms)%ntod )
		enddo
		
		call map2tod(scans(ms),maps)
		call invCw_mult_tod(scans(ms))


	
		do diode=1,QUIET_NDIODES_MAX
			if (scans(ms)%flags(diode)==1) then
				call subtractTimestreamFromOffset(scans(ms)%timestreams(diode),scans(ms)%offsets_2(diode))
				call destroyTimestream(scans(ms)%timestreams(diode))
				
!				if (rank==0) then
!					write(*,*) "Looking at offsets", ms, diode, scans(ms)%offsets_2(diode)%na
!					write(*,*) scans(ms)%offsets_2(diode)%values
!					do p=1,scans(ms)%offsets_2(diode)%na
!						if (isnan(scans(ms)%offsets_2(diode)%values(p))) write(*,*) "NAN:", ms, diode, p
!					enddo
!				endif
			endif
		enddo

	enddo

end subroutine generate_naive_and_subtract


subroutine report_covariance_cut(covariance,covariance_limit)
	type(ds_covariance) :: covariance
	real(dp) :: covariance_limit
	character(256) :: message
	integer ncut, nnan, p
	real(dp) :: c
	
	nnan=0
	ncut=0
	
	if (covariance%has_p) then
		do p=1,covariance%npix
			c = covariance%QQ(p)*covariance%UU(p) - covariance%QU(p)**2
			if (.not. (c .lt. covariance_limit)) ncut=ncut+1
			if (isnan(c)) nnan=nnan+1
		enddo
	endif
	write(message,'("Cut ",I0," pixels, of which ",I0, " were NaN")') ncut,nnan
	call ds_log(message,ds_feedback_noisy)

end subroutine

function best_covariance(cov)
	!We want the smallest covariance.
	type(ds_covariance) :: cov
	real(dp) :: best_covariance, c
	integer p
	
	call ds_assert(cov%has_p .and. (.not. cov%has_t), "Have not coded best_covariance for T maps")
	best_covariance = 1.0e30
	do p=1,cov%npix
		c = cov%QQ(p)*cov%UU(p)-cov%QU(p)**2
		if (c .lt. best_covariance) best_covariance = c
	enddo
	
end function best_covariance

subroutine remove_missed_pixels(scan, covariance, covariance_limit)
	type(ds_modulescan) :: scan
	type(ds_covariance) :: covariance
	real(dp) :: covariance_limit
	integer :: t, p
	real(dp) :: c
	
	
	if (covariance%has_t) then
		do t=1,scan%ntod
			p = scan%pointing(t)
			if (p == bad_pixel) cycle
			c = covariance%TT(p)
			if (isnan(c)) scan%pointing(t)=bad_pixel
		enddo
	endif
	
	if (covariance%has_p) then
		do t=1,scan%ntod
			p = scan%pointing(t)
			if (p == bad_pixel) cycle
			c = covariance%QQ(p)*covariance%UU(p) - covariance%QU(p)**2
			if (.not. (c .lt. covariance_limit)) scan%pointing(t)=bad_pixel
		enddo
	endif
	
end subroutine remove_missed_pixels

subroutine allocate_scans_compute_stats(correlator, nside, npix, scanList, ces_statistics, diode_statistics, originalIndices, my_nscan)

	type(ds_correlator) :: correlator
	integer nscan, nproc
	integer rank
	integer(i8b), allocatable, dimension(:) :: originalIndices
	integer npix
	integer my_nscan
	type(ds_modulescanList) :: scanList
	integer, allocatable, dimension(:) :: nFile_proc, file_owner
	real(dp), allocatable, dimension(:,:) :: ces_statistics
	real(dp), allocatable, dimension(:,:,:) :: diode_statistics
	character(256), allocatable, dimension(:) :: fileList
	character(256) :: current_filename
	integer nextra 
	integer nside
	integer nmin, i, j, p, current_filenumber
	type(ds_modulescaninfo), pointer :: scanInfo
	
	
	integer nFile
	nproc = correlator%nproc
	nFile = count_all_filenames(scanList)
	call collect_all_filenames(scanList, fileList)

	allocate(nFile_proc(0:nproc-1))
	nmin=floor((1.0*nFile)/nproc)
	nextra = nFile - nmin*nproc
	nFile_proc = nmin
	call ds_assert(nextra>=0,"Broken nextra - this is Joe's fault.")
	if (nextra>0) nFile_proc(0:nextra-1) = nmin + 1
	!Determine the number of moduleScans each process is responsible for.
	!Each process knows this number for each process.
	!nmin is the minimum number of moduleScans held by a process.
	!Some have more to make up the difference.

	!Some useful checks and feedback about the allocation of moduleScans
	call ds_assert(sum(nFile_proc) == nFile, "Mistake in assignment of dets - not all assigned.")

	allocate(file_owner(nFile))
	p=1
	do i=0,nproc-1
		do j=1,nFile_proc(i)
			file_owner(p)=i
			p=p+1
		enddo
	enddo

	
	scanInfo=>scanList%first
	current_filenumber=0
	current_filename=""
	my_nscan = 0
	do
		if (current_filename .ne. scanInfo%filename) then
			current_filenumber = current_filenumber + 1
			current_filename=fileList(current_filenumber)
		endif
		scanInfo%owner = file_owner(current_filenumber)
		if (scanInfo%owner==correlator%proc) my_nscan = my_nscan+1
		if (.not. associated(scanInfo%next)) exit
		scanInfo=>scanInfo%next
	enddo
	if (correlator%proc==0) call ds_log("Computing hit map and statistic list",ds_feedback_quiet)
	call compute_hit_map_and_stats(fileList, file_owner, ces_statistics, diode_statistics, nside, originalIndices, npix, correlator%comm, correlator%proc)  !This computes the number of pixels npix in the map.
	
	
	deallocate(nFile_proc)
	deallocate(fileList)

end subroutine allocate_scans_compute_stats


subroutine process_scans(full_file_data, current_filename, scanInfoList, scans, originalIndices, correlator, npix, opt)
	type(ds_modulescaninfo) :: scanInfo
	type(ds_correlator) :: correlator
	type(ds_modulescanlist) :: scanInfoList
	type(ds_modulescan), dimension(0:) :: scans
	type(Lx_struct) :: full_file_data
	character(256) :: current_filename
	type(ds_oslo_options) :: opt
	integer(i8b), allocatable, dimension(:) :: originalIndices
	integer npix
	

	scanInfo = scanInfoList%first
	do
		if(scanInfo%owner==rank .and. scanInfo%filename==current_filename) then
			call process_scan(scanInfo,scans(scanInfo%modscan),full_file_data, originalIndices, correlator, opt, npix)
		endif
		if(.not. associated(scanInfo%next)) exit
		scanInfo=scanInfo%next
	enddo
end subroutine process_scans


subroutine process_scan(scanInfo, scan, full_file_data, originalIndices, correlator, opt, npix)
	type(ds_modulescan) :: scan
	integer(i8b), allocatable, dimension(:) :: originalIndices
	type(ds_modulescaninfo) :: scanInfo
	type(ds_oslo_options) :: opt
	type(ds_correlator) :: correlator
	type(Lx_struct) :: full_file_data
	integer ntod, na, i
	integer npix
	
	call setup_moduleScan_ndiodes(scan,scanInfo%ndiodes_t,scanInfo%ndiodes_p)
	ntod = size(full_file_data%tod,1)
	if (opt%should_cut_timestreams) ntod = ntod * opt%cutfraction
	na = ntod/opt%offsetLength
	ntod = na * opt%offsetLength
	scan%ntod=ntod
	scan%owns_az=.false.
!	subroutine get_scan_from_fits(info,full_data,gain,moduleScan,originalIndices, opt,needsAz, do_T, do_P)

	call get_scan_from_fits(scanInfo,full_file_data,scan,originalIndices, opt, .false., opt%do_temperature, opt%do_polarization)
!	if (opt%simulate) call simulate_moduleScan(scan,simulation, opt)
	
	
	scan%inv_Cw=0.0

	!This is a good place to get the priors if we want to do it directly from the data.
	if (opt%data_prior) call prepare_one_data_prior(correlator,opt%offsetLength,scan)

	!build inverse covariance matrix from variances and correlations
	!we already filled in the sigmas above
	call make_inv_Cw(scan,opt%do_correlations)
	!apply the inverse covariance matrix in moduleScan to the timestreams in moduleScan
	call invCw_mult_tod(scan)

	call prepareTriMap(scan%maps, npix, opt%do_temperature, opt%do_polarization)
	!add module scan to accumulated healpix map
	call add2rhs(scan,scan%maps)
	!deproject moduleScan%timestream -> moduleScan offsets
	!When we have done this we no longer need the timestream so we deallocate it		

!	write(*,*) "Q SCAN MAP MAX", maxval(scan%maps%Q%map)

	do i=1,QUIET_NDIODES_MAX
		if (scan%flags(i)==1) then
			call deprojectTimestreamOntoOffset(scan%timestreams(i),scan%offsets_2(i))
			call destroyTimestream(scan%timestreams(i))
		endif
	enddo

end subroutine process_scan
#if 0
#pragma mark -
#pragma mark LOADING FITS DATA
#pragma mark -
#endif


subroutine load_fits_maps(sim, do_T, do_P, filename)
	use fitstools
	type(ds_trimap) :: sim
	logical :: do_T, do_P
	character(*) :: filename
	real(dp), dimension(:,:), allocatable :: map_data
	integer, parameter :: nmaps = 3
	integer nside, npix
	real(dp) :: null_value
	logical any_null
	integer unit, status, hdutype
	character(80) :: comment

	call FTGIOU(unit, status)					
	call FTNOPN(unit,trim(filename),FITS_READ_ONLY,status)
	call FTMAHD(unit,2,hdutype,status)
	call FTGKYK(unit,"NSIDE",nside,comment,status)
	call FTFIOU(unit, status)			
	
	call ds_assert(status==0, "Could not read NSIDE parameter from simulated map: "//trim(filename) )
	npix = nside2npix(nside)
	allocate(map_data(0:npix-1,nmaps))
	call read_bintab(filename, map_data, npix, nmaps, null_value, any_null)
	call prepareTriMap(sim,npix, do_T, do_P, zero_based = .true.)  !Healpix maps are zero-based.
	
	if (do_T) sim%T%map = map_data(:,1)
	if (do_P) sim%Q%map = map_data(:,2)
	if (do_P) sim%U%map = map_data(:,3)
		
	deallocate(map_data)

	
end subroutine load_fits_maps



	function find_hdu_with_name(unit,name) result(n)
		integer :: unit
		character(*) :: name
		integer n
		integer :: original_hdu
		integer number_of_hdu,hdutype, status
		character(80) :: extname, comment
		integer hdu
		n=-100
		status=0
		!Record the original hdu number so we can return to it.
		call FTGHDT(unit,original_hdu,status)
		call FTTHDU(unit, number_of_hdu, status)
		do hdu = 1,number_of_hdu
			call FTMAHD(unit,hdu,hdutype,status)
			call FTGKYS(unit,"EXTNAME",extname,comment,status)
			if (status/=0) then
				status=0
			else
				if (trim(extname)==trim(name)) then
					n=hdu
					exit
				endif
			endif
		enddo
		call FTMAHD(unit,original_hdu,hdutype,status)
		
	end function find_hdu_with_name




subroutine get_scan_from_fits(info,full_data,moduleScan,originalIndices, opt,needsAz, do_T, do_P)
  use quiet_assembly_mod
	use quiet_gain_mod

	type(ds_oslo_options) :: opt
	real(dp),dimension(QUIET_NDIODES_MAX) :: gain
	type(ds_moduleScanInfo) :: info
	integer(i8b), allocatable, dimension(:) :: originalIndices
	logical :: do_T, do_P
	type(ds_moduleScan) :: moduleScan
	type(Lx_struct) :: full_data
	real(dp) :: mu
	logical :: needsAz
	real(dp), dimension(2)::linearFit
	integer, parameter :: psi_field = 1
	integer, parameter :: az_field = 2
	logical, save :: firstTime=.true.
	real(dp) :: delta_psi
	integer start,finish
	integer pixel
	integer k
	real(dp) :: theta,phi
	integer ntod,na,i,t,diode
	real(dp) :: plane_x,plane_y, plane_half_size, plane_center_ra_radians, plane_center_dec_radians
	integer :: status, pixel_x, pixel_y
	type(quiet_assembly) :: assembly
	integer :: coord_in, coord_out
	real(dp), dimension(:), allocatable :: gains
	integer :: overall_diode_index
	real(dp) :: mean_gain
!	type(common_info) :: shared_info
!	logical :: do_T, do_P
	integer :: module_number, ndiode, healpix_pixel
	
	do_T = opt%do_temperature
	do_P = opt%do_polarization
!	coord_in = parse_coord_name(coord_name_in)
!	coord_out = parse_coord_name(coord_name_out)
!	call ds_generate_common_info(info,shared_info)
!	call reorganize_data_assembly(do_T, opt%coord_in, opt%coord_out, opt%nside, full_data, assembly, shared_info)
	
	
	module_number = info%id3
	ndiode = size(full_data%sigma0)
	
	!Set up the ID numbers in the moduleScans so we can keep track of them later
	moduleScan%id(1) = info%id1
	moduleScan%id(2) = info%id2
	moduleScan%id(3) = info%id3
	
	moduleScan%has_T = do_T
	moduleScan%has_P = do_P
	
	!Determine the size of this timestream
	
	! start =  info%first_index
	! finish = info%last_index
	! ntod = finish - start + 1

	ntod = size(full_data%tod,1)
	moduleScan%flags = info%diode_flag

	if (opt%should_cut_timestreams) then
		ntod = ntod * opt%cutfraction
	endif
	
	
	call ds_assert(ntod>0,"Zero-sized (or negative) TOD passed to get_scan_from_fits")
	
	na = ntod/opt%offsetLength
	ntod = na * opt%offsetLength
!	finish = start+ntod-1
	
	!Set up the offsets
	do i=1,QUIET_NDIODES_MAX
		if (moduleScan%flags(i)==1) then
			call prepareOffsets(moduleScan%offsets_2(i),na,opt%offsetLength,azimuth_flag=.false.)
		endif
	enddo
	
	
	!Set up the pixel pointing
	allocate(moduleScan%pointing(ntod))
	
	do t=1,ntod
		theta = full_data%point(2,t,module_number)
		phi = full_data%point(1,t,module_number)
		call ang2pix_nest(opt%nside,theta,phi,healpix_pixel)  !NOW USING NEST BY DEFAULT
		moduleScan%pointing(t) = originalIndices(healpix_pixel)
#warning now doing nested pixels by default - fix elsewhere
		!assembly%pix(1,t,1)  !Indices: (main/partner),sample,diode  !NB We ignore partner and assume all diodes point at the same thing for now.
	enddo

	
	moduleScan%inv_Cw=0.0

	if (do_P) then
		allocate(moduleScan%theta(ntod))
		do t=1,ntod
			moduleScan%theta(t) = full_data%point(3,t,module_number)  !assembly%point(3,t,1,1)  !Indices : (phi/theta/psi),sample,diode,(main/partner)
		enddo
!		if (info%owner==0) write(*,*) module_number, "MIN,MAX THETA = ",minval(moduleScan%theta), maxval(moduleScan%theta)
		call read_dpsi(modulescan)
!		if (info%owner==0) write(*,*) module_number, "DPSI = ", modulescan%dpsi
	endif

	do diode=1,QUIET_NDIODES_MAX
		overall_diode_index = diode + QUIET_NDIODES_MAX*(module_number-1)  !From 1 ..
		if ((.not. do_T) .and. (diode .le. info%ndiodes_T)) cycle
		if ((.not. do_P) .and. (diode .gt. info%ndiodes_T)) cycle
		if (moduleScan%flags(diode)==0) cycle
		call prepareTimestream(moduleScan%timestreams(diode),ntod)
		
		if (opt%simulate) cycle

		if (moduleScan%ndiodes_P .gt. 0)  then !ie This is a polarization module for QUIET
			do t=1,ntod
				moduleScan%timestreams(diode)%timestream(t) = full_data%tod(t,overall_diode_index)
			enddo
		else  !This is a temperature module so we must apply the partner splitting
			call ds_assert(.false.,"Cannot do temperature map-making yet since there is a differencing scheme")
			do t=1,ntod
				moduleScan%timestreams(diode)%timestream(t) = full_data%tod(t,overall_diode_index)
			enddo
			
		endif

!#warning REMOVED GAIN DIVISION
		call apply_gain(moduleScan%timestreams(diode),full_data%gain(:,overall_diode_index), full_data%time(1:ntod), full_data%time_gain)

		mean_gain = sum(full_data%gain(:,overall_diode_index))/size(full_data%gain(:,overall_diode_index))
!		allocate(gains(ntod))
!		call get_gains(full_data%time(1:ntod), overall_diode_index, gains)
!		moduleScan%timestreams(diode)%timestream = moduleScan%timestreams(diode)%timestream/gains
!		deallocate(gains)

!			moduleScan%timestreams(diode)%timestream(t) = moduleScan%timestreams(diode)%timestream(t) / full_data%gain(overall_diode_index,t)

		!Subtract a linear trend from the data
		if (opt%subtract_linear) then
			linearFit=least_square_fit(moduleScan%timestreams(diode)%timestream,ntod)

			call ds_log(message, ds_feedback_debug)
			do t=1,ntod
				moduleScan%timestreams(diode)%timestream(t)=moduleScan%timestreams(diode)%timestream(t)-linearFit(1)*t-linearFit(2)
			enddo
		endif

		!Subtract the mean of the data
		if(opt%subtractMeans) then			
			mu = sum(moduleScan%timestreams(diode)%timestream)/ntod
			if (.not. mu==mu) write(*,*) "BAD DIODE:", module_number, diode
			do t=1,ntod
				moduleScan%timestreams(diode)%timestream(t) = &
				moduleScan%timestreams(diode)%timestream(t) - mu						
			enddo
		endif


		
!#warning sigma set to 1		
		moduleScan%noise%sigma(diode) = full_data%sigma0(overall_diode_index) / mean_gain
		moduleScan%noise%alpha(diode) = full_data%alpha(overall_diode_index)
		moduleScan%noise%fknee(diode) = full_data%fknee(overall_diode_index)
		
		
!		write(*,*) "STDEV, SIGMA0 = ", stdev(moduleScan%timestreams(diode)%timestream), moduleScan%noise%sigma(diode)
		
		do t=1,ntod
			if (isnan(moduleScan%timestreams(diode)%timestream(t))) then
				moduleScan%flags(diode)=0
				write(message,'(I0,": Removed diode (NaN): ",A," ",I0," ",I0)') rank, trim(info%filename), module_number, diode
				call ds_log(message,ds_feedback_quiet)
				exit
			endif
		enddo
		
		
	enddo

	
end subroutine get_scan_from_fits


function stdev(x)
	real(dp), dimension(:) :: x
	real(sp) :: stdev
	integer i
	integer n
	real(dp) xsum, x2sum
	
	xsum = 0.0
	x2sum = 0.0
	
	do i=lbound(x,1), ubound(x,1)
		xsum = xsum + x(i)
		x2sum = x2sum + x(i)**2
	enddo
	
	n = size(x)
	stdev = sqrt(x2sum/n + (xsum/n)**2)
end function stdev

subroutine read_applied_cuts(unit, cards)
	use fl_lists
	integer unit
	type(fl_string80_list) :: cards
	character(80) :: card
	character(8) :: keyword
	integer i, status
	i=0
	call FTPMRK
	keyword=""
	do
		write(keyword,'("DSCUT",I0.3)') i
		call FTGCRD(unit,keyword,card,status)
		if (status .ne. 0) exit
		call fl_append_string80_list(cards,card)
		i=i+1
	enddo
	call FTCMRK
end subroutine read_applied_cuts



subroutine apply_gain(timestream,sparse_gain, dense_time, sparse_time)
	use math_tools
	type(ds_timestream) :: timestream
	real(sp), dimension(:) :: sparse_gain
	real(dp), dimension(timestream%nt) :: dense_time
	real(dp), dimension(:) :: sparse_time
	real(dp), dimension(:), allocatable :: sparse_gain_dp, dense_gain
	integer t
	call ds_assert(size(sparse_time)==size(sparse_gain),"UNEXPECTED SHAPE OF GAIN MODULES")
	
	allocate(sparse_gain_dp(size(sparse_gain)))
	allocate(dense_gain(timestream%nt))
	sparse_gain_dp = sparse_gain
	call lin_interpol(sparse_time, sparse_gain_dp, dense_time, dense_gain)
	do t=1,timestream%nt
		timestream%timestream(t) = timestream%timestream(t)/dense_gain(t)
	enddo
	deallocate(dense_gain)
	deallocate(sparse_gain_dp)

end subroutine apply_gain

subroutine read_dpsi(scan)
	use quiet_module_mod
	type(ds_moduleScan) :: scan
	integer diode
	integer modnum
	
	modnum=scan%id(3)
	
	do diode=1,scan%ndiodes_p
		scan%dpsi(diode) = get_diode_angle(modnum-1, diode-1)
	enddo

end	subroutine read_dpsi


#if 0
#pragma mark -
#pragma mark LISTS OF SCANS TO PROCESS
#pragma mark -
#endif
	subroutine buildModuleScanList(modScanList, opt, cards)
	! The job of this function is to build up a list of scans to be read.
	! These are not yet read at this stage but will be later.
	! The modScanList object is a linked-list of modScanInfo objects.
	
	
		use fits_helper
		use fl_lists
		use quiet_ces_mod
		use quiet_assembly_mod
		use quiet_lx_mod
		use quiet_module_mod
		use quiet_utils
		use quiet_fileutils
		use quiet_task_mod
		use quiet_acceptlist_mod
		use quiet_filter_mod
		use quiet_target_mod

!!!
		! Read the headers of all the input files and find out how many scans are in each and how to find them.
		!!!

		type(ds_moduleScanList) :: modScanList
		type(ds_oslo_options) :: opt
		type(fl_string80_list), optional :: cards

!		integer, dimension(QUIET_NDIODES_MAX) :: diode_flag
!		character(256), dimension(:) :: filename_list
		character(256) :: filename
		integer id1,id2, id3
		integer unit,hdutype
		integer nelements
		character(80) :: comment
		integer f, status, diode
		integer :: ndiodes_t, ndiodes_p
		character(8) :: keyword
		character(8) :: value
		logical :: do_T, do_P
		integer nscan,s
		integer i
!		type(quiet_ces) :: ces
		type(quiet_target)   :: target
		type(acceptlist)     :: alist
		character(len=512)   :: target_name, acceptfile
		integer(i4b),          dimension(:),     allocatable :: cnums, map2mask
		real(dp),              dimension(:,:),   allocatable :: stats
	  type(quiet_ces_info) :: ces_info
		integer ces_index, ces, cnum
		integer accepted_diodes_this_ces


		do_T = opt%do_temperature

		do_P = opt%do_polarization
		
		target_name="patch_gc"
		
!QUIET STUFF		
!  call initialize_modules
		call get_parameter(0, opt%quiet_parameter_file, 'ACCEPTLIST',            par_string=acceptfile)

		call initialize_accept_list(acceptfile, alist)

		call init_target(target, alist)

		call filter_object(target, target_name)
		if (rank==0) call ds_log(message, ds_feedback_quiet)
		!  call collect_ces_info(target, stats, map2mask)
		call get_accepted_ceses(target%alist, cnums)
		write(message,'(A,I0,A)') "There are", get_num_ces(), " accepted CES."
		call init_moduleScanList(modScanList)


		do ces=1,size(cnums)
			cnum = cnums(ces)
			ces_index = lookup_ces(cnums(ces))
!			cnums(ces) is the index into alist%status
!			not sure what ces_index means.  Index into something else?

!			call ds_assert(is_accepted(target%alist,cnums(ces)), "Something has gone wrong in buildModuleScanList - accept failure.")
			call get_ces_info(cnum, ces_info)
			!For now, just read one CES file, for testing.
			id1 = cnum
			id2 = 0
			ndiodes_t = QUIET_NDIODES_T
			ndiodes_p = QUIET_NDIODES_P
			filename=trim(ces_info%l3file)
			accepted_diodes_this_ces=0
			do i=1,QUIET_N_MODULES
				if (.not. any(target%alist%status(:,i-1,cnum)==REJECTED_NONE)) cycle
				if ((.not. do_T) .and. is_temperature_module(i-1)) cycle
				if ((.not. do_P) .and. is_polarization_module(i-1)) cycle

				
				if (is_polarization_module(i-1)) then
					ndiodes_p = 4
					ndiodes_t = 0
				endif
				if (is_temperature_module(i-1)) then
					ndiodes_t = 4
					ndiodes_p = 0
				endif

				id3 = i

				call appendModuleScanInfo(modScanList,id1,id2,id3,ndiodes_t, ndiodes_p,.false.,filename)
				modScanList%last%diode_flag=0
				do diode=1,QUIET_NDIODES_MAX
					if (target%alist%status(diode-1,i-1,cnum)==REJECTED_NONE) then
						modScanList%last%diode_flag(diode) = 1
						accepted_diodes_this_ces = accepted_diodes_this_ces+1
					endif
				enddo
				modScanList%last%ndiodes = QUIET_NDIODES_MAX
			end do
			
		enddo

	end subroutine buildModuleScanList



subroutine init_moduleScanList(L)
	type(ds_moduleScanList) :: L
	nullify(L%first)
	nullify(L%last)
	L%length=0
end subroutine init_moduleScanList

function makeModuleScanInfo(id1,id2,id3,ndiodes_T, ndiodes_P,needs_azfilter,filename) result(output)
	type(ds_moduleScanInfo), pointer :: output
	integer id1,id2,id3,ndiodes_T, ndiodes_P
	logical needs_azfilter
	character(len=256) filename

	allocate(output)

	output%id1=id1
	output%id2=id2
	output%id3=id3
	output%ndiodes_T = ndiodes_T
	output%ndiodes_P = ndiodes_P
	output%filename=filename
	
	nullify(output%next)
	
end function makeModuleScanInfo

subroutine appendModuleScanInfo(list,id1,id2,id3,ndiodes_t, ndiodes_p,needs_azfilter,filename)
	type(ds_moduleScanList) :: list
	integer id1,id2,id3,ndiodes_t, ndiodes_p
	logical needs_azfilter
	character(len=256) filename
	type(ds_moduleScanInfo), pointer :: info
	info => makeModuleScanInfo(id1,id2,id3,ndiodes_t,ndiodes_p,needs_azfilter,filename)
	if (list%length==0) then
		list%first=>info
		list%last=>info
		list%length=1
	else
		list%last%next=>info
		list%last=>info
		list%length=list%length+1
	endif

end subroutine appendModuleScanInfo

subroutine destroyModuleScanList(list)
	type(ds_moduleScanList) :: list
	type(ds_moduleScanInfo), pointer :: info,next


	info => list%first
	if (.not. associated(info)) return
	do
		next=>info%next
		deallocate(info)
		if (.not. associated(next)) exit
		info=>next
	enddo
	nullify(list%first)
	nullify(list%last)
	list%length=0

end subroutine destroyModuleScanList

#if 0
#pragma mark -
#pragma mark PREPARING DATA
#pragma mark -
#endif

	subroutine repixelizeData(correlator,moduleScans,maxIndex,originalIndices,maps,hitMapFilename)
		!This subroutine resets the pointing of a ds_pointing instance so that the 
		!pixel indices it uses run from 1..npix (and can therefore be used as array indices)
		!The maxIndex argument specifies the number of pixels in the input class,
		!for example if you are converting a set of healpix indices then maxIndex should be 12*Nside^2
		type(ds_correlator) :: correlator
		type(ds_modulescan), dimension(0:correlator%my_nmodules-1) :: moduleScans
!		type(ds_covariance) :: covariance
		integer(i8b) :: maxIndex
		integer(i8b), allocatable, dimension(:) :: originalIndices
		type(ds_trimap) :: maps
		character(len=*), optional :: hitMapFilename
		integer (i4b), allocatable, dimension(:) :: hitCount
		integer i,d,n,p,q,t,nbad
		integer ierror

		!Each process builds its own hit count map in healpix format
		allocate(hitCount(0:maxIndex-1))
		
		hitCount=0
		do i=0,correlator%my_nmodules-1
			do d=1,moduleScans(i)%ndiodes
				if (moduleScans(i)%flags(d)==1) then
					do t=1,moduleScans(i)%ntod
						p=moduleScans(i)%pointing(t)
						if (p>=0) hitCount(p)=hitCount(p)+1
					enddo
				endif
			enddo
		enddo
		
		write(message,*) "Computed hit counts: ", rank
		if (rank==0) call ds_log(message,ds_feedback_debug)
		
		!The processes all sum their hit count map to get the global hit count map.
		ierror=0
		write(message,*) "Reached barrier: ", rank
		if (rank==0) call ds_log(message,ds_feedback_debug)
		
		call MPI_Barrier(correlator%comm,ierror)
		write(message,*) "Crossed barrier: ", rank
		if (rank==0) call ds_log(message,ds_feedback_debug)

		call MPI_AllReduce(MPI_IN_PLACE,hitCount,maxIndex,MPI_INTEGER,MPI_SUM,correlator%comm,ierror)
		call ds_assert(ierror==0,"Error in MPI_AllReduce in repixeizeData")

		write(message,*) "Shared hit counts among processors"
		if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)


		!If requested, the hit count map is saved to file
		if (present(hitMapFilename) .and. correlator%proc==0) then
			!Save hit map
           open(file=hitMapFilename,unit=25)
           do p=0,maxIndex-1
              write(25,*) p,hitCount(p)
           enddo
			write(message,*) "Saved hit counts to file"
			if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
		endif



		!Count the number of hit pixels to get the total npix
		n=0
		do p=0,maxIndex-1
			if (hitCount(p)>0) 	n=n+1
		enddo

		!added by DWPS - code to handle degenerate pixels
		!if the pixels have only a single hit, they are bad. So set pointing to bad_pixel flag for those
		!pixels (ds_mapping knows how to handle them).  Also remove from the map, by removing number of 
		!bad pixels from npix and setting the hitcount for those pixels to 0.
		nbad=0
		do p=0,maxIndex-1	
			if(hitCount(p)==1) nbad = nbad+1
		enddo

		do i=0,correlator%my_nmodules-1
			do d=1,moduleScans(i)%ndiodes
				if (moduleScans(i)%flags(d)==1) then
					do t=1,moduleScans(i)%ntod
						if(hitCount(moduleScans(i)%pointing(t))==1) then
							moduleScans(i)%pointing(t) = bad_pixel
						endif
					enddo
				endif
			enddo
		enddo
		
		write(message,*) "Handled degenerate pixels"
		if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
		

		!remove bad pixels from the map
		forall(p=0:maxIndex-1 , hitCount(p)==1) hitCount(p)=0

		!remove bad pixels from npix	
		n=n-nbad

		!end of added by DWPS

		write(message,*) "Built global hit map: npix = ",n
		if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)


		!The root process should preserve the mapping that takes the new indices back to the old ones.
		!The originalIndices array stores that data.
		if (correlator%proc==0) then
			allocate(originalIndices(n))
			q=1
			do p=0,maxIndex-1
				if (hitCount(p)>0) then
					originalIndices(q)=p
					q=q+1
				endif
			enddo
			
			write(message,*) "Computed new indices"
			if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
			

			!		open(unit=14,file="pixels.dat")
			!		do p=1,n
			!			write(14,*) originalIndices(p)
			!		enddo
			!		close(14)
			!		open(unit=14,file="hitcount.dat")
			!		do p=0,maxIndex-1
			!			write(14,*) hitCount(p)
			!		enddo
			!		close(14)

		endif

		!Relabel the hit count map so that it contains the new index for each healpix pixel
		!instead of the hit count.  This will be useful for re-pixelizing the pointings and maps
		q=1
		do p=0,maxIndex-1
			if (hitCount(p)>0) then
				hitCount(p)=q
				q=q+1
			else
				hitCount(p)=-1
			endif
		enddo
		
		write(message,*) "Relabelled hit count"
		if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
		

		!Re-pixelize the pointings in the moduleScans
		do i=0,correlator%my_nmodules-1
			do t=1,moduleScans(i)%ntod
				p=moduleScans(i)%pointing(t)
				if (p>=0) then
					moduleScans(i)%pointing(t) = hitCount(p)
				else 
					moduleScans(i)%pointing(t) = bad_pixel
				endif
			enddo
		enddo

	write(message,*) "Repixelized pointings"
	if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
	

		!Repixelize the maps, reducing them in size from the full healpix sized ones to
		!ones containing only the pixels hit somewhere in one of the timestreams
		if (maps%has_T) call repixelizeMap(maps%T,hitCount,n)
		if (maps%has_P) call repixelizeMap(maps%Q,hitCount,n)
		if (maps%has_P) call repixelizeMap(maps%U,hitCount,n)
		
		write(message,*) "Repixelized maps"
		if (correlator%proc==0) call ds_log(message,ds_feedback_quiet)
		
		
		deallocate(hitCount)

	end subroutine repixelizeData


	subroutine repixelizeMap(map,new_pix,npix)
		type(ds_map) :: map
		integer, dimension(0:map%npix-1) :: new_pix
		integer npix
		integer p,i
		integer, dimension(:), allocatable :: new_indices
		real(dp), dimension(:), allocatable :: new_map
		
		allocate(new_map(npix))
		allocate(new_indices(npix))
		do i=0,map%npix-1
			p = new_pix(i)
			if (p .ge. 1) then
				new_map(p) = map%map(i)
				new_indices(p) = i
			endif
		enddo
		call destroyMap(map)
		call prepareMap(map,npix)
		map%map = new_map
		map%indices = new_indices
		deallocate(new_map)
		deallocate(new_indices)
	end subroutine
	
	

	subroutine setupScans(scans,scanList,rank)
		!Perform misc setup
		type(ds_moduleScanList) :: scanList
		type(ds_moduleScan), pointer, dimension(:) :: scans
		integer :: rank

		type(ds_moduleScan), pointer :: scan
		type(ds_moduleScanInfo), pointer :: scanInfo
		integer ms,m,run,seg,diode

		scanInfo=>scanList%first
		ms=0
		do 	
			if (scanInfo%owner .ne. rank) then
				if (.not. associated(scanInfo%next)) exit
				scanInfo=>scanInfo%next
				cycle		
			endif
			scan=>scans(ms)
			scanInfo%modscan=ms
			allocate(scan%flags(scaninfo%ndiodes_t+scaninfo%ndiodes_p))

			scan%flags = scanInfo%diode_flag
			if (.not. associated(scanInfo%next)) exit
			scanInfo=>scanInfo%next
			ms=ms+1				
		enddo

	end subroutine setupScans




#if 0
#pragma mark -
#pragma mark SIMULATION FUNCTIONS
#pragma mark -
#endif

subroutine simulate_moduleScan(moduleScan,sim,opt)


 	type(ds_moduleScan) :: moduleScan
 	type(ds_oslo_options) :: opt
 	type(ds_trimap) :: sim
! !	real(dp), dimension(:) :: az
! !	integer rank
! 	real(dp), dimension(QUIET_NDIODES_MAX) :: gain
! 	integer diode
! 	integer nbin
! 	integer unit
! 	integer t
! 	logical is_temperature
! 	
! 	logical, save :: first_sim = .false.
stop 'Have not coded simulate_moduleScan - use library.'
! 	
! 	do diode=1,oslo_NDIODES_MAX
! 		if (moduleScan%flags(diode) .ne. 1) cycle
! 		
! 		do t=1,moduleScan%timestreams(diode)%nt
! 			moduleScan%timestreams(diode)%timestream(t) = 0.0
! 		enddo
! 	
! 		call ds_simulate_white(moduleScan%timestreams(diode), scan%noise%sigma(diode) )
! 		call ds_simulate_signals(moduleScan,sim, gain, opt%do_temperature, opt%do_polarization)
! 
! 	enddo

end subroutine simulate_moduleScan



subroutine ds_simulate_white(timestream,sigma)
	use random
	type(ds_timestream) :: timestream
	real(dp), intent(in) :: sigma
	integer :: t
	
	do t=1,timestream%nt
		timestream%timestream(t) = timestream%timestream(t) + random_normal() * sigma !I tried to put the sigma afterwards, multiplying the whole lot, but it crashed the compiler
	enddo
end subroutine ds_simulate_white


subroutine ds_simulate_signals(scan,sim, gain,do_T, do_P)
 	type(ds_moduleScan) :: scan
 	type(ds_trimap) :: sim
 	real(dp), dimension(QUIET_NDIODES_MAX) :: gain
 	logical do_T, do_P
! 	integer diode, t, p
! 	integer npix_sim, nside_sim
! 	
! 	call ds_assert(do_T .or. do_P, "Must simulate either temperature or polarization")
! 	
! 	if (do_T) call ds_assert(sim%has_T, "No temperature simulation map provided in ds_simulate_signal")
! 	if (do_P) call ds_assert(sim%has_P, "No polarization simulation map provided in ds_simulate_signal")
! 
! 	if (do_T) then
! 		npix_sim = sim%T%npix
! 	else
! 		npix_sim = sim%Q%npix
! 	endif
! 	
! 	nside_sim = npix2nside(npix_sim)
! 	
! 	
! ! TODO: Add support for simulating from different resoluation maps	
! 	
! 	if (do_T) then
! 		do diode=1,oslo_NDIODES_T
! 			if (scan%flags(diode) .ne. 1) cycle
! 			do t=1,scan%timestreams(diode)%nt
! 				p=scan%pointing(t)
! 				scan%timestreams(diode)%timestream(t) = scan%timestreams(diode)%timestream(t) + sim%T%map(p) * gain(diode)
! 			enddo
! 		enddo
! 	endif
! 
! 	if (do_P) then
! 		do diode=oslo_NDIODES_T+1,oslo_NDIODES_MAX
! 			if (scan%flags(diode) .ne. 1) cycle
! 			do t=1,scan%timestreams(diode)%nt
! 				p=scan%pointing(t)
! 				stop 'Have not coded simulated polarizations.  Not too hard though'
! 				scan%timestreams(diode)%timestream(t) = scan%timestreams(diode)%timestream(t) + sim%T%map(p) * cos(2*scan%theta(t)) * gain(diode)
! 			enddo
! 		enddo
! 	endif
! 	
	
end subroutine ds_simulate_signals

#if 0
#pragma mark -
#pragma mark MISC UTILS
#pragma mark -
#endif
subroutine testValueForEvil(value,variable_name,zero_is_bad)
	real(dp) :: value
	character(*) :: variable_name
	character(256) ::  message
	logical, optional :: zero_is_bad
	
#ifdef USE_IEEE
	write(message,*) variable_name, " is NaN :", value
	call ds_assert(.not. ieee_is_nan(value),message)
	write(message,*) variable_name, " is infinite: ", value
	call ds_assert(ieee_is_finite(value),message)
#else
write(message,*) variable_name, " is NaN/infinite :", value
	call ds_assert(ds_isfinite(value),message)

#endif
	if (present(zero_is_bad)) then
		if (zero_is_bad) then
			write(message,*) variable_name, " is zero"
			call ds_assert(value/=0,message)
		endif
	endif

end subroutine testValueForEvil



function least_square_fit(y,n) result(mc)
	real(dp), dimension(2) :: mc
	real(dp), dimension(1:):: y
	real(dp) :: xysum,ysum,x2sum,xsum,m,c
	integer :: i,n


	xsum=n*(n+1.0)/2.0
	x2sum=n*(n+1.0)*(2*n+1.0)/6.0


	xysum=0.0
	ysum=0.0
	do i=1,n
		ysum=ysum+y(i)
		xysum=xysum+y(i)*i
	enddo

	m=(n*xysum - ysum*xsum) / (n*x2sum - xsum**2)
	c=ysum/n-m*xsum/n

	mc(1)=m
	mc(2)=c

end function least_square_fit





subroutine get_parallactic_angle(n,lat, dec, ra, par)
	!All arguments should be in radians
	integer n,i
	real(dp) :: lat
	external sla_pa
	double precision sla_pa
	real(dp), dimension(1:n) :: dec, ra, par
	do i=1,n
		par(i) = sla_PA (ra(i), dec(i), lat)
	enddo

end subroutine

subroutine save_name_mapping(list,filename)
	type(ds_moduleScanList) :: list
	type(ds_moduleScanInfo), pointer :: info
	character(*) :: filename
	integer unit
	
	unit=ds_get_lun()
	open(unit=unit,file=filename,status='replace')
	info => list%first
	do
!		write(unit,*) info%id1,info%id2,trim(info%filename),info%first_index,info%last_index
		
		if (associated(info%next)) then
			info=>info%next
		else
			exit
		endif
	enddo
	close(unit)
end subroutine save_name_mapping


	subroutine log_time_memory(name,debug)
		character(*) :: name
		logical debug
		character(256) :: message
		real(dp) :: mem
		if (.not. debug) return
		if (rank .ne. 0) return
		mem = get_mem_use()/1024d0**3
		write(message,fmt="(a,' [MEM:',f8.4,']')") trim(name), mem
		call ds_log_milestone(trim(message))
	end subroutine

  subroutine initialize_libquiet_modules(opt)
	  use quiet_mpi_mod
!	  use quiet_calib_mod
!	  use tod2map_mapmaker
!	  use tod2map_utils
	use quiet_pointing_mod
	  use quiet_target_mod
	  use quiet_i2qu_mod
	  use quiet_fileutils
	  use quiet_todsim_mod
	  use quiet_task_mod
	  use rngmod
	use quiet_ces_mod
	  use quiet_utils
	  use quiet_filter_mod
	  use quiet_assembly_mod
	  use quiet_mpi_mod
	  use quiet_patch_mod
	  use quiet_system_mod
	  use quiet_status_mod
	use quiet_target_mod
!	use quiet_validation_mod
!	use quiet_target_mod

	  use quiet_lx_mod
	
	  use ziggurat
!	use quiet_ces_mod
	
!	type(ces_index_set) :: targets
	type(ds_oslo_options) :: opt
	integer :: unit
    integer(i4b)     :: base_seed, i
    type(planck_rng) :: rng_handle
    real(dp)         :: mem
	logical :: debug
	debug = opt%debug

	unit = ds_get_lun()

    ! Set up auxilliary modules
    i = 14716
!    call get_fft3_magic_number(i, opt%fft3_magic_file)

!	call log_time_memory("noise estimation module",debug)
!    call initialize_noise_estimation_mod

	call log_time_memory('assembly module',debug)
    call initialize_quiet_assembly_mod(unit, opt%quiet_parameter_file)

	call log_time_memory('pointing module',debug)
    call initialize_quiet_pointing_mod(opt%quiet_parameter_file, apparent_correction = .true.)

	call initialize_quiet_hdf_mod()
	
	call log_time_memory('ces mod',debug)
	call initialize_ces_mod(opt%quiet_parameter_file)
	
	call log_time_memory('target module',debug)
    call initialize_target_mod(opt%quiet_parameter_file)

	call log_time_memory('pointing module',debug)
    call initialize_quiet_pointing_mod(opt%quiet_parameter_file, apparent_correction = .false.)

! 	call log_time_memory('mapmaker module',debug)
!     call initialize_map_maker(opt%quiet_parameter_file)

	call log_time_memory('gain module',debug)
    call initialize_gain_mod(opt%quiet_parameter_file)

!	call log_time_memory('scan database module',debug)
!    call initialize_target_mod(opt%quiet_parameter_file)

	call log_time_memory('filter module',debug)
    call initialize_filter_mod(opt%quiet_parameter_file)


	call log_time_memory('patch module',debug)
    call initialize_quiet_patch_mod(opt%quiet_parameter_file)



	call log_time_memory('module module (hee hee)',debug)
	call initialize_module_mod(opt%quiet_parameter_file)
!     if (opt%apply_i2qu) then
! 		call log_time_memory('i2qu module',debug)
!        call initialize_i2qu(unit, opt%quiet_parameter_file)
!     end if
!     if (analyze_simulated_data) then 
! 		call log_time_memory('tod sim module',debug)
!        call get_parameter(unit, opt%quiet_parameter_file, 'BASE_SEED', par_int=base_seed)
!        call initialize_random_seeds(MPI_COMM_WORLD, base_seed, rng_handle)
!        base_seed = nint(rand_uni(rng_handle)*1.d7)
!        call initialize_todsim_mod(unit, base_seed, opt%quiet_parameter_file)
!        call zigset(base_seed+1)
!     end if

    ! Read filelist for processing
! 	call log_time_memory('target list',debug)
!     call get_targets(opt%categories, opt%objects, opt%splits, targets)

  end subroutine initialize_libquiet_modules

subroutine get_fft3_magic_number(k, filename)
  implicit none

  integer(i4b),     intent(inout)          :: k
  character(len=*), intent(in),   optional :: filename

  integer(i4b)       :: i, unit, pos
  integer(i4b), save :: n
  integer(i4b), allocatable, dimension(:), save :: list

  if (present(filename)) then
     unit = getlun()
     open(unit, file=trim(filename))
     read(unit,*) n
     if (allocated(list)) deallocate(list)
     allocate(list(n))
     do i = 1, n
        read(unit,*) list(i)
     end do
  end if

  if (.not. allocated(list)) return

  pos = locate(list, k)
  if (pos > 1 .and. pos < n) then
     if (k - list(pos) < 50 .and. k-list(pos) >= 0) then
        k = list(pos)
     end if
  end if
end subroutine get_fft3_magic_number




!This will not correctly count the pixels, because some will be counted twice, I think.
!But I do not care as I only want to see which pixels are hit, and that will be right.
subroutine compute_hit_map_and_stats(file_list, file_owner, stats, diode_stats, nside, new_indices, total_hit_pixels, comm, rank)
	use quiet_target_mod
	use quiet_ces_mod
	use quiet_hdf_mod
	implicit none
	character(256), dimension(:) :: file_list
	integer, dimension(:) :: file_owner
	integer nside, rank
	integer comm
	integer(i4b), dimension(:),   allocatable :: hits, my_hits
	real(dp), dimension(:,:),   allocatable :: stats, my_stats
	real(dp),     dimension(:,:,:), allocatable :: diode_stats, my_diode_stats
    
	integer(i8b), dimension(:),   allocatable :: new_indices
	integer total_hit_pixels, i, j, p


	type(ds_moduleScanInfo), pointer	:: info
	type(hdf_file)                          :: file_pointer
	integer(i4b), dimension(:),   allocatable :: file_pixel_size, pixels
	integer :: degrade_factor, degrade_steps
	integer :: ierror, file_nside, npix, lowres_index
	integer number_files
	integer ndiodes
	
	
	ndiodes = get_num_modules() * get_num_diodes()
	
	number_files = size(file_list)
	npix = 12*nside**2
	
	
	allocate(hits(0:npix-1))
	allocate(my_hits(0:npix-1))
	hits = 0
	my_hits = 0

	allocate(stats(number_files,STAT_NUM))
	allocate(my_stats(number_files,STAT_NUM))
    stats = 0.0
	my_stats = 0.0
	
	allocate(diode_stats(ndiodes,number_files,NUM_MOD_STATS))
	allocate(my_diode_stats(ndiodes,number_files,NUM_MOD_STATS))
	
	diode_stats = 0.0
	my_diode_stats = 0.0
	
	do i=1,size(file_list)
		if (file_owner(i) .ne. rank) cycle 
		call open_hdf_file(file_list(i), file_pointer, "r")
		call ds_assert(file_pointer%status==0,"Failed to open file: "//trim(file_list(i)))
		call get_size_hdf(file_pointer, "pixels", file_pixel_size)
		allocate(pixels(file_pixel_size(1)))
		deallocate(file_pixel_size)

		call read_hdf(file_pointer, "pixels", pixels)
		call read_hdf(file_pointer, "nside", file_nside)
		degrade_steps = nint(log(real(file_nside/nside,dp))/log(2d0))*2
		degrade_factor = 2**degrade_steps
		do p = 1, size(pixels)
			lowres_index = pixels(p)/degrade_factor
			my_hits(lowres_index) = my_hits(lowres_index) + 1
	     end do
		deallocate(pixels)
		call read_hdf(file_pointer, "stats", my_stats(i,:))
		call read_hdf(file_pointer, "diode_stats", my_diode_stats (:,i,:))
       
		call close_hdf_file(file_pointer)
	enddo


	call MPI_Allreduce(my_hits, hits, npix, MPI_INTEGER, MPI_SUM, comm, ierror)
    call MPI_Allreduce(stats, my_stats, size(stats), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
    call MPI_Allreduce(diode_stats, my_diode_stats, size(diode_stats), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)

	deallocate(my_hits)
	deallocate(my_stats)
	deallocate(my_diode_stats)

	allocate(new_indices(0:npix-1))
	new_indices=-1
	j = 1
	do i = 0, npix-1
		if(hits(i) .lt. 2) cycle
		new_indices(i) = j
		j = j+1
  end do
  deallocate(hits)
  total_hit_pixels = j-1
end subroutine compute_hit_map_and_stats




function count_my_filenames(scanList, rank) result(n)
	type(ds_moduleScanList) :: scanList
	type(ds_moduleScanInfo), pointer :: info
	integer rank
	integer i,n
	character(len=256) :: filename
	info => scanList%first
	filename=""
	n=0
	do i=1,scanList%length
		if (filename .ne. info%filename) then
			filename=info%filename
		 	if (info%owner==rank) n=n+1
		endif
		if (associated(info%next)) info => info%next
	enddo
end function 


subroutine collect_my_filenames(scanlist, rank, my_files)
	type(ds_modulescanlist) :: scanlist
	integer rank
	character(256), dimension(:), allocatable :: my_files
	character(256) :: filename
	integer n, i, j
	type(ds_modulescaninfo), pointer :: info
	
	n=count_my_filenames(scanList, rank)
	allocate(my_files(n))
	filename=""
	j=1
	info=>scanList%first
	do i=1,scanList%length
		if (filename .ne. info%filename) then
			filename=info%filename
			if (info%owner==rank) then
				my_files(j)=filename
				j=j+1
			endif
			
		endif
		if (associated(info%next)) info => info%next
	enddo
	
end subroutine collect_my_filenames


function count_all_filenames(scanList) result(n)
	type(ds_moduleScanList) :: scanList
	type(ds_moduleScanInfo), pointer :: info
	integer rank
	integer i,n
	character(len=256) :: filename
	info => scanList%first
	filename=""
	n=0
	do i=1,scanList%length
		if (filename .ne. info%filename) then
			n=n+1
			filename=info%filename
		endif
		if (associated(info%next)) info => info%next
	enddo
end function 


subroutine collect_all_filenames(scanlist, files)
	type(ds_modulescanlist) :: scanlist
	integer rank
	character(256), dimension(:), allocatable :: files
	character(256) :: filename
	integer n, i, j
	type(ds_modulescaninfo), pointer :: info
	
	n=count_all_filenames(scanList)
	allocate(files(n))
	filename=""
	j=1
	info=>scanlist%first
	do i=1,scanList%length
		if (filename .ne. info%filename) then
			filename=info%filename
			files(j)=filename
			j=j+1			
		endif
		if (associated(info%next)) info => info%next
	enddo
	
end subroutine collect_all_filenames


end module ds_oslotools

