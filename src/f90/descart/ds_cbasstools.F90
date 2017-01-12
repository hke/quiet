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
	use ds_cbass_option_utils
#ifndef __GFORTRAN__
	use ieee_arithmetic
#endif
	use fl_lists
	implicit none

	real(dp), parameter :: DEGRA = 0.0174532925
! 	real(dp), parameter :: CBASS_NORTHERN_LATITUDE_DEG = 37.233888  !From the OVRO website
! 	real(dp), parameter :: CBASS_NORTHERN_LATITUDE_RAD = CBASS_NORTHERN_LATITUDE_DEG * DEGRA
	integer, parameter :: FITS_READ_ONLY = 0
	character(8), parameter :: NAXIS2_KEYWORD = "NAXIS2  "

	
#if 0
#pragma mark -
#pragma mark DATA TYPES
#endif	
	type ds_moduleScanInfo
		integer id1,id2,id3
		integer ndiodes
		integer owner
		logical needs_azfilter
		character(len=256) filename
		integer modScan !The local index into moduleScans
		type(ds_moduleScanInfo),pointer :: next
		integer n
		type(ds_noiseInfo) :: noise
	end type ds_moduleScanInfo

	type ds_file_data
	     real(dp), dimension(:), allocatable :: I1,I2,Q1,Q2,U1,U2
	     real(dp), dimension(:), allocatable :: MJD,Ra,Dec,Az,El,theta
	     integer n
	end type ds_file_data



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
	subroutine readDataAssignWork(fileList,scans,noiseInfo,correlator,originalIndices,maps,covariance,opt,output_fits_cards)
#if 0	
#pragma mark -
#endif
	! fileList - input, array of strings of length 256
	! scans - output, array of (module,scan) pairs (ds_modulescan)
	! noiseInfo - output, array of noise info structures (ds_noiseinfo)
	! correlator - input/output.  object containing data on which processor owns what, and noise correlations.  On input, should have %comm set.  Other stuff set on output
	! originalIndices - output.  contains map from our pixel numbering to healpix.
	! maps - output.  (I,Q,U map.)  Contains naive map on output.
	! covariance - output.  the naive map covariance. (ds_covariance)
	! opt - input.  The options file.
	! output_fits_cards - output.  List of cards to propagate to the FITS header of the output map.
	

	!This is the only subroutine that you need to re-write to use 


		implicit none
		!Subroutine arguments
		type(ds_moduleScan), pointer, dimension(:) :: scans
		type(ds_noiseInfo), pointer, dimension(:) :: noiseInfo
		type(ds_correlator) :: correlator
		integer(i8b), allocatable, dimension(:) :: originalIndices
		type(ds_trimap) :: maps, buffer_maps
		type(ds_covariance) :: covariance
		type(ds_cbass_options) :: opt
		integer(i8b) :: full_map_pixels
		character(256), dimension(:) :: fileList
		type(fl_string80_list) :: output_fits_cards

		logical, parameter ::  do_correlations = .false.
		integer progress

		!Information about the data to be loaded.
		integer nmodules
		integer(i4b) npix

		!Work assignment information
		integer rank,nproc
		integer nmin,nextra
		integer nScan,my_nScan
		integer, allocatable, dimension(:) :: nScan_proc

		!Things needed to load the data from file
		type(module_struct), dimension(:), allocatable :: full_file_data
		logical haveBadScans, thisScanBad, needsAz
		integer status
		integer ordering
		integer scan_number_in_file
		character(256) :: current_filename

		!Things needed to process the data and put it in the moduleScan structure.
		real(dp) :: gain(CBASS_NDIODES_MAX)
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


		!Set up the size and rank information and save it in the correlator too.
		!This requires that correlator%comm has already been set in the driver.
		unit = ds_get_lun()
		call MPI_Comm_size(correlator%comm, nproc,ierr)
		call MPI_Comm_rank(correlator%comm, rank, ierr)
		correlator%nproc = nproc
		correlator%proc = rank

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
		call ds_assert(.not. associated(noiseInfo),"Noise informatio should not already be allocated in readDataAssignWork")
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
		
		if (rank==0 .and. opt%cut_extension_name/="") call ds_log("Reading scans from extension: "//trim(opt%cut_extension_name), ds_feedback_quiet)
		
		call buildModuleScanList(scanList, fileList, opt, output_fits_cards)
		if (opt%save_name_mapping) call save_name_mapping(scanList,opt%name_mapping_file)

		!This is the total number of moduleScans for all the processors

		nScan = scanList%length
		if(rank==0) then
			write(message,*) "Total number of (module,scans) = ", nScan
			call ds_log(message,ds_feedback_quiet)
		endif

		if (nScan==0) then
		    write(*,*) "NO SCANS!  CANNOT MAKE A MAP WITHOUT SCANS!"
		    call MPI_Finalize(ierr)
		    stop
		endif		


		!Determine the number of moduleScans each process is responsible for.
		!Each process knows this number for each process.
		!nmin is the minimum number of moduleScans held by a process.
		!Some have more to make up the difference.
		allocate(nScan_proc(0:nproc-1))
		nmin=floor((1.0*nScan)/nproc)
		nextra = nScan - nmin*nproc
		nScan_proc = nmin
		call ds_assert(nextra>=0,"Broken nextra - this is Joe's fault.")
		if (nextra>0) nScan_proc(0:nextra-1) = nmin +1

		!Some useful checks and feedback about the allocation of moduleScans
		call ds_assert(sum(nScan_proc) == nScan, "Mistake in assignment of dets - not all assigned.")

		scanInfo=>scanList%first
		do i=0,nproc-1
			do j=0,nScan_proc(i)-1
				scanInfo%owner=i
				scanInfo=>scanInfo%next
			enddo
		enddo
		
		my_nScan = nScan_proc(rank)
		correlator%my_nmodules = my_nScan
		if(rank==0) then
			write(message,*) "Approx number of modules per proc:", my_nScan
			call ds_log(message,ds_feedback_quiet)
		endif
		

#if 0
#pragma mark -  ALLOCATE SPACE
#endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Setup space for this proc
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!Now we have that info we check how many moduleScans this process has and allocate space
		!and initialize the pointers inside
		allocate(scans(0:my_nScan-1))
		allocate(noiseInfo(0:my_nScan-1))
		do ms=0,my_nScan-1
			scan=>scans(ms)
			call setup_moduleScan_ndiodes(scan,CBASS_NDIODES_T,CBASS_NDIODES_P)
			nullify(scan%theta)
			nullify(scan%pointing) !DWPS: always nullify pointers when you make them
			allocate(noiseInfo(ms)%sigma(CBASS_NDIODES_MAX))
			allocate(noiseInfo(ms)%fknee(CBASS_NDIODES_MAX))
			allocate(noiseInfo(ms)%alpha(CBASS_NDIODES_MAX))
			allocate(noiseInfo(ms)%corrs(CBASS_NDIODES_MAX,CBASS_NDIODES_MAX))
		enddo


		if ((.not. do_correlations) .and. correlator%proc==0) write(*,*) "Cross-correlation deactivated"

		!Allocate the maps
		call prepareTriMap(maps,npix, opt%do_temperature, opt%do_polarization, zero_based = .true.)  !Healpix maps are zero-based.  For now this is a healpix map.
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
		scanInfo=>scanList%first
		haveBadScans=.false.
		progress=0
		
		!Loop through the moduleScanList.  Check if this proc is responsible for the corresponding moduleScan.
		!If so, read it in and process it.
		do 
			progress = progress +1
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! End loop if finished
			! Cycle loop if not owner of this scan
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			if (.not. scanInfo%owner==rank) then
				if (.not. associated(scanInfo%next)) exit
				scanInfo=>scanInfo%next
				cycle
			endif

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Load new FITS file
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!Load the new file data in if this is a new filename.
			if (.not. current_filename==scanInfo%filename  )     then
				
				if (current_filename .ne. "") then
					stop 'dealloc l2 data full_file_data'
				endif
				scan_number_in_file = 1

        		ces     = get_ces_info(targets%list(k)%ces(i))
				call l2_read(scanInfo%filename, data=full_file_data,status=ierr)
	           
				write(message,'("Loaded file: ",A)') ,trim(scanInfo%filename)
				call ds_log(message,ds_feedback_noisy)
				if (rank==0) then
					write(message,'(A,F7.1,A)') "Approximate load progress: ", (progress*100.0)/my_nScan,"%"
					call ds_log(message,ds_feedback_noisy)
				endif

				thisScanBad=.false.
				needsAz = opt%use_azimuth_offsets				
				current_filename = scanInfo%filename
			endif
			
			!Find the correct location in the tod_data for the current diode data, 
			!and the correct moduleScan object to save that data in.

#if 0
#pragma mark - - BUILD STRUCTURE
#endif

			ms = scanInfo%modScan
			scan=>scans(ms)
			write(message, '(A,I3)') " - Loaded scan ", scan_number_in_file
			call ds_log(message,ds_feedback_debug)
			scan_number_in_file = scan_number_in_file + 1
			ntod = scanInfo%last_index - scanInfo%first_index + 1

			if (opt%should_cut_timestreams) ntod = ntod * opt%cutfraction
			na = ntod/opt%offsetLength
			ntod = na * opt%offsetLength
			
			scan%ntod=ntod
			call buildNoiseInfo(scanInfo,NoiseInfo(ms),gain,mjd,opt)
			scan%owns_az=.false.
			stop 'Extract.'
			call get_scan_from_fits(scanInfo,full_file_data,noiseinfo(ms),gain,scan,opt,needsAz, opt%do_temperature, opt%do_polarization)

			if (opt%simulate) call simulate_moduleScan(scan,noiseinfo(ms),simulation, gain, opt)




			if (needsAz .and. scan%owns_az) then
				needsAz=.false.
!				azOwner=>scan
			endif
!			call share_azimuth(scan,azOwner)
			scan%inv_Cw=0.0

			!This is a good place to get the priors if we want to do it directly from the data.
			if (opt%data_prior) call prepare_one_data_prior(correlator,opt%offsetLength,noiseinfo(ms),scan)

			
			!build inverse covariance matrix from variances and correlations
			!we already filled in the sigmas above
			call make_inv_Cw(noiseInfo(ms),scan,do_correlations)
			!apply the inverse covariance matrix in moduleScan to the timestreams in moduleScan
			call invCw_mult_tod(scan)



#if 0
#pragma mark - -  ACCUMULATE MAP
#endif
			!add module scan to accumulated healpix map
			call add2rhs(scan,maps)
			!deproject moduleScan%timestream -> moduleScan offsets
			!When we have done this we no longer need the timestream so we deallocate it		
			do i=1,CBASS_NDIODES_MAX
				if (scan%flags(i)==1) then
					call deprojectTimestreamOntoOffset(scan%timestreams(i),scan%offsets_2(i))
					call destroyTimestream(scan%timestreams(i))
				endif
			enddo
			if (.not. associated(scanInfo%next)) exit
			scanInfo=>scanInfo%next
		enddo

		deallocate(nScan_proc)

#if 0
#pragma mark - FINALIZE MAP
#endif
		write(message,*) "Rank ", rank," loaded all data."
		call ds_log(message,ds_feedback_debug)

		!If any of the scans were bad and did not load, quit here
		!The user should supply a better accepted list!
		!We may want to reconsider this behaviour.
		call ds_assert(.not. haveBadScans, "Bad scans reported")

		if (rank==0) call ds_log_milestone("DATA_LOADED")



		!Repixelize both the pointings and the accumulated maps
		if (opt%plane_map) then
			full_map_pixels = opt%plane_size**2
		else
			full_map_pixels = nside2npix(opt%nside)
		endif
		
        if (opt%save_hits) then
           call repixelizeData(correlator,scans,full_map_pixels,originalIndices,maps,hitMapFilename=opt%hits_filename)           
        else
           call repixelizeData(correlator,scans,full_map_pixels,originalIndices,maps)
        endif


		
		if (opt%do_temperature) then
			npix = maps%T%npix
		else 
			npix = maps%Q%npix
		endif
		
		
		!Sum the accumulated maps accross all the processes, thereby summing for all modules for all scans
		if (opt%do_temperature)  call MPI_Allreduce(MPI_IN_PLACE, maps%T%map, maps%T%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
		if (opt%do_polarization) call MPI_Allreduce(MPI_IN_PLACE, maps%Q%map, maps%Q%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
		if (opt%do_polarization) call MPI_Allreduce(MPI_IN_PLACE, maps%U%map, maps%U%npix, MPI_REAL8,MPI_SUM,correlator%comm,ierr)

		call make_naivecov(npix,covariance,correlator,scans, opt%do_temperature, opt%do_polarization)
		call cov_mult(maps,covariance)



		!Subtract the naive map from the offsets?
		do ms=0,correlator%my_nmodules-1
			scan=>scans(ms)
			do diode=1,CBASS_NDIODES_MAX
				if (scan%flags(diode)==1) call prepareTimestream(scan%timestreams(diode),scan%ntod )
			enddo
			
			call map2tod(scan,maps)
			call invCw_mult_tod(scan)

			do diode=1,CBASS_NDIODES_MAX
				if (scan%flags(diode)==1 ) then
					call subtractTimestreamFromOffset(scan%timestreams(diode),scan%offsets_2(diode))
					call destroyTimestream(scan%timestreams(diode))
				endif
			enddo

		enddo
		
		call destroyModuleScanList(scanList)



	end subroutine readDataAssignWork
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




subroutine get_scan_from_fits(info,full_data,noise,gain,moduleScan,opt,needsAz, do_T, do_P)
  use quiet_assembly_mod

	type(ds_cbass_options) :: opt
	real(dp),dimension(CBASS_NDIODES_MAX) :: gain
	type(ds_moduleScanInfo) :: info
	logical :: do_T, do_P
	type(ds_noiseInfo) :: noise
	type(ds_moduleScan) :: moduleScan
	type(module_struct), dimension(:) :: full_data
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
	type(common_info) :: shared_info
	logical :: do_T, do_P
	
	do_T = opt%do_temperature
	do_P = opt%do_polarization
	coord_in = parse_coord_name(coord_name_in)
	coord_out = parse_coord_name(coord_name_out)
	call ds_generate_common_info(info,shared_info)
	call reorganize_data_assembly(do_T, opt%coord_in, opt%coord_out, opt%nside, full_data, assembly, shared_info)
	
	
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

	ntod = assembly%num_samp


	if (opt%should_cut_timestreams) then
		ntod = ntod * opt%cutfraction
	endif
	
	
	call ds_assert(ntod>0,"Zero-sized (or negative) TOD passed to get_scan_from_fits")
	
	na = ntod/opt%offsetLength
	ntod = na * opt%offsetLength
	finish = start+ntod-1
	
	!Set up the offsets
	do i=1,CBASS_NDIODES_MAX
		if (moduleScan%flags(i)==1) then
			call prepareOffsets(moduleScan%offsets_2(i),na,opt%offsetLength,azimuth_flag=.false.)
		endif
	enddo
	
	
	!Set up the pixel pointing
	allocate(moduleScan%pointing(ntod))
	stop 'HERE: Get pointing from assembly'
	
	do t=1,ntod,
		moduleScan%pointing(t) = assembly%pix(1,t,1)  !Indices: (main/partner),sample,diode  !NB We ignore partner and assume all diodes point at the same thing for now.
	enddo

	
	moduleScan%inv_Cw=0.0

	if (do_P) then
		allocate(moduleScan%theta(ntod))
		do t=1,ntod
			moduleScan%theta(t) = assembly%point(3,t,1,1)  !Indices : (phi/theta/psi),sample,diode,(main/partner)
		enddo
		call read_dpsi(modulescan%dpsi)
	endif
	
	do diode=1,assembly%num_diodes
		! if ((.not. do_T) .and. (diode .le. ndiodes_T)) cycle
		! if ((.not. do_P) .and. (diode .gt. ndiodes_T)) cycle
		
		call prepareTimestream(moduleScan%timestreams(diode),ntod)
		
		if (opt%simulate) cycle

		do t=1,ntod
			moduleScan%timestreams(diode) = assembly%tod(t,diode)
		enddo
		
		do t=1,ntod
			moduleScan%timestreams(diode)%timestream(t) = moduleScan%timestreams(diode)%timestream(t) / gain(diode)
		enddo

		!Subtract a linear trend from the data
		if (opt%subtract_linear) then
			linearFit=least_square_fit(moduleScan%timestreams(diode)%timestream,ntod)
			write(message,*) "Subtracting linear:",linearFit
			call ds_log(message, ds_feedback_debug)
			do t=1,ntod
				moduleScan%timestreams(diode)%timestream(t)=moduleScan%timestreams(diode)%timestream(t)-linearFit(1)*t-linearFit(2)
			enddo
		endif

		!Subtract the mean of the data
		if(opt%subtractMeans) then				
			mu = sum(moduleScan%timestreams(diode)%timestream)/ntod
			do t=1,ntod
				moduleScan%timestreams(diode)%timestream(t) = &
				moduleScan%timestreams(diode)%timestream(t) - mu						
			enddo
		endif
		
		
	enddo
	stop 'HERE: Deallocate assembly'	
	
end subroutine get_scan_from_fits


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





subroutine read_dpsi(dpsi)
	real(dp), dimension(CBASS_NDIODES_P) :: dpsi
	integer diode

	stop 'read dpsi'

end	subroutine read_dpsi


#if 0
#pragma mark -
#pragma mark LISTS OF SCANS TO PROCESS
#pragma mark -
#endif
	subroutine buildModuleScanList(modScanList,filename_list, opt, cards)
	! The job of this function is to build up a list of scans to be read.
	! These are not yet read at this stage but will be later.
	! The modScanList object is a linked-list of modScanInfo objects.
	
	! The CBASS data format has data in the first extension header of FITS files and
	! in the other headers it has cut information that defines which chunks 
	! of the primary to use.
	
	   use fits_helper
		use fl_lists
		!!!
		! Read the headers of all the input files and find out how many scans are in each and how to find them.
		!!!

		type(ds_moduleScanList) :: modScanList
		type(ds_cbass_options) :: opt
		type(fl_string80_list), optional :: cards

		integer, dimension(CBASS_NDIODES_MAX) :: diode_flag
		character(256), dimension(:) :: filename_list
		character(256) :: filename
		integer id1,id2
		integer unit,hdutype
		integer nelements
		character(80) :: comment
		integer f, status
		character(8) :: keyword
		character(8) :: value
		logical :: do_T, do_P
		integer nscan,s
		integer i
		integer number_of_hdu, cuts_hdu
		integer, parameter :: number_noise_parameters = 3
		real(dp), allocatable, dimension(:,:,:) :: noise_data
		

		do_T = opt%do_temperature
		do_P = opt%do_polarization
		
		call init_moduleScanList(modScanList)
		

		do  !Loop through tasks.
		
			id1 = 
			id2 = 
			id3 = 
			ndiodes = 
			filename = 
			
			call appendModuleScanInfo(modScanList,id1,id2,id3,ndiodes,.false.,filename)

			modScanList%last%noise%sigma =
			modScanList%last%noise%fknee =
			modScanList%last%noise%alpha =
			modScanList%last%noise%corrs = 
			do i=1,ndiodes
				modScanList%last%noise%corrs(i,i)=
			enddo
			
			
		enddo


	end subroutine buildModuleScanList



subroutine init_moduleScanList(L)
	type(ds_moduleScanList) :: L
	nullify(L%first)
	nullify(L%last)
	L%length=0
end subroutine init_moduleScanList

function makeModuleScanInfo(id1,id2,id3,ndiodes,needs_azfilter,filename) result(output)
	type(ds_moduleScanInfo), pointer :: output
	integer id1,id2,id3,ndiodes
	logical needs_azfilter
	character(len=256) filename
	integer, dimension(CBASS_NDIODES_MAX) :: diode_flag

	allocate(output)

	output%id1=id1
	output%id2=id2
	output%id3=id3
	output%ndiodes = ndiodes
	output%filename=filename
	
	allocate(output%noise%sigma(ndiodes))
	allocate(output%noise%fknee(ndiodes))
	allocate(output%noise%alpha(ndiodes))
	allocate(output%noise%corrs(ndiodes_max,ndiodes_max))
	
	nullify(output%next)
	
end function makeModuleScanInfo

subroutine appendModuleScanInfo(list,id1,id2,id3,ndiode,needs_azfilter,filename)
	type(ds_moduleScanList) :: list
	integer id1,id2,id3,ndiode
	logical needs_azfilter
	character(len=256) filename
	type(ds_moduleScanInfo), pointer :: info
	info => makeModuleScanInfo(id1,id2,id3,ndiode,needs_azfilter,filename)
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
			do d=1,CBASS_NDIODES_MAX
				if (moduleScans(i)%flags(d)==1) then
					do t=1,moduleScans(i)%ntod
						p=moduleScans(i)%pointing(t)
						if (p>=0) hitCount(p)=hitCount(p)+1
					enddo
				endif
			enddo
		enddo

		!The processes all sum their hit count map to get the global hit count map.
		ierror=0
		call MPI_AllReduce(MPI_IN_PLACE,hitCount,maxIndex,MPI_INTEGER,MPI_SUM,correlator%comm,ierror)
		call ds_assert(ierror==0,"Error in MPI_AllReduce in repixeizeData")

		!If requested, the hit count map is saved to file
		if (present(hitMapFilename) .and. correlator%proc==0) then
			!Save hit map
           open(file=hitMapFilename,unit=25)
           do p=0,maxIndex-1
              write(25,*) p,hitCount(p)
           enddo
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
			do d=1,CBASS_NDIODES_MAX
				if (moduleScans(i)%flags(d)==1) then
					do t=1,moduleScans(i)%ntod
						if(hitCount(moduleScans(i)%pointing(t))==1) then
							moduleScans(i)%pointing(t) = bad_pixel
						endif
					enddo
				endif
			enddo
		enddo

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


		!Repixelize the maps, reducing them in size from the full healpix sized ones to
		!ones containing only the pixels hit somewhere in one of the timestreams
		if (maps%has_T) call repixelizeMap(maps%T,hitCount,n)
		if (maps%has_P) call repixelizeMap(maps%Q,hitCount,n)
		if (maps%has_P) call repixelizeMap(maps%U,hitCount,n)
		
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

			scan%flags = scanInfo%diode_flag
			if (.not. associated(scanInfo%next)) exit
			scanInfo=>scanInfo%next
			ms=ms+1				
		enddo

	end subroutine setupScans





subroutine buildNoiseInfo(info,noise,gain,mjd,opt)
	type(ds_moduleScanInfo) :: info
	type(ds_noiseInfo) :: noise
	type(ds_cbass_options) :: opt
	real(dp) :: gain(CBASS_NDIODES_MAX)
	integer, parameter :: number_of_noise_parameters = 3
	real(dp) :: mjd
	integer diode,diode2
	
	call ds_assert(associated(noise%sigma),"Noise array (sigma) not associated in buildNoiseInfo")
	call ds_assert(associated(noise%fknee),"Noise array (fknee) not associated in buildNoiseInfo")
	call ds_assert(associated(noise%alpha),"Noise array (alpha) not associated in buildNoiseInfo")


	!noise_info
	do diode=1,CBASS_NDIODES_MAX

		if (info%diode_flag(diode)<=0) cycle
#warning Gain set to one. Should be read from file
		gain(diode) = 1.0
		noise%sigma(diode) = info%noise%sigma(diode)
		noise%alpha(diode) = info%noise%alpha(diode)
		noise%fknee(diode) = info%noise%fknee(diode)

		call testValueForEvil(gain(diode),"gain",zero_is_bad=.true.)
		call testValueForEvil(noise%sigma(diode),"sigma",zero_is_bad=.true.)
		call testValueForEvil(noise%fknee(diode),"fknee",zero_is_bad=.false.)
		call testValueForEvil(noise%alpha(diode),"alpha",zero_is_bad=.false.)
	enddo
	noise%corrs = info%noise%corrs
	
	deallocate(info%noise%sigma)
	deallocate(info%noise%alpha)
	deallocate(info%noise%fknee)
	deallocate(info%noise%corrs)
	

end subroutine buildNoiseInfo



#if 0
#pragma mark -
#pragma mark SIMULATION FUNCTIONS
#pragma mark -
#endif

subroutine simulate_moduleScan(moduleScan,noise,sim,gain,opt)
	type(ds_moduleScan) :: moduleScan
	type(ds_noiseInfo) :: noise
	type(ds_cbass_options) :: opt
	type(ds_trimap) :: sim
!	real(dp), dimension(:) :: az
!	integer rank
	real(dp), dimension(CBASS_NDIODES_MAX) :: gain
	integer diode
	integer nbin
	integer unit
	integer t
	logical is_temperature
	
	logical, save :: first_sim = .false.
	
	do diode=1,CBASS_NDIODES_MAX
		if (moduleScan%flags(diode) .ne. 1) cycle
		
		do t=1,moduleScan%timestreams(diode)%nt
			moduleScan%timestreams(diode)%timestream(t) = 0.0
		enddo
	
		call ds_simulate_white(moduleScan%timestreams(diode), noise%sigma(diode) )
		call ds_simulate_signals(moduleScan,sim, gain, opt%do_temperature, opt%do_polarization)

	enddo

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
	real(dp), dimension(CBASS_NDIODES_MAX) :: gain
	logical do_T, do_P
	integer diode, t, p
	integer npix_sim, nside_sim
	
	call ds_assert(do_T .or. do_P, "Must simulate either temperature or polarization")
	
	if (do_T) call ds_assert(sim%has_T, "No temperature simulation map provided in ds_simulate_signal")
	if (do_P) call ds_assert(sim%has_P, "No polarization simulation map provided in ds_simulate_signal")

	if (do_T) then
		npix_sim = sim%T%npix
	else
		npix_sim = sim%Q%npix
	endif
	
	nside_sim = npix2nside(npix_sim)
	
	
! TODO: Add support for simulating from different resoluation maps	
	
	if (do_T) then
		do diode=1,CBASS_NDIODES_T
			if (scan%flags(diode) .ne. 1) cycle
			do t=1,scan%timestreams(diode)%nt
				p=scan%pointing(t)
				scan%timestreams(diode)%timestream(t) = scan%timestreams(diode)%timestream(t) + sim%T%map(p) * gain(diode)
			enddo
		enddo
	endif

	if (do_P) then
		do diode=CBASS_NDIODES_T+1,CBASS_NDIODES_MAX
			if (scan%flags(diode) .ne. 1) cycle
			do t=1,scan%timestreams(diode)%nt
				p=scan%pointing(t)
				stop 'Have not coded simulated polarizations.  Not too hard though'
				scan%timestreams(diode)%timestream(t) = scan%timestreams(diode)%timestream(t) + sim%T%map(p) * cos(2*scan%theta(t)) * gain(diode)
			enddo
		enddo
	endif
	
	
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
		write(unit,*) info%id1,info%id2,trim(info%filename),info%first_index,info%last_index
		
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
		mem = get_mem_use()/1024d0**3
		write(message,fmt="(a,f8.4,a)") ' mem: ', mem, 'noise estimation module'
	    call ds_log_milestone(trim(message))
	end subroutine

  subroutine initialize_libquiet_modules(parfile,opt,targets)
	  use quiet_mpi_mod
	  use quiet_calib_mod
	  use tod2map_mapmaker
	  use tod2map_utils
	  use quiet_target_mod
	  use quiet_i2qu_mod
	  use quiet_fileutils
	  use quiet_assembly_mod
	  use quiet_todsim_mod
	  use quiet_task_mod
	  use rngmod
	  use ziggurat

    integer(i4b)     :: base_seed, i
    type(planck_rng) :: rng_handle
    real(dp)         :: mem
	logical :: debug
	debug = opt%debug

    ! Set up auxilliary modules
    i = 14716
    call get_fft3_magic_number(i, opt%fft3_magic_file)

!	call log_time_memory("noise estimation module",debug)
!    call initialize_noise_estimation_mod

	call log_time_memory('assembly module',debug)
    call initialize_quiet_assembly_mod(unit, parfile)
	
	call log_time_memory('validation module',debug)
    call initialize_quiet_validation_mod(unit, parfile)

	call log_time_memory('pointing module',debug)
    call initialize_quiet_pointing_mod(unit, parfile, apparent_correction = .not. approximate_hor)

	call log_time_memory('mapmaker module',debug)
    call initialize_map_maker(parfile)

	call log_time_memory('calibration module',debug)
    call initialize_quiet_calib_mod(unit, parfile)

	call log_time_memory('gain module',debug)
    call initialize_gain_mod(unit, parfile)

	call log_time_memory('scan database module',debug)
    call initialize_target_mod(parfile)

	call log_time_memory('filter module',debug)
    call initialize_filter_mod(unit, parfile)

    if (apply_i2qu) then
		call log_time_memory('i2qu module',debug)
       call initialize_i2qu(unit, parfile)
    end if
    if (analyze_simulated_data) then 
		call log_time_memory('tod sim module',debug)
       call get_parameter(unit, parfile, 'BASE_SEED', par_int=base_seed)
       call initialize_random_seeds(MPI_COMM_WORLD, base_seed, rng_handle)
       base_seed = nint(rand_uni(rng_handle)*1.d7)
       call initialize_todsim_mod(unit, base_seed, parfile)
       call zigset(base_seed+1)
    end if

    ! Read filelist for processing
	call log_time_memory('target list',debug)
    call get_targets(opt%categories, opt%objects, opt%splits, targets)

  end subroutine initialize_libquiet_modules



end module ds_oslotools

