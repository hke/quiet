program master_filter
	use ds_fitstools
	use ds_multidetector
	
	implicit none
	integer narg,ierror,rank,arg,d,nside,ndet,n
	character(256) :: todname,todroot,parameter_filename,noiseFile,pixelFile,pointroot,pointfile,outputName, filtered_tod_name
	type(ds_timestream), allocatable, dimension(:) :: timestream	
	type(ds_correlator) :: correlator
	type(ds_detpointing), dimension(:), allocatable :: pointing
	type(ds_map) :: qmap, umap
	real(dp), allocatable, dimension(:) :: whiteNoise, fknee, filter
	real(dp), allocatable, dimension(:,:) :: rho
	real(dp) :: nyquist, delta_t, datarate
	integer,allocatable,dimension(:) :: originalIndices
	integer npix,p
	integer(dp) :: totpix
	
	narg = command_argument_count()
	call ds_assert(narg > 0, "Program Usage: descarts_quiet parameter_file [more parameter files ... ]")

	call MPI_Init(ierror)
	call mpi_comm_rank(MPI_COMM_WORLD,rank,ierror)
	call ds_assert(ierror==0, 'Failed to init mpi in descarts_quiet')


	do arg=1,narg


	call get_command_argument(arg,parameter_filename)

	!read parameters
	call read_params
	allocate(whiteNoise(0:ndet-1))
	allocate(fknee(0:ndet-1))
	allocate(rho(0:ndet-1,0:ndet-1))

	call read_noisefile(noiseFile,ndet,nyquist,whiteNoise,fknee,rho)
	delta_t = 1.0_dp / (2*nyquist)
	correlator%comm = MPI_COMM_WORLD
	call assign_detectors(correlator,ndet)

	allocate(timestream(0:correlator%my_ndet-1))
	allocate(pointing(0:correlator%my_ndet/2-1))
	
	do n=0,correlator%my_ndet-1
	
		d = correlator%my_det(n)
		todname = trim(todroot) // trim(intToStr(d)) // ".fits"
		call read_simple_tod(todname, timestream(n))
  write(*,*) 'loaded ',trim(todname)
		if (.not. allocated(filter) ) allocate(filter(0:timestream(n)%nt-1))
  
        
		call applyFilter(timestream(n),nyquist,fknee(d))!,whiteNoise(d))  !DWPS: filter should have whiteNoise of unity
!		timestream(n)%timestream = filter * timestream(n)%timestream
		
		if (.not. (filtered_tod_name == '')) then
			todname = trim(filtered_tod_name) // trim(intToStr(d)) // ".fits"
			call write_simple_tod(todname,delta_t,nside,timestream(n))
		endif
	
	enddo
	if (correlator%proc==0) call ds_log("TOD read and filtered",ds_feedback_quiet)

! if(correlator%proc==0) then
!    print*,'writing filtered tod to file...'
!    open(unit=1234,file="filtered_tod.bin",form="binary",status="replace")
!    write(1234) timestream(0)%timestream
!    close(1234)
! endif
! call mpi_barrier(correlator%comm,ierror)
! stop


	do d=0,correlator%my_ndet/2-1
		p = correlator%my_det(det_qindex(d)) / 2
		pointfile = trim(pointroot) // "pair_" // trim(intToStr(p)) // ".fits"
		call readDetPointing(pointfile,nside,pointing(d),datarate)
	enddo
	if (correlator%proc==0) call ds_log("Pointing read",ds_feedback_quiet)


	totPix = 12*nside*nside
	call makeSharedPointingContinuous(correlator,pointing,totPix,pixelFile,originalIndices)
	npix = pointing(0)%npix
	if (correlator%proc==0) call ds_log("Shared pointing complete.",ds_feedback_quiet)



	call prepareMap(qmap,npix)
	call prepareMap(umap,npix)
	if(correlator%proc==0) then
		call ds_assert(size(qmap%indices)==size(originalIndices),"map size not equal to hit pixels size")
		qmap%indices=originalIndices
		umap%indices=originalIndices
		deallocate(originalIndices)
	endif


call makeNaiveMap(timestream,qmap,umap,pointing,correlator)
if (correlator%proc==0) call ds_log("Made naive map",ds_feedback_quiet)

if (correlator%proc == 0) then
	call savePolMaps(outputName, qmap,umap,nside)
endif

        deallocate(whiteNoise)
        deallocate(fknee)
        deallocate(rho)
        deallocate(filter)
!        do d= 0, correlator%my_ndet
!           call destroyTimestream(timestream(d))
!           call destroyPointing(pointing(d))
!        enddo
        deallocate(timestream)
        deallocate(pointing)
        deallocate(correlator%owner)
        deallocate(correlator%my_det)

	enddo

	call MPI_Finalize(ierror)


	contains

	subroutine read_params
		use inifile
		integer param_lun
		logical param_status
		type(TIniFile) :: param_file
		logical, parameter :: FAIL_IF_NO_PARAM = .true.
		logical, parameter :: OKAY_IF_NO_PARAM = .false.
		
	
		param_lun = ds_get_lun()
		call Ini_Open_File(param_file, parameter_filename, param_lun, param_status, .false.)
		call ds_assert(.not. param_status, "Could not open parameter file " // parameter_filename)
		todroot = trim( Ini_Read_String_File(param_file, 'tod_filename', FAIL_IF_NO_PARAM) )
		filtered_tod_name = trim( Ini_Read_String_File(param_file, 'filtered_todroot', OKAY_IF_NO_PARAM) )
		noiseFile = trim( Ini_Read_String_File(param_file, 'noise_filename', FAIL_IF_NO_PARAM) )
		ndet = Ini_Read_Int_File(param_file,'number_detectors')
		pointroot = trim(ini_read_string_file(param_file, 'pointing_root', FAIL_IF_NO_PARAM) )
		pixelFile = trim( Ini_Read_String_File(param_file, 'pixel_filename', FAIL_IF_NO_PARAM) )
		outputName = trim( Ini_Read_String_File(param_file, 'filtered_map', FAIL_IF_NO_PARAM) )
		nside = Ini_Read_Int_File(param_file,'nside')
		ds_global_feedback = Ini_Read_Int_File(param_file,'feedback',2)
		if (rank==0) call ds_log("Done reading parameters",ds_feedback_quiet)
		call Ini_Close_File(param_file)

	end subroutine read_params


	subroutine applyFilter(stream,nyq,fk)!,v)
		type(ds_timestream) :: stream
		real(dp) :: nyq,fk!,v
		complex(dpc), allocatable, dimension(:) :: power, stream_f
		real(dp) dk
		integer k
		allocate(power(0:stream%nt/2))
		allocate(stream_f(0:stream%nt/2))
		!Build the power spectrum
		power(0)=cmplx(0.0_dp,0.0_dp,kind=dp)
		dk = 2.0_dp * nyq / stream%nt
		do k=1,stream%nt/2
			power(k) = dcmplx(1.0_dp / k / dk, 0.0_dp)
		enddo
		power = (1.0_dp + fk * power)
		power = 1.0_dp / sqrt(power)
		call ds_fft(stream%timestream,stream_f)
		stream_f = stream_f * power
		call ds_ifft(stream_f, stream%timestream)
		deallocate(power,stream_f)
			
	end subroutine applyFilter


end program master_filter
