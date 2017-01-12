module ds_filters
	use ds_types
	!Project specific filters
	use quiet_tod_filter_mod
	implicit none
	

	type ds_filter_options
		real(dp) :: 	DEFAULT_NU_HIGHPASS,DEFAULT_ALPHA_HIGHPASS
		real(dp) :: 	DEFAULT_NU_LOWPASS,DEFAULT_ALPHA_LOWPASS
		logical :: APPLY_HIGHPASS_FILTER, APPLY_LOWPASS_FILTER, APPLY_SCANFREQ_FILTER
		real(dp) :: SCANFREQ_WIDTH
		integer :: SCANFREQ_HARMONICS
		real(dp) :: SAMPLING_RATE
	
	end type ds_filter_options

	type(ds_filter_options) :: default_filters
	character(256), private :: message
	real(dp), parameter :: QUIET_SAMPLING_RATE = 50.0

contains

subroutine read_filter_options(parfile,report,filter_options)
		use quiet_utils
		character(*) :: parfile
		type(ds_filter_options), optional :: filter_options
		type(ds_filter_options) :: opt
		integer unit
		logical, optional :: report
		logical :: do_report
		
		
		do_report=.true.
		if (present(report)) do_report=report
		unit=ds_get_lun()
		
		opt%SAMPLING_RATE=QUIET_SAMPLING_RATE
		
		call get_parameter(unit, parfile, 'NU_HIGHPASS',                   par_dp=opt%DEFAULT_NU_HIGHPASS)    
		call get_parameter(unit, parfile, 'ALPHA_HIGHPASS',                par_dp=opt%DEFAULT_ALPHA_HIGHPASS)    
		call get_parameter(unit, parfile, 'NU_LOWPASS',                    par_dp=opt%DEFAULT_NU_LOWPASS)    
		call get_parameter(unit, parfile, 'ALPHA_LOWPASS',                 par_dp=opt%DEFAULT_ALPHA_LOWPASS)    
		call get_parameter(unit, parfile, 'APPLY_HIGHPASS_FILTER',         par_lgt=opt%APPLY_HIGHPASS_FILTER)    
		call get_parameter(unit, parfile, 'APPLY_LOWPASS_FILTER',          par_lgt=opt%APPLY_LOWPASS_FILTER)    
		call get_parameter(unit, parfile, 'APPLY_SCANFREQ_FILTER',         par_lgt=opt%APPLY_SCANFREQ_FILTER)    			
		call get_parameter(unit, parfile, 'SCANFREQ_FILTER_WIDTH',         par_dp=opt%SCANFREQ_WIDTH)    
		call get_parameter(unit, parfile, 'SCANFREQ_FILTER_NUM_HARMONICS', par_int=opt%SCANFREQ_HARMONICS)
		if (do_report) then
			call ds_log("Filtering has been activated.",ds_feedback_quiet)
			if (opt%APPLY_LOWPASS_FILTER) then
				write(message,'(A,F7.3,A,F8.3)') "Low pass filter: nu=",opt%DEFAULT_NU_LOWPASS," alpha=",opt%DEFAULT_ALPHA_LOWPASS
				call ds_log(message,ds_feedback_quiet)
			else
				call ds_log("Low pass filter is OFF.",ds_feedback_quiet)
			endif
			if (opt%APPLY_HIGHPASS_FILTER) then
				write(message,'(A,F5.3,A,F6.3)') "High pass filter: nu=",opt%DEFAULT_NU_HIGHPASS," alpha=",opt%DEFAULT_ALPHA_HIGHPASS
				call ds_log(message,ds_feedback_quiet)
			else
				call ds_log("High pass filter is OFF.",ds_feedback_quiet)
			endif
			if (opt%APPLY_SCANFREQ_FILTER) then
				write(message,'(A,F8.3,A,I)') "Scan-frequency filter: width=",opt%SCANFREQ_WIDTH," harmonics=",opt%SCANFREQ_HARMONICS
				call ds_log(message,ds_feedback_quiet)
			else
				call ds_log("Scan-frequency filter is OFF.",ds_feedback_quiet)
			endif
		endif
		
		if (present(filter_options)) then
			filter_options=opt
		else
			default_filters=opt
		endif
		
end subroutine read_filter_options

	subroutine ds_build_filter(run,segment,modnum,diode, ntod, filter, scanfreq, filter_options)
		!uses descart convention for diode number
		integer ntod
		integer :: run,segment,modnum,diode
		real(dp), dimension(0:ntod/2), intent(out) :: filter
		real(dp) :: lowpass_freq,lowpass_alpha
		real(dp) :: highpass_freq,highpass_alpha
		real(dp), optional :: scanfreq
		integer status
		type(ds_filter_options), optional :: filter_options
		type(ds_filter_options) :: opt
		
		filter = 1.0d0
		
		if (present(filter_options)) then
			opt=filter_options
		else
			opt=default_filters
		endif
		

		if (opt%APPLY_HIGHPASS_FILTER) then
			call get_calibration(run=run, segment=segment, module=modnum, diode=diode-1, highpass_freq=highpass_freq, lookup_status=status)
			if (status /= 0) then
				call get_calibration(run=run, segment=segment, module=modnum, diode=diode-1,highpass_alpha=highpass_alpha)
			else
				highpass_freq  = opt%DEFAULT_NU_HIGHPASS
				highpass_alpha = opt%DEFAULT_ALPHA_HIGHPASS
			end if
			call apodize_noise_filter_fft(.true., ntod, opt%SAMPLING_RATE, highpass_freq, highpass_alpha, .true., filter)
		end if

		if (opt%APPLY_LOWPASS_FILTER) then
			call get_calibration(run=run, segment=segment, module=modnum, diode=diode-1, lowpass_freq=lowpass_freq, lookup_status=status)
			if (status /= 0) then
				call get_calibration(run=run, segment=segment, module=modnum, diode=diode-1,lowpass_alpha=lowpass_alpha)
			else
				lowpass_freq  = opt%DEFAULT_NU_HIGHPASS
				lowpass_alpha = opt%DEFAULT_ALPHA_HIGHPASS
			end if
			call apodize_noise_filter_fft(.true., ntod, opt%SAMPLING_RATE, lowpass_freq, lowpass_alpha, .true., filter)
		end if
		
		if (opt%APPLY_SCANFREQ_FILTER) then
			call apply_scanfreq_filter_fft(.true., ntod, opt%SAMPLING_RATE, scanfreq, opt%SCANFREQ_WIDTH,opt%SCANFREQ_HARMONICS, filter)
		endif
!		apply_scanfreq_filter_fft
	end subroutine ds_build_filter


subroutine applyTod2mapFilter(run,segment,modnum,moduleScan, filter_options)
	use quiet_noise_filter_mod
	use quiet_noise_mod

	integer run,segment,modnum
	integer ntod
	integer diode
	real(dp), allocatable, dimension(:) :: filter, tod
	type(ds_moduleScan) :: moduleScan
	type(ds_filter_options), optional :: filter_options
	type(ds_filter_options)  :: opt

	if (present(filter_options)) then
		opt=filter_options
	else
		opt=default_filters
	endif


	ntod = size(moduleScan%pointing)
	allocate(filter(0:ntod/2))
	allocate(tod(1:ntod))

	do diode=1,ndiodes_max
		if (.not. moduleScan%diode_flag(diode)==1) cycle
		if (opt%APPLY_SCANFREQ_FILTER) then
			call ds_build_filter(run,segment,modnum,diode, ntod, filter,moduleScan%scanfreq,filter_options=opt)
		else
			call ds_build_filter(run,segment,modnum,diode, ntod, filter,filter_options=opt)
		endif
		tod = moduleScan%timestream(diode)%timestream(1:ntod)
		call apply_noise_filter_fft(ntod, filter_options%SAMPLING_RATE, filter, tod)
		moduleScan%timestream(diode)%timestream(1:ntod) = tod
	enddo

	deallocate(tod)
	deallocate(filter)

end subroutine applyTod2mapFilter



end module ds_filters