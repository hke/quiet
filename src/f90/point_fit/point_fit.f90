program point_fit
  use point_fit_mod
  use point_fit_params
  use quiet_pointing_mod
  use quiet_mpi_mod
  implicit none
  character(len=512)    :: arg, parfile, pointfile
  type(point_data)      :: data
  logical(lgt)          :: only_apply_mount
  integer(i4b)          :: mod, npar, err
  real(dp),     allocatable :: params(:), samps(:,:), step(:)

  call getarg(1, parfile)
  call getarg(2, pointfile)
  call read_point_data(data, pointfile)
 
  call get_parameter(0, parfile, 'ONLY_APPLY_MOUNT_MODEL', par_lgt=only_apply_mount)

  call init_point_fit_mod(parfile)
  call init_point_fit_params(parfile)

  ! Set up the parameter space
  npar = get_npar()
  allocate(params(npar))!, samps(npar,100000), step(npar))

!  params = [ &
!   -6.231061688296636E-004, -5.131266559350742E-004,  &
!    7.11943638418153-real(2*pi), -1.258521721056886E-004, -6.423962615057877E-005, &
!    5.146644510698710E-003,  6.74885924754226-real(2*pi), 0.212620855079277,      &
!    4.43187037518670,        0.203628387523333,      -3.941870931406408E-004, &
!    -8.037428355525987E-005 ]
!  step   = [10.0, 3.0, 16000.0, 6.0, 8.0, 5.0, 600.0, 800.0, 4000.0, 10000.0, 10.0, 5.0]*0.000002
!
!  call get_samples_mc(data, params, step, samps)
!  stop
  
  if (.not. only_apply_mount) then

     ! Output fit with all data points
     params = get_best_multi(data, 4)
     call plot_fit(data, params,   .false., "after.txt")
     call plot_fit(data, params*0, .false., "before.txt")

     write(*,*) 'params = ', params*RAD2DEG
     
     write(*,*) 'All data points:'
     call print_header
     call print_params(params)
     
     ! Cut bad data
     call clean_data(data, 0.d0*params)
     
     ! Fit with smaller data set
     params = get_best_multi(data)
     call plot_fit(data, params,   .false., "after_clean.txt")

     write(*,*) 'Good data points:'
     call print_header
     call print_params(params)
     
  end if

  ! Output fit with default mount model
  call plot_fit(data, params,   .true.,  "after_mount.txt")

  deallocate(params, t_horns, saved)

end program
