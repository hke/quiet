! Analyse all the diodes together, looking for a point source.
! In some cases the source is strong enough to determine properties
! for each single diode. In these cases, the free parameters
! are diode amplitudes, fwhms and source position. In other
! cases, the signal is too weak, and the diode properties need
! be fixed at fiducial values, while the free parameters are
! source amplitude and source position.

! TODO: Clean up the map stuff. It is horribly messy right now.

program point_scan
  use quiet_task_mod
  use quiet_detector_mod
  use quiet_pointing_mod
  use quiet_pixspace_mod
  use quiet_lx_mod
  use quiet_ces_mod
  use powell_mod
  use quiet_mpi_utils
  use quiet_patch_detect_mod
  use quiet_shared_output_mod
use rngmod
  implicit none

  type diode_sample_data
     integer(i4b), allocatable :: inds(:)
     real(dp),     allocatable :: sim(:)
  end type

  type horn_sample_data
     integer(i4b), allocatable :: inds(:)
  end type

  type sample_data
     type(diode_sample_data), allocatable :: diodes(:)
     type(horn_sample_data),  allocatable :: horns(:)
  end type

  type diode_info
     real(dp), allocatable :: amps(:)     ! (ngroup)
     real(dp), allocatable :: stokes(:,:) ! (nstok=3,ngroup)
  end type

  type validator
     real(dp)                  :: sigma_sig_tot, sigma_null_tot, sn_tot
     real(dp),     allocatable :: sigma_sig(:), sigma_null(:), sn(:)
     logical(lgt)              :: accepted_tot
     logical(lgt), allocatable :: accepted(:)
  end type

  ! Everything needed to 'paint' the source on the sky of each diode.
  ! stokes is the stokes *fraction*, where stokes(1) == 1, and
  ! stokes(2:3) are the polarization fractions. Multiply by amp
  ! to get the actual amplitude in each component
  !
  ! I am not satisfied with the way things are organized here: The variables
  ! do not have consistent roles: Sometimes amp is the true source amplitude,
  ! and diodes%amps are fit, while sometimes diodes%amps are given, and
  ! the source amplitude is fit, for example.
  type source_model
     character(len=512)            :: name
     real(dp)                      :: point(3), point0(3), amp, stokes(3), chisq, dist
     real(dp)                      :: delta_dk
     real(dp),         allocatable :: fwhm(:)       ! (nhorn)
     real(dp),         allocatable :: poss(:,:)     ! (2,nhorn)
     real(dp),         allocatable :: angs(:)       ! (ndi)
     type(diode_info), allocatable :: diodes(:)     ! (ndi)
     logical(lgt),     allocatable :: active(:)     ! (ndi)
     integer(i4b),     allocatable :: comps(:)      ! (num_active_comps)
     integer(i4b),     allocatable :: compsIQUV(:)  ! (num_active_IQUV)

     real(dp),         allocatable :: template(:,:) ! (nmask,3)
     integer(i4b),     allocatable :: map2mask(:)   ! (0:npix-1)
     integer(i4b),     allocatable :: pixels(:)     ! (nmask)

     ! cache stuff
     real(dp),         allocatable :: nhit(:,:)     ! (2,ndi)

     integer(i4b)                  :: fit_pos, fit_amp, fit_fwhm, fit_ang, fit_stokes
     integer(i4b)                  :: fit_dk, nside, order
     logical(lgt)                  :: point_source
  end type

  real(dp)                  :: powell_dk_real(4)
  integer(i4b)              :: powell_dk_int(3)

  character(len=512)        :: parfile, odir, source_pos, source_amp, source_stokes
  character(len=512)        :: fit_what, fit_diodes, components, mapname
  character(len=512)        :: map_dir(4), tod_dir
  character(len=512)        :: dimapdir, modmapdir, source_template
  real(dp)                  :: comp_lim, rlim, nsamp_ratio, tmp, noise_tol, foo
  real(dp)                  :: rlim_raw, plot_rad, sigma_lim, null_lim, sn_lim, vlim
  integer(i4b)              :: myid, nproc, err, nces, di, cnum, nparam, debug, dfile
  integer(i4b)              :: plot_n, i, n, fit_step
  logical(lgt)              :: strong_source
  type(task_list)           :: tasks
  type(lx_struct)           :: data
  type(quiet_ces_info)      :: ces
  type(source_model)        :: source, psource
  type(sample_data)         :: samps
  type(validator)           :: val
  type(shared_ofile)        :: full_file(4)
  type(shared_ofile)        :: point_file(4)
  real(dp),     allocatable :: params(:)
  character(len=4), parameter :: otypes(4) = ["good","half","bad ","full"]

  call getarg(1, parfile)

  call initialize_ces_mod(parfile)
  call initialize_quiet_pointing_mod(parfile, apply_mount_model=.false.)
  call initialize_patch_detect_mod(parfile)
  call get_parameter(0, parfile, 'OUTPUT_DIR',      par_string=odir)
  call get_parameter(0, parfile, 'DEBUG',           par_int=debug)
  call get_parameter(0, parfile, 'SOURCE_POS',      par_string=source_pos)
  call get_parameter(0, parfile, 'SOURCE_TEMPLATE', par_string=source_template)
  call get_parameter(0, parfile, 'SOURCE_AMP',      par_string=source_amp)
  call get_parameter(0, parfile, 'SOURCE_STOKES',   par_string=source_stokes)
  call get_parameter(0, parfile, 'FIT_WHAT',        par_string=fit_what)
  call get_parameter(0, parfile, 'FIT_DIODES',      par_string=fit_diodes)
  call get_parameter(0, parfile, 'COMPONENTS',      par_string=components)
  call get_parameter(0, parfile, 'COMPONENT_THRESHOLD', par_dp=comp_lim)
  call get_parameter(0, parfile, 'FIT_RANGE',       par_dp=rlim)
  call get_parameter(0, parfile, 'FIT_PLOT_RAD',    par_dp=plot_rad)
  call get_parameter(0, parfile, 'STRONG_SOURCE',   par_lgt=strong_source)
  nsamp_ratio = 0.6
  rlim        = rlim*DEG2RAD
  rlim_raw    = rlim*4
  plot_rad    = plot_rad*DEG2RAD
  vlim        = 0.1
  plot_n      = 200
  sigma_lim   = 50d0
  null_lim    = 10d0
  sn_lim      = 5d0
  noise_tol   = 1.5d0 ! Factor we think the noise might be off by

  call mpi_init(err)
  call mpi_comm_rank(MPI_COMM_WORLD, myid,  err)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, err)
  call dset(id=myid, level=debug)

  call print_host_mapping

  map_dir   = trim(odir) // "/maps_" // otypes
  tod_dir   = trim(odir) // "/tods"
  dimapdir  = trim(odir) // "/dimaps"
  modmapdir = trim(odir) // "/modmaps"
  do i = 1, size(otypes)
     call open_shared_ofile(full_file(i), trim(odir) // "/full_"//trim(otypes(i))//".txt", MPI_COMM_WORLD)
     call open_shared_ofile(point_file(i),trim(odir) // "/point_"//trim(otypes(i))//".txt",MPI_COMM_WORLD)
     call mkdirs(map_dir(i), .false.)
  end do
  call mkdirs(dimapdir, .false.)
  call mkdirs(modmapdir, .false.)

  nces = get_num_ces()
  call init_task_list(tasks, trim(odir)//"/lock.dat", nces, MPI_COMM_WORLD)
  do while(get_next_task(tasks, cnum))
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " processing ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(nces)) // ")"
     fit_step = 1
     call dmem("ces start")
     call read_l3_file(ces%l3file, data);              call dmem("l3 read")
     call setup_model(source, data);                   call dmem("setup model")
     call setup_psource(psource)

     call filter_tod(data, source)
     if(source%point_source) then
        call calc_point(data, source, coord_tele, coord_gal); call dmem("calc point")
        if(strong_source) then
           call get_relevant_samples_point(source, data, rlim_raw, samps)
           call find_start_point(source, data, samps)
        end if
        call get_relevant_samples_point(source, data, rlim, samps)
     else
        call calc_point(data, source, coord_tele, coord_gal); call dmem("calc point")
        call get_relevant_samples_template(source, data, rlim, vlim, samps)
     end if
     !call measure_noise(source, data, rlim)

     call prune_unhit(source, samps);                  call dmem("prune")
     if(all(.not. source%active)) goto 1 ! nothing to do

     ! Powell needs a flattened model. This is expanded
     ! back into source inside calc_lnl
     nparam = model_size(source)
     allocate(params(nparam))
     call transform_model(source, params, nparam, 1);  call dmem("transform")
     call powell(params, calc_lnl, err);               call dmem("solve")
     call validate_model(source, data, samps, val)
     call output_model(ces%cid, source, data, val);    call dmem("output")
     call output_maps (ces%cid, source, data, samps, val)

!call get_relevant_samples_point(source, data, 1d5, samps)
!call simulate_tod(source, data, samps, foo)
!call output_tod  (ces%cid, data, samps)

     deallocate(params)
1    call free_lx_struct(data)
     call free_sample_data(samps)
     call free_source_model(source)
     call free_source_model(psource)
  end do
  do i = 1, size(otypes)
     call close_shared_ofile(full_file(i))
     call close_shared_ofile(point_file(i))
  end do
  call mpi_finalize(err)
  write(*,*) myid, "finished"
contains

  ! Uses globals source_pos, source_amp, source_stokes, components,
  ! fit_what, fit_diodes, comp_lim
  subroutine setup_model(source, data)
    implicit none
    type(source_model), intent(inout) :: source
    type(lx_struct),    intent(in)    :: data
    character(len=64),  allocatable   :: toks(:)
    logical(lgt)                      :: comp_active(3)
    integer(i4b)       :: ndi, nmod, di, mod, i, n, nh
    ndi  = size(quiet_diodes)
    nmod = size(quiet_horns)
    call parse_source_point(source_pos, data%time(1), source%point0, source%dist, source%name)
    source%point_source = .true.
    if(source_template /= "" .and. source_template /= "none") then
       call read_map(source%template, source%pixels, source%nside, source%order, &
        & source_template)
       allocate(source%map2mask(0:12*source%nside**2-1))
       source%map2mask = 0
       do i = 1, size(source%pixels)
          source%map2mask(source%pixels(i)) = i
       end do
       source%point_source = .false.
    end if

    source%point = source%point0
    call parse_source_stokes(source_stokes, source%stokes)
    read(source_amp,*) source%amp

    allocate(source%diodes(ndi),  source%fwhm(nmod), source%active(ndi))
    allocate(source%poss(2,nmod), source%angs(ndi))
    ! activate the relevant diodes
    source%active    = .false.
    n = num_tokens(fit_diodes, ", ")
    allocate(toks(n))
    call get_tokens(fit_diodes, ", ", toks)
    call assert(n > 0, "No diodes to fit for!")
    comp_active(1) = index(components,"T") > 0
    comp_active(2) = index(components,"Q") > 0 .or. index(components,"U") > 0
    comp_active(3) = index(components,"V") > 0
    call wherei(comp_active, source%comps)
    call wherei(comp_active([1,2,2,3]), source%compsIQUV)

    if(toks(1) == "*") then
       ! Choose as many diodes as possible
       do i = 1, ndi
          if(any(comp_active .and. abs(quiet_diodes(i)%stokes) > comp_lim)) source%active(i) = .true.
       end do
       source%active = source%active .and. quiet_diodes%ok
    else
       do i = 1, n
          read(toks(i),*) di
          call assert(di > 0 .and. di <= size(quiet_diodes), "Invalid diode " // trim(itoa(di)))
          source%active(di) = .true.
       end do
    end if
    deallocate(toks)

    ! What to fit for?
    source%fit_pos = 0; source%fit_amp    = 0; source%fit_fwhm = 0
    source%fit_ang = 0; source%fit_stokes = 0; source%fit_dk   = 0
    n = num_tokens(fit_what, ", ")
    allocate(toks(n))
    call get_tokens(fit_what, ", ", toks)
    if(any(toks == "pos"))    source%fit_pos    = 1
    if(any(toks == "poss"))   source%fit_pos    = 2
    if(any(toks == "dk"))     source%fit_dk     = 1
    if(any(toks == "dk2"))    source%fit_dk     = 2
    if(any(toks == "ang"))    source%fit_ang    = 1
    if(any(toks == "angs"))   source%fit_ang    = 2
    if(any(toks == "amp"))    source%fit_amp    = 1
    if(any(toks == "amps"))   source%fit_amp    = 2
    if(any(toks == "fwhm"))   source%fit_fwhm   = 1
    if(any(toks == "stokes")) source%fit_stokes = 1
    if(any(toks == "stokess"))source%fit_amp    = 3
    deallocate(toks)

    source%fwhm      = quiet_horns%fwhm
    source%angs      = quiet_diodes%psi
    source%delta_dk  = 0
    source%poss      = spread(source%point(1:2),2,nmod)
    do di = 1, ndi
       mod = quiet_diodes(di)%horn
       nh  = size(quiet_horns(mod)%amps)
       allocate(source%diodes(di)%amps(nh))
       allocate(source%diodes(di)%stokes(3,nh))
       source%diodes(di)%amps = quiet_horns(mod)%amps * data%gain(1, di)
       source%diodes(di)%stokes = spread(quiet_diodes(di)%stokes,2,nh)
    end do
    allocate(source%nhit(2,ndi))
  end subroutine

  subroutine setup_psource(source)
    implicit none
    type(source_model), intent(inout) :: source
    integer(i4b)       :: ndi, nmod, di, mod, i, n, nh
    ndi  = size(quiet_diodes)
    nmod = size(quiet_horns)
    allocate(source%fwhm(nmod))
    allocate(source%poss(2,nmod), source%angs(ndi))
    allocate(source%diodes(size(quiet_diodes)))
    allocate(source%nhit(2,ndi))
    source%point     = NaN
    source%amp       = NaN
    source%delta_dk  = NaN
    source%stokes    = NaN
    source%fwhm      = NaN
    source%angs      = NaN
    source%poss      = NaN
    do di = 1, ndi
       mod = quiet_diodes(di)%horn
       nh  = size(quiet_horns(mod)%amps)
       allocate(source%diodes(di)%amps(nh))
       allocate(source%diodes(di)%stokes(3,nh))
       source%diodes(di)%amps   = NaN
       source%diodes(di)%stokes = NaN
    end do
  end subroutine

  subroutine find_start_point(source, data, samps)
    implicit none
    type(source_model), intent(inout) :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(in)    :: samps
    source%point(1:2) = pos_from_map(source, data, samps)
    source%poss       = spread(source%point(1:2),2,size(quiet_horns))
  end subroutine

  ! Measure the overall noise level of the samples *not* hit.
  ! We have to measure the noise because l3gen doesn't do this correctly
  ! for strong sources. Also, we do not filter away the 1/f noise here,
  ! so using only the white noise floor is inconsistent. This function
  ! has a large overlap with get_relevant_samples.
  subroutine measure_noise(source, data, rlim)
    implicit none
    type(source_model), intent(in)    :: source
    type(lx_struct),    intent(inout) :: data
    real(dp)                          :: rlim, r, sigma
    integer(i4b)                      :: di, ndi, horn, nh, mod, nsamp, i
    integer(i4b),       allocatable   :: inds(:)
    logical(lgt),       allocatable   :: mask(:)
    ndi   = size(quiet_diodes)
    nsamp = size(data%point,2)
    allocate(mask(nsamp))
    do di = 1, ndi
       if(.not. source%active(di)) cycle
       nh   = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       mask = .false.
       do horn = 1, nh
          mod = quiet_horns(quiet_diodes(di)%horn)%groups(horn)
          do i = 1, nsamp
             r = polangdist(real(data%point([2,1],i,mod),dp), source%point([2,1]))
             mask(i) = mask(i) .or. r < rlim
          end do
       end do
       call wherei(.not. mask, inds)
       sigma = sqrt(sum(data%tod(inds,di)**2)/size(inds))
       data%sigma0(di) = sigma
       deallocate(inds)
    end do
    deallocate(mask)
  end subroutine

  function calc_lnl(p)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    integer(i4b),                    allocatable :: diodes(:)
    real(dp)                                     :: calc_lnl, chisq, mval
    integer(i4b)                                 :: i, ndof
    character(len=512)                           :: fname
    calc_lnl = 1d100
    call transform_model(source, p, i, -1);            !call dmem("transform2")
    if(source%point(2) < 0 .or. source%point(2) > pi) return
    if(any(source%poss(2,:) < 0 .or. source%poss(2,:) > pi)) return
    if(any(source%fwhm <= 0)) return
    if(any(source%stokes(2:2) < 0)) return
    if(any(source%fwhm > 300./60*pi/180)) return
    call simulate_tod(source, data, samps, mval);      !call dmem("simulate")
    if(mval < 0.5 .or. mval /= mval) return ! Require significant overlap
    call calc_chisq(source, data, samps, chisq, ndof); !call dmem("calc chisq")
    source%chisq = chisq/ndof
    calc_lnl = source%chisq
    if(mod(fit_step-1,10)==0) then
       if(dtest(2)) write(*,*) fit_step-1, source%chisq, p
       call wherei(source%active, diodes)
       if(dtest(10)) then
          fname = trim(odir) // "/dumps/ces" // trim(itoa(ces%cid,4)) // "/" // trim(itoa(fit_step-1,5))//".dat"
          call mkdirs(fname, .true.)
          call dump_map(source, data, samps, fname, diodes)
       end if
       if(dtest(11)) then
          do i = 1, size(diodes)
             fname = trim(odir) // "/didumps/ces" // trim(itoa(ces%cid,4)) // &
                  & "/di" // trim(itoa(diodes(i),3)) // "/" // trim(itoa(fit_step-1,5))//".dat"
             call mkdirs(fname, .true.)
             call dump_map(source, data, samps, fname, diodes(i:i))
          end do
       end if
       deallocate(diodes)
    end if
    call copy_psource(source, psource);                !call dmem("copy")
    fit_step = fit_step+1
  end function

  ! Find the chisquare and number of degrees of freedom for the
  ! simulated samples. We use a white noise model here, since
  ! we do not have continuous samples.
  subroutine calc_chisq(source, data, samps, chisq, ndof)
    implicit none
    type(source_model),intent(in)  :: source
    type(lx_struct),   intent(in)  :: data
    type(sample_data), intent(in)  :: samps
    real(dp),          intent(out) :: chisq
    integer(i4b),      intent(out) :: ndof
    integer(i4b)                   :: di
    ndof  = 0
    chisq = 0
    do di = 1, size(samps%diodes)
       if(.not. source%active(di)) cycle
       ndof  = ndof  + size(samps%diodes(di)%inds)
       chisq = chisq + sum((samps%diodes(di)%sim - &
        & data%tod(samps%diodes(di)%inds,di))**2) / data%sigma0(di)**2
    end do
  end subroutine

  subroutine output_model(cid, source, data, val)
    implicit none
    integer(i4b),       intent(in) :: cid
    type(source_model), intent(in) :: source
    type(lx_struct),    intent(in) :: data
    type(validator),    intent(in) :: val
    ! Output in four classes:
    !  Fully ok, individually ok, bad and everything
    ! We do this for full and point
    call output_full(cid, source, data, val, &
     & val%accepted .and. val%accepted_tot, full_file(1))
    call output_full(cid, source, data, val, &
     & val%accepted .and. .not. val%accepted_tot, full_file(2))
    call output_full(cid, source, data, val, .not. val%accepted, full_file(3))
    call output_full(cid, source, data, val, val%accepted .or. .not. val%accepted, full_file(4))
    call output_point(cid, source, data, val, &
     & val%accepted .and. val%accepted_tot, point_file(1))
    call output_point(cid, source, data, val, &
     & val%accepted .and. .not. val%accepted_tot, point_file(2))
    call output_point(cid, source, data, val, .not. val%accepted, point_file(3))
    call output_point(cid, source, data, val, val%accepted .or. .not. val%accepted, point_file(4))
  end subroutine

  subroutine output_full(cid, source, data, val, mask, fd)
    implicit none
    integer(i4b),       intent(in) :: cid
    type(source_model), intent(in) :: source
    type(lx_struct),    intent(in) :: data
    type(validator),    intent(in) :: val
    logical(lgt),       intent(in) :: mask(:)
    type(shared_ofile)             :: fd
    character(len=1024)            :: line
    real(dp)                       :: p0(3), p(3), psi, sn_ratio, T0, gain
    integer(i4b)                   :: mod, di, horn, nh, h

    2 format(i5,i4,i3,f14.7,6f10.4,5e15.7,e15.7,f10.4,6e15.7,f10.4)
    do di = 1, size(quiet_diodes)
       if(.not. (source%active(di) .and. mask(di))) cycle
       mod  = quiet_diodes(di)%horn
       nh   = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       do horn = 1, nh
          h  = quiet_horns(mod)%groups(horn)
          p0 = source%point0; p(1:2) = source%poss(:,mod+1); p(3) = source%angs(di)
          T0 = get_observed_object_temperature(source%name, source%fwhm(h+1)*RAD2DEG*60, data%time(1), quiet_diodes(di)%freq)
          if(T0 == 0) T0 = 1
          ! amps is in uK/K(ant), amp is K(thermo). T0 is in K(ant),
          ! must convert source%amp to ant and divide by 1000 to get milli
          ! from micro.
          gain = 1d-3/ant2thermo(get_diode_freq(di))*source%diodes(di)%amps(horn)*source%amp/T0
          call swap_coordinate_convention(p(1),p(2),p(3),coord_gal)
          call swap_coordinate_convention(p0(1),p0(2),p0(3),coord_gal)
          write(line,2) ces%cid, di, horn, data%time(1), p0*RAD2DEG, &
           & p*RAD2DEG, source%amp, source%diodes(di)%amps(horn), &
           & source%diodes(di)%amps(horn)/source%diodes(di)%amps(1), source%dist, &
           & T0, gain, &
           & source%fwhm(h+1)*RAD2DEG, val%sigma_null_tot, val%sigma_sig_tot, &
           & val%sn_tot, val%sigma_null(di), val%sigma_sig(di), val%sn(di), &
           & source%delta_dk*RAD2DEG
          call write_shared_ofile(fd, trim(line))
       end do
    end do
    call flush_shared_ofile(fd)
  end subroutine

  ! Old format was: mjd az_tobs el_tobs az el dk obj mod cid
  ! This format depended on the focalplane layout, which means that
  ! point fit needed to be rerun every time we changed this.
  ! To avoid this, we will instead report the actual horizontal
  ! coordinates each horn is pointing in when observing the source.
  subroutine output_point(cid, source, data, val, mask, fd)
    implicit none
    integer(i4b),       intent(in) :: cid
    type(source_model), intent(in) :: source
    type(lx_struct),    intent(in) :: data
    type(validator),    intent(in) :: val
    logical(lgt),       intent(in) :: mask(:)
    type(shared_ofile)             :: fd
    character(len=1024)            :: line
    real(dp)                       :: p0(3), p(3), psi, sn_ratio
    integer(i4b)                   :: mod, di, horn, nh, h
    1 format(f14.7, 5f10.4, a10, i3, i8)
    call coord_convert(coord_gal, source%point0(1), source%point0(2), source%point0(3), &
         & coord_hor, p0(1), p0(2), p0(3), mjd=data%time(1))
    call swap_coordinate_convention(p0(1), p0(2), p0(3), coord_hor)
    p(3) = data%orig_point(3,1)
    if(.false. .and. source%fit_pos == 1) then
       if(any(mask)) then
          call cc_invdk(coord_gal, source%point(1), source%point(2), psi, coord_tele, &
               & p(1), p(2), p(3), data%time(1), -1)
          call swap_coordinate_convention(p(1), p(2), p(3), coord_tele)
          write(line,1) data%time(1), p0(1:2)*RAD2DEG, p*RAD2DEG, trim(ces%object), -1, cid
          call write_shared_ofile(fd, trim(line))
       end if
    elseif(.true. .or. source%fit_pos == 2) then
       do mod = 0, size(quiet_horns)-1
          if(.not. any(quiet_diodes%horn == mod .and. source%active .and. mask)) cycle
          p(3) = data%orig_point(3,1)
          call cc_invdk(coord_gal, source%poss(1,mod+1), source%poss(2,mod+1), psi, &
               & coord_tele, p(1), p(2), p(3), data%time(1), mod)
          call swap_coordinate_convention(p(1), p(2), p(3), coord_tele)
          write(line,1) data%time(1), p0(1:2)*RAD2DEG, p*RAD2DEG, trim(ces%object), mod, cid
          call write_shared_ofile(fd, trim(line))
       end do
    end if
    call flush_shared_ofile(fd)
  end subroutine

  subroutine output_maps(cid, source, data, samps, val)
    implicit none
    integer(i4b),       intent(in) :: cid
    type(source_model), intent(in) :: source
    type(lx_struct),    intent(in) :: data
    type(sample_data),  intent(in) :: samps
    type(validator),    intent(in) :: val
    integer(i4b),       allocatable:: diodes(:)
    call wherei(source%active .and. val%accepted .and. val%accepted_tot, diodes)
    mapname = trim(map_dir(1)) // "/map" // trim(itoa(cid,4)) // ".dat"
    call dump_map(source, data, samps, mapname, diodes)

    call wherei(source%active .and. val%accepted .and. .not. val%accepted_tot, diodes)
    mapname = trim(map_dir(2)) // "/map" // trim(itoa(cid,4)) // ".dat"
    call dump_map(source, data, samps, mapname, diodes)

    call wherei(source%active .and. .not. val%accepted, diodes)
    mapname = trim(map_dir(3)) // "/map" // trim(itoa(cid,4)) // ".dat"
    call dump_map(source, data, samps, mapname, diodes)

    call wherei(source%active, diodes)
    mapname = trim(map_dir(4)) // "/map" // trim(itoa(cid,4)) // ".dat"
    call dump_map(source, data, samps, mapname, diodes)

    call dump_dimaps (source, data, samps, dimapdir,  cid)
    call dump_modmaps(source, data, samps, modmapdir, cid)
    deallocate(diodes)
  end subroutine

  ! Build a binned map of both data and simulation. Could use some
  ! refactoring. Only include the components we are interested in.
  subroutine make_map(source, data, samps, width, height, res, diodes)
    implicit none
    type(source_model),intent(in)  :: source
    type(lx_struct),   intent(in)  :: data
    type(sample_data), intent(in)  :: samps
    real(dp),          intent(out) :: res(:,:,:,:)
    integer(i4b),      intent(in)  :: diodes(:)
    real(dp),          allocatable :: rhs(:,:,:,:), div(:,:,:,:), stupid(:)
    integer(i4b),      allocatable :: inds(:)
    integer(i4b)                   :: n1,n2,di,idi,ndi, mod, i, j, h, x, y, status, ncomp
    integer(i4b)                   :: nh, horn
    real(dp)                       :: theta0, phi0, psi, phase(4,1), w, gain, lon, lat
    real(dp)                       :: width, height, mask(4)
    mask = 0; mask(source%compsIQUV) = 1
    n1 = size(res,2); n2 = size(res,3);
    ncomp = 4; ndi = size(quiet_diodes)
    theta0 = source%point0(2); phi0 = source%point0(1)
    ! Build equation set
    allocate(rhs(ncomp,n1,n2,2), div(ncomp,ncomp,n1,n2))
    rhs = 0; div = 0
    do idi = 1, size(diodes)
       di = diodes(idi)
       nh   = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       do horn = 1, 1!nh
          mod   = quiet_horns(quiet_diodes(di)%horn)%groups(horn)
          do i = 1, size(samps%diodes(di)%inds)
             j = samps%diodes(di)%inds(i)
             y = floor((data%point(2,j,mod)-theta0)/height*n1/2) + n1/2
             x = floor((ang_diff(real(data%point(1,j,mod),dp),phi0))/width*n2/2)+n2/2
             if(x < 1 .or. x > n2 .or. y < 1 .or. y > n1) cycle
             psi   = data%point(3,j,mod) + source%angs(di)
             gain  = 1d-6 / ant2thermo(get_diode_freq(di))
             phase(:,1) = [ 1d0, cos(2*psi), sin(2*psi), 1d0 ] * &
                  & source%diodes(di)%stokes([1,2,2,3],horn) * mask * &
                  & source%diodes(di)%amps(horn) * gain
             if(data%sigma0(di) == 0 .or. data%sigma0(di) /= data%sigma0(di)) cycle
             w     = data%sigma0(di)**(-2)
             rhs(:,y,x,1) = rhs(:,y,x,1) + data%tod(j,di)          * phase(:,1) * w
             rhs(:,y,x,2) = rhs(:,y,x,2) + samps%diodes(di)%sim(i) * phase(:,1) * w
             div(:,:,y,x) = div(:,:,y,x) + matmul(phase,transpose(phase)) * w
          end do
       end do
    end do

    ! Solve
    res = 0.d0 !NaN
    do x = 1, n2
       do y = 1, n1
          call wherei(get_diag(div(:,:,y,x)) /= 0, inds)
          if(size(inds) > 0) then
             allocate(stupid(size(inds)))
             call solve_system_eiglim(div(inds,inds,y,x), rhs(inds,y,x,1), stupid, 1d-12,status)
             res(inds,y,x,1) = stupid
             call solve_system_eiglim(div(inds,inds,y,x), rhs(inds,y,x,2), stupid, 1d-12,status)
             res(inds,y,x,2) = stupid
             deallocate(stupid)
          end if
          deallocate(inds)
       end do
    end do
    deallocate(rhs,div)
  end subroutine

  subroutine dump_map(source, data, samps, file, diodes)
    implicit none
    type(source_model),intent(in)  :: source
    type(lx_struct),   intent(in)  :: data
    type(sample_data), intent(in)  :: samps
    character(len=*),  intent(in)  :: file
    integer(i4b),      intent(in)  :: diodes(:)
    real(dp),          allocatable :: res(:,:,:,:)
    integer(i4b)                   :: n, di, ndi, mod, i, j, h, x, y, status, ncomp
    integer(i4b)                   :: nh, horn, unit
    real(dp)                       :: theta0, phi0, rad, psi, phase(4,1), w, gain, lon,lat
    logical(lgt)                   :: ok
    if(size(diodes) == 0) return

    n = plot_n; ncomp = 4; rad = plot_rad
    allocate(res(ncomp,n,n,2))
    call make_map(source, data, samps, rad, rad, res, diodes)

    where(abs(res) < 1d-99) res = 0 ! stupid ifort

    ! And output in row-major order
    unit = getlun()
    open(unit,file=file)
    do y = 1, n
       do x = 1, n
          lon =  (x-n/2)*rad/(n/2) * RAD2DEG
          lat =  (y-n/2)*rad/(n/2) * RAD2DEG
          write(unit,'(2e15.7)',advance="no") lon, lat
          do i = 1, ncomp
             write(unit,'(2e15.7)',advance="no") res(i,y,x,1:2)
          end do
          write(unit,*)
       end do
       write(unit,*)
    end do
    close(unit)

    open(unit,file=trim(file)//'.true')
    do i = 1, size(source%poss,2)
       ok = .false.
       do j = 1, size(diodes)
          horn = quiet_diodes(diodes(j))%horn
          if(any(quiet_horns(horn)%groups == i-1)) ok = .true.
       end do
       if(ok) write(unit,*) ((source%poss(:,i)-source%point0(1:2))) * RAD2DEG
    end do
    close(unit)

  end subroutine

  subroutine dump_dimaps(source, data, samps, dir, cid)
    implicit none
    type(source_model),intent(in)  :: source
    type(lx_struct),   intent(in)  :: data
    type(sample_data), intent(in)  :: samps
    character(len=*),  intent(in)  :: dir
    integer(i4b),      intent(in)  :: cid
    character(len=512)             :: mapname
    integer(i4b)                   :: di
    do di = 1, size(source%active)
       if(.not. source%active(di)) cycle
       mapname = trim(dir) // "/map" // trim(itoa(cid,4)) // "_" // &
        & trim(itoa(di,3)) // ".dat"
       call dump_map(source, data, samps, mapname, [di])
    end do
  end subroutine

  subroutine dump_modmaps(source, data, samps, dir, cid)
    implicit none
    type(source_model),intent(in)  :: source
    type(lx_struct),   intent(in)  :: data
    type(sample_data), intent(in)  :: samps
    character(len=*),  intent(in)  :: dir
    integer(i4b),      intent(in)  :: cid
    character(len=512)             :: mapname
    integer(i4b)                   :: mod
    integer(i4b),      allocatable :: dis(:)
    do mod = 0, size(quiet_horns)-1
       call wherei(quiet_diodes%horn == mod .and. source%active, dis)
       if(size(dis) > 0) then
          mapname = trim(dir) // "/map" // trim(itoa(cid,4)) // "_" // &
           & trim(itoa(mod,3)) // ".dat"
          call dump_map(source, data, samps, mapname, dis)
       end if
       deallocate(dis)
    end do
  end subroutine

  subroutine validate_model(source, data, samps, val)
    implicit none
    type(source_model), intent(inout) :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(in)    :: samps
    type(validator),    intent(inout) :: val
    integer(i4b)                      :: di, i, ndi, nsamp, ntot, nok
    real(dp)                          :: sigsum, nullsum, mean, var, sigma

    call free_validator(val)
    val%accepted_tot = .true.
    ndi              = size(quiet_diodes)
    allocate(val%sigma_sig(ndi),val%sigma_null(ndi),val%accepted(ndi),val%sn(ndi))
    val%accepted     = source%active

    ! For each active diode, calculate the significance of detection
    ! with and without subtracting the model.
    val%sigma_sig  = 0
    val%sigma_null = 0
    val%sn         = 0
    val%sigma_sig_tot  = 0
    val%sigma_null_tot = 0
    val%sn_tot         = 0

    ntot = 0
    nok  = 0
    do di = 1, ndi
       if(.not. source%active(di)) cycle
       nsamp = size(samps%diodes(di)%inds)

       ! If sigma_null is high, it is either because of resudials or
       ! because our noise model is wrong. I often see that sigma_null
       ! is correlated with sigma_sig, which is what one would expect if
       ! sigma0 is fluctuating incorrectly.
       !
       ! We can make ourselves immune to sigma0 by instead measuring
       ! the sig/null ratio. The rationale for this is that we are
       ! willing to accept bad weather and contamination etc. as long
       ! as our signal is much stronger than it. The distribution used
       ! for sig/null is mathematically dubious here. But it work very
       ! well as a number for distinguishing bad scans from good ones
       !
       ! Conversely, this means that we will reject diodes with too
       ! low sig/null. Which means that weak sources like patch_gc
       ! must be treated differently, since they will have sig ~ null
       ! for a single diode.

       sigma   = data%sigma0(di)*noise_tol ! Noise is uncertain, so add a margin.
       nullsum = sum((data%tod(samps%diodes(di)%inds,di) - samps%diodes(di)%sim)**2)
       sigsum  = sum(data%tod(samps%diodes(di)%inds,di)**2)
       val%sigma_null(di) = (nullsum/sigma**2-nsamp)/sqrt(2d0*nsamp)
       val%sigma_sig(di)  = (sigsum/ sigma**2-nsamp)/sqrt(2d0*nsamp)

       !val%sn(di) = (sigsum/nullsum-1)/(2d0*nsamp**(-0.5))
       mean = nsamp/(nsamp-2)
       var  = 4d0*nsamp*(nsamp-1d0)/((nsamp-2d0)**2*(nsamp-4d0))
       val%sn(di) = (sigsum/nullsum-mean)/sqrt(var)

       ! We wish to combine a set of these sn values to an sn_tot. To do
       ! this, we note that sn is approximately gaussian, with mean 0 and
       ! variance 1, so the sum of N of these has variance N.

       ! For strong sources, require detection in every diode
       !if(val%sn(di) < sn_lim .and. strong) then
       if(.not. any(val%accepted)) val%accepted_tot = .false.
       if(val%sn(di) < sn_lim .or. val%sigma_null(di) > null_lim) then
          val%accepted(di) = .false.
       else
          val%sigma_sig_tot  = val%sigma_sig_tot  + val%sigma_sig (di)*sqrt(2d0*nsamp)
          val%sigma_null_tot = val%sigma_null_tot + val%sigma_null(di)*sqrt(2d0*nsamp)
          val%sn_tot         = val%sn_tot         + val%sn(di)
          ntot = ntot + nsamp
          nok  = nok  + 1
       end if
    end do
    if(ntot == 0) then
       val%accepted_tot = .false.
    else
       val%sigma_sig_tot  = val%sigma_sig_tot /sqrt(2d0*ntot)
       val%sigma_null_tot = val%sigma_null_tot/sqrt(2d0*ntot)
       val%sn_tot         = val%sn_tot/sqrt(1d0*nok)
       if(val%sigma_null_tot > sigma_lim .or. val%sn_tot < sn_lim) &
        & val%accepted_tot = .false.
    end if
  end subroutine

  ! point is the galactic pointing for each sample-diode,
  ! based on a coord_convert with no mount model.
  ! source is the initial source model. The initial source model
  ! shouldn't be too incorrect.
  subroutine get_relevant_samples_point(source, data, rlim, samps)
    implicit none
    type(source_model), intent(in)    :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(inout) :: samps
    integer(i4b)                      :: di, ndi, horn, nh, mod, nsamp, i, dii
    real(dp)                          :: r, rlim, mpoint(3)
    logical(lgt),       allocatable   :: mask(:)
    ndi   = size(quiet_diodes)
    nsamp = size(data%point,2)
    call free_sample_data(samps)
    allocate(samps%diodes(ndi), samps%horns(size(quiet_horns)), mask(nsamp))
    do di = 1, ndi
       if(.not. source%active(di)) then
          allocate(samps%diodes(di)%inds(0))
          cycle
       end if
       nh   = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       mask = .false.
       do horn = 1, nh
          mod = quiet_horns(quiet_diodes(di)%horn)%groups(horn)
          do i = 1, nsamp
             r = polangdist(real(data%point([2,1],i,mod),dp), source%point([2,1]))
             mask(i) = mask(i) .or. r < rlim
          end do
       end do
       call wherei(mask, samps%diodes(di)%inds)
       allocate(samps%diodes(di)%sim(size(samps%diodes(di)%inds)))
       samps%diodes(di)%sim = 0
    end do
    call build_horn_inds(nsamp, samps)
    deallocate(mask)
  end subroutine

  ! For each pixel in the template, calculate v = sqrt((vals(1:4)*stokes([1,2,2,3]))**2)
  ! Make a mask based on v > vmax*vlim. Accept any samples closer than rlim to these.
  subroutine get_relevant_samples_template(source, data, rlim, vlim, samps)
    implicit none
    type(source_model), intent(in)    :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(inout) :: samps
    real(dp),           intent(in)    :: rlim, vlim
    integer(i4b)                      :: ndi, nsamp, npix_full, npix
    integer(i4b)                      :: di, n, i, j, mod, horn, nh
    integer(i4b)                      :: pix, npixlist, ncomp, nside
    integer(i4b),       allocatable   :: pixlist(:), map2mask(:)
    real(dp)                          :: vals(4), vec(3), pixsize
    real(dp),           allocatable   :: v(:),     map(:,:)
    logical(lgt),       allocatable   :: vmask(:), rmask(:), tmask(:)
    type(degrade_info)                :: deg
    ndi    = size(quiet_diodes)
    nsamp  = size(data%point,2)
    ncomp  = size(source%template,2)
    call free_sample_data(samps)

    ! We will do the relevant sample analysis on a degraded map in
    ! order to both reduce noise and speed things up.
    nside = 256
    call prepare_degrade(deg, source%pixels, source%order, source%nside, nside)
    allocate(map(size(deg%opix),ncomp))
    call degrade_map(source%template, deg, map)
    npix       = size(deg%opix)
    npix_full  = 12*nside**2
    allocate(map2mask(0:npix_full-1))
    map2mask = 0
    do i = 1, npix
       map2mask(deg%opix(i)) = i
    end do
    ! Subtract median in case the map has a large offset
    do i = 1, ncomp
       map(:,i) = map(:,i) - median(map(:,i))
    end do

    ! For each active horn we must calculate a v map
    allocate(samps%diodes(ndi), samps%horns(size(quiet_horns)))
    allocate(v(npix), vmask(npix), rmask(npix), tmask(nsamp))
    pixsize  = sqrt(4*pi/12/nside**2)
    npixlist = nint(pi*(rlim+2*pixsize)**2/(4*pi)*12*nside**2)+10
    allocate(pixlist(npixlist))
    do di = 1, ndi
       if(.not. source%active(di)) then
          allocate(samps%diodes(di)%inds(0))
          cycle
       end if
       do i = 1, npix
          vals = 0
          do j = 1, ncomp
             if(healok(map(i,j))) vals(j) = map(i,j)
          end do
          v(i) = sqrt(sum((vals*quiet_diodes(di)%stokes([1,2,2,3]))**2))
       end do
       vmask = v/maxval(v) > vlim
       rmask = .false.
       ! This step is pretty heavy. Find all nearby pixels.
       do i = 1, npix
          if(.not. vmask(i)) cycle
          vec = pix2vec(nside, source%order, deg%opix(i))
          call query_disc(nside, vec, rlim, pixlist, n, &
           & nest=source%order-1, inclusive=1)
          do j = 1, n
             if(map2mask(pixlist(j)) <= 0) cycle
             rmask(map2mask(pixlist(j))) = .true.
          end do
       end do
       ! Then see which samples hit the mask
       nh   = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       tmask = .false.
       do horn = 1, nh
          mod = quiet_horns(quiet_diodes(di)%horn)%groups(horn)
          do i = 1, nsamp
             pix = ang2pix(nside, source%order, &
              & real(data%point(2,i,mod),dp), real(data%point(1,i,mod),dp))
             if(map2mask(pix) <= 0) cycle
             tmask(i) = tmask(i) .or. rmask(map2mask(pix))
          end do
       end do
       call wherei(tmask, samps%diodes(di)%inds)
       allocate(samps%diodes(di)%sim(size(samps%diodes(di)%inds)))
       samps%diodes(di)%sim = 0
    end do
    call free_degrade_info(deg)
    deallocate(vmask, rmask, tmask, pixlist, map2mask, map)
    call build_horn_inds(nsamp, samps)
  end subroutine

  subroutine build_horn_inds(nsamp, samps)
    implicit none
    type(sample_data)  :: samps
    integer(i4b)       :: mod, di, dii, nsamp
    logical(lgt), allocatable :: mask(:)
    allocate(mask(nsamp))
    do mod = 0, size(quiet_horns)-1
       mask = .false.
       do dii = 0, size(quiet_horns(mod)%diodes)-1
          di = quiet_horns(mod)%diodes(dii)
          mask(samps%diodes(di)%inds) = .true.
       end do
       call wherei(mask, samps%horns(mod+1)%inds)
    end do
    deallocate(mask)
  end subroutine

  subroutine prune_unhit(source, samps)
    implicit none
    type(source_model), intent(inout) :: source
    type(sample_data),  intent(in)    :: samps
    integer(i4b)                      :: di, nmax
    ! Horns with significantly fewer samples than the max horn
    ! will be dropped, as they make the fit much slower, and produce
    ! poor fits in the end.
    nmax = 0
    do di = 1, size(quiet_diodes)
       nmax = max(nmax,size(samps%diodes(di)%inds))
    end do
    do di = 1, size(quiet_diodes)
       if(size(samps%diodes(di)%inds) <= nmax*nsamp_ratio) source%active(di) = .false.
    end do
  end subroutine

  subroutine simulate_tod(source, data, samps, mval)
    implicit none
    type(source_model), intent(inout) :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(inout) :: samps
    real(dp),           intent(out), optional :: mval
    if(source%point_source) then
       call simulate_tod_point(source, data, samps, mval)
    else
       call simulate_tod_template(source, data, samps, mval)
    end if
  end subroutine

  ! Produce a simulated time stream in V based on the source
  ! and pointing data
  subroutine simulate_tod_point(source, data, samps, mval)
    implicit none
    type(source_model), intent(in)    :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(inout) :: samps
    real(dp),           intent(out), optional :: mval
    integer(i4b)                      :: di, nh, horn, mod, i, j
    real(dp)                          :: r, val, sigma, psi, pamp, gain, gauss
    if(present(mval)) mval = 0
    do di = 1, size(samps%diodes)
       if(.not. source%active(di)) cycle
       nh   = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       samps%diodes(di)%sim = 0
       do horn = 1, nh
          mod   = quiet_horns(quiet_diodes(di)%horn)%groups(horn)
          sigma = source%fwhm(mod+1)*fwhm2sigma
          do i = 1, size(samps%diodes(di)%inds)
             j    = samps%diodes(di)%inds(i)
             r    = polangdist(real(data%point([2,1],j,mod),dp), source%poss([2,1],mod+1))
             psi  = data%point(3,j,mod) + source%angs(di)
             pamp = cos(2*(psi-source%point(3)))
             ! source%amp is in K (thermo). Our effective gain is
             ! amps*gain, and is in units V/K(thermo). "gain" has
             ! unit 1d-6 K(ant)/K(thermo).
             ! This means that amps has unit V/K(thermo) *1d6 * K(thermo)/K(ant).
             ! = uV/K(ant)
             gain = 1d-6 / ant2thermo(get_diode_freq(di))
             gauss= exp(-0.5*(r/sigma)**2)
             if(present(mval)) mval = max(mval, gauss)
             val  = source%amp * source%diodes(di)%amps(horn) * gauss * gain
             samps%diodes(di)%sim(i) = samps%diodes(di)%sim(i) + & 
              & val * sum(source%stokes*source%diodes(di)%stokes(:,horn)*[1d0,pamp,1d0])
          end do
       end do
    end do
  end subroutine

  ! source%poss(:,horn) is the apparent position of the center of the
  ! galaxy. We need to apply a rotation to all our samples which moves
  ! the center from source%point to source%poss(:,horn) without rotating
  ! the plane. This rotation should be
  ! euler(phi_to, -(theta_to-theta_from), -phi_from)
  subroutine simulate_tod_template(source, data, samps, mval)
    implicit none
    type(source_model), intent(inout) :: source
    type(lx_struct),    intent(in)    :: data
    type(sample_data),  intent(inout) :: samps
    real(dp),           intent(out), optional :: mval
    real(dp),           allocatable   :: point(:,:,:)
    integer(i4b)                      :: horn_map(0:size(quiet_horns)-1)
    integer(i4b),       allocatable   :: map_horn(:)
    logical(lgt),       allocatable   :: done(:,:)
    logical(lgt)                      :: horn_active(0:size(quiet_horns)-1)
    integer(i4b)                      :: di, nh, horn, mod, i, j, nhorn, h, d, pix, ind
    integer(i4b)                      :: nhit(2)
    real(dp)                          :: rots(3,3,0:size(quiet_horns)-1)
    real(dp)                          :: r, val(3), sigma, psi, pamp, gain, gauss, foo
    real(dp)                          :: mat(3,3), p(3)
    if(present(mval)) mval = 0
    ! Create our shifted coordinate system. In this system, the galaxy
    ! has been moved from source%point to source%poss(:,horn).
    ! First get the horns we need to consider
    horn_active = .false.
    do di = 1, size(samps%diodes)
       if(.not. source%active(di)) cycle
       horn_active(quiet_horns(quiet_diodes(di)%horn)%groups) = .true.
    end do
    nhorn = count(horn_active)
    allocate(map_horn(nhorn),done(size(data%time),nhorn))
    i = 0
    horn_map = -1
    do mod = 0, size(quiet_horns)-1
       if(.not. horn_active(mod)) cycle
       i = i+1
       horn_map(mod) = i
       map_horn(i)   = mod
    end do

    ! The rotation matrices are static per horn
    do h = 1, nhorn
       horn = map_horn(h)
       call compute_euler_matrix_zyz(source%point0(1), &
        & source%point0(2)-source%poss(2,horn+1), -source%poss(1,horn+1), rots(:,:,horn))
    end do

    ! Keep track of what coordinates we already have
    allocate(point(3, size(data%time), nhorn))
    done = .false.

    ! And use this to draw the map. No beam is applied here, so don't
    ! fit fwhm with a template.
    source%nhit = psource%nhit
    do di = 1, size(samps%diodes)
       if(.not. source%active(di)) cycle
       nh = size(quiet_horns(quiet_diodes(di)%horn)%groups)
       ! We do not need to recalculate the signal for a horn-diode if
       ! nothing for that horn-diode has changed.
       if(source%amp == psource%amp .and. &
        & all(source%poss(:,quiet_horns(quiet_diodes(di)%horn)%groups+1) == &
        &  psource%poss(:,quiet_horns(quiet_diodes(di)%horn)%groups+1)) .and. &
        & all(source%diodes(di)%amps == psource%diodes(di)%amps) .and. &
        & all(source%diodes(di)%stokes == psource%diodes(di)%stokes)) then
          cycle
       end if
       source%nhit(:,di)    = 0
       samps%diodes(di)%sim = 0
       do horn = 1, nh
          mod   = quiet_horns(quiet_diodes(di)%horn)%groups(horn)
          h     = horn_map(mod)
          do i = 1, size(samps%diodes(di)%inds)
             j = samps%diodes(di)%inds(i)
             ! If we don't already know the pointing, compute it.
             ! This could be cached for a significant speedup.
             if(.not. done(j, h)) then
                p = data%point(:,j,mod)
                call compute_euler_matrix_zyz(p(1),p(2),p(3),mat)
                call convert_euler_matrix_to_angles_zyz(matmul(rots(:,:,mod),mat), &
                 & point(1,j,h), point(2,j,h), point(3,j,h))
                done(j,h) = .true.
             end if
             ! Look up this position in the template
             pix = ang2pix(source%nside, source%order, point(2,j,h), point(1,j,h))
             ind = source%map2mask(pix)
             source%nhit(2,di) = source%nhit(2,di)+1
             if(ind <= 0) cycle
             source%nhit(1,di) = source%nhit(1,di)+1
             val = source%template(ind,1:3)
             ! Rotate polarization part
             psi  = point(3,j,h)
             pamp = sum(val(2:3)*[cos(2*psi),sin(2*psi)])
             gain = 1d-6 / ant2thermo(get_diode_freq(di)) * source%amp * source%diodes(di)%amps(horn)

             samps%diodes(di)%sim(i) = samps%diodes(di)%sim(i) + &
              & gain * sum(source%diodes(di)%stokes(1:2,horn) * [val(1),pamp])
          end do
       end do
    end do
    deallocate(point,map_horn,done)
    if(present(mval)) mval = real(sum(source%nhit(1,:)))/sum(source%nhit(2,:))
  end subroutine

  subroutine calc_point(data, source, isys, osys)
    implicit none
    type(lx_struct) :: data
    type(source_model) :: source
    integer(i4b) :: i, nsamp, mod, nmod, isys, osys
    real(dp)     :: op(3), np(3), mat(3,3)
    logical(lgt) :: relevant_horn(0:size(quiet_horns)-1)
    nsamp = size(data%tod,1)
    nmod  = get_num_modules()
    relevant_horn = .false.
    do i = 1, size(source%active)
       if(.not. source%active(i)) cycle
       relevant_horn(quiet_horns(quiet_diodes(i)%horn)%groups) = .true.
    end do
    if(allocated(data%point)) deallocate(data%point)
    allocate(data%point(3,nsamp,0:nmod-1))
    do i = 1, nsamp
      op = data%orig_point(:,i)
      call swap_coordinate_convention(op(1), op(2), op(3), isys)
      ! Slow but accurate
      !do mod = 0, nmod-1
      !   call coord_convert(isys, op(1), op(2), op(3), osys, np(1), np(2), np(3), &
      !    & mjd=data%time(i), mod=mod)
      !   data%point(:,i,mod) = np
      !end do
      !! Fast but slightly inaccurate
      call coord_convert(isys, op(1), op(2), op(3), osys, np(1), np(2), np(3), &
       & mjd=data%time(i), euler=mat)
      do mod = 0, nmod-1
         if(.not. relevant_horn(mod)) cycle
         call rot2angles(matmul(mat, rot_module2boresight(mod,-1)), np(1), np(2), np(3))
         data%point(:,i,mod) = np
      end do
    end do
    data%coord_sys = osys
  end subroutine

  subroutine filter_tod(data, source)
    implicit none
    type(lx_struct)          :: data
    type(source_model)       :: source
    real(sp),    allocatable :: tod(:,:)
    complex(sp), allocatable :: ft(:,:)
    integer(i4b),allocatable :: inds(:)
    integer(i4b)             :: ndi, n, nf, m
    n   = size(data%tod,1)
    nf  = n/2+1
    call wherei(source%active, inds)
    m   = size(inds)
    allocate(ft(nf,m),tod(n,m))
    tod = data%tod(:,inds)
    call fft_multi(tod, ft,  1)
    ft(1:100,:) = 0
    call fft_multi(tod, ft, -1)
    data%tod(:,inds) = tod
    deallocate(ft,tod,inds)
  end subroutine

  function model_size(source) result(n)
    implicit none
    type(source_model), intent(inout) :: source
    integer(i4b)                      :: n
    real(dp)                          :: foo(1)
    call transform_model(source, foo, n, 0)
  end function

  subroutine transform_model(source, arr, i, dir)
    implicit none
    type(source_model)  :: source
    real(dp)            :: arr(:), foo(1), stokamp(3), stoktmp(size(source%comps)), tmp
    real(dp)            :: mount_param(num_mount_par), op(3)
    real(dp)            :: Pgal(3,3), Pmodel(3,3), F(3,3), D(3,3), psi, bar(3)
    logical(lgt)        :: active_horns(size(quiet_horns))
    integer(i4b)        :: dir, i, mod, h, j, k
    i = 0
    if(source%fit_pos == 1) then
       ! With this model, the horns are not fit individually, so copy
       ! a single pointing to all
       call transform(source%point(1:2),  arr, i, dir)
       source%poss = spread(source%point(1:2), 2, size(quiet_horns))
    end if
    foo(1) = source%amp
    if(source%fit_amp == 1)    call transform(foo,                arr, i, dir)
    source%amp = foo(1)
    if(source%fit_ang == 1)    call transform(source%point(3:3),  arr, i, dir)
    if(source%fit_stokes == 1) call transform(source%stokes(2:2), arr, i, dir)
    do di = 1, size(source%diodes)
       if(.not. source%active(di)) cycle
       if(source%fit_amp == 2) call transform(source%diodes(di)%amps, arr, i, dir)
       if(source%fit_amp == 3) then
          ! To make the fit as linear as possible, we transform the combined amplitude
          do h = 1, size(source%diodes(di)%amps)
             stokamp = source%diodes(di)%amps(h) * source%diodes(di)%stokes(:,h)
             stoktmp = stokamp(source%comps)
             call transform(stoktmp, arr, i, dir)
             stokamp(source%comps) = stoktmp
             tmp = stokamp(maxloc(abs(stokamp),1))
             if(abs(tmp) < 1d-30) tmp = 1
             source%diodes(di)%amps(h) = tmp
             source%diodes(di)%stokes(:,h) = stokamp / tmp
          end do
       end if
       if(source%fit_ang == 2) call transform(source%angs(di:di), arr, i, dir)
    end do
    active_horns = .false.
    do j = 1, size(quiet_diodes)
       if(.not. source%active(j)) cycle
       active_horns(quiet_horns(quiet_diodes(j)%horn)%groups+1) = .true.
    end do
    do mod = 1, size(quiet_horns)
       if(.not. active_horns(mod)) cycle
       if(source%fit_fwhm == 1) call transform(source%fwhm(mod:mod), arr, i, dir)
       if(source%fit_pos  == 2) call transform(source%poss(1:2,mod), arr, i, dir)
    end do
    ! If we fit delta_dk, we need to rebuild the pointing. There are two
    ! ways of doing this. One is to actually reprocess the pointings themselves.
    ! This is the most straightforward one, but is slower. The other possibility
    ! is to account for this in the poss array and source psi. Start with the
    ! position pos, apply the focalplane with deck offset, and then unapply
    ! the focalplane without deck offset. This is then how jupiter will appear.
    ! I'll use the straightforward method here.
    !
    ! If we fit a template instead of a point source we will also
    ! need to rebuild the pointing.
    if(source%fit_dk == 1) then
       foo(1) = source%delta_dk
       call transform(foo, arr, i, dir)
       source%delta_dk = foo(1)
       if(dir < 0) then
          mount_param = 0; mount_param(1) = source%delta_dk
          call set_mount_override(.true., mount_param)
          do mod = 1, size(quiet_horns)
             if(.not. active_horns(mod)) cycle
             do j = 1, size(samps%horns(mod)%inds)
                k = samps%horns(mod)%inds(j)
                op = data%orig_point(:,k)
                call swap_coordinate_convention(op(1), op(2), op(3), coord_tele)
                call coord_convert(coord_tele, op(1), op(2), op(3), &
                 & coord_gal, op(1), op(2), op(3), mjd=data%time(k), mod=mod-1)
                data%point(:,k,mod-1) = op
             end do
          end do
          call set_mount_override(.false.)
       end if
    end if
    if(source%fit_dk == 2) then
       foo(1) = source%delta_dk
       call transform(foo, arr, i, dir)
       source%delta_dk = foo(1)
       ! Ignoring the dynamic mount corrections,
       ! Pgal = GBDCF, while our model Pmodel GBF, so
       ! Pmodel = Pgal F'C'D'F. For the time being, I will
       ! ignore collimation here, so Pmodel = Pgal F'D'F
       do mod = 0, size(quiet_horns)-1
          call cc_invdk(coord_gal, source%point(1), source%point(2), psi, &
           & coord_tele, bar(1), bar(2), real(data%orig_point(3,1),dp), mod=mod, &
           & mjd=data%time(1))
          Pgal = angles2rot(source%point(1), source%point(2), psi)
          F    = rot_module2boresight(mod, -1)
          call compute_euler_matrix_zyz(0d0, 0d0, source%delta_dk, D)
          Pmodel = matmul(Pgal,matmul(transpose(F),matmul(transpose(D),F)))
          call rot2angles(Pmodel, source%poss(1,mod+1), source%poss(2,mod+1), psi)
       end do
    end if
  end subroutine

  subroutine free_sample_data(ind)
    implicit none
    type(sample_data) :: ind
    integer(i4b)       :: i
    if(allocated(ind%diodes)) then
       do i = 1, size(ind%diodes)
          if(allocated(ind%diodes(i)%inds)) deallocate(ind%diodes(i)%inds)
          if(allocated(ind%diodes(i)%sim))  deallocate(ind%diodes(i)%sim)
       end do
       deallocate(ind%diodes)
    end if
    if(allocated(ind%horns)) then
       do i = 1, size(ind%horns)
          if(allocated(ind%horns(i)%inds)) deallocate(ind%horns(i)%inds)
       end do
       deallocate(ind%horns)
    end if
  end subroutine

  ! This copies the variable part of source1 into source2.
  ! Its purpose is to be able to compare source and source2 later,
  ! and see if any parameters have changed.
  subroutine copy_psource(source1, source2)
    implicit none
    type(source_model) :: source1, source2
    integer(i4b)       :: di
    source2%point    = source1%point
    source2%amp      = source1%amp
    source2%stokes   = source1%stokes
    source2%delta_dk = source1%delta_dk
    source2%fwhm     = source1%fwhm
    source2%poss     = source1%poss
    source2%angs     = source1%angs
    do di = 1, size(quiet_diodes)
       source2%diodes(di)%amps   = source1%diodes(di)%amps
       source2%diodes(di)%stokes = source1%diodes(di)%stokes
    end do
    source2%nhit     = source1%nhit
  end subroutine

  subroutine free_source_model(source)
    implicit none
    type(source_model) :: source
    integer(i4b)       :: i
    if(allocated(source%fwhm)) deallocate(source%fwhm)
    if(allocated(source%poss)) deallocate(source%poss)
    if(allocated(source%angs)) deallocate(source%angs)
    if(allocated(source%diodes)) then
       do i = 1, size(source%diodes)
          if(allocated(source%diodes(i)%amps))   deallocate(source%diodes(i)%amps)
          if(allocated(source%diodes(i)%stokes)) deallocate(source%diodes(i)%stokes)
       end do
       deallocate(source%diodes)
    end if
    if(allocated(source%active)) deallocate(source%active)
    if(allocated(source%comps))  deallocate(source%comps)
    if(allocated(source%compsIQUV))  deallocate(source%compsIQUV)
    if(allocated(source%template)) deallocate(source%template)
    if(allocated(source%map2mask)) deallocate(source%map2mask)
    if(allocated(source%pixels)) deallocate(source%pixels)
    if(allocated(source%nhit)) deallocate(source%nhit)
  end subroutine

  subroutine free_validator(v)
    implicit none
    type(validator) :: v
    if(allocated(v%sigma_sig))  deallocate(v%sigma_sig)
    if(allocated(v%sigma_null)) deallocate(v%sigma_null)
    if(allocated(v%sn))         deallocate(v%sn)
    if(allocated(v%accepted))   deallocate(v%accepted)
  end subroutine

  subroutine transform(b, a, i, dir)
    implicit none
    real(dp)     :: a(:), b(:)
    integer(i4b) :: i, dir
    if(dir > 0) then
       a(1+i:size(b)+i) = b
    elseif(dir < 0) then
       b = a(1+i:size(b)+i)
    end if
    i = i + size(b)
  end subroutine

  ! Two possible formats: lon, lat, [psi]
  ! or obj, [psi]
  subroutine parse_source_point(desc, mjd, point, dist, name)
    implicit none
    character(len=*), intent(in)  :: desc, name
    real(dp),         intent(in)  :: mjd
    real(dp),         intent(out) :: point(:), dist
    character(len=64),allocatable :: toks(:)
    integer(i4b)                  :: i, n, id
    point = 0
    n     = num_tokens(desc, ",")
    call assert(n >= 1, "Source pointing too short: " // trim(desc))
    allocate(toks(n))
    call get_tokens(desc, ",", toks)
    do i = 1, n; toks(i) = adjustl(toks(i)); end do
    ! Is desc an object?
    source%name = toks(1)
    id = lookup_patch(toks(1), patches)
    if(id > 0) then
       point(1:2) = get_patch_pos_single(patches(id), mjd, coord_gal)
       if(n >= 2) read(toks(2),*) point(3)
       point(3) = point(3)*DEG2RAD
       dist     = get_patch_dist(patches(id), mjd)
    else
        do i = 1, min(n,3); read(toks(i),*) point(i); end do
        point = point * DEG2RAD
        dist  = infinity
    end if
    deallocate(toks)
  end subroutine

  subroutine parse_source_stokes(desc, stokes)
    implicit none
    character(len=*), intent(in)  :: desc
    real(dp),         intent(out) :: stokes(:)
    character(len=64)             :: toks(3)
    integer(i4b)                  :: i, n
    n     = num_tokens(desc, ",")
    call assert(n == 3, "Invalid stokes: " // trim(desc))
    call get_tokens(desc, ",", toks)
    do i = 1, 3; read(toks(i),*) stokes(i); end do
  end subroutine

  ! Like coord convert, but the sys2 psi angle is known, while the sys1 psi angle isn't
  subroutine cc_invdk(sys1, phi1, theta1, psi1, sys2, phi2, theta2, psi2, mjd, mod)
    implicit none
    integer(i4b) :: sys1, sys2, mod, err
    real(dp)     :: phi1, phi2, theta1, theta2, psi1, psi2, mjd, p(1), foo
    p = 0
    powell_dk_real = [ phi1, theta1, psi2, mjd ]
    powell_dk_int  = [ sys1, sys2, mod ]
    call powell(p, powell_calc_dk, err)
    psi1 = p(1)
    call coord_convert(sys1, phi1, theta1, psi1, sys2, phi2, theta2, foo, mjd, mod)
  end subroutine

  function powell_calc_dk(p) result(res)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: p
    real(dp)                                     :: foo(2), psi2, res
    call coord_convert(powell_dk_int(1), powell_dk_real(1), powell_dk_real(2), &
     & p(1), powell_dk_int(2), foo(1), foo(2), psi2, powell_dk_real(4), powell_dk_int(3))
    res = ang_diff(psi2, powell_dk_real(3))**2
  end function

  function pos_from_map(source, data, samps) result(pos)
    implicit none
    type(source_model) :: source
    type(lx_struct)    :: data
    type(sample_data)  :: samps
    real(dp)           :: pos(2), rad
    integer(i4b)       :: nbin, i
    real(dp)           :: rms, mean, vals(size(source%compsIQUV)), locs(2,size(source%compsIQUV))
    real(dp)           :: avgloc(2)
    real(dp),     allocatable, dimension(:,:,:,:) :: map
    integer(i4b), allocatable                     :: inds(:)

    ! Build map
    nbin = 100
    rad  = rlim*3
    allocate(map(4,nbin,nbin,2))
    call wherei(source%active, inds)
    call make_map(source, data, samps, rad, rad, map, inds)
    deallocate(inds)
    where(map /= map) map = 0
    mean = sum(map(source%compsIQUV,:,:,1))/size(map(source%compsIQUV,:,:,1))
    rms  = sqrt(sum((map(source%compsIQUV,:,:,1)-mean)**2)/size(map(source%compsIQUV,:,:,1)))
    if(rms == 0) rms = 1
    do i = 1, size(source%compsIQUV)
       if(source%comps(i) == 1) then
          vals(i)   = maxval(map(i,:,:,1))
          locs(:,i) = maxloc(map(i,:,:,1))
       else
          vals(i)   = maxval(abs(map(i,:,:,1)))
          locs(:,i) = maxloc(abs(map(i,:,:,1)))
       end if
    end do
    if(any(vals/rms > 4)) then
       avgloc = sum(locs * spread(vals,1,2),2) / sum(vals) ! (row,col) -> (-lat,lon)
       pos(1) = source%point(1) + (avgloc(2)-nbin/2-0.5)*rad/(nbin/2)
       pos(2) = source%point(2) + (avgloc(1)-nbin/2-0.5)*rad/(nbin/2)
    else
       ! Should we give up in this case?
       pos = source%point(1:2)
    end if
    deallocate(map)
  end function

  function get_observed_object_temperature(object, fwhm, mjd, nu0)
    implicit none

    character(len=*), intent(in) :: object
    real(dp),         intent(in) :: mjd, fwhm, nu0
    real(dp)                     :: get_observed_object_temperature

    integer(i4b)       :: id
    real(dp)           :: omega_b, omega_p, d_fid, omega_ref, f_A, f_d, d, sigma
    real(dp)           :: D_w, R_pole, R_equ
    real(dp), dimension(2) :: pos
    real(dp), dimension(2) :: r_jup
    real(dp), dimension(5) :: nu, T_p, T_p2

    id       = name2eph(trim(object))
    d        = ephem_dist(id, mjd)
    sigma    = fwhm / 60.d0 / sqrt(8.d0*log(2.d0)) * DTOR
    omega_b  = 2.d0*pi*sigma**2

    if (trim(object) == 'jupiter') then
       nu        = [22.85d0, 33.11d0, 40.92d0, 60.41d0, 93.25d0]
       T_p       = [136.2d0, 147.2d0, 154.4d0, 165.d0,  173.d0]
       call spline(nu, T_p, 1.d30, 1.d30, T_p2)
       omega_ref = 2.481d-8
       d_fid     = 5.2d0
       D_w       = 0.d0
       R_pole    = 66854.d0 ! in km; Weiland et al. 2011
       R_equ     = 71492.d0 ! in km
    else
       get_observed_object_temperature = 0.d0
       return
    end if
    
    f_d = (d/d_fid)**2
    f_A = 1.d0                 ! D_w = 0.d0 for now; sub-0.1% error

!    get_observed_object_temperature = T_p * omega_ref/omega_b * f_A / f_d
    get_observed_object_temperature = splint(nu, T_p, T_p2, nu0) * omega_ref/omega_b * f_A / f_d

  end function get_observed_object_temperature

  subroutine output_tod(cid, data, samps)
    implicit none
    integer(i4b),       intent(in) :: cid
    type(lx_struct),    intent(in) :: data
    type(sample_data),  intent(in) :: samps
    integer(i4b)                   :: i, j, k, n, ndi
    character(len=512)             :: mapname
    type(hdf_file)                 :: hfile
    real(dp),  allocatable         :: sims(:,:)
    mapname = trim(tod_dir) // "/tod" // trim(itoa(cid,4)) // ".dat"
    n       = size(data%time)
    ndi     = size(quiet_diodes)
    allocate(sims(n,ndi))
    sims = NaN
    do i = 1, ndi
       if(.not. allocated(samps%diodes)) cycle
       do j = 1, size(samps%diodes(i)%sim)
          k = samps%diodes(i)%inds(j)
          sims(k,i) = samps%diodes(i)%sim(j)
       end do
    end do
    call mkdirs(mapname, .true.)
    call open_hdf_file(mapname, hfile, "w")
    call write_hdf(hfile, "time", data%time)
    call write_hdf(hfile, "tod",  data%tod)
    call write_hdf(hfile, "sim",  sims)
    call close_hdf_file(hfile)
    deallocate(sims)
  end subroutine

end program
