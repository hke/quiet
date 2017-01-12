! Given an input mount model and a source specification,
! generate level3-files with simulated data.
program point_sim
  use quiet_task_mod
  use quiet_detector_mod
  use quiet_pointing_mod
  use quiet_lx_mod
  use quiet_ces_mod
  use powell_mod
  use quiet_mpi_utils
  use quiet_patch_detect_mod
  use quiet_shared_output_mod
  implicit none

  type source_model
     real(dp)                      :: point(3), amp, stokes(3)
  end type

  character(len=512)    :: parfile, odir, source_pos, source_amp, source_stokes, fname
  integer(i4b)          :: myid, nproc, err, nces, di, cnum, nparam, debug
  type(task_list)       :: tasks
  type(lx_struct)       :: data
  type(quiet_ces_info)  :: ces
  type(source_model)    :: source
  type(planck_rng)      :: rng

  call getarg(1, parfile)

  call initialize_ces_mod(parfile)
  call initialize_quiet_pointing_mod(parfile)
  call initialize_patch_detect_mod(parfile)
  call get_parameter(0, parfile, 'OUTPUT_DIR',    par_string=odir)
  call get_parameter(0, parfile, 'DEBUG',         par_int=debug)
  call get_parameter(0, parfile, 'SOURCE_POS',    par_string=source_pos)
  call get_parameter(0, parfile, 'SOURCE_AMP',    par_string=source_amp)
  call get_parameter(0, parfile, 'SOURCE_STOKES', par_string=source_stokes)

  call mpi_init(err)
  call mpi_comm_rank(MPI_COMM_WORLD, myid,  err)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, err)
  call dset(id=myid, level=debug)

  call print_host_mapping
  call mkdirs(odir, .false.)

  nces = get_num_ces()
  call init_task_list(tasks, trim(odir)//"/lock.dat", nces, MPI_COMM_WORLD)
  do while(get_next_task(tasks, cnum))
     call get_ces_info(cnum, ces)
     write(*,fmt="(i3,a,i4,a)") myid, " processing ces ", ces%cid, &
      & " (" // trim(itoa(cnum)) // "/" // trim(itoa(nces)) // ")"
     call dmem("ces start")
     call rand_init(rng, ces%cid)
     call read_l3_file(ces%l3file, data);              call dmem("l3 read")
     call calc_point(data, coord_tele, coord_gal);     call dmem("calc point")
     call setup_model(source, data);                   call dmem("setup model")
     call simulate_tod(source, data);                  call dmem("simulate tod")
     fname = trim(odir)//"/level3/"//trim(ces%object)//"/"//trim(ces%object)//"_"//trim(itoa(ces%cid))//".hdf"
     call mkdirs(fname, .true.)
     call write_l3_file(fname, data)
     call free_lx_struct(data)
  end do
  call mpi_finalize(err)

contains

  subroutine setup_model(source, data)
    implicit none
    type(source_model), intent(inout) :: source
    type(lx_struct),    intent(in)    :: data
    call parse_source_point(source_pos, data%time(1), source%point)
    call parse_source_stokes(source_stokes, source%stokes)
    read(source_amp,*) source%amp
  end subroutine

  subroutine simulate_tod(source, data)
    implicit none
    type(source_model), intent(in)    :: source
    type(lx_struct),    intent(inout) :: data
    integer(i4b)                      :: di, nh, horn, mod1, mod2, i, j
    real(dp)                          :: r, val, sigma, psi, pamp, gain, gauss
    data%tod = 0
    do di = 1, size(data%tod,2)
       mod1 = quiet_diodes(di)%horn
       nh   = size(quiet_horns(mod1)%groups)
       do horn = 1, nh
          mod2  = quiet_horns(mod1)%groups(horn)
          sigma = quiet_horns(mod2)%fwhm*fwhm2sigma
          do i = 1, size(data%tod,1)
             r    = polangdist(real(data%point([2,1],i,mod2),dp), source%point([2,1]))
             psi  = data%point(3,i,mod2) + quiet_diodes(di)%psi
             pamp = cos(2*(psi-source%point(3)))
             gain = 1d-6 / ant2thermo(get_diode_freq(di)) * data%gain(1, di)
             gauss= exp(-0.5*(r/sigma)**2)
             val  = source%amp * quiet_horns(mod1)%amps(horn) * gauss * gain
             data%tod(i,di) = data%tod(i,di) + & 
              & val * sum(source%stokes*quiet_diodes(di)%stokes*[1d0,pamp,1d0])
          end do
       end do
       ! Add white noise
       do i = 1, size(data%tod,1)
          data%tod(i,di) = data%tod(i,di) + data%sigma0(di)*rand_gauss(rng)
       end do
    end do
  end subroutine

  subroutine calc_point(data, isys, osys)
    implicit none
    type(lx_struct) :: data
    integer(i4b) :: i, nsamp, mod, nmod, isys, osys
    real(dp)     :: op(3), np(3), mat(3,3)
    nsamp = size(data%tod,1)
    nmod  = get_num_modules()
    if(allocated(data%point)) deallocate(data%point)
    allocate(data%point(3,nsamp,0:nmod-1))
    do i = 1, nsamp
      op = data%orig_point(:,i)
      call swap_coordinate_convention(op(1), op(2), op(3), isys)
      ! Slow but accurate
      do mod = 0, nmod-1
         call coord_convert(isys, op(1), op(2), op(3), osys, np(1), np(2), np(3), &
          & mjd=data%time(i), mod=mod)
         data%point(:,i,mod) = np
      end do
      !! Fast but slightly inaccurate
      !call coord_convert(isys, op(1), op(2), op(3), osys, np(1), np(2), np(3), &
      ! & mjd=data%time(i), euler=mat)
      !do mod = 0, nmod-1
      !   call rot2angles(matmul(mat, rot_module2boresight(mod,-1)), np(1), np(2), np(3))
      !   data%point(:,i,mod) = np
      !end do
    end do
    data%coord_sys = osys
  end subroutine

  ! Two possible formats: lon, lat, [psi]
  ! or obj, [psi]
  subroutine parse_source_point(desc, mjd, point)
    implicit none
    character(len=*), intent(in)  :: desc
    real(dp),         intent(in)  :: mjd
    real(dp),         intent(out) :: point(:)
    character(len=64),allocatable :: toks(:)
    integer(i4b)                  :: i, n, id
    point = 0
    n     = num_tokens(desc, ",")
    call assert(n >= 1, "Source pointing too short: " // trim(desc))
    allocate(toks(n))
    call get_tokens(desc, ",", toks)
    do i = 1, n; toks(i) = adjustl(toks(i)); end do
    ! Is desc an object?
    id = lookup_patch(toks(1), patches)
    if(id > 0) then
       point(1:2) = get_patch_pos_single(patches(id), mjd, coord_gal)
       if(n >= 2) read(toks(2),*) point(3)
       point(3) = point(3)*DEG2RAD
    else
        do i = 1, min(n,3); read(toks(i),*) point(i); end do
        point = point * DEG2RAD
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

end program
