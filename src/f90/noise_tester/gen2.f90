program gen
  use quiet_utils
  use quiet_lx_mod
  use quiet_fileutils
  use quiet_module_mod
  use quiet_pointing_mod
  use rngmod
  implicit none
  character(len=512)        :: arg, l3file, cmdfile, ofile, type, cmd, fname, parfile, str
  character(len=32)         :: posstr(2)
  real(dp)                  :: sigma, fknee, alpha, nu, center(2), amp, gain, pos(2)
  real(dp)                  :: stokes(3), psi, r, rad, fwhm, pol, p(3)
  integer(i4b)              :: i, j, k, m, n, nsim, seed, maxn, unit, ntemp
  integer(i4b)              :: nside, order, pix, di, mod
  type(lx_struct)           :: data
  type(hdf_file)            :: hfile
  type(planck_rng)          :: rng
  real(dp),     allocatable :: tod(:,:), templates(:,:), mask(:), point(:,:), map(:,:)
  real(dp),     allocatable :: sig(:), ps(:)
  complex(dpc), allocatable :: ft(:)
  integer(i4b), allocatable :: comps(:)

  call getarg(1, parfile)
  call getarg(2, l3file)
  call getarg(3, cmdfile)
  call getarg(4, ofile)

  call initialize_module_mod(parfile)
  call initialize_quiet_pointing_mod(parfile)

  seed = 1
  di   = 1
  mod  = quiet_diodes(di)%horn
  call rand_init(rng, seed)
  allocate(comps(2)); comps = [ 2, 3 ]

  ! Read a level3-file, which will act as a template
  call read_l3_file(l3file, data)
  ! Scan through command file, in order to determine the number
  ! of templates (these are added to the tod after an amplitude).
  n     = min(10000,size(data%time))
  nsim  = 100
  ntemp = 0
  unit  = getlun()
  open(unit,file=cmdfile,action="read",status="old")
  do
     read(unit,*,end=1) type
     if(type == "tod") ntemp = ntemp + 1
  end do
1 continue
  allocate(tod(n,nsim), templates(n,ntemp),mask(n),point(3,n),sig(n))
  mask      = 1
  tod       = 0
  templates = 0
  gain      = data%gain(1,di) * 1d-9 * ant2thermo(quiet_diodes(di)%freq)
  point     = data%point(:,1:n,1)
  center    = average_pointing(point)

  ! Process the command file properly
  ntemp = 0
  rewind(unit)
  do
     read(unit,*,end=2) type
     backspace(unit)
     read(unit,'(a)') str
     write(stderr,'(a)') trim(str)
     backspace(unit)
     select case(type)
        case("tod")
           read(unit,*) type, cmd
           backspace(unit)
           select case(cmd)
              case("map")
                 read(unit,*) type, cmd, fname, amp
                 call read_map(map, order, fname, nside=nside)
                 where(abs((map-hpx_dbadval)/hpx_dbadval) < 1e-5) map = 0
                 do i = 1, n
                    pix    = ang2pix(nside, order, point(2,i), point(1,i))
                    psi    = point(3,i)
                    stokes = [ 1d0, 2*cos(psi), sin(2*psi) ]
                    sig(i) = sum(map(pix,comps)*stokes(comps)) * gain
                 end do
                 deallocate(map)
              case("hor")
                 read(unit,*) type, cmd, fname, amp
                 call read_map(map, order, fname, nside=nside)
                 where(abs((map-hpx_dbadval)/hpx_dbadval) < 1e-5) map = 0
                 do i = 1, n
                    p      = data%orig_point(:,i)
                    call swap_coordinate_convention(p(1),p(2),p(3),coord_tele)
                    call coord_convert(coord_tele, p(1), p(2), p(3), coord_hor, &
                     & p(1), p(2), p(3), data%time(i), mod, quiet_diodes(di)%sub)
                    pix    = ang2pix(nside, order, p(2), p(1))
                    stokes = [ 1d0, 2*cos(p(3)), sin(2*p(3)) ]
                    sig(i) = sum(map(pix,comps)*stokes(comps)) * gain
                 end do
                 sig = sig*amp
                 deallocate(map)
              case("point") ! Point source
                 read(unit,*) type, cmd, posstr, fwhm, pol, amp
                 sigma = fwhm/sqrt(8*log(2d0))*DEG2RAD
                 call parse_pos(posstr, center, pos)
                 do i = 1, n
                    r   = polangdist(point([2,1],i), pos([2,1]))
                    psi = point(3,i)
                    stokes = [ 1d0, 2*cos(psi)*pol, 2*sin(psi)*pol ]
                    sig(i) = sum(stokes(comps))*exp(-0.5*r**2/sigma**2) * gain
                 end do
                 sig = sig*amp
              case default
                 call assert(.false., "Unknown command " // trim(cmd))
           end select
           ntemp = ntemp + 1
           templates(:,ntemp) = sig
           tod = tod + spread(sig, 2, nsim)
        case("noise")
           read(unit,*) type, sigma, fknee, alpha
           ! Set up power spectrum
           m = n/2+1
           allocate(ps(m), ft(m))
           ps(1) = 0
           do i = 2, m
              nu    = ind2freq(i, data%samprate, m)
              ps(i) = (sigma**2 * (1.d0 + (nu/fknee)**alpha))
           end do
           ! Create realizations
           call rand_init(rng, seed)
           do k = 1, nsim
              do i = 1, m
                 ft(i) = cmplx(rand_gauss(rng),rand_gauss(rng))/sqrt(2.0)*ps(i)**0.5
              end do
              call fft(sig, ft, -1)
              tod(:,k) = tod(:,k) + sig
           end do
           deallocate(ft, ps)
        case("mask")
           read(unit,*) type, cmd
           backspace(unit)
           select case(cmd)
              case("map")
                 read(unit,*) type, cmd, fname
                 call read_map(map, order, fname, nside=nside)
                 where(map == hpx_dbadval) map = 0
                 do i = 1, n
                    pix = ang2pix(nside, order, point(2,i), point(1,i))
                    sig(i) = map(pix,1)
                 end do
                 deallocate(map)
              case("point")
                 read(unit,*) type, cmd, posstr, rad
                 rad = rad*DEG2RAD
                 call parse_pos(posstr, center, pos)
                 sig = 1
                 do i = 1, n
                    r = polangdist(point([2,1],i), pos([2,1]))
                    if(rad > 0) then
                       if(r <  rad) then; sig(i) = 0; else; sig(i) = 1; end if
                    else
                       if(r < -rad) then; sig(i) = 1; else; sig(i) = 0; end if
                    end if
                 end do
              case default
                 call assert(.false., "Unknown mask " // trim(cmd))
           end select
           mask = mask * sig
     end select
  end do
2 continue

  ! And finally output the result
  call open_hdf_file(ofile, hfile, "w")
  call write_hdf(hfile, "data",       tod)
  call write_hdf(Hfile, "templates", templates)
  call write_hdf(hfile, "mask",      mask)
  call write_hdf(hfile, "samprate",  data%samprate)
  call close_hdf_file(hfile)

contains

  subroutine parse_pos(strpos, defpos, pos)
    implicit none
    character(len=*) :: strpos(:)
    real(dp)         :: defpos(:), pos(:)
    integer(i4b)     :: i
    pos = defpos
    do i = 1, size(pos)
       read(strpos(i),*,err=1) pos(i)
       pos(i) = pos(i) * DEG2RAD
       1 continue
    end do
  end subroutine

end program
