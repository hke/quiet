program gen
  use quiet_utils
  use quiet_lx_mod
  use quiet_fileutils
  use rngmod
  implicit none
  character(len=512)        :: arg, l3file, mapfile, maskfile, ofile
  real(dp)                  :: sigma, fknee, alpha, nu, samprate, gain
  integer(i4b)              :: i, j, k, m, n, nsim, nside_mask, order_mask, seed, maxn
  integer(i4b)              :: nside_map, order_map
  type(lx_struct)           :: data
  type(hdf_file)            :: hfile
  type(planck_rng)          :: rng
  real(dp),     allocatable :: point(:,:), mask(:,:), map(:,:), todmap(:)
  integer(i4b), allocatable :: hits(:,:), todmask(:)
  real(sp),     allocatable :: ps(:), vals(:,:), noise(:), time(:)
  complex(spc), allocatable :: ft(:)

  ! Read a level1-file, from which we will get the pointing for the first horn.
  ! Create an automatic mask based on this.
  ! Create a cmb realization
  ! Create a 1/f realization

  seed  = 1
  gain  = 1e-6
!  maxn  = 1000

  if(iargc() < 8) then
     write(stderr,'(a)') "Syntax: gen l3file maskfile mapfile nsim sigma fknee alpha ofile"
     stop
  end if
  call getarg(1, l3file)
  call getarg(2, maskfile)
  call getarg(3, mapfile)
  call getarg(4, arg); read(arg,*) nsim
  call getarg(5, arg); read(arg,*) sigma
  call getarg(6, arg); read(arg,*) fknee
  call getarg(7, arg); read(arg,*) alpha
  call getarg(8, ofile)

  call read_l3_file(l3file, data)
  n = size(data%point,2)
  allocate(point(3,n),vals(n,nsim),todmask(n),todmap(n))
  point    = data%point(:,:n,1)
  samprate = data%samprate
  call free_lx_struct(data)

  call read_map(mask, order_mask, maskfile, nside=nside_mask)
  if(mapfile == '') then
     nside_map = 1
     order_map = ring
     allocate(map(0:12*nside_map**2-1,1))
     map = 0
  else
     call read_map(map,  order_map,  mapfile, nside=nside_map)
  end if

  ! Build our tod and mask
  todmask = 0
  do i = 1, n
     j = ang2pix(nside_mask, order_mask, point(2,i), point(1,i))
     if(all(mask(j,:)>0.1)) todmask(i) = 1
     j = ang2pix(nside_map,  order_map,  point(2,i), point(1,i))
     if(size(map,2) == 1) then
        todmap(i) = map(j,1)
     elseif(size(map,2) == 3) then
        todmap(i) = sum(map(j,2:3)*[cos(2*point(3,i)),sin(2*point(3,i))])
     end if
  end do
  todmap = todmap * gain

  ! Set up power spectrum
  m = n/2+1
  allocate(ps(m), ft(m), noise(n))
  ps(1) = 0
  do i = 2, m
     nu    = ind2freq(i, samprate, m)
     ps(i) = (sigma**2 * (1.d0 + (nu/fknee)**alpha))
  end do

  ! Create realizations
  call rand_init(rng, seed)
  do k = 1, nsim
     do i = 1, m
        ft(i) = cmplx(rand_gauss(rng),rand_gauss(rng))/sqrt(2.0)*ps(i)**0.5
     end do
     call fft(noise, ft, -1)
     vals(:,k) = noise + todmap
     ! Add strange stuff to masked region
     nu = rand_gauss(rng)
     do i = 1, n
        if(todmask(i)==0) vals(i,k) = vals(i,k)/(0.1+sin(nu*i/10)**2)
     end do
  end do

  allocate(time(n))
  do i = 1, n
     time(i) = (i-1)/samprate
  end do

  ! Write to hdf
  call open_hdf_file(ofile, hfile, "w")
  call write_hdf(hfile, "mask", todmask)
  call write_hdf(hfile, "data", vals)
  call write_hdf(hfile, "samprate", samprate)
  call write_hdf(hfile, "time", time)
  call close_hdf_file(hfile)

end program
