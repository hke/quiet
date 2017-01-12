! We work throughout with healpix-convention coordinates.
module point_fit_mod
  use quiet_utils
  use quiet_mpi_mod
  use quiet_pointing_mod
  use powell_mod
  use quiet_ephem_mod
  use point_fit_params
  implicit none

  integer(i4b), private ::  npar = 0

  type point_data
     integer(i4b)                                 :: n
     real(dp),          allocatable, dimension(:) :: mjd, az0, el0, az, el, dk
     real(dp),          allocatable, dimension(:) :: sigma0, sigma
     integer(i4b),      allocatable, dimension(:) :: mod, cid
     logical(lgt),      allocatable, dimension(:) :: radio
     logical(lgt),      allocatable, dimension(:) :: active
     character(len=64), allocatable, dimension(:) :: object
  end type

  type(point_data), private :: data

contains

  subroutine init_point_fit_mod(parfile)
    implicit none
    character(len=*), intent(in) :: parfile
    logical(lgt),     save       :: initialized = .false.
    if(initialized) return
    call init_point_fit_params(parfile)
    npar = get_npar()
    call initialize_quiet_pointing_mod(parfile)
    initialized = .true.
  end subroutine

  function get_best_multi(data_in, attempts) result(params)
    implicit none
    type(point_data), intent(in)           :: data_in
    integer(i4b),     intent(in), optional :: attempts
    integer(i4b)                           :: natt, i, j
    real(dp)                               :: p(npar), chisq, params(npar), best_chisq, dist
    natt = 3; if(present(attempts)) natt = attempts
    best_chisq = infinity
    call init_point_data(data, data_in%n)
    data   = data_in
    do i = 1, natt
       call init_params(p)
       p = get_best_fit(data_in, p)
       p([3,7]) = modulo(p([3,7]),2.d0*pi)
       chisq = log_likelihood(p)
       dist  = stddev_distance(p)
       if(i == 1) call print_header
       call print_params(p)
       write(*,*) "scatter: ", dist
       if(chisq < best_chisq) then
          best_chisq = chisq
          params     = p
       end if
    end do
    call free_point_data(data)
  end function

  function get_best_fit(data_in, p) result(params)
    implicit none
    type(point_data), intent(in) :: data_in
    real(dp)                     :: params(npar)
    real(dp), optional           :: p(npar)
    integer(i4b)                 :: error
    call assert(npar > 0, "npar must be set before using point_fit_mod!")
    params = 0; if(present(p)) params = p
    call powell(params, log_likelihood, error, niter=10000)
    if(error/=0) write(stderr,*) "Error in powell fit: " // trim(itoa(error))
  end function

  subroutine get_samples_mc(data_in, p0, step, samps)
    implicit none
    type(point_data) :: data_in
    integer(i4b)     :: nstep, myid, nproc, i, j, k, n, err, unit, decimate
    real(dp)         :: p0(:), samps(:,:), step(:)
    real(dp)         :: p(size(p0),2), mllik(2), aprob, nok, ntot
    call init_point_data(data, data_in%n)
    data = data_in
    n    = size(p0)
    nstep= size(samps,2)
    decimate = 4
    unit = getlun()
    !call mpi_comm_rank(mpi_comm_world, myid,  err)
    !call mpi_comm_size(mpi_comm_world, nproc, err)
    myid = 0; nproc = 1
    samps = 0
    p(:,1) = p0
    mllik(1) = log_likelihood(p0)
    nok = 0; ntot = 0
    open(unit,file="samps.txt")
    do i = 1+myid, nstep, nproc
       do k = 1, decimate
          do j = 1, n
             p(j,2) = p(j,1) + rand_gauss(rng)*step(j)
          end do
          mllik(2) = log_likelihood(p(:,2))
          aprob = exp(mllik(1)-mllik(2))
          ntot = ntot+1
          if(rand_uni(rng) < aprob) then ! accept
             p(:,1)   = p(:,2)
             mllik(1) = mllik(2)
             nok      = nok+1
          end if
       end do
       write(unit,'(i8,'//trim(itoa(n+3))//'e15.7)'), i, p(:,1), nok/ntot, mllik(1), stddev_distance(p(:,1))
       !samps(:,i) = p(:,1)
    end do
    close(unit)
    !call mpi_allreduce...
    call free_point_data(data)
  end subroutine

  subroutine clean_data(data, params)
    implicit none
    type(point_data), intent(inout) :: data
    real(dp),         intent(in)    :: params(npar)
    integer(i4b)                    :: unit, i
    real(dp)                        :: phi, theta, psi, r
    call apply_params(params)
    do i = 1, data%n
       if (data%active(i)) then
          call coord_convert(coord_tele, data%az(i), data%el(i), data%dk(i), &
               & coord_hor, phi, theta, psi, mod=data%mod(i), mjd=data%mjd(i))
!          write(*,*) theta, phi
!          write(*,*) data%el0(i), data%az0(i)
          r = polangdist([theta, phi], [data%el0(i), data%az0(i)])
          write(*,*) 'r = ', r*RAD2DEG*60 
          if (r*RAD2DEG*60 > 40.d0) then
             write(*,*) 'Observation ', i, ' rejected, r = ', r*RAD2DEG*60 
             data%active(i) = .false.
          end if
       end if
    end do
  end subroutine clean_data

  subroutine plot_fit(data, params, default_mount, fname, unit)
    implicit none
    type(point_data), intent(in) :: data
    real(dp),         intent(in) :: params(npar)
    logical(lgt),     intent(in)  :: default_mount
    character(len=*), intent(in), optional :: fname
    integer(i4b),     intent(in), optional :: unit
    integer(i4b)                 :: unit_, i, j, k
    real(dp)                     :: phi, theta, psi
    real(dp),              dimension(2)   :: mu
    real(dp),              dimension(2,2) :: sigma
    real(dp), allocatable, dimension(:,:) :: p
    if (default_mount) then
       call set_mount_override(.false.)
    else
       call apply_params(params)
    end if
    unit_ = getlun(); if(present(unit)) unit_ = unit
    if(present(fname)) open(unit_,file=fname)
    allocate(p(data%n,2))
    k = 0
    do i = 1, data%n
       if (.not. data%active(i)) cycle
       k = k+1
       call coord_convert(coord_tele, data%az(i), data%el(i), data%dk(i), &
        & coord_hor, phi, theta, psi, mod=data%mod(i), mjd=data%mjd(i))
       phi = modulo(phi, 2.d0*pi)
       p(k,:) = [theta-data%el0(i), (phi-data%az0(i))*sin(theta)] * RAD2DEG*60
!       write(*,*) theta, phi
!       write(*,*) data%el0(i), data%az0(i)
!       stop
       write(unit_,'(f13.7,5f9.5,3e15.7,a10,2i8)') data%mjd(i), phi, theta, &
        & data%az0(i), data%el0(i), data%dk(i), (phi-data%az0(i))*RAD2DEG*60, &
        & p(k,:), &
        & " "//data%object(i), data%cid(i), data%mod(i)
    end do
    if(present(fname)) close(unit_)
    ! Find best-fit 2D gaussian
    do i = 1, 2
       mu(i) = mean(p(1:k,i))
    end do
    mu = 0.d0
    do i = 1, 2
       do j = 1, 2
          sigma(i,j) = sum((p(1:k,i)-mu(i))*(p(1:k,j)-mu(j))) / real(k-1,dp)
       end do
    end do
    call get_eigenvalues(sigma, mu)
    mu = sqrt(mu)
    write(*,fmt='(a,3f8.3)') '(Major, minor, mean) semiaxis = ', mu, 0.5*sum(mu)
    deallocate(p)
  end subroutine

  function log_likelihood(params) result(lnl)
    use healpix_types
    implicit none
    real(dp), optional, intent(in) :: params(:)
    real(dp)                       :: lnl, phi, theta, psi
    real(dp), allocatable          :: r(:)
    integer(i4b)                   :: i
!if(any(params<-pi.or.params>2*pi)) then; lnl = 1e20; return; end if
    call apply_params(params)
    allocate(r(data%n))
    r = 0.d0
    do i = 1, data%n
       if (data%active(i)) then
          call coord_convert(coord_tele, data%az(i), data%el(i), data%dk(i), &
               & coord_hor, phi, theta, psi, mod=data%mod(i), mjd=data%mjd(i))
          r(i) = polangdist([theta, phi], [data%el0(i), data%az0(i)])
       end if
    end do
    lnl = 0.5*dot_product(r,r/(data%sigma**2+data%sigma0**2))
    deallocate(r)
  end function

  function stddev_distance(params) result(dist)
    use healpix_types
    implicit none
    real(dp), optional, intent(in) :: params(:)
    real(dp)                       :: dist, phi, theta, psi, r
    integer(i4b)                   :: i, n
    call apply_params(params)
    dist = 0.d0
    n    = 0
    do i = 1, data%n
       if (data%active(i)) then
          call coord_convert(coord_tele, data%az(i), data%el(i), data%dk(i), &
               & coord_hor, phi, theta, psi, mod=data%mod(i), mjd=data%mjd(i))
          r    = polangdist([theta, phi], [data%el0(i), data%az0(i)])*RAD2DEG*60.
          dist = dist + r**2
          n    = n+1
       end if
    end do
    dist = sqrt(dist/n)
  end function

  subroutine init_point_data(d, n)
    implicit none
    type(point_data), intent(inout) :: d
    integer(i4b),      intent(in)    :: n
    call free_point_data(d)
    d%n = n
    allocate(d%mjd(n), d%az0(n), d%el0(n), d%az(n), d%el(n), d%dk(n))
    allocate(d%mod(n), d%radio(n), d%object(n), d%sigma0(n), d%sigma(n), d%active(n), d%cid(n))
  end subroutine

  ! If anything is allocated, all is assumed to be, for simplicity
  subroutine free_point_data(d)
    implicit none
    type(point_data), intent(inout) :: d
    if(.not. allocated(d%mjd)) return
    deallocate(d%mjd, d%az0, d%el0, d%az, d%el, d%dk, d%mod, d%radio, d%object, d%active)
    deallocate(d%sigma0, d%sigma, d%cid)
  end subroutine

  subroutine read_point_data(d, fname)
    implicit none
    type(point_data), intent(inout) :: d
    character(len=*),  intent(in)    :: fname
    character(len=1)                 :: sym
    integer(i4b)                     :: unit, i, j, k, m, n
    real(dp)                         :: mjd, psi, theta, phi, ang(2)
    unit = getlun()
    n    = 0
    psi  = 0
    open(unit,file=fname,action="read",status="old")
    do
       read(unit,*,end=1) mjd
       n = n+1
    end do
1   call init_point_data(d, n)
    rewind(unit)
    do i = 1, n
       read(unit,*) d%mjd(i), d%az0(i), d%el0(i), d%az(i), d%el(i), d%dk(i), &
        & d%object(i), d%mod(i), d%cid(i)
       ! To radians and healpix
       d%az0(i) = d%az0(i)*DEG2RAD; d%el0(i) = d%el0(i)*DEG2RAD
       d%az(i)  = (d%az(i)) *DEG2RAD; d%el(i)  = d%el(i) *DEG2RAD; d%dk(i) = d%dk(i)*DEG2RAD
       call swap_coordinate_convention(d%az(i),  d%el(i),  d%dk(i), coord_tele)
       call swap_coordinate_convention(d%az0(i), d%el0(i), psi,     coord_hor)
       
       d%az(i) = modulo(d%az(i),2*pi)
       d%az0(i) = modulo(d%az0(i),2*pi)

    end do
    d%sigma  = 2d0/60/180*pi
    d%sigma0 = 0.1d0/60/180*pi
    d%active = .true.
    close(unit)
  end subroutine

end module
