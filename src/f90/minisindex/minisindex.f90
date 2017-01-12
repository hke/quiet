program minisindex
  use quiet_mapfile_mod
  use quiet_fileutils
  use rngmod
  implicit none

  integer(i4b)          :: i, j, nfreq, nside, npix, nmaps, n, ntot, unit, polarization, k, numbin, ordering
  integer(i8b)          :: m, seed
  real(dp)              :: beta_min, beta_max, sigma_reg, nu0_synch
  logical(lgt)          :: inv
  character(len=512)    :: maskfile, filename, parfile, outfile
  integer(i4b)          :: mod, npar
  character(len=1)      :: itext
  logical(lgt), dimension(3) :: inc_pol
  type(planck_rng)      :: rng_handle
  real(dp),     allocatable, dimension(:)     :: freqs, a2t, beta, lnL, offset, amps
  real(dp),     allocatable, dimension(:,:)   :: mask, map, map_in, T, fg, cov_in
  real(dp),     allocatable, dimension(:,:,:) :: invcov
  integer(i4b), pointer,     dimension(:)     :: map2mask, mask2map

  call getarg(1, parfile)

  ! Read parameters
  unit = getlun()
  call get_parameter(unit, parfile, 'MASKFILE',     par_string=maskfile)
  call get_parameter(unit, parfile, 'NUMFREQ',      par_int=nfreq)
  call get_parameter(unit, parfile, 'NSIDE',        par_int=nside)
  call get_parameter(unit, parfile, 'NMAPS',        par_int=nmaps)
  call get_parameter(unit, parfile, 'INCLUDE_T',    par_lgt=inc_pol(1))
  call get_parameter(unit, parfile, 'INCLUDE_Q',    par_lgt=inc_pol(2))
  call get_parameter(unit, parfile, 'INCLUDE_U',    par_lgt=inc_pol(3))
  call get_parameter(unit, parfile, 'OUTFILE',      par_string=outfile)
  call get_parameter(unit, parfile, 'BETA_MIN',     par_dp=beta_min)
  call get_parameter(unit, parfile, 'BETA_MAX',     par_dp=beta_max)
  call get_parameter(unit, parfile, 'BETA_NUMBIN',  par_int=numbin)
  call get_parameter(unit, parfile, 'SIGMA_REG',    par_dp=sigma_reg)
  call get_parameter(unit, parfile, 'NU0_SYNCH',     par_dp=nu0_synch)
  !call get_parameter(unit, parfile, 'SEED',         par_int=seed)
  seed = 182394
  npix = 12*nside**2
  call rand_init(rng_handle, seed)
  
  ! Read mask
  call read_map(mask, ordering, maskfile)
  do i = 0, npix-1
     if (mask(i,1) == -1.6375d30) then
        mask(i,1) = 0.d0
     else
        mask(i,1) = 1.d0
     end if
  end do
  if (ordering == 1) call convert_ring2nest(nside, mask(:,1))
  call get_map2mask_from_mask(mask(:,1), map2mask, mask2map)
  n = sum(mask(:,1))
  ntot = 0
  do i = 1, 3
     if (inc_pol(i)) ntot = ntot + n
  end do
  write(*,*) 'Total number of pixels = ', ntot

  ! Read maps and covariance matrices
  allocate(freqs(nfreq), a2t(nfreq))
  allocate(map(ntot,nfreq), fg(ntot,nfreq))
  allocate(invcov(ntot,ntot,nfreq))
  do i = 1, nfreq
     call int2string(i, itext)
     call get_parameter(unit, parfile, 'FREQUENCY'//itext, par_dp=freqs(i))
     a2t(i) = ant2thermo(freqs(i))
     call get_parameter(unit, parfile, 'MAPNAME'//itext, par_string=filename)
     call read_map(map_in, ordering, filename)
     k = 1
     do j = 1, 3
        if (inc_pol(j)) then
           map(k:k+n-1,i) = map_in(mask2map,j)
           k = k+n
        end if
     end do
     deallocate(map_in)

     call get_parameter(unit, parfile, 'COVNAME'//itext, par_string=filename)
     call read_covmatrix(unit, filename, ordering, polarization, cov_in, inv, m)
     if (all(inc_pol(2:3))) then
        invcov(:,:,i) = cov_in
     else if (inc_pol(2)) then
        invcov(:,:,i) = cov_in(1:n,1:n)
     else if (inc_pol(3)) then
        invcov(:,:,i) = cov_in(n+1:2*n,n+1:2*n)
     end if
     deallocate(cov_in)
     if (i == 2) then
        open(58,file='data.dat')
        do j = 1, ntot
           write(58,*) real(map(j,:),sp), sqrt(real(invcov(j,j,:),sp))
        end do
        close(58)
     end if

     if (.not. inv) then
        do j = 1, ntot
           invcov(j,j,i) = invcov(j,j,i) + sigma_reg**2           
           map(j,i)      = map(j,i)      + sigma_reg * rand_gauss(rng_handle)
        end do
        call invert_matrix(invcov(:,:,i))
     end if
  end do

  ! Compute grid in beta
  allocate(beta(numbin), lnL(numbin))
  allocate(offset(nfreq-1), amps(ntot))
  lnL = 0.d0
  do i = 1, numbin
     beta(i) = beta_min + (i+0.5d0)*(beta_max-beta_min)/numbin
     call compute_fg(a2t, freqs, map, invcov, beta(i), fg)
     do j = 1, nfreq
        lnL(i) = lnL(i) - 0.5d0 * sum((map(:,j)-fg(:,j)) * matmul(invcov(:,:,j), map(:,j)-fg(:,j)))
     end do
     write(*,*) 'beta = ', beta(i), ', chisq = ', -2.d0*lnL(i)
  end do

  ! Normalize and output likelihood
  lnL = exp(lnL - maxval(lnL))
  open(unit,file=trim(outfile))
  do i = 1, numbin
     write(unit,*) beta(i), lnL(i)
  end do
  close(unit)

  deallocate(map, invcov)

contains

  subroutine compute_fg(a2t, freqs, map, invcov, beta, fg)
    implicit none

    real(dp),                   intent(in)  :: beta
    real(dp), dimension(:),     intent(in)  :: a2t, freqs
    real(dp), dimension(:,:),   intent(in)  :: map
    real(dp), dimension(:,:,:), intent(in)  :: invcov
    real(dp), dimension(:,:),   intent(out) :: fg

    integer(i4b) :: i, j, n
    real(dp), allocatable, dimension(:)   :: b, x
    real(dp), allocatable, dimension(:,:) :: A

    n = size(map,1)

    allocate(A(n,n), b(n), x(n))

    ! Set up equation set
    A = 0.d0
    b = 0.d0
    do k = 1, nfreq
       do i = 1, n
          do j = 1, n
             A(j,i) = A(j,i) + a2t(k)**2 * (freqs(k)/nu0_synch)**(2*beta) * invcov(j,i,k)
          end do
       end do
       b = b + a2t(k) * (freqs(k)/nu0_synch)**beta * matmul(invcov(:,:,k),map(:,k))
    end do
    
    ! Solve system
    call solve_linear_system(A, x, b)

    ! Project to frequencies
    do j = 1, nfreq
       do i = 1, n
          fg(i,j) = x(i) * a2t(j) * (freqs(j)/nu0_synch)**beta
       end do
    end do

    deallocate(A, b, x)

  end subroutine compute_fg


end program minisindex
