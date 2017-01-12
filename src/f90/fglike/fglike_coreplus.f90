program fglike
  use healpix_types
  use rngmod
  use quiet_utils
  implicit none

  real(dp), parameter :: k_B     = 1.3806503d-23
  real(dp), parameter :: h       = 6.626068d-34
  real(dp), parameter :: c       = 2.99792458d8
  real(dp), parameter :: T_0     = 2.726d0

  integer(i4b) :: i, j, k, l, m, n, q, p, r, numband, numcomp, numsim, seed
  real(dp)     :: A_cmb, A_dust, min_cmb, max_cmb, min, fg_min, A_range(2,3), dA(3), T_d, alpha, rms0, beta_s
  real(dp)     :: A_in(3), beta(2), dbeta(2), s, mu, sigma, deg_fac, pixsize, polfactor, npix, sky_area
  real(dp), allocatable, dimension(:,:)  :: x, summary, lnL_sub
  real(dp), allocatable, dimension(:)    :: marg, d, a2t, rms, nu, NET, NET_array, ndet, obstime, beam
  real(dp), allocatable, dimension(:,:,:)  :: lnL, dust_spec
  type(planck_rng) :: handle

  numband = 21
  numcomp = 3
  numsim  = 100
  seed    = 8481
  
  call rand_init(handle, seed)
  
  allocate(d(numband), rms(numband), a2t(numband), nu(numband), beam(numband))
  allocate(NET(numband), NET_array(numband), ndet(numband), obstime(numband))

  !obstime   = [2.5d0, 2.5d0, 2.5d0] * 24.d0*3600.d0      ! Observation time in seconds
  !NET       = [53.368632d0, 194.34570d0, 534.35253d0] ! uK sqrt(s)
  !ndet      = [3000.d0, 2167.d0, 833.d0]              ! Number of detectors
  !polfactor = sqrt(2.d0)                              ! 1 for radiometers, sqrt(2) for bolometers/MKIDs
  !NET_array = NET / sqrt(ndet) * polfactor
  !sky_area  = 41000.d0                                 ! Square degrees
  pixsize   = 60.d0                                   ! arcmin
  !npix      = sky_area / (pixsize/60.d0)**2

  !rms = [23.2*14.d0, 11.d0*9.5d0, 6.d0*7.2d0, 12.d0*5.d0, 43.d0*5.d0]/pixsize*sqrt(2.d0) !NET_array / sqrt(obstime) * sqrt(npix)
  !rms(1) = 0.2d0
!  rms = [23.2d0*14.d0, 11.d0*9.5d0, 6.d0*7.2d0, 12.d0*5.d0, 43.d0*5.d0]/pixsize !NET_array / sqrt(obstime) * sqrt(npix)
  nu  = [60.d0, 70.d0, 80.d0, 90.d0, 100.d0, 115.d0, 130.d0, 145.d0, 160.d0, 175.d0, 195.d0, 220.d0, &
       & 255.d0, 295.d0, 340.d0, 390.d0, 450.d0, 520.d0, 600.d0, 700.d0, 800.d0]*1d9
  rms = [13.8d0, 12.9d0, 8.6d0, 6.7d0, 6.1d0, 4.8d0, 4.2d0, 3.6d0, 3.2d0, 3.3d0, 3.5d0, 4.5d0, 5.8d0, &
       & 11.4d0, 20.7d0, 47.6d0, 100.9d0, 256.1d0, 779.2d0, 4657.8d0, 20505.9d0]
  beam = [14.d0, 12.d0, 10.5d0, 9.33d0, 8.4d0, 7.3d0, 6.46d0, 5.79d0, 5.25d0, 4.8d0, 4.31d0, 3.82d0, &
       & 3.29d0, 2.85d0, 2.47d0, 2.15d0, 1.87d0, 1.62d0, 1.4d0, 1.2d0, 1.05d0]
  rms = rms / pixsize
!  rms(2) = 2.23
!  rms(3) = 6.7
!  nu  = [150.d0, 217.d0, 353.d0]*1.d9
  do i = 1, numband
     a2t(i) = compute_ant2thermo_single(nu(i))
  end do
  write(*,*) 'a2t  = ', real(a2t,sp)
  write(*,*) 'rms  = ', real(rms,sp)
  write(*,*) 'rms0 = ', 1.d0/sqrt(sum(1.d0/rms**2)) * pixsize
  write(*,*)

  A_in = [1.d0, 100.d0, 1.d0]
  beta = [1.6d0, 18.d0]
  dbeta = [0.1d0, 1d0]*1 ! Set to zero if you don't want to marginalize over spectral parameters
  beta_s = -3.d0

  n = 50
  m = 10
  A_range(1,1) = 0.5d0
  A_range(2,1) = 1.5d0
  A_range(1,2) = 100.d0-1.d0
  A_range(2,2) = 100.d0+1.d0
  A_range(1,3) = 1.d0-0.01d0
  A_range(2,3) = 1.d0+0.01d0

  allocate(lnL(n,n,n), x(n, numcomp), marg(n), lnL_sub(m,m), summary(numsim,3))
  allocate(dust_spec(numband,m,m))
  do j = 1, numcomp
     dA(j)        = (A_range(2,j)-A_range(1,j)) / (n-1.d0)
     do i = 1, n
        x(i,j) = A_range(1,j) + (i-1)*dA(j)
     end do
  end do

  do j = 1, numband
     do q = 1, m
        do p = 1, m
           alpha = (beta(1)-dbeta(1)) + (q-1) * 2*dbeta(1)/(m-1)
           T_d = (beta(2)-dbeta(2)) + (p-1) * 2*dbeta(2)/(m-1)
           dust_spec(j,q,p) = compute_one_comp_dust_spectrum(nu(j), nu(15), alpha, T_d)
        end do
     end do
  end do


  open(59,file='cmb_allsims.dat')
  do i = 1, numsim
     ! Generate simulation
     do j = 1,  numband
        d(j) = A_in(1) + A_in(2) * a2t(j) * compute_one_comp_dust_spectrum(nu(j), nu(15), beta(1), beta(2)) + &
             & A_in(3) * a2t(j) * compute_synch_spectrum(nu(j), nu(1), beta_s) +  rms(j) * rand_gauss(handle)
     end do

     lnL = 0.d0
     do j = 1, n
        !write(*,*) '    j = ', j
        do k = 1, n
           do r = 1, n

!              lnL_sub = 0.d0
!              do q = 1, m
!                 alpha = (beta(1)-dbeta(1)) + (q-1) * 2*dbeta(1)/(m-1)
!                 do p = 1, m
!                    T_d = (beta(2)-dbeta(2)) + (p-1) * 2*dbeta(2)/(m-1)
!                    do l = 1, numband
                       
                       !                    s        = x(j,1) + x(k,2) * a2t(l) * compute_one_comp_dust_spectrum(nu(l), nu(1), beta(1), beta(2))
                       !s        = x(j,1) + x(k,2) * a2t(l) * compute_one_comp_dust_spectrum(nu(l), nu(1), alpha, T_d) + x(r,3) * a2t(l) * compute_synch_spectrum(nu(l), nu(1), beta_s) 
!                       s        = x(j,1) + x(k,2) * a2t(l) * dust_spec(l,q,p) + x(r,3) * a2t(l) * compute_synch_spectrum(nu(l), nu(1), beta_s) 
!                       lnL_sub(q,p) = lnL_sub(q,p) - 0.5d0*((d(l)-s)/rms(l))**2
!                    end do
!                 end do
!              end do
              !lnL(j,k,r) = sum(exp(lnL_sub))

              do l = 1, numband
                 s        = x(j,1) + x(k,2) * a2t(l) * compute_one_comp_dust_spectrum(nu(l), nu(15), beta(1), beta(2)) + x(r,3) * a2t(l) * compute_synch_spectrum(nu(l), nu(1), beta_s) 
                 lnL(j,k,r) = lnL(j,k,r) -0.5d0*((d(l)-s)/rms(l))**2
              end do
              lnL(j,k,r) = exp(lnL(j,k,r))

           end do           
        end do
     end do

     marg = 0.d0
     do j = 1, n
        marg(j) = sum(lnL(j,:,:))
     end do
     marg = marg / (sum(marg) * dA(1))

     mu      = sum(marg * x(:,1) * dA(1))
     sigma   = sqrt(sum(marg * (x(:,1)-mu)**2 * dA(1)))
     deg_fac = sigma / rms(1)
     summary(i,:) = [mu, sigma, sigma*pixsize]
     write(*,*) i, real(mu,sp), real(sigma,sp), real(sigma*pixsize,sp)
     write(59,*) i, real(mu,sp), real(sigma,sp), real(sigma*pixsize,sp)

     if (i == 1) then
        open(58,file='cmb_amp.dat')
        do j = 1, n
           write(58,*) x(j,1), marg(j)
        end do
        close(58)
     end if
     

  end do
  close(59)

  write(*,*) 'CMB amplitude mean     = ', mean(summary(:,1))
  write(*,*) 'CMB amplitude RMS      = ', mean(summary(:,2))
  write(*,*) 'RMS degradation factor = ', mean(summary(:,3))


contains

  function compute_ant2thermo_single(frequency)
    implicit none

    real(dp), intent(in)  :: frequency
    real(dp)              :: compute_ant2thermo_single

    real(dp)     :: x

    x = h*frequency / (k_B*T_0)
    compute_ant2thermo_single = (exp(x)-1.d0)**2 / (x**2 * exp(x))
    
  end function compute_ant2thermo_single

  function compute_one_comp_dust_spectrum(nu, nu_ref, alpha, T_d)
    implicit none
    
    real(dp),     intent(in)  :: nu, nu_ref, alpha, T_d
    real(dp)                  :: compute_one_comp_dust_spectrum

    real(dp) :: x

    x = h / (k_B*T_d)
    compute_one_comp_dust_spectrum = &
         & (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu / nu_ref)**(alpha+1.d0)

  end function compute_one_comp_dust_spectrum

  function compute_synch_spectrum(nu, nu_ref, beta)
    implicit none
    
    real(dp),     intent(in)  :: nu, nu_ref, beta
    real(dp)                  :: compute_synch_spectrum

    real(dp) :: x

    compute_synch_spectrum = (nu/nu_ref)**beta

  end function compute_synch_spectrum

end program fglike
