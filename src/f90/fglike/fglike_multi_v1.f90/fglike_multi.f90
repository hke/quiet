program fglike
  use healpix_types
  use rngmod
  use quiet_utils
  implicit none

  include "mpif.h"


  real(dp), parameter :: k_B     = 1.3806503d-23
  real(dp), parameter :: h       = 6.626068d-34
  real(dp), parameter :: c       = 2.99792458d8
  real(dp), parameter :: T_0     = 2.726d0

  integer(i4b) :: i, j, k, l, m, n, q, numband, numcomp, numsim, seed, num_setup, numsamp, burnin, numpar, root
  integer(i4b) :: num_delta_nu, accept, status
  logical(lgt) :: output_chain, log_delta_nu
  real(dp)     :: A_cmb, A_dust, min_cmb, max_cmb, min, fg_min, T_d, alpha, rms0, nu_rms, nu0, lnL0, lnL_prop
  real(dp)     :: A_in(2), beta(2), dbeta(2), s, mu, sigma, deg_fac, pixsize, polfactor, npix, sky_area, nu_cut_low, nu_cut_high
  real(dp)     :: nu_info(3), delta_nu_info(3), delta_nu, t1, t2, ndet, fpa, fpa_ref, rms_per_det, sig, sig_ref
  real(dp)     :: prior_synch_bias, prior_dust_bias, sigma1, sigma2, norm, norm_ref
  real(dp), allocatable, dimension(:,:)  :: lnL, x, lnL_sub, summary, samples, L_prop, prior, F, cov
  real(dp), allocatable, dimension(:)    :: marg, d, a2t, rms, nu, parmean
  real(dp), allocatable, dimension(:)    :: A_ref, beta_ref, nu_ref, degfac, p0, p_prop, eta
  character(len=256) :: parfile, noise_model, outfile
  type(planck_rng) :: handle, handle_sim
  integer(i4b) :: ierr, myid, numprocs

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0

  call getarg(1, parfile)
  call get_parameter(58, parfile, 'NUMBAND',    par_int=numband)
  call get_parameter(58, parfile, 'MIN_DELTA_NU', par_dp=delta_nu_info(1))
  call get_parameter(58, parfile, 'MAX_DELTA_NU', par_dp=delta_nu_info(2))
  call get_parameter(58, parfile, 'NUM_DELTA_NU', par_int=num_delta_nu)
  call get_parameter(58, parfile, 'LOG_DELTA_NU', par_lgt=log_delta_nu)
  call get_parameter(58, parfile, 'NUMSIM',     par_int=numsim)
  call get_parameter(58, parfile, 'SEED',       par_int=seed)
  call get_parameter(58, parfile, 'NU_MIN',     par_dp=nu_info(1))
  call get_parameter(58, parfile, 'NU_MAX',     par_dp=nu_info(2))
  call get_parameter(58, parfile, 'NUM_SETUP',  par_int=num_setup)
  call rand_init(handle, seed+100*myid)

  numcomp = 3
  allocate(A_ref(numcomp), beta_ref(3), nu_ref(numcomp), prior(2,2))
  call get_parameter(58, parfile, 'REFERENCE_CMB_AMP',    par_dp=A_ref(1))
  call get_parameter(58, parfile, 'REFERENCE_SYNCH_AMP',  par_dp=A_ref(2))
  call get_parameter(58, parfile, 'REFERENCE_DUST_AMP',   par_dp=A_ref(3))
  call get_parameter(58, parfile, 'REFERENCE_SYNCH_BETA', par_dp=beta_ref(1))
  call get_parameter(58, parfile, 'REFERENCE_DUST_BETA',  par_dp=beta_ref(2))
  call get_parameter(58, parfile, 'REFERENCE_DUST_T',     par_dp=beta_ref(3))
  nu_ref(1)  = 100.d0
  call get_parameter(58, parfile, 'REFERENCE_SYNCH_NU0',  par_dp=nu_ref(2))
  call get_parameter(58, parfile, 'REFERENCE_DUST_NU0',   par_dp=nu_ref(3))
  call get_parameter(58, parfile, 'PRIOR_SYNCH_BETA_BIAS',  par_dp=prior_synch_bias)
  call get_parameter(58, parfile, 'PRIOR_SYNCH_BETA_RMS', par_dp=prior(1,2))
  call get_parameter(58, parfile, 'PRIOR_DUST_BETA_RMS',  par_dp=prior(2,2))
  call get_parameter(58, parfile, 'PRIOR_DUST_BETA_BIAS',  par_dp=prior_dust_bias)
  call get_parameter(58, parfile, 'REFERENCE_NOISE_RMS',  par_dp=rms0)
  call get_parameter(58, parfile, 'NOISE_MODEL',          par_string=noise_model)
  call get_parameter(58, parfile, 'NUM_MCMC_STEP',        par_int=numsamp)
  call get_parameter(58, parfile, 'OUTFILE',              par_string=outfile)
  nu_ref     = nu_ref * 1d9
  prior(1,1) = beta_ref(1)+prior_synch_bias
  prior(2,1) = beta_ref(2)+prior_dust_bias
  burnin = 150000
  nu_rms = 100.d9   ! Arbitrary
  output_chain = .false.
  nu_cut_high  = 2000.d9 ! No information above 2 THz
  !if (myid == 0) output_chain = .true.
  
  numpar = 5
  allocate(d(numband), rms(numband), a2t(numband), nu(numband), parmean(numpar))
  allocate(degfac(numsim), p0(numpar), p_prop(numpar), eta(numpar), L_prop(numpar,numpar), cov(numpar,numpar))
  allocate(samples(0:numpar, numsamp), F(numband,numcomp))
  allocate(summary(numsim,3))

  if (myid == 0) open(60,file=trim(outfile),recl=1024)
  do q = 1, num_delta_nu
     if (log_delta_nu) then
        delta_nu = delta_nu_info(1) * (delta_nu_info(2)/delta_nu_info(1))**((q-1.d0)/(num_delta_nu-1.d0))
     else
        delta_nu = delta_nu_info(1) + real(q-1,dp)/real(num_delta_nu-1,dp) * (delta_nu_info(2)-delta_nu_info(1))
     end if
     do i = 1, num_setup

        ! Set up frequencies
        nu0 = nu_info(1) * (nu_info(2)/nu_info(1))**((i-1.d0)/(num_setup-1.d0))
        do j = 1, numband
           nu(j)  = nu0 * delta_nu**(j-1) * 1d9
           a2t(j) = compute_ant2thermo_single(nu(j))
           
           if (nu(j) < nu_cut_high) then
              F(j,1) = 1.d0
              F(j,2) = a2t(j) * compute_synch_spectrum(nu(j), nu_ref(2), beta_ref(1))
              F(j,3) = a2t(j) * compute_one_comp_dust_spectrum(nu(j), nu_ref(3), beta_ref(2), beta_ref(3))
           else
              F(j,:) = 0.d0
           end if
        end do
        
        ! Set up noise rms
        norm     = 0.d0
        fpa      = 0.d0
        fpa_ref  = 0.d0
        do j = 1, numband
           if (trim(noise_model) == 'uniform') then
              rms(j) = 1.d0
              norm = norm + 1.d0/rms(j)**2
           else if (trim(noise_model) == 'prop_to_inv_nu') then
              rms(j) = (nu_rms/nu(j)) 
              norm = norm + 1.d0/rms(j)**2
           else if (trim(noise_model) == 'signal_to_noise') then
              rms(j) = A_ref(1) / a2t(j) + &    ! CMB
                   & A_ref(2) * compute_synch_spectrum(nu(j), nu_ref(2), beta_ref(1)) + &
                   & A_ref(3) * compute_one_comp_dust_spectrum(nu(j), nu_ref(3), beta_ref(2), beta_ref(3)) 
           else if (trim(noise_model) == 'radiometer+focalplane') then
              if (nu(j) > nu_cut_high) then
                 rms(j)   = 1.d30
              else
                 ndet     = 1
                 rms(j)   = rms0 * (nu_rms/nu(j))**0.5d0
                 norm     = norm + ndet * (nu_rms/nu(j))**2
              end if
           else if (trim(noise_model) == 'signal_to_noise_v2') then
              if (nu(j) > nu_cut_high) then
                 rms(j) = 1.d30
              else
                 sig = A_ref(1) / a2t(j) + &    ! CMB
                      & A_ref(2) * compute_synch_spectrum(nu(j), nu_ref(2), beta_ref(1)) + &
                      & A_ref(3) * compute_one_comp_dust_spectrum(nu(j), nu_ref(3), beta_ref(2), beta_ref(3)) 
                 sig_ref = 0.d0*A_ref(1) / compute_ant2thermo_single(nu_rms) + &    ! CMB
                      & A_ref(2) * compute_synch_spectrum(nu_rms, nu_ref(2), beta_ref(1)) + &
                      & A_ref(3) * compute_one_comp_dust_spectrum(nu_rms, nu_ref(3), beta_ref(2), beta_ref(3)) 
                 rms(j)   = rms0 * (nu_rms/nu(j))**0.5d0
                 ndet     = ((sig_ref/rms0) / (sig/rms(j)))**2
                 norm     = norm + ndet * (nu_rms/nu(j))**2
                 rms(j)   = rms(j) / sqrt(ndet)
              end if
           else
              write(*,*) 'Unsupported noise model = ', trim(noise_model)
              stop
           end if
        end do
        if (trim(noise_model) == 'radiometer+focalplane' .or. trim(noise_model) == 'signal_to_noise_v2') then
           rms = rms * sqrt(norm)
           if (myid == 0) then
              write(*,*) 'nu      = ', real(nu*1e-9,sp)
              write(*,*) 'rms_RJ  = ', real(rms,sp)
              write(*,*) 'norm    = ', real(sqrt(norm),sp)
           end if
           where (nu > nu_cut_high) 
              rms = 1.d30
           elsewhere
              rms = rms * a2t
           end where
           if (myid == 0) then
              write(*,*) 'rms_CMB = ', real(rms,sp)
           end if
        else
           norm = sqrt(1.d0/norm) 
           if (myid == 0) then
              write(*,*) 'nu      = ', real(nu,sp)
              write(*,*) 'rms_RJ  = ', real(rms0*rms/norm,sp)
              write(*,*) 'rms_CMB = ', real(rms0*rms/norm*a2t,sp)
           end if
           rms  = rms0 * rms / norm * a2t
        end if
        
        if (.true. .or. i == 1) then
           ! Set up proposal matrix -> Fisher matrix
           L_prop = 0.d0
           do j = 1, numband
              do k = 1, 3
                 do l = 1, 3
                    L_prop(k,l) = L_prop(k,l) + F(j,k) * F(j,l) / rms(j)**2
                 end do
              end do
           end do
           L_prop(4,4) = 1.d0/0.02d0**2
           L_prop(5,5) = 1.d0/0.02d0**2     
           cov = L_prop
           call invert_matrix(cov)
           call cholesky_decompose(cov, L_prop, status)
           L_prop = L_prop * 0.3d0
        end if
        
        degfac  = 0.d0
        summary = 0.d0
        if (status == 0) then
           call wall_time(t1)
           do j = 1+myid, numsim, numprocs
              
              ! Simulate data set
              do k = 1, numband
                 d(k) = A_ref(1) + &    ! CMB
                      & A_ref(2) * a2t(k) * compute_synch_spectrum(nu(k), nu_ref(2), beta_ref(1)) + &
                      & A_ref(3) * a2t(k) * compute_one_comp_dust_spectrum(nu(k), nu_ref(3), beta_ref(2), beta_ref(3)+0.d0) + &
                      & rms(k) * rand_gauss(handle)
              end do
              
              ! Run Markov chain until convergence
              p0     = [A_ref(1), A_ref(2), A_ref(3), beta_ref(1), beta_ref(2)]
              lnL0   = compute_log_like(p0, d, nu, rms, prior)
              accept = 0
              if (output_chain) open(58,file='samples.dat',recl=1024)
              do k = 1, numsamp
                 ! Propose new sample
                 do l = 1, numpar
                    eta(l) = rand_gauss(handle)
                 end do
                 p_prop   = p0 + matmul(L_prop,eta)
                 lnL_prop = compute_log_like(p_prop, d, nu, rms, prior)
                 
                 if (exp(lnL_prop-lnL0) > rand_uni(handle)) then
                    p0   = p_prop
                    lnL0 = lnL_prop
                    accept = accept+1
                 end if
                 
                 samples(0,k)        = lnL0
                 samples(1:numpar,k) = p0
                 
                 if (output_chain) then
                    if (mod(k,100) == 0) write(58,*) k, real(samples(:,k),sp)
                 end if
                 
                 if (k == 50000 .or. k == 100000 .or. k == 150000) then
                    do l = 1, numpar
                       parmean(l) = mean(samples(l,1:k))
                    end do
                    do l = 1, numpar
                       do m = 1, numpar
                          cov(l,m) = mean((samples(l,1:k)-parmean(l))*(samples(m,1:k)-parmean(m)))
                       end do
                    end do
                    call cholesky_decompose(cov, L_prop, status)
                    if (status /= 0) exit
                    L_prop = 0.5d0*L_prop
                 end if
                 
                 if (mod(k,100000) == 0 .and. k > 400000) then
                    ! Check for convergence
                    sigma1 = sqrt(variance(samples(1,burnin:k/2)))
                    sigma2 = sqrt(variance(samples(1,k/2+1:k)))
                    !write(*,*) sigma1, sigma2, abs(sigma1-sigma2)/abs(sigma1+sigma2)
                    if (abs(sigma1-sigma2)/abs(sigma1+sigma2) < 1.d-3) then
                       exit
                    end if
                 end if
              end do
              !write(*,*) 'chain length = ', k, ', accept rate = ', real(accept,sp)/real(k,sp)
              if (output_chain) then
                 close(58)
                 call mpi_finalize(ierr)
                 stop
              end if
              
              mu      = mean(samples(1,burnin:k-1))
              sigma   = sqrt(variance(samples(1,burnin:k-1)))
              deg_fac = sigma / (1.d0/(sqrt(sum(1.d0/rms**2))))
              summary(j,:) = [mu, sigma, deg_fac]
              
           end do
           call mpi_allreduce(MPI_IN_PLACE, summary, size(summary), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
           call mpi_allreduce(MPI_IN_PLACE, status, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
           
           call wall_time(t2)
        end if

        if (status /= 0) then
           summary(:,1) =   0.d0
           summary(:,2) = 100.d0
           summary(:,3) = 100.d0
        end if
           !write(*,*) 'CPU time = ', t2-t1
           
        if (myid == 0) then
           write(*,fmt='(a,f8.3,a,f8.3,a,f8.3,a,f8.3)') 'delta_nu = ', delta_nu, ', nu = ', nu(1)*1d-9, ', cpu = ', t2-t1, &
                & ', rms = ', mean(summary(:,2))
           write(*,*)
           write(60,fmt='(f8.3,f8.3,3f16.8)') delta_nu, nu(1)*1d-9, &
                & mean(summary(:,1)), mean(summary(:,2)), mean(summary(:,3))
           
           
        end if
        !stop

     end do
  end do
  if (myid == 0) close(60)

  !write(*,*) 'CMB amplitude mean     = ', mean(summary(:,1))
  !write(*,*) 'CMB amplitude RMS      = ', mean(summary(:,2))
  !write(*,*) 'RMS degradation factor = ', mean(summary(:,3))

  call mpi_finalize(ierr)

contains

  function compute_log_like(p, d, nu, rms, prior)
    implicit none
    real(dp), dimension(1:),    intent(in)  :: p, d, nu, rms
    real(dp), dimension(1:,1:), intent(in)  :: prior
    real(dp)                                :: compute_log_like

    integer(i4b) :: i, j
    real(dp)     :: s

    ! Compute chi-square term
    compute_log_like = 0.d0
    do i = 1, size(d)
       if (nu(i) < nu_cut_high) then
          s = p(1) + &    ! CMB
               & p(2) * a2t(i) * compute_synch_spectrum(nu(i), nu_ref(2), p(4)) + &
               & p(3) * a2t(i) * compute_one_comp_dust_spectrum(nu(i), nu_ref(3), p(5), beta_ref(3))
          compute_log_like = compute_log_like - 0.5d0 * ((d(i)-s)/rms(i))**2
       end if
    end do

    ! Add prior term
    compute_log_like = compute_log_like - 0.5d0 * ((p(4)-prior(1,1))/prior(1,2))**2 - &
         & 0.5d0 * ((p(5)-prior(2,1))/prior(2,2))**2

  end function compute_log_like

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

    if (nu > nu_cut_high) then
       compute_one_comp_dust_spectrum = 0.d0
    else
       x = h / (k_B*T_d)
       compute_one_comp_dust_spectrum = &
            & (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu / nu_ref)**(alpha+1.d0)
    end if

  end function compute_one_comp_dust_spectrum


  function compute_synch_spectrum(nu, nu_ref, beta)
    implicit none
    
    real(dp),     intent(in)  :: nu, nu_ref, beta
    real(dp)                  :: compute_synch_spectrum

    if (nu > nu_cut_high) then
       compute_synch_spectrum = 0.d0
    else
       compute_synch_spectrum = (nu/nu_ref)**beta
    end if


  end function compute_synch_spectrum

end program fglike
