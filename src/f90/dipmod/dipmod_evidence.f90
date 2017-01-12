module dipmod_evidence
  use dipmod_utils
  implicit none

contains

  subroutine compute_evidence(chain, unit, rng_handle, map, templates_in, beam, & 
       & mask, sigma_n, numsamp, fid_spectrum, l_pivot, prior_f, prior_q, &
       & prior_n, pls, evidence)
    implicit none

    integer(i4b),                      intent(in)    :: chain, numsamp, unit, l_pivot
    real(dp),                          intent(in)    :: sigma_n, prior_f, prior_q, prior_n
    real(sp),     dimension(0:,1:),    intent(in)    :: map, templates_in
    real(dp),     dimension(0:),       intent(in)    :: beam, fid_spectrum
    logical(lgt), dimension(0:),       intent(in)    :: mask
    real(dp),     dimension(0:,1:),    intent(in)    :: pls
    type(planck_rng),                  intent(inout) :: rng_handle
    real(dp),                          intent(out)   :: evidence

    integer(i4b)     :: i, j, k, l, m, ii, jj, nside, npix, numval, nbands, numdata, lmax, info
    integer(i4b)     :: numpar, numaccept, numstep, jm1
    integer(i4b)     :: ldb, nrhs, n, lda, numband, i2, numtemp
    logical(lgt)     :: exist, sample_diagonally, accept, reject_sample, outside_prior
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, ln_detC, scale_factor
    real(dp)         :: X_j, X_jm1, L_j, L_jm1, lnL_0, avg_L, X_old, t, fg_marginalization
    real(dp)         :: abs_eta, n_prop, q_prop
    character(len=1) :: uplo
    character(len=5) :: chain_text
    character(len=128) :: filename

    integer(i4b), allocatable, dimension(:)       :: mask2map
    real(dp),     allocatable, dimension(:)       :: f_n, f_par, f_par_prop, params, params_prop
    real(dp),     allocatable, dimension(:)       :: cls, cb, cb_prop, eta, L_x, centroid
    real(dp),     allocatable, dimension(:,:)     :: data, covar_temp, templates
    real(dp),     allocatable, dimension(:,:)     :: covar_mat, L_param
    real(dp),     allocatable, dimension(:,:)     :: vectors, trig_angles, priors, samples, dE

    npix    = size(map)
    nside   = nint(sqrt(real(npix,sp)/12.))
    lmax    = size(fid_spectrum)-1
    numband = 1
    numtemp = size(templates_in(0,:))

    fg_marginalization = 1.d3

    n            = 30
    scale_factor = 1.5d0
    numstep      = 1000
    t            = real(n,dp) / real(n+1,dp)

    numval = 0
    do i = 0, npix-1
       if (mask(i)) numval = numval+1
    end do
    numdata = numval * numband

    allocate(mask2map(numval))
    j = 1
    do i = 0, npix-1
       if (mask(i)) then
          mask2map(j) = i
          j           = j+1
       end if
    end do

    ! Pre-compute cos(theta) etc.
    allocate(vectors(3,numval))
    allocate(trig_angles(4,numval))
    do i = 1, numval
       call pix2vec_ring(nside, mask2map(i), vectors(:,i))

       call pix2ang_ring(nside, mask2map(i), theta, phi)
       trig_angles(1,i) = cos(theta)
       trig_angles(2,i) = sin(theta)
       trig_angles(3,i) = cos(phi)
       trig_angles(4,i) = sin(phi)
   end do

    ! Organize and initialize data
    allocate(data(numdata,1))
    allocate(covar_mat(numdata,numdata))
    allocate(templates(numdata,numtemp))

    allocate(cls(0:lmax))
    do j = 1, numval
       data(j,1) = map(mask2map(j),1)
       do i = 1, numtemp
          templates(j,i) = templates_in(mask2map(j),i)
       end do
    end do

    allocate(covar_temp(numdata, numdata))
    call compute_covar_temp(trig_angles, numtemp, templates, fg_marginalization, covar_temp)


    cls = fid_spectrum

    allocate(f_n(numdata))
    allocate(f_par(3))
    allocate(f_par_prop(3))

    f_par = 0.d0
    f_n  = 0.d0

    numpar = 5

    allocate(priors(numpar,2))

    priors(1,1)   =  -1.d0
    priors(1,2)   =  1.d0
    priors(2,1)   =  0.d0
    priors(2,2)   =  2.d0*pi
    priors(3,1)   =  0.d0
    priors(3,2)   =  prior_f
    priors(4,1)   = 1.d0 - prior_q
    priors(4,2)   = 1.d0 + prior_q
    priors(5,1)   = -prior_n
    priors(5,2)   =  prior_n

    ! Initialize the sample array
    allocate(samples(n, numpar))
    allocate(L_x(n))
    allocate(L_param(numpar,numpar))
    allocate(params_prop(numpar))
    allocate(centroid(numpar))
    allocate(eta(numpar))

    do i = 1, n

       write(*,*) 'Initializing sample no. ', i, ' of ', n

       do j = 1, numpar
          samples(i,j) = priors(j,1) + rand_uni(rng_handle) * (priors(j,2)-priors(j,1))
       end do

       f_par_prop = samples(i,1:3)

       q_prop = samples(i,4)
       n_prop = samples(i,5)

       cls = 0.d0
       do l = 2, lmax
          cls(l) = q_prop * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n_prop * beam(l)**2
       end do

       call f_par2f_n(vectors, f_par_prop, f_n)

       call compute_invC_detC(mask, f_n, sigma_n, cls, pls, covar_temp, covar_mat, ln_detC)
       
       ! Compute likelihood
       call compute_likelihood(data, covar_mat, ln_detC, L_x(i))

       write(*,*) real(f_par_prop,sp), real(q_prop,sp), real(n_prop,sp)
       write(*,*) i, L_x(i)

    end do

    L_j = 1.d30
    do i = 1, n
       if (L_x(i) < L_j) then
          L_j = L_x(i)
          j   = i
       end if
    end do

    ! Compute evidence 
    allocate(dE(numstep,2))
    dE    = 0.d0
    X_old = 1.d0


    call int2string(chain, chain_text)
    filename = 'chain_no' // chain_text // '.dat'
    open(unit,file=trim(filename), recl=1024)
    do i = 2, numstep

       do k = 1, numpar
          centroid(k) = sum(samples(:,k)) / real(n,dp)
       end do

       call compute_covar_mat(samples, L_param)
       abs_eta = 0.d0
       do k = 1, n
          call cholesky_solve(L_param, samples(k,:)-centroid, eta)
          abs_eta = max(abs_eta, sqrt(sum(eta*eta)))
       end do
       L_param = L_param * abs_eta
       L_param = L_param * scale_factor


!!$       do k = 1, numpar
!!$          priors(k,1) = minval(samples(:,k))
!!$          priors(k,2) = maxval(samples(:,k))
!!$       end do
!!$       priors = priors * scale_factor


       ! Generate new sample
       reject_sample = .true.
       do while (reject_sample)

          eta = 1000.d0
          do while (sqrt(sum(eta*eta)) > 1.d0) 
             do k = 1, numpar
                if (priors(k,1) < priors(k,2)) then
                   eta(k) = 2.d0*rand_uni(rng_handle) - 1.d0
                else
                   eta(k) = 0.d0
                end if
             end do
          end do

          params_prop = centroid + matmul(L_param, eta)

          outside_prior = .false.
          do k = 1, numpar
             if (priors(k,1) < priors(k,2)) then
                if (params_prop(k) < priors(k,1)) outside_prior = .true.
                if (params_prop(k) > priors(k,2)) outside_prior = .true.
             else
                params_prop(k) = priors(k,1)
             end if
          end do

          if (outside_prior) cycle

          f_par_prop     = params_prop(1:3)

          q_prop = params_prop(4)
          n_prop = params_prop(5)
          
          cls = 0.d0
          do l = 2, lmax
             cls(l) = q_prop * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n_prop * &
                  & beam(l)**2
          end do
          
          
          call f_par2f_n(vectors, f_par_prop, f_n)

          call compute_invC_detC(mask, f_n, sigma_n, cls, pls, covar_temp, covar_mat, ln_detC)
          
          ! Compute likelihood
          call compute_likelihood(data, covar_mat, ln_detC, lnL_prop)

          if (lnL_prop > L_x(j)) then

             if (chain == 1) write(*,*) 'Accepted -- lnL_prop = ', real(lnL_prop,sp)

             reject_sample = .false.

          else

             if (chain == 1) write(*,*) 'Rejected -- lnL_prop = ', real(lnL_prop,sp)

!             write(*,*) 'Rejected proposal -- lnL_prop = ', real(lnL_prop,sp), &
!                  & ', lnL_min = ', real(L_x(j))

!             write(*,*) L_x

          end if

       end do

       ! Store results
       dE(i,1) = L_x(j)
       dE(i,2) = X_old * t

       write(unit,*) i, real(samples(j,:),sp), real(dE(i,:),sp)

       ! Replace lowest point with proposed point
       samples(j,:) = params_prop
       L_x(j) = lnL_prop
       X_old  = X_old * t

       ! Find the new lowest point
       L_j = 1.d30
       do k = 1, n
          if (L_x(k) < L_j) then
             L_j = L_x(k)
             j   = k
          end if
       end do

    end do
    close(unit)

    lnL_0 = maxval(dE(:,1))

    dE(:,1) = dE(:,1) - lnL_0

    evidence = 0.d0
    do i = 1, numstep
       if (dE(i,2) > 0.d0) then
          evidence = evidence + exp(dE(i,1)) * dE(i,2)
       end if
    end do

    ! Add the remainder
    L_x = L_x - lnL_0
    L_x = exp(L_x)
    avg_L = sum(L_x) / real(n,dp)

    evidence = evidence + avg_L * X_old

    evidence = log(evidence) + lnL_0

    write(*,*) 'Final log(E) = ', evidence

    deallocate(mask2map)
    deallocate(cls)
    deallocate(data)
    deallocate(covar_mat)

  end subroutine compute_evidence

  subroutine compute_covar_mat(samples, L)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: samples
    real(dp), dimension(1:,1:), intent(out) :: L

    integer(i4b) :: i, j, numpar, n
    real(dp), allocatable, dimension(:)   :: avg
    real(dp), allocatable, dimension(:,:) :: covar

    numpar = size(L(:,1))
    n      = size(samples(:,1))

    ! Compute covariance matrix
    allocate(avg(numpar))
    allocate(covar(numpar, numpar))
  
    do i = 1, numpar
       avg(i) = sum(samples(:,i)) / real(n,dp)
    end do

    do i = 1, numpar
       do j = 1, numpar
          covar(i,j) = sum((samples(:,i)-avg(i)) * (samples(:,j)-avg(j)))
          covar(i,j) = covar(i,j) / real(n-1,dp)
          covar(j,i) = covar(i,j)
       end do
       
       if (covar(i,i) == 0.d0) then
          covar(:,i) = 0.d0
          covar(i,:) = 0.d0
          covar(i,i) = 1.d0
       end if
!       write(*,*) i, real(covar(i,:),sp)
    end do

    call cholesky_decompose(covar, L)
  
    deallocate(avg)
    deallocate(covar)

  end subroutine compute_covar_mat

end module dipmod_evidence
