module dipmod_mcmc
  use healpix_types
  use dipmod_utils
  implicit none

contains

  subroutine generate_chain(chain, unit, rng_handle, map, templates_in, beam, mask, sigma_n, &
       & numsamp, fid_spectrum, l_pivot, prior_f, f_rms, prior_q, q_rms, prior_n, n_rms, &
       & theta_rms, pls)
    implicit none

    integer(i4b),                      intent(in)    :: chain, numsamp, unit, l_pivot
    real(dp),                          intent(in)    :: sigma_n, prior_f, f_rms, theta_rms
    real(dp),                          intent(in)    :: prior_q, q_rms, prior_n, n_rms
    real(sp),     dimension(0:,1:),    intent(in)    :: map, templates_in
    real(dp),     dimension(0:),       intent(in)    :: beam, fid_spectrum
    logical(lgt), dimension(0:),       intent(in)    :: mask
    real(dp),     dimension(0:,1:),    intent(in)    :: pls
    type(planck_rng),                  intent(inout) :: rng_handle

    integer(i4b)     :: i, j, k, l, m, ii, jj, nside, npix, numval, nbands, numdata, lmax, info
    integer(i4b)     :: numbin, numpar, numaccept, numband, numtemp
    logical(lgt)     :: exist, sample_diagonally, accept
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, ln_detC, scale_factor
    real(dp)         :: fg_marginalization, q_prop, n_prop, rms_costheta, rms_phi, q, n
    real(dp)         :: cl_min, cl_max, d_cl, num_cl_bin
    character(len=1) :: uplo
    character(len=2) :: chain_text
    character(len=128) :: filename

    integer(i4b), allocatable, dimension(:)       :: mask2map
    real(dp),     allocatable, dimension(:)       :: f_n, f_par, f_par_prop, params, params_prop
    real(dp),     allocatable, dimension(:)       :: cls, cb, cb_prop, eta
    real(dp),     allocatable, dimension(:,:)     :: data, x, templates
    real(dp),     allocatable, dimension(:,:)     :: covar_mat, L_param, covar_temp, covar_noise
    real(dp),     allocatable, dimension(:,:)     :: trig_angles, vectors
    real(sp),     allocatable, dimension(:,:)     :: f_map
    real(dp),                  dimension(3)       :: vector
    real(dp),                  dimension(3,3)     :: R

    npix    = size(map)
    nside   = nint(sqrt(real(npix,sp)/12.))
    lmax    = size(fid_spectrum)-1
    numband = 1
    numtemp = size(templates_in(1,:))

    scale_factor = 0.2d0
    fg_marginalization = 1.d3

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
    allocate(templates(numdata,numtemp))
    allocate(covar_mat(numdata,numdata))

    allocate(cls(0:lmax))
    do j = 1, numval
       data(j,1) = map(mask2map(j),1)
       do i = 1, numtemp
          templates(j,i) = templates_in(mask2map(j),i)
       end do
    end do

    allocate(covar_temp(numdata, numdata))
    call compute_covar_temp(trig_angles, numtemp, templates, fg_marginalization, covar_temp)

!    allocate(covar_noise(0:npix-1, 0:npix-1))
!    open(58,file='smoothed_WMAP3_w1to4_rms.unf', form='unformatted')
!    read(58) covar_noise
!    write(*,*) minval(covar_noise), maxval(covar_noise)
!    write(*,*) minval(covar_temp), maxval(covar_temp)
!    stop
!    covar_noise = 0.d0
!    close(58)

!    do i = 1, numdata
!       do j = 1, numdata
!          covar_temp(i,j) = covar_temp(i,j) + covar_noise(mask2map(i),mask2map(j))
!       end do
!    end do
!    deallocate(covar_noise)
    

    q = 0.9d0
    n = 0.1d0
!    q = 1.d0
!    n = 0.d0

    cls = 0.d0
    do l = 2, lmax
       cls(l) = q * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n * beam(l)**2
    end do
    

    allocate(f_n(numdata))
    allocate(f_par(3))
    allocate(f_par_prop(3))

    f_par(1) = -0.5379478 !2.d0 * rand_uni(rng_handle) - 1.d0
    f_par(2) =  3.149815 !2.d0*pi * rand_uni(rng_handle)
    f_par(3) =  0.1397533 !prior_f * rand_uni(rng_handle)
    call f_par2f_n(vectors, f_par, f_n)

    call compute_invC_detC(mask, f_n, sigma_n, cls, pls, covar_temp, covar_mat, ln_detC)
    call compute_likelihood(data, covar_mat, ln_detC, lnL)

    inquire(file='param_covar_mat.unf',exist=exist)
    if (exist) then

       numpar = 5
       
       allocate(params(numpar))
       allocate(params_prop(numpar))
       allocate(eta(numpar))

       sample_diagonally = .false.
       allocate(L_param(numpar, numpar))
       open(unit,file='param_covar_mat.unf', form='unformatted')
       read(unit) L_param
       close(unit)

       params(1) = 2.d0 * rand_uni(rng_handle) - 1.d0
       params(2) = 2.d0*pi * rand_uni(rng_handle)
       params(3) = prior_f * rand_uni(rng_handle)
       params(4) = 1.d0
       params(5) = 0.d0

       sample_diagonally = .false.
       
    else

       sample_diagonally = .true.

    end if


    ! Run the chains
    call int2string(chain, chain_text)
    filename = 'chain_no' // chain_text // '.dat'
    open(unit,file=trim(filename), recl=1024)
    numaccept = 0
    do i = 1, numsamp

       ! Draw proposal
       if (sample_diagonally) then

          vector(1) = 0.d0
          vector(2) = 0.d0
          vector(3) = 1.d0

          theta = rand_uni(rng_handle) * theta_rms
          phi   = rand_uni(rng_handle) * 2.d0*pi
          call compute_euler_matrix(0.d0, theta, phi, R)

          vector = matmul(R, vector)

          call compute_euler_matrix_zyz(f_par(2), acos(f_par(1)), 0.d0, R)
          vector = matmul(R, vector)          
          
          call vec2ang(vector, f_par_prop(1), f_par_prop(2))
          f_par_prop(1) = cos(f_par_prop(1))

          f_par_prop(3) = f_par(3) + rand_gauss(rng_handle) * f_rms
          
          q_prop = q + rand_gauss(rng_handle) * q_rms
          n_prop = n + rand_gauss(rng_handle) * n_rms
          
       else

          do j = 1, numpar
             eta(j) = rand_gauss(rng_handle)
          end do

          params_prop = params + scale_factor * matmul(L_param, eta)

          f_par_prop(1:3) = params_prop(1:3)
          q_prop          = params_prop(4)
          n_prop          = params_prop(5)

!          q_prop = q + rand_gauss(rng_handle) * q_rms
!          n_prop = n + rand_gauss(rng_handle) * n_rms

       end if
       
       
       if (f_par_prop(1) > -1.d0 .and. f_par_prop(1) < 1.d0 .and. f_par_prop(3) > 0.d0) then
          
          cls = 0.d0
          do l = 2, lmax
             cls(l) = q_prop * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n_prop * beam(l)**2
          end do
                    
          call f_par2f_n(vectors, f_par_prop, f_n)

          call compute_invC_detC(mask, f_n, sigma_n, cls, pls, covar_temp, covar_mat, ln_detC)
          

          ! Compute likelihood
          call compute_likelihood(data, covar_mat, ln_detC, lnL_prop)
          
          
          ! Apply Metropolis-Hastings rule
          if (lnL_prop > lnL) then
             accept = .true.
          else
             if (rand_uni(rng_handle) < exp(lnL_prop-lnL)) then
                accept = .true.
             else
                accept = .false.
             end if
          end if
          
          if (accept) then
             numaccept = numaccept + 1
             
             q      = q_prop
             n      = n_prop
             f_par  = f_par_prop
             lnL    = lnL_prop

             if (.not. sample_diagonally) then
                params(1:3) = f_par_prop
                params(4)   = q_prop
                params(5)   = n_prop
             end if
             
          end if
          
       end if
       
       ! Output sample
       write(unit,*) i, real(f_par,sp), real(q,sp), real(n,sp), lnL
       
       if (mod(i,10) == 0 .and. i > 0) &
            & write(*,*) 'Acceptance ratio = ', real(numaccept,sp)/real(i,sp)
       
    end do
    close(unit)

    deallocate(mask2map)
    deallocate(cls)
    deallocate(data)
    deallocate(covar_mat)
    deallocate(covar_temp)

  end subroutine generate_chain



end module dipmod_mcmc
