module modulation_mcmc
  use healpix_types
  use modulation_utils
  implicit none

contains

  subroutine evaluate_single_point(map, beam, sigma_n, fid_spectrum, l_pivot, l_mod, theta_s, phi_s, A_s, q_s, n_s, pls)
    implicit none

    integer(i4b),                      intent(in)    :: l_pivot, l_mod
    real(dp),                          intent(in)    :: sigma_n
    real(dp),                          intent(in)    :: theta_s, phi_s, A_s, q_s, n_s
    real(dp),     dimension(0:,1:),    intent(in)    :: map
    real(dp),     dimension(0:),       intent(in)    :: beam, fid_spectrum
    real(dp),     dimension(0:,1:),    intent(in)    :: pls

    integer(i4b)     :: i, j, k, l, m, ii, jj, lmax, info
    integer(i4b)     :: numbin, numpar, numaccept, numtemp
    logical(lgt)     :: exist, sample_diagonally, accept
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, ln_detC, scale_factor
    real(dp)         :: fg_marginalization, q_prop, n_prop, rms_costheta, rms_phi, q, n
    real(dp)         :: cl_min, cl_max, d_cl, num_cl_bin, chisq, chisq_prop
    character(len=1) :: uplo
    character(len=2) :: chain_text
    character(len=128) :: filename

    real(dp),     allocatable, dimension(:)       :: f_n, f_par, f_par_prop, params, params_prop
    real(dp),     allocatable, dimension(:)       :: cls, cb, cb_prop, eta
    real(dp),     allocatable, dimension(:,:)     :: data, x
    real(dp),     allocatable, dimension(:,:)     :: covar_mat, L_param, L_buffer, covar_temp, covar_noise
    real(dp),     allocatable, dimension(:,:)     :: vectors
    real(dp),     allocatable, dimension(:,:)     :: f_map
    real(dp),                  dimension(3)       :: vector
    real(dp),                  dimension(3,3)     :: R

    lmax    = size(fid_spectrum)-1

    ! Pre-compute cos(theta) etc.
    allocate(vectors(3,numdata))
    do i = 1, numdata
       call pix2vec_ring(nside, mask2map(i), vectors(:,i))
    end do


    ! Organize and initialize data
    allocate(data(numdata,1))
    allocate(covar_mat(numdata,numdata))

    allocate(cls(0:lmax))
    do j = 1, numdata
       data(j,1) = map(mask2map(j),1)
    end do

    allocate(covar_temp(numdata, numdata))
    covar_temp = 0.d0

    cls = 0.d0
    do l = 2, lmax
       cls(l) = q_s * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n_s * beam(l)**2
    end do
    
    call compute_high_l_covar(l_mod, cls, pls, covar_temp)

    allocate(f_n(numdata))
    allocate(f_par(3))

    f_par(1) = cos(max(min(theta_s,pi),0.d0))
    f_par(2) = phi_s
    f_par(3) = A_s
    call f_par2f_n(vectors, f_par, f_n)

    call compute_invC_detC(f_n, sigma_n, l_mod, cls, pls, covar_temp, covar_mat, ln_detC)
    call compute_likelihood(data, covar_mat, ln_detC, lnL, chisq)

    write(*,*) '     '
    write(*,*) '     Amplitude      = ', A_s
    write(*,*) '     Latitude       = ', real((0.5d0*pi-theta_s)*180.d0/pi,sp)
    write(*,*) '     Longitude      = ', real(phi_s*180.d0/pi,sp)
    write(*,*) '     q              = ', real(q_s,sp)
    write(*,*) '     n              = ', real(n_s,sp)
    write(*,*) '     '
    write(*,*) '     Log-likelihood = ', lnL
    write(*,*) '     Chi-square     = ', chisq
    write(*,*) '     Npix           = ', numdata
    write(*,*) '     Reduced chisq  = ', (chisq-numdata) / sqrt(2.d0*numdata)
    write(*,*) ''

    deallocate(data, covar_mat, f_n, f_par)

  end subroutine evaluate_single_point

  subroutine generate_chain(prefix, chain, unit, rng_handle, map, beam, sigma_n, &
       & numsamp, fid_spectrum, l_pivot, l_mod, prior_f, f_rms, prior_q, q_rms, prior_n, n_rms, &
       & theta_rms, pls)
    implicit none

    character(len=*),                  intent(in)    :: prefix
    integer(i4b),                      intent(in)    :: chain, numsamp, unit, l_pivot, l_mod
    real(dp),                          intent(in)    :: sigma_n, prior_f, f_rms, theta_rms
    real(dp),                          intent(in)    :: prior_q, q_rms, prior_n, n_rms
    real(dp),     dimension(0:,1:),    intent(in)    :: map
    real(dp),     dimension(0:),       intent(in)    :: beam, fid_spectrum
    real(dp),     dimension(0:,1:),    intent(in)    :: pls
    type(planck_rng),                  intent(inout) :: rng_handle

    integer(i4b)     :: i, j, k, l, m, ii, jj, lmax, info
    integer(i4b)     :: numbin, numpar, numaccept, numtemp
    logical(lgt)     :: exist, sample_diagonally, accept
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, ln_detC, scale_factor
    real(dp)         :: fg_marginalization, q_prop, n_prop, rms_costheta, rms_phi, q, n
    real(dp)         :: cl_min, cl_max, d_cl, num_cl_bin, chisq, chisq_prop
    character(len=1) :: uplo
    character(len=2) :: chain_text
    character(len=128) :: filename

    real(dp),     allocatable, dimension(:)       :: f_n, f_par, f_par_prop, params, params_prop
    real(dp),     allocatable, dimension(:)       :: cls, cb, cb_prop, eta
    real(dp),     allocatable, dimension(:,:)     :: data, x
    real(dp),     allocatable, dimension(:,:)     :: covar_mat, L_param, L_buffer, covar_temp, covar_noise
    real(dp),     allocatable, dimension(:,:)     :: vectors
    real(dp),     allocatable, dimension(:,:)     :: f_map
    real(dp),                  dimension(3)       :: vector
    real(dp),                  dimension(3,3)     :: R

    lmax    = size(fid_spectrum)-1

    scale_factor = 0.5d0
    fg_marginalization = 1.d3

    ! Pre-compute cos(theta) etc.
    allocate(vectors(3,numdata))
    do i = 1, numdata
       call pix2vec_ring(nside, mask2map(i), vectors(:,i))
    end do


    ! Organize and initialize data
    allocate(data(numdata,1))
    allocate(covar_mat(numdata,numdata))

    allocate(cls(0:lmax))
    do j = 1, numdata
       data(j,1) = map(mask2map(j),1)
    end do

    allocate(covar_temp(numdata, numdata))
    covar_temp = 0.d0

    q = 1.d0
    n = 0.d0

    cls = 0.d0
    do l = 2, lmax
       cls(l) = q * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n * beam(l)**2
    end do

    write(*,*) real(cls,sp)
    
    call compute_high_l_covar(l_mod, cls, pls, covar_temp)

    allocate(f_n(numdata))
    allocate(f_par(3))
    allocate(f_par_prop(3))

    f_par(1) = 2.d0 * rand_uni(rng_handle) - 1.d0
    f_par(2) = 2.d0*pi * rand_uni(rng_handle)
    f_par(3) = prior_f * rand_uni(rng_handle)
    call f_par2f_n(vectors, f_par, f_n)

    call compute_invC_detC(f_n, sigma_n, l_mod, cls, pls, covar_temp, covar_mat, ln_detC)
    call compute_likelihood(data, covar_mat, ln_detC, lnL, chisq)

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

       allocate(L_buffer(5,5))
       L_buffer = matmul(L_param, transpose(L_param))
       L_buffer(1:2,:) = 0.d0
       L_buffer(:,1:2) = 0.d0
       L_buffer(1,1)   = 1.d0
       L_buffer(2,2)   = 1.d0
       call cholesky_decompose(L_buffer, L_param)
       L_param(1:2,:) = 0.d0
       L_param(:,1:2) = 0.d0
       deallocate(L_buffer)

       params(1) = f_par(1)
       params(2) = f_par(2)
       params(3) = f_par(3)
       params(4) = 1.d0
       params(5) = 0.d0

       sample_diagonally = .false.
       
    else

       sample_diagonally = .true.

    end if


    ! Run the chains
    call int2string(chain, chain_text)
    filename = 'chain_'//trim(prefix)//'_no' // chain_text // '.dat'
    open(unit,file=trim(filename), recl=1024)
    numaccept = 0
    do i = 1, numsamp

       write(*,*) 'Chain = ', chain, ' -- Generating sample no. ', i, ' of ', numsamp

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

          ! Draw (A,n,alpha)
          do j = 1, numpar
             eta(j) = rand_gauss(rng_handle)
          end do
          params_prop = params + scale_factor * matmul(L_param, eta)

          ! Draw direction
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


          f_par_prop(3)   = params_prop(3)
          q_prop          = params_prop(4)
          n_prop          = params_prop(5)

       end if
       
       
       if (f_par_prop(1) > -1.d0 .and. f_par_prop(1) < 1.d0 .and. f_par_prop(3) > 0.d0) then
          
          cls = 0.d0
          do l = 2, lmax
             cls(l) = q_prop * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n_prop * beam(l)**2
          end do
                    
          call f_par2f_n(vectors, f_par_prop, f_n)

          call compute_invC_detC(f_n, sigma_n, l_mod, cls, pls, covar_temp, covar_mat, ln_detC)
          

          ! Compute likelihood
          call compute_likelihood(data, covar_mat, ln_detC, lnL_prop, chisq_prop)
          
          
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
             chisq  = chisq_prop

             if (.not. sample_diagonally) then
                params(1:3) = f_par_prop
                params(4)   = q_prop
                params(5)   = n_prop
             end if
             
          end if
          
       end if
       
       ! Output sample
       write(unit,*) i, real(f_par,sp), real(q,sp), real(n,sp), lnL, (chisq-numdata)/sqrt(2.d0*numdata)
       
       if (mod(i,10) == 0 .and. i > 0) &
            & write(*,*) 'Acceptance ratio = ', real(numaccept,sp)/real(i,sp)
       
    end do
    close(unit)

    deallocate(cls)
    deallocate(data)
    deallocate(covar_mat)
    deallocate(covar_temp)

  end subroutine generate_chain



end module modulation_mcmc
