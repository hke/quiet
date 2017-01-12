module modulation_chisq
  use modulation_utils
  use modulation_evidence
  implicit none

contains

  subroutine compute_chisq_map(unit, map, templates_in, beam, & 
       & sigma_n, fid_spectrum, l_pivot, pls)
    implicit none

    integer(i4b),                      intent(in)    :: unit, l_pivot
    real(dp),                          intent(in)    :: sigma_n
    real(dp),     dimension(0:,1:),    intent(in)    :: map, templates_in
    real(dp),     dimension(0:),       intent(in)    :: beam, fid_spectrum
    real(dp),     dimension(0:,1:),    intent(in)    :: pls

    integer(i4b)     :: i, j, k, l, m, ii, jj, lmax, info
    integer(i4b)     :: numpar, numaccept, numstep, jm1
    integer(i4b)     :: ldb, nrhs, n, lda, i2
    logical(lgt)     :: exist, sample_diagonally, accept, reject_sample, outside_prior
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, ln_detC, scale_factor
    real(dp)         :: X_j, X_jm1, L_j, L_jm1, lnL_0, avg_L, X_old, t, fg_marginalization
    real(dp)         :: abs_eta, n_prop, q_prop
    character(len=1) :: uplo
    character(len=2) :: chain_text
    character(len=128) :: filename

    real(dp),     allocatable, dimension(:)       :: f_n, f_par, f_par_prop, params, params_prop
    real(dp),     allocatable, dimension(:)       :: cls, cb, cb_prop, eta, L_x, centroid
    real(dp),     allocatable, dimension(:,:)     :: data, data2, covar_temp, templates
    real(dp),     allocatable, dimension(:,:)     :: covar_mat, L_param
    real(dp),     allocatable, dimension(:,:)     :: vectors, trig_angles, priors, samples, dE
    real(dp),     allocatable, dimension(:,:)     :: output_map

    lmax    = size(fid_spectrum)-1

    fg_marginalization = 1.d3

    ! Pre-compute cos(theta) etc.
    allocate(vectors(3,numdata))
    allocate(trig_angles(4,numdata))
    do i = 1, numdata
       call pix2vec_ring(nside, mask2map(i), vectors(:,i))
   end do

    ! Organize and initialize data
    allocate(data(numdata,1))
    allocate(data2(numdata,1))
    allocate(covar_mat(numdata,numdata))
    allocate(templates(numdata,numtemp))

    allocate(cls(0:lmax))
    do j = 1, numdata
       data(j,1) = map(mask2map(j),1)
       do i = 1, numtemp
          templates(j,i) = templates_in(mask2map(j),i)
       end do
    end do

    allocate(covar_temp(numdata, numdata))
    covar_temp = 0.d0

    cls      = fid_spectrum
    q_prop   = 1.d0
    n_prop   = 0.d0

    cls = 0.d0
    do l = 2, lmax
       cls(l) = q_prop * fid_spectrum(l) * (real(l,dp)/real(l_pivot,dp))**n_prop * beam(l)**2
    end do

!    call compute_invC_detC(sigma_n, cls, pls, covar_temp, covar_mat, ln_detC)


    ! Solve for C^-1 d
    ldb     = numdata
    nrhs    = 1
    uplo    = 'L'
    call dpotrs(uplo, numdata, nrhs, covar_mat, numdata, data, ldb, info)

    data = data**2

    write(*,*) 'Total chisq = ', 2.*sum(data)

    allocate(output_map(1,0:npix-1))
    output_map = -1.6375e30
    do i = 1, numdata
       output_map(1,mask2map(i)) = data(i,1)
    end do

    call write_map('chisq_map.fits', output_map)

    deallocate(data)
    deallocate(data2)
    deallocate(covar_mat)
    deallocate(covar_temp)
    deallocate(output_map)
              
  end subroutine compute_chisq_map

end module modulation_chisq
