module dipmod_utils
  use quiet_fileutils
  use sort_utils
  use quiet_mapfile_mod
  use quiet_utils
  use rngmod
  use alm_tools
  implicit none


contains

  subroutine compute_likelihood(data, covar_mat, ln_detC, lnL)
    implicit none

    real(dp),                   intent(in)  :: ln_detC
    real(dp), dimension(1:,1:), intent(in)  :: data, covar_mat
    real(dp),                   intent(out) :: lnL
    
    integer(i4b)     :: i, numdata, ldb, nrhs, info
    character(len=1) :: uplo
    real(dp), allocatable, dimension(:,:) :: x, x_t

    numdata = size(data)

    ! Solve for C^-1 d
    allocate(x(numdata,1))
    allocate(x_t(numdata,1))
    x(:,1)  = data(:,1)
    x_t     = x
    ldb     = numdata
    nrhs    = 1
    uplo    = 'L'
    call dpotrs(uplo, numdata, nrhs, covar_mat, numdata, x, ldb, info)

    ! Return result
    lnL = -0.5d0 * (sum(x_t(:,1)*x(:,1)) + ln_detC)
!    write(*,*) 'chisq = ', sum(x_t(:,1)*x(:,1))
!    write(*,*) 'lnL = ', lnL
!    stop

    deallocate(x)
    deallocate(x_t)

  end subroutine compute_likelihood


  subroutine compute_invC_detC(mask, f_n, sigma_n, cls, pls, covar_temp, covar_mat, ln_detC)
    implicit none

    logical(lgt), dimension(0:),    intent(in)  :: mask
    real(dp),                       intent(in)  :: sigma_n
    real(dp),     dimension(0:),    intent(in)  :: cls
    real(dp),     dimension(1:),    intent(in)  :: f_n
    real(dp),     dimension(0:,1:), intent(in)  :: pls
    real(dp),     dimension(1:,1:), intent(in)  :: covar_temp
    real(dp),     dimension(1:,1:), intent(out) :: covar_mat
    real(dp),                       intent(out) :: ln_detC

    integer(i4b)     :: i, j, ii, jj, numval, lmax, l, npix, nside
    integer(i4b)     :: ldb, nrhs, n, lda, info
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, t1, t2
    character(len=1) :: uplo
    type(planck_rng) :: rng_handle

    real(dp), allocatable, dimension(:) :: pl
    real(dp),              dimension(3) :: v1, v2

    numval = size(covar_mat(:,1))
    lmax   = size(cls)-1
    npix   = size(mask)
    nside  = nint(sqrt(real(npix,sp)/12.))

    ! Construct covariance matrix
    if (size(pls) > 1) then

       ii = 1
       do i = 1, numval
          
          do j = i, numval
             
             val = 0.d0
             do l = 2, lmax
                val = val + pls(l,ii) * cls(l)
             end do
             
             covar_mat(i,j) = val
             covar_mat(j,i) = val
             
             ii = ii+1
             
          end do
       end do

    else

       call cpu_time(t1)
       allocate(pl(0:lmax))

       ii = 0
       do i = 0, npix-1
          if (mask(i)) then

             ii = ii+1
             
             call pix2vec_ring(nside, i, v1)
             
             jj = 0
             do j = 0, i
                if (mask(j)) then
                   
                   jj = jj+1

                   call pix2vec_ring(nside, j, v2)
                   
                   theta = sum(v1*v2)
                   if (theta > 1.d0) then
                      theta = 0.d0
                   else if (theta < -1.d0) then
                      theta = pi
                   else
                      theta = acos(theta)
                   end if
                
                   ! Compute the associated Legendre polynomials
                   call comp_normalised_Plm(lmax, 0, theta, pl)
       
                   ! Convert Spherical harmonics to Pl(cos(theta))
                   do l = 0, lmax
                      pl(l) = real(2*l+1,dp) * pl(l) * sqrt(4.d0*pi/(real(2*l+1,dp)))
                   end do
                   pl = pl / (4.d0*pi)


                   val = 0.d0
                   do l = 2, lmax
                      val = val + pl(l) * cls(l)
                   end do
                   
                   covar_mat(ii,jj) = val
                   covar_mat(jj,ii) = val

                end if
             end do
             
          end if
       end do

       call cpu_time(t2)

    end if

    ! Modulate field
    do i = 1, numval
       do j = 1, numval
          covar_mat(i,j) = (1.d0+f_n(i)) * covar_mat(i,j) * (1.d0+f_n(j))
       end do
    end do

    ! Add noise 
    do j = 1, numval
       covar_mat(j,j) = covar_mat(j,j) + sigma_n**2
    end do

    ! Add template marginalization
    covar_mat = covar_mat + covar_temp


    ! Cholesky factorize matrix
    uplo = 'L'
    n    = numval
    lda  = numval
    call cpu_time(t1)
    call dpotrf(uplo, n, covar_mat, lda, info)
    call cpu_time(t2)


    if (info /= 0) then
       write(*,*) 'Singular covariance matrix. info = ', info
    end if

    ! Compute determinant
    ln_detC = 0.d0
    do i = 1, numval
       ln_detC = ln_detC + log(covar_mat(i,i))
    end do
    ln_detC = 2.d0 * ln_detC

  end subroutine compute_invC_detC


  subroutine flms2f_n(trig_angles, flms, f_n)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: trig_angles
    real(dp), dimension(1:),    intent(in)  :: flms
    real(dp), dimension(1:),    intent(out) :: f_n

    integer(i4b) :: i, numval
    real(dp)     :: sqrt_3over4pi

    sqrt_3over4pi = sqrt(3.d0/(4.d0*pi))

    numval = size(f_n)
    
    f_n = 0.d0
    do i = 1, numval
       f_n(i)  =        - sqrt_3over4pi  * trig_angles(2,i) * trig_angles(4,i) * flms(2)
       f_n(i)  = f_n(i) + sqrt_3over4pi  * trig_angles(1,i)                    * flms(3)
       f_n(i)  = f_n(i) + sqrt_3over4pi  * trig_angles(2,i) * trig_angles(3,i) * flms(4)
    end do

  end subroutine flms2f_n


  subroutine f_par2f_n(vectors, f_par, f_n)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: vectors
    real(dp), dimension(1:),    intent(in)  :: f_par
    real(dp), dimension(1:),    intent(out) :: f_n

    integer(i4b) :: i, numval
    real(dp)     :: theta, phi, amp
    real(dp), dimension(3) :: vector_0

    numval = size(f_n)

    theta       = acos(f_par(1))
    phi         = f_par(2)
    amp         = f_par(3)

    vector_0(1) = sin(theta) * cos(phi)
    vector_0(2) = sin(theta) * sin(phi)
    vector_0(3) = cos(theta)

    f_n = 0.d0
    do i = 1, numval
       f_n(i)  = amp * sum(vectors(:,i) * vector_0)
    end do

  end subroutine f_par2f_n


  subroutine md_lms2md_n(trig_angles, templates, md_lms, md_n)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: trig_angles, templates
    real(dp), dimension(1:),    intent(in)  :: md_lms
    real(dp), dimension(1:),    intent(out) :: md_n

    integer(i4b) :: i, j, numval, numtemp
    real(dp)     :: one_over_sqrt4pi, sqrt_3over4pi

    sqrt_3over4pi = sqrt(3.d0/(4.d0*pi))
    one_over_sqrt4pi = sqrt(1.d0/(4.d0*pi))

    numval  = size(md_n)
    numtemp = size(md_lms)-4

    md_n = 0.d0
    do i = 1, numval
       md_n(i)  = one_over_sqrt4pi * md_lms(1)
       md_n(i)  = md_n(i) - sqrt_3over4pi  * trig_angles(2,i) * trig_angles(4,i) * md_lms(2)
       md_n(i)  = md_n(i) + sqrt_3over4pi  * trig_angles(1,i)                    * md_lms(3)
       md_n(i)  = md_n(i) + sqrt_3over4pi  * trig_angles(2,i) * trig_angles(3,i) * md_lms(4)

       do j = 1, numtemp
          md_n(i)  = md_n(i) + templates(i,j) * md_lms(4+j)
       end do
    end do

  end subroutine md_lms2md_n



  subroutine precompute_pls(mask, lmax, pls)
    implicit none

    integer(i4b),                           intent(in)  :: lmax
    logical(lgt),          dimension(0:),   intent(in)  :: mask
    real(dp),     pointer, dimension(:,:)              :: pls
    
    real(dp)     :: theta
    integer(i4b) :: i, ii, j, jj, l, m, numval, npix, nside, numpair

    real(dp),                  dimension(3) :: v1, v2

    npix   = size(mask)
    nside  = nint(sqrt(real(npix,sp)/12.))

    numval = 0
    do i = 0, npix-1
       if (mask(i)) numval = numval+1
    end do
    
    numpair = numval*(numval+1)/2

    allocate(pls(0:lmax, numpair))

    ii = 1
    do i = 0, npix-1 
       if (mask(i)) then

          call pix2vec_ring(nside, i, v1)
          
          do j = i, npix-1
             if (mask(j)) then

                call pix2vec_ring(nside, j, v2)

                theta = sum(v1*v2)
                if (theta > 1.d0) then
                   theta = 0.d0
                else if (theta < -1.d0) then
                   theta = pi
                else
                   theta = acos(theta)
                end if
                
                ! Compute the associated Legendre polynomials
                call comp_normalised_Plm(lmax, 0, theta, pls(:,ii))
       
                ! Convert Spherical harmonics to Pl(cos(theta))
                do l = 0, lmax
                   pls(l,ii) = real(2*l+1,dp) * pls(l,ii) * sqrt(4.d0*pi/(real(2*l+1,dp)))
                end do
                pls(:,ii) = pls(:,ii) / (4.d0*pi)

                ii = ii+1
             end if
          end do

       end if
    end do


  end subroutine precompute_pls


  subroutine compute_covar_temp(trig_angles, num_fg_temp, templates, fg_fact, covar_temp)
    implicit none

    integer(i4b),                   intent(in)  :: num_fg_temp
    real(dp),                       intent(in)  :: fg_fact
    real(dp),     dimension(1:,1:), intent(in)  :: trig_angles, templates
    real(dp),     dimension(1:,1:), intent(out) :: covar_temp

    integer(i4b) :: i, j, k, numval
    real(dp), allocatable, dimension(:,:) :: md_temp

    numval = size(trig_angles(1,:))

    allocate(md_temp(numval,4))

    do i = 1, numval
       md_temp(i,1) =   sqrt(1.d0/(4.d0*pi))
       md_temp(i,2) = - sqrt(3.d0/(4.d0*pi)) * trig_angles(2,i) * trig_angles(4,i) 
       md_temp(i,3) =   sqrt(3.d0/(4.d0*pi)) * trig_angles(1,i)                    
       md_temp(i,4) =   sqrt(3.d0/(4.d0*pi)) * trig_angles(2,i) * trig_angles(3,i) 
    end do

    covar_temp = 0.d0

    do j = 1, numval
       do k = 1, numval
          do i = 1, 4
             covar_temp(j,k) = covar_temp(j,k) + md_temp(j,i) * md_temp(k,i)
          end do

          do i = 1, num_fg_temp
             covar_temp(j,k) = covar_temp(j,k) + templates(j,i) * templates(k,i)
          end do
       end do
    end do

    covar_temp = covar_temp * fg_fact

    deallocate(md_temp)

  end subroutine compute_covar_temp


  subroutine generate_simulation(rng_handle, beam, sigma_n, cls, cmbmap)
    implicit none

    type(planck_rng),                intent(inout) :: rng_handle
    real(dp),                        intent(in)    :: sigma_n
    real(sp),         dimension(0:), intent(in)    :: beam, cls
    real(sp),         dimension(0:), intent(out)   :: cmbmap

    real(dp)     :: cl
    integer(i4b) :: i, j, l, m, nside, npix, lmax
    complex(spc), allocatable, dimension(:,:,:) :: alms

    lmax = size(beam)-1
    npix = size(cmbmap)
    nside = nint(sqrt(real(npix,sp)/12.))

    allocate(alms(1,0:lmax,0:lmax))
    
    alms = cmplx(0.d0,0.d0)
    do l = 2, lmax
       alms(1,l,0) = cmplx(rand_gauss(rng_handle),0.d0) * sqrt(cls(l)) * beam(l)
       do m = 1, l
          alms(1,l,m) = cmplx(rand_gauss(rng_handle),rand_gauss(rng_handle)) * &
               & sqrt(cls(l)) * beam(l) / sqrt(2.d0)
       end do
    end do

    call alm2map(nside, lmax, lmax, alms, cmbmap)    
    
    do i = 0, npix-1
       cmbmap(i) = cmbmap(i) + sigma_n * rand_gauss(rng_handle)
    end do

    deallocate(alms)

  end subroutine generate_simulation


end module dipmod_utils
