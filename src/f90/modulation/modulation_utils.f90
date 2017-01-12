module modulation_utils
  use healpix_types
  use pix_tools
  use rngmod
  use fitstools
  use math_tools
  use head_fits
  implicit none

  integer(i4b)                              :: nside, npix, numtemp, numdata
  integer(i4b), allocatable, dimension(:)   :: mask2map
  real(dp),     allocatable, dimension(:,:) :: T

contains

  subroutine initialize_util_mod(parfile)
    implicit none

    character(len=*), intent(in) :: parfile

    integer(i4b)       :: i, j, ordering, unit
    character(len=2)   :: i_text
    character(len=256) :: filename
    real(dp)           :: nullval
    logical(lgt)       :: anynull
    real(dp), allocatable, dimension(:,:) :: map

    unit = 30
    call get_parameter(unit, parfile, 'MASKFILE',     par_string=filename)
    call get_parameter(unit, parfile, 'NUMTEMPLATES', par_int=numtemp)
    numtemp = numtemp + 4

    i = getsize_fits(filename, nside=nside, ordering=ordering)
    npix = 12*nside**2
    allocate(map(0:npix-1,1))
    call read_bintab(trim(filename), map, npix, 1, nullval, anynull)
    if (ordering == 2) call convert_nest2ring(nside,map(:,1))
    numdata = count(map == 1.d0)

    allocate(mask2map(numdata))
    j = 1
    do i = 0, npix-1
       if (map(i,1) == 1.d0) then
          mask2map(j) = i
          j           = j+1
       end if
    end do

    allocate(T(numdata,numtemp))
    do i = 1, numdata
       T(i,1) = 1.d0
       call pix2vec_ring(nside, mask2map(i), T(i,2:4))
    end do

    if (numtemp > 4) then
       do i = 1, numtemp-4
          call int2string(i, i_text)
          call get_parameter(unit, parfile, 'TEMPLATE'//i_text, par_string=filename)
          j = getsize_fits(filename, nside=nside, ordering=ordering)
          npix = 12*nside**2
          call read_bintab(trim(filename), map, npix, 1, nullval, anynull)
          if (ordering == 2) call convert_nest2ring(nside,map(:,1))
          T(:,4+i) = map(mask2map,1)
       end do
    end if

    deallocate(map)

  end subroutine initialize_util_mod

  subroutine compute_likelihood(data, covar_mat, ln_detC, lnL, chisq)
    implicit none

    real(dp),                   intent(in)  :: ln_detC
    real(dp), dimension(1:,1:), intent(in)  :: data, covar_mat
    real(dp),                   intent(out) :: lnL, chisq
    
    integer(i4b)     :: i, ldb, nrhs, info
    character(len=1) :: uplo
    real(dp), allocatable, dimension(:,:) :: x, x_t

    ! Solve for C^-1 d
    allocate(x(numdata,1))
    allocate(x_t(numdata,1))
    x(:,1)  = data(:,1)
    x_t     = x
!    ldb     = numdata
!    nrhs    = 1
!    uplo    = 'L'
!    call dpotrs(uplo, numdata, nrhs, covar_mat, numdata, x, ldb, info)
    x = matmul(covar_mat, x)
    

    ! Return result
    chisq = sum(x_t(:,1)*x(:,1))
    lnL = -0.5d0 * (chisq + ln_detC)
!    write(*,*) 'chisq = ', sum(x_t(:,1)*x(:,1)), numdata, (sum(x_t(:,1)*x(:,1))-numdata) / sqrt(2.d0*numdata)
!    write(*,*) 'ln_det_C = ', ln_detC
!    write(*,*) 'lnL = ', lnL
!    write(*,*) 'sum_data = ', sum(data)
!    stop

    deallocate(x)
    deallocate(x_t)

  end subroutine compute_likelihood


  subroutine compute_invC_detC(f_n, sigma_n, l_mod, cls, pls, covar_temp, covar_mat, ln_detC)
    implicit none

    integer(i4b),                   intent(in)  :: l_mod
    real(dp),                       intent(in)  :: sigma_n
    real(dp),     dimension(0:),    intent(in)  :: cls
    real(dp),     dimension(1:),    intent(in)  :: f_n
    real(dp),     dimension(0:,1:), intent(in)  :: pls
    real(dp),     dimension(1:,1:), intent(in)  :: covar_temp
    real(dp),     dimension(1:,1:), intent(out) :: covar_mat
    real(dp),                       intent(out) :: ln_detC

    integer(i4b)     :: i, j, ii, jj, lmax, l
    integer(i4b)     :: ldb, nrhs, n, lda, info
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, t1, t2
    character(len=1) :: uplo
    type(planck_rng) :: rng_handle

    real(dp), allocatable, dimension(:)   :: pl, eigenvals
    real(dp), allocatable, dimension(:,:) :: invC_T, invC_T_W, W
    real(dp),              dimension(3)   :: v1, v2

    lmax   = size(cls)-1

    ! Construct covariance matrix
    if (size(pls) > 1) then
       ii = 1
       do i = 1, numdata
          do j = i, numdata
             val = 0.d0
             do l = 2, l_mod
                val = val + pls(l,ii) * cls(l)
             end do
             covar_mat(i,j) = val
             covar_mat(j,i) = val
             ii = ii+1
          end do
       end do
    else

       call wall_time(t1)
       allocate(pl(0:lmax))
       do i = 1, numdata
          call pix2vec_ring(nside, mask2map(i), v1)
          do j = 1, i
             call pix2vec_ring(nside, mask2map(j), v2)
             theta = acos(max(min(sum(v1*v2),1.d0),-1.d0))
                
             ! Compute the associated Legendre polynomials
             call comp_normalised_Plm(lmax, 0, theta, pl)
       
             ! Convert Spherical harmonics to Pl(cos(theta))
             do l = 0, l_mod
                pl(l) = real(2*l+1,dp) * pl(l) * sqrt(4.d0*pi/(real(2*l+1,dp)))
             end do
             pl = pl / (4.d0*pi)

             val = 0.d0
             do l = 2, l_mod
                val = val + pl(l) * cls(l)
             end do
                   
             covar_mat(i,j) = val
             covar_mat(j,i) = val

          end do
       end do
       call wall_time(t2)
    end if

    ! Modulate field
    do i = 1, numdata
       do j = 1, numdata
          covar_mat(i,j) = (1.d0+f_n(i)) * covar_mat(i,j) * (1.d0+f_n(j))
       end do
    end do

    ! Add noise 
    do j = 1, numdata
       covar_mat(j,j) = covar_mat(j,j) + sigma_n**2
    end do

    ! Add template marginalization
    covar_mat = covar_mat + covar_temp

    ! Cholesky factorize matrix
    uplo = 'L'
    n    = numdata
    lda  = numdata
    call cpu_time(t1)
    call dpotrf(uplo, n, covar_mat, lda, info)
    call cpu_time(t2)


    if (info /= 0) then
       write(*,*) 'Singular covariance matrix. info = ', info
    end if

    ! Compute determinant
    ln_detC = 0.d0
    do i = 1, numdata
       ln_detC = ln_detC + log(covar_mat(i,i))
    end do
    ln_detC = 2.d0 * ln_detC

    ! Invert matrix 
    call dpotri(uplo, n, covar_mat, lda, info)
    
    ! Symmetrize
    do i = 1, numdata
       do j = i+1, numdata
          covar_mat(i,j) = covar_mat(j,i)
       end do
    end do

    ! Add Woodbury correction
    allocate(invC_T(numdata, numtemp), W(numtemp, numtemp), invC_T_W(numdata,numtemp))
    invC_T = matmul(covar_mat, T)
    W      = matmul(transpose(T), invC_T)
    
    call invert_matrix(W)
    invC_T_W = matmul(invC_T, W)
    covar_mat = covar_mat - matmul(invC_T, transpose(invC_T_W))
    deallocate(invC_T, W, invC_T_W)
    
    ! Compute log determinant
    allocate(eigenvals(numdata))
    call get_eigenvalues(covar_mat, eigenvals)
    ln_detC = 0.d0
    do i = numtemp+1, numdata
       ln_detC = ln_detC - log(eigenvals(i))
    end do
    deallocate(eigenvals)

  end subroutine compute_invC_detC




  subroutine compute_high_l_covar(l_mod, cls, pls, covar_temp)
    implicit none

    integer(i4b),                   intent(in)    :: l_mod
    real(dp),     dimension(0:),    intent(in)    :: cls
    real(dp),     dimension(0:,1:), intent(in)    :: pls
    real(dp),     dimension(1:,1:), intent(inout) :: covar_temp

    integer(i4b)     :: i, j, ii, jj, lmax, l
    integer(i4b)     :: ldb, nrhs, n, lda, info
    real(dp)         :: val, ln_det, lnL, lnL_prop, theta, phi, t1, t2
    character(len=1) :: uplo
    type(planck_rng) :: rng_handle

    real(dp), allocatable, dimension(:) :: pl
    real(dp),              dimension(3) :: v1, v2

    lmax   = size(cls)-1

    ! Construct covariance matrix
    if (size(pls) > 1) then

       ii = 1
       do i = 1, numdata
          do j = i, numdata
             val = 0.d0
             do l = l_mod+1, lmax
                val = val + pls(l,ii) * cls(l)
             end do
             covar_temp(i,j) = covar_temp(i,j) + val
             if (i /= j) covar_temp(j,i) = covar_temp(j,i) + val
             ii = ii+1
          end do
       end do

    else

       call cpu_time(t1)
       allocate(pl(0:lmax))

       do i = 1, numdata
          call pix2vec_ring(nside, mask2map(i), v1)
          do j = 1, i
             call pix2vec_ring(nside, mask2map(j), v2)
             theta = acos(max(min(sum(v1*v2),1.d0),-1.d0))
             
             ! Compute the associated Legendre polynomials
             call comp_normalised_Plm(lmax, 0, theta, pl)
             
             ! Convert Spherical harmonics to Pl(cos(theta))
             do l = l_mod+1, lmax
                pl(l) = real(2*l+1,dp) * pl(l) * sqrt(4.d0*pi/(real(2*l+1,dp)))
             end do
             pl = pl / (4.d0*pi)
             
             val = 0.d0
             do l = l_mod+1, lmax
                val = val + pl(l) * cls(l)
             end do
             
             covar_temp(i,j) = covar_temp(i,j) + val
             if (i /= j) covar_temp(j,i) = covar_temp(j,i) + val
          end do
       end do

       call cpu_time(t2)
       deallocate(pl)
       
    end if

  end subroutine compute_high_l_covar



  subroutine flms2f_n(trig_angles, flms, f_n)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: trig_angles
    real(dp), dimension(1:),    intent(in)  :: flms
    real(dp), dimension(1:),    intent(out) :: f_n

    integer(i4b) :: i
    real(dp)     :: sqrt_3over4pi

    sqrt_3over4pi = sqrt(3.d0/(4.d0*pi))

    f_n = 0.d0
    do i = 1, numdata
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

    integer(i4b) :: i
    real(dp)     :: theta, phi, amp
    real(dp), dimension(3) :: vector_0

    theta       = acos(f_par(1))
    phi         = f_par(2)
    amp         = f_par(3)

    vector_0(1) = sin(theta) * cos(phi)
    vector_0(2) = sin(theta) * sin(phi)
    vector_0(3) = cos(theta)

    f_n = 0.d0
    do i = 1, numdata
       f_n(i)  = amp * sum(vectors(:,i) * vector_0)
    end do

  end subroutine f_par2f_n


  subroutine md_lms2md_n(trig_angles, templates, md_lms, md_n)
    implicit none

    real(dp), dimension(1:,1:), intent(in)  :: trig_angles, templates
    real(dp), dimension(1:),    intent(in)  :: md_lms
    real(dp), dimension(1:),    intent(out) :: md_n

    integer(i4b) :: i, j, numtemp
    real(dp)     :: one_over_sqrt4pi, sqrt_3over4pi

    sqrt_3over4pi    = sqrt(3.d0/(4.d0*pi))
    one_over_sqrt4pi = sqrt(1.d0/(4.d0*pi))

    numtemp = size(md_lms)-4

    md_n = 0.d0
    do i = 1, numdata
       md_n(i)  = one_over_sqrt4pi * md_lms(1)
       md_n(i)  = md_n(i) - sqrt_3over4pi  * trig_angles(2,i) * trig_angles(4,i) * md_lms(2)
       md_n(i)  = md_n(i) + sqrt_3over4pi  * trig_angles(1,i)                    * md_lms(3)
       md_n(i)  = md_n(i) + sqrt_3over4pi  * trig_angles(2,i) * trig_angles(3,i) * md_lms(4)

       do j = 1, numtemp
          md_n(i)  = md_n(i) + templates(i,j) * md_lms(4+j)
       end do
    end do

  end subroutine md_lms2md_n



  subroutine precompute_pls(lmax, pls)
    implicit none

    integer(i4b),                           intent(in)  :: lmax
    real(dp),     pointer, dimension(:,:)              :: pls
    
    real(dp)     :: theta
    integer(i4b) :: i, ii, j, jj, l, m, numpair

    real(dp),                  dimension(3) :: v1, v2

    numpair = numdata*(numdata+1)/2

    allocate(pls(0:lmax, numpair))

    ii = 1
    do i = 1, numdata
       call pix2vec_ring(nside, mask2map(i), v1)
       do j = i, numdata
          call pix2vec_ring(nside, mask2map(j), v2)
          theta = acos(max(min(sum(v1*v2),1.d0),-1.d0))
          call comp_normalised_Plm(lmax, 0, theta, pls(:,ii))

          ! Convert Spherical harmonics to Pl(cos(theta))
          do l = 0, lmax
             pls(l,ii) = real(2*l+1,dp) * pls(l,ii) * sqrt(4.d0*pi/(real(2*l+1,dp)))
          end do
          pls(:,ii) = pls(:,ii) / (4.d0*pi)
          
          ii = ii+1
       end do
    end do

  end subroutine precompute_pls


  ! Convention: First phi around z, then theta around x, then psi around z
  subroutine compute_euler_matrix(phi, theta, psi, euler_matrix)
    implicit none
                                                                                
    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix
 
    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi
     
    sphi = sin(phi)
    cphi = cos(phi)
     
    sth  = sin(theta)
    cth  = cos(theta)
 
    spsi = sin(psi)
    cpsi = cos(psi)
 
    euler_matrix(1,1) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(1,2) =  sphi * cpsi + cth * cphi * spsi
    euler_matrix(1,3) =                sth * spsi
    euler_matrix(2,1) = -cphi * spsi - cth * sphi * cpsi
    euler_matrix(2,2) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(2,3) =                sth * cpsi
    euler_matrix(3,1) =                sth * sphi
    euler_matrix(3,2) =              - sth * cphi
    euler_matrix(3,3) =                cth
 
  end subroutine compute_euler_matrix

  ! Convention: First psi around z, then theta around y, then phi around z
  subroutine compute_euler_matrix_zyz(phi, theta, psi, euler_matrix)
    implicit none
                                                                                
    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix
 
    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi
     
    sphi = sin(phi)
    cphi = cos(phi)
     
    sth  = sin(theta)
    cth  = cos(theta)
 
    spsi = sin(psi)
    cpsi = cos(psi)
 
    euler_matrix(1,1) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(1,2) = -sphi * cpsi - cth * cphi * spsi
    euler_matrix(1,3) =                sth * cphi
    euler_matrix(2,1) =  cphi * spsi + cth * sphi * cpsi
    euler_matrix(2,2) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(2,3) =                sth * sphi
    euler_matrix(3,1) =              - sth * cpsi
    euler_matrix(3,2) =                sth * spsi
    euler_matrix(3,3) =                cth
 
  end subroutine compute_euler_matrix_zyz


  ! Small utility for converting an integer to a string
  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string
  

  ! *****************************************************************************************
  !
  ! Routine for reading one parameter from a parameter file
  !     Example usage:   call get_parameter(21, "parfile.txt", "NSIDE", par_int=nside)
  !
  ! *****************************************************************************************

  subroutine get_parameter(unit, parfile, parname, par_int, par_char, par_string, par_sp, par_dp, par_lgt)
    implicit none

    integer(i4b),      intent(in)  :: unit
    character(len=*),  intent(in)  :: parfile, parname
    integer(i4b),      intent(out), optional :: par_int
    logical(lgt),      intent(out), optional :: par_lgt
    character(len=1),  intent(out), optional :: par_char
    character(len=*),  intent(out), optional :: par_string
    real(sp),          intent(out), optional :: par_sp
    real(dp),          intent(out), optional :: par_dp


    integer(i4b)        :: i
    character(len=128)  :: string, variable, value
    character(len=1)    :: equals


    open(unit, file=trim(parfile))

    do while (.true.)
       read(unit,*,end=1) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*) variable, equals, value

       if (trim(variable) == trim(parname)) then

          if (present(par_int)) then
             read(value,*) par_int
          else if (present(par_char)) then
             read(value,*) par_char
          else if (present(par_string)) then
             read(value,*) par_string
          else if (present(par_sp)) then
             read(value,*) par_sp
          else if (present(par_dp)) then
             read(value,*) par_dp
          else if (present(par_lgt)) then
             read(value,*) par_lgt
          end if

          close(unit)
          return

       end if

    end do

1   write(*,*) 'GET_PARAMETER:    Critical error -- parameter not found'
    write(*,*) 'GET_PARAMETER:       Parameter file = ', trim(parfile)
    write(*,*) 'GET_PARAMETER:       Parameter name = ', trim(parname)

    close(unit)
    stop

  end subroutine get_parameter

  subroutine read_beam(beamfile, lmax, beam)
    implicit none

    character(len=128),                   intent(in)  :: beamfile
    integer(i4b),                         intent(in)  :: lmax
    real(dp),           dimension(0:,1:), intent(out) :: beam

    integer(i4b) :: i, l, nlheader, nmaps
    real(sp)     :: sigma_sq
    real(sp),          allocatable, dimension(:,:)   :: inbeam
    character(len=80),              dimension(1:180) :: header
    
    nmaps = size(beam(0,:))

    ! Seem to remember there is something weird going on with the WMAP beams when reading only 
    ! one component at a time. Remove this wrapper if you feel more comfortable with that...
    allocate(inbeam(0:lmax,4))

    nlheader = size(header)
    call fits2cl(beamfile, inbeam, lmax, 4, header)

    do i = 1, nmaps
       beam(:,i) = inbeam(:,i)
    end do

    deallocate(inbeam)

  end subroutine read_beam

  subroutine read_ringweights(nside, nmaps, weights)
    implicit none

    integer(i4b),                          intent(in)  :: nside, nmaps
    real(dp),     pointer, dimension(:,:)              :: weights

    character(len=128)  :: weight_file
    character(len=5)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    allocate(weights(1:2*nside,nmaps))
    
    call int2string(nside, nside_text)
    weight_file = 'weight_ring_n' // nside_text // '.fits'
    inquire(file=weight_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(weight_file, weights, 2*nside, nmaps, nullval, anynull)
       weights = 1.d0 + weights
    else
       weights = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Weight file ', trim(weight_file), ' not found. '
       write(*,*) 'Using unity weights in the spherical harmonic transforms.'
       write(*,*) ''
    end if

  end subroutine read_ringweights

  subroutine read_pixwin(nside, nmaps, pixwin)
    implicit none

    integer(i4b),                          intent(in)  :: nside, nmaps
    real(dp),     pointer, dimension(:,:)              :: pixwin

    integer(i4b)        :: nc
    character(len=128)  :: pixwin_file
    character(len=4)    :: nside_text
    logical(lgt)        :: exist, anynull
    real(dp)            :: nullval

    allocate(pixwin(0:4*nside,nmaps))
    
    if (nmaps == 3) then
       nc = 2
    else
       nc = 1
    end if

    call int2string(nside, nside_text)
    pixwin_file = 'pixel_window_n' // nside_text // '.fits'
    inquire(file=pixwin_file, exist=exist)
    if (exist) then
       nullval = 0.d0
       call read_dbintab(pixwin_file, pixwin(0:4*nside,1:nc), 4*nside+1, nc, nullval, anynull)
       if (nmaps == 3) pixwin(:,3) = pixwin(:,2)
    else
       pixwin = 1.d0
       write(*,*) ''
       write(*,*) 'Warning! Pixel window file ', trim(pixwin_file), ' not found. '
       write(*,*) 'Using unity weights.'
       write(*,*) ''
    end if

  end subroutine read_pixwin



  subroutine write_map(filename, map)
    implicit none

    character(len=*),                   intent(in)  :: filename
    real(dp),         dimension(0:,1:), intent(in)  :: map

    integer(i4b)   :: npix, nlheader, nmaps, i, nside
    logical(lgt)   :: exist, polarization

    character(len=80), dimension(1:120)    :: header

    npix         = size(map(:,1))
    nside        = nint(sqrt(real(npix,sp)/12.))
    nmaps        = size(map(0,:))
    polarization = (nmaps == 3)


    inquire(file=trim(filename),exist=exist)
    if (exist) call system('rm ' // trim(filename))
    

    !-----------------------------------------------------------------------
    !                      write the map to FITS file
    !  This is copied from the synfast.f90 file in the Healpix package
    !-----------------------------------------------------------------------
    
    nlheader = SIZE(header)
    do i=1,nlheader
       header(i) = ""
    enddo

    ! start putting information relative to this code and run
    call add_card(header)
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Sky Map Pixelisation Specific Keywords    ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"PIXTYPE","HEALPIX","HEALPIX Pixelisation")
    call add_card(header,"ORDERING","RING",  "Pixel ordering scheme, either RING or NESTED")
    call add_card(header,"NSIDE"   ,nside,   "Resolution parameter for HEALPIX")
    call add_card(header,"FIRSTPIX",0,"First pixel # (0 based)")
    call add_card(header,"LASTPIX",npix-1,"Last pixel # (0 based)")
    call add_card(header,"BAD_DATA",  HPX_DBADVAL ,"Sentinel value given to bad pixels")
    call add_card(header) ! blank line

    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Planck Simulation Specific Keywords      ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"EXTNAME","Commander Gibbs sample")
    call add_card(header,"POLCCONV","COSMO"," Coord. convention for polarisation (COSMO/IAU)")

    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"COMMENT","     Data Description Specific Keywords       ")
    call add_card(header,"COMMENT","-----------------------------------------------")
    call add_card(header,"INDXSCHM","IMPLICIT"," Indexing : IMPLICIT or EXPLICIT")
    call add_card(header,"GRAIN", 0, " Grain of pixel indexing")
    call add_card(header,"COMMENT","GRAIN=0 : no indexing of pixel data                         (IMPLICIT)")
    call add_card(header,"COMMENT","GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)")
    call add_card(header,"COMMENT","GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)")
    call add_card(header) ! blank line
    call add_card(header,"POLAR",polarization," Polarisation included (True/False)")
    call add_card(header,"DERIV",0," Derivative included (0, 1 or 2)")

    call add_card(header) ! blank line
    call add_card(header,"TTYPE1", "TEMPERATURE","Temperature map")
    call add_card(header,"TUNIT1", 'muK',"map unit")
    call add_card(header)

    if (polarization) then
       call add_card(header,"TTYPE2", "Q-POLARISATION","Q Polarisation map")
       call add_card(header,"TUNIT2", 'muK',"map unit")
       call add_card(header)
       
       call add_card(header,"TTYPE3", "U-POLARISATION","U Polarisation map")
       call add_card(header,"TUNIT3", 'muK',"map unit")
       call add_card(header)
    endif
    call add_card(header,"COMMENT","*************************************")

    call write_bintab(map, npix, nmaps, header, nlheader, filename)

  end subroutine write_map


end module modulation_utils
