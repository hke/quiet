!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a module to the CMB Gibbs sampler gibcl      !
! June 2011, Unni Fuskeland                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module gibcl_noise_mod
  use healpix_types
  use alm_tools
  use gibcl_utils
  use quiet_utils
  use quiet_fileutils
  implicit none

  private
  public :: multiply_by_invN, initialize_noise_mod, get_N_l, remove_monopol_dipole

  real(dp)     :: sigma_noise
  integer(i4b) :: npix, nside, lmax, nmaps
  real(dp), allocatable, dimension(:)   :: diag_invN, diag_invN_sqroot
  real(dp), allocatable, dimension(:,:) :: w8ring
    real(dp), pointer, dimension(:,:) :: mask


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_noise_mod(parfile)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    character(len=*), intent(in) :: parfile

    integer(i4b) :: unit, ordering, i
    character(len=256) :: maskfile
!    real(dp), pointer, dimension(:,:) :: mask

    unit = getlun()

    call get_parameter(unit, parfile,  'SIGMA_NOISE', par_dp=sigma_noise)
    call get_parameter(unit, parfile,  'LMAX',        par_int=lmax)
    call get_parameter(unit, parfile,  'NSIDE',       par_int=nside)
    call get_parameter(unit, parfile,  'MASKFILE',    par_string=maskfile)
    call get_parameter(unit, parfile,  'SIGMA_NOISE', par_dp=sigma_noise)
    npix  = 12*nside**2
    nmaps = 1

    !Read in mask file (galaxy mask)
    write(*,*) "Reading mask from file ",trim(maskfile)
    call read_map(mask,ordering,maskfile)
    if (ordering == 2) call convert_nest2ring(nside, mask(:,1))

    call write_map_with_mask(parfile,mask(:,1))

    !Make rms noise, covariance matrix
    allocate(diag_invN(0:npix-1))
    allocate(diag_invN_sqroot(0:npix-1))
    diag_invN = 0.d0
    diag_invN_sqroot = 0.d0
    do i = 0, npix-1
       if(mask(i,1) .gt. 0.5) then
          diag_invN(i) = 1.d0/ (sigma_noise**2)
          diag_invN_sqroot(i) = 1.d0/ (sigma_noise)
       endif
    end do
!    deallocate(mask)

    ! Allocate pixel weights
    allocate(w8ring(1:2*nside,1))
    w8ring = npix/(4.d0*pi)

  end subroutine initialize_noise_mod
  


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine multiply_by_invN(sqrt, alms)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !take a real array and multiply by the inverse noise covariance matrix
    !if sqrt=.true. it multiplies by the square root off the cov. matrix
    implicit none
    
    logical(lgt),            intent(in)    :: sqrt
    real(dp), dimension(0:), intent(inout) :: alms
    
    integer(i4b) :: i
    complex(dpc), allocatable, dimension(:,:,:) :: alms_cmplx
    real(dp),     allocatable, dimension(:,:)   :: map
    real(dp),     allocatable, dimension(:)     :: buffer

    allocate(alms_cmplx(1,0:lmax,0:lmax), map(0:npix-1,1))
    call alm_real2complex(alms, alms_cmplx)
    call alm2map(nside,lmax,lmax, alms_cmplx, map(:,1))

    if (sqrt) then
       map(:,1) = diag_invN_sqroot * map(:,1)
       call map2alm(nside,lmax,lmax,map(:,1), alms_cmplx,[0.d0,0.d0],sqrt(w8ring))
    else
       map(:,1) = diag_invN        * map(:,1)
       call map2alm(nside,lmax,lmax,map(:,1), alms_cmplx,[0.d0,0.d0],w8ring)
    end if

    call alm_complex2real(alms_cmplx, alms)

    deallocate(alms_cmplx)
    
  end subroutine multiply_by_invN


  !!!!!!!!!!!!!!!!!!
  function get_N_l()
  !!!!!!!!!!!!!!!!!!
    implicit none

    real(dp) :: get_N_l

    get_N_l = sigma_noise**2 * 4.d0*pi / real(npix,dp)

  end function get_N_l


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_map_with_mask(parfile,mask)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !apply a mask to a map and write it out to file
    implicit none

    character(len=*), intent(in)          :: parfile
    real(dp), dimension(0:),intent(in)    :: mask

    character(len=256) :: wmapfile,outfile
    integer(i4b) :: unit, ordering, i,nmaps=1
    real(dp),     pointer, dimension(:,:)   :: map
    real(dp)                       :: healnan=-1.6375d30

    unit = getlun()
    call get_parameter(unit, parfile,  'WMAPFILE',    par_string=wmapfile)
    allocate(map(0:npix-1,1:3))
    call read_map(map,ordering,wmapfile)
    nside = npix2nside(size(map,1))
    if (ordering == 2) call convert_nest2ring(nside, map(:,1))

    do i = 0, npix-1
       if(mask(i) .lt. 0.5) then
          map(i,1) = healnan
       endif
    end do

    outfile='withmask_'//wmapfile
    call write_map(map, ordering=ring, outfile)
  end subroutine write_map_with_mask


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine remove_monopol_dipole(map)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    real(dp)   ,         intent(inout), dimension(0:)  :: map
    integer(i4b)   :: i, counter
    real(dp)       :: sum_unmasked,mean_map

    sum_unmasked = 0.d0
    counter = 0
    do i = 0,npix-1
       if(mask(i,1) .gt. 0.5 ) then
          sum_unmasked = sum_unmasked + map(i)
          counter = counter + 1
       end if
    end do
    mean_map = sum_unmasked / real(counter,dp)
    map = map - mean_map

  end subroutine remove_monopol_dipole



end module gibcl_noise_mod
