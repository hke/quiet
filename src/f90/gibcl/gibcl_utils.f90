!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some useful utilities for the CMB gibbs sampler gibcl  !
! June 2011, Unni Fuskeland                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module gibcl_utils
  use healpix_types
  use alm_tools
  use quiet_fileutils
  implicit none

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine alm_complex2real(alm_complex,alm_real)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !converts the complex alm_Tlm(T,l,m) to the real alm_real(0:)
     implicit none

     complex(dpc), dimension(1:,0:,0:), intent(in)     :: alm_complex
     real(dp),     dimension(0:),    intent(out)    :: alm_real

     integer(i4b)                :: ii, l, m, lmax
     real(dp)                    :: sqrt_of_two

     lmax = size(alm_complex,2)-1
     sqrt_of_two = sqrt(2.d0)

     alm_real= 0.d0
     do l = 0,lmax
        ii= l*(l+1)
        alm_real(ii)= REAL(alm_complex(1,l,0))
        do m = 1,l
           alm_real(ii-m) = sqrt_of_two * AIMAG(alm_complex(1,l,m))
           alm_real(ii+m) = sqrt_of_two * REAL(alm_complex(1,l,m))
        end do
     end do

   end subroutine alm_complex2real


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine alm_real2complex(alm_real,alm_complex)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !converts real alm_real(0:) to the complex alm_Tlm(T,l,m)
     implicit none

     complex(dpc), dimension(1:,0:,0:), intent(out)     :: alm_complex
     real(dp),     dimension(0:),    intent(in)    :: alm_real

     integer(i4b)                :: ii, l, m, lmax
     real(dp)                    :: one_over_sqrt_of_two

     lmax = size(alm_complex,2)-1
     one_over_sqrt_of_two = 1.d0 / sqrt(2.d0)

     alm_complex= cmplx(0.d0,0.d0)
     do l = 0,lmax
        ii= l*(l+1)
        alm_complex(1,l,0) = cmplx(alm_real(ii),0)
        do m = 1,l
           alm_complex(1,l,m) = one_over_sqrt_of_two * cmplx(alm_real(ii+m),alm_real(ii-m))
        end do
     end do

   end subroutine alm_real2complex

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine map2alm_real(map, alms)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     !converts a map to a real alm array 
     implicit none

     real(dp), dimension(0:),    intent(in)  :: map
     real(dp), dimension(0:),    intent(out) :: alms

     integer(i4b) :: nside, npix, lmax
     real(dp),     allocatable, dimension(:,:)   :: w8ring
     complex(dpc), allocatable, dimension(:,:,:) :: alms_tmp

     npix  = size(map,1)
     nside = nint(sqrt(real(npix,dp)/12.d0))
     lmax  = sqrt(real(size(alms),dp))-1
     allocate(w8ring(1:2*nside,1), alms_tmp(1,0:lmax,0:lmax))
     w8ring = 1.d0

     call map2alm(nside, lmax, lmax, map, alms_tmp,[0.d0,0.d0],w8ring)
     call alm_complex2real(alms_tmp, alms)

     deallocate(w8ring, alms_tmp)

   end subroutine map2alm_real

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine real_multiplication(cls,alms)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !multiplies cls(0:lmax) with alms(0:lmax(lmax+2))
     implicit none
     real(dp), dimension(0:), intent(in)               :: cls
     real(dp), dimension(0:), intent(inout)            :: alms

     integer(i4b)                                      :: l, i, m, lmax

     lmax = size(cls)-1

     i = 0
     do l = 0, lmax
        do m = -l, l
           alms(i) = cls(l) * alms(i)
           i    = i+1
        end do
     end do
    
   end subroutine real_multiplication



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine lag_beam(beam_fwhm, beam)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !generates a beam from a Full With Half Maximum value
     implicit none

     real(dp)   , intent(in)                     :: beam_fwhm
     real(dp), dimension(0:), intent(inout)      :: beam

     integer(i4b)                                :: l
     real(dp)                                    :: sigma_beam_squared
     
     !beam_fwhm is in degrees, must convert to radians
     sigma_beam_squared = (beam_fwhm * pi / 180. )**2 / (8. * log(2.) ) 

     do l = 0 , size(beam)-1
        beam(l) = exp( -0.5 * l * (l+1) * sigma_beam_squared )
     end do

   end subroutine lag_beam


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine plot_cls(cls,outfile)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! plot Cl in the standard way (multiplied by l(l+1)/4pi)
     implicit none

     real(dp), dimension(0:,1:), intent(in)     :: cls
     character(len=256),    intent(in)          :: outfile
     integer(i4b)                               :: l
     real(dp)                                   :: y_cls
     logical(lgt)  :: chatty

     chatty = .false.   
     if (chatty) write(*,*) "writing to ",trim(outfile)," , plot_cls_Nl routine"
     call mkdirs(outfile,.true.)  !this makes the directory if it doesnt exist
     open(41, file=trim(outfile),action="write",status="unknown")

     do l = 0,size(cls(:,1))-1
        y_cls = cls(l,1)  * l *(l+1)/(2*pi)
        write(41,*) y_cls
     end do
        
     close(41)

   end subroutine plot_cls



end module gibcl_utils
