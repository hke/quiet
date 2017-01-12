module likelihood_models
  use healpix_types
  implicit none


contains

  subroutine compute_power_law_spectrum(lmin, lmax, Q_rms_ps, n, cls)
    implicit none

    integer(i4b),                       intent(in)  :: lmin, lmax
    real(dp),                           intent(in)  :: Q_rms_ps, n
    real(dp),     dimension(lmin:lmax), intent(out) :: cls

    integer(i4b) :: l
    real(dp), allocatable, dimension(:) :: cls_int

    allocate(cls_int(0:lmax))

    cls_int(2)   = 4.d0*pi/5.d0 * Q_rms_ps**2
    cls_int(0:1) = cls_int(2)
    do l = 3, lmax
       cls_int(l) = cls_int(l-1) * real(2*(l-1)+n-1,dp)/real(2*(l-1)+5-n)
    end do

    cls(lmin:lmax) = cls_int(lmin:lmax)

    deallocate(cls_int)

  end subroutine compute_power_law_spectrum

  subroutine compute_fiducial_spectrum(cls_theory, lmin, lmax, l_pivot, q, n, cls)
    implicit none

    integer(i4b),                       intent(in)  :: lmin, lmax, l_pivot
    real(dp),                           intent(in)  :: q, n
    real(sp),     dimension(0:,1:),     intent(in)  :: cls_theory   
    real(dp),     dimension(lmin:lmax), intent(out) :: cls

    integer(i4b) :: l
    real(dp), allocatable, dimension(:) :: cls_int

    do l = lmin, lmax
       cls(l) = q * (real(l,dp)/real(l_pivot,dp))**n * real(cls_theory(l,1),dp)
    end do

  end subroutine compute_fiducial_spectrum

end module likelihood_models
