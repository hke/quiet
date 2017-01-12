module comap_utils
  use quiet_utils
  use quiet_fft_mod
  implicit none

contains

  subroutine read_comap_tod
    implicit none

  end subroutine read_comap_tod
  
  
  subroutine apply_filter_fft(nu_cut, tod)
    implicit none
    
    integer(i4b),                intent(in)    :: nu_cut
    real(dp),     dimension(1:), intent(inout) :: tod
    
    integer(i4b) :: n
    complex(dpc), allocatable, dimension(:) :: tod_fft

    n = size(tod)/2+1
    allocate(tod_fft(n))

    call fft(tod, tod_fft, 1)
    tod_fft(1:nu_cut) = 0.d0
    call fft(tod, tod_fft, -1)
    deallocate(tod_fft)
  end subroutine apply_filter_fft


end module comap_utils
