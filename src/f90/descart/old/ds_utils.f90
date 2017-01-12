module ds_utils
implicit none


integer ds_global_feedback
integer, parameter :: ds_feedback_silent = 0, ds_feedback_quiet = 1, ds_feedback_noisy = 2, ds_feedback_debug = 3

interface ds_checkAllocate
	module procedure ds_checkAllocate_I4, ds_checkAllocate_I8, ds_checkAllocate_R4, ds_checkAllocate_R8
end interface


contains

function det_qindex(d)
integer det_qindex, d
det_qindex = 2*d
end function det_qindex

function det_uindex(d)
integer det_uindex, d
det_uindex = 2*d+1
end function det_uindex


subroutine ds_logTime(message)
	integer count, rate, cmax
	real(4) seconds
	character(*) message
	call system_clock(count, rate,cmax)
	seconds = 1.0*count
	seconds = seconds/rate
	write(*,*) message," time = ", seconds
end subroutine ds_logTime

subroutine ds_log(message,threshold)
character(*) :: message
integer :: threshold
if (ds_global_feedback .ge. threshold) write(*,*) trim(adjustl(message))
end subroutine ds_log

subroutine ds_assert(predicate,errorMsg)
!test an assertion.  used to validate inputs to functions.
	logical :: predicate
	character(len=*) :: errorMsg
	if (.not. predicate) then
		write(*,*) 'Error'
		write(*,*) errorMsg
		stop
	endif
end subroutine ds_assert


function ds_get_lun() result(f)
!Get an unused logical unit number for a file.
integer :: f
logical :: inUse
f=50
do
inquire(UNIT=f, opened=inUse)
if (.not. inUse) exit
f=f+1
enddo
end function ds_get_lun

  function IntToStr(I)
   integer , intent(in) :: I
   character(LEN=30) IntToStr

   write (IntToStr,*) i
   IntToStr = adjustl(IntToStr)

  end function IntToStr





subroutine ds_checkAllocate_I4(field,n)
integer(4), dimension(:), allocatable :: field
integer n
    if (allocated(field)) then
      if (size(field) .ne. n) then
        deallocate(field)
        allocate(field(n))
		field= 0
      endif
    else
      allocate(field(n))
	  field= 0
    endif
end subroutine ds_checkAllocate_I4

subroutine ds_checkAllocate_I8(field,n)
integer(8), dimension(:), allocatable :: field
integer n
    if (allocated(field)) then
      if (size(field) .ne. n) then
        deallocate(field)
        allocate(field(n))
		field= 0
      endif
    else
      allocate(field(n))
	  field= 0
    endif
end subroutine ds_checkAllocate_I8


subroutine ds_checkAllocate_R4(field,n)
real(4), dimension(:), allocatable :: field
integer n
    if (allocated(field)) then
      if (size(field) .ne. n) then
        deallocate(field)
        allocate(field(n))
		field= 0.0
      endif
    else
      allocate(field(n))
	  field= 0.0
    endif
end subroutine ds_checkAllocate_R4




subroutine ds_checkAllocate_R8(field,n)
real(8), dimension(:), allocatable :: field
integer n
    if (allocated(field)) then
      if (size(field) .ne. n) then
        deallocate(field)
        allocate(field(n))
		field= 0.0_8
      endif
    else
      allocate(field(n))
	  field= 0.0_8
    endif
end subroutine ds_checkAllocate_R8



end module ds_utils
