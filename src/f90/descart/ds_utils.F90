module ds_utils
implicit none


integer ds_global_feedback
integer, parameter :: ds_feedback_silent = 0, ds_feedback_quiet = 1, ds_feedback_noisy = 2, ds_feedback_debug = 3

interface ds_checkAllocate
	module procedure ds_checkAllocate_I4, ds_checkAllocate_I8, ds_checkAllocate_R4, ds_checkAllocate_R8
end interface

interface ds_isfinite
	module procedure ds_isfinite_R4, ds_isfinite_R8
end interface

character(256), private :: last_milestone
real(4), private :: last_time
contains

subroutine ds_init_milestone()
	last_milestone = "START"
	last_time = ds_get_time()
end subroutine ds_init_milestone

subroutine ds_log_milestone(milestone)
	character(*) :: milestone
	real(4) current_time, delta_time
	current_time = ds_get_time()
	delta_time = current_time - last_time
	!Deal with time overflow
	if (delta_time<0) delta_time = delta_time + ds_get_max_time()
	write(*,"('TIMELOG time taken from ',A,' to ',A,':',F7.3)") trim(last_milestone), trim(milestone), delta_time
	last_time = current_time
	last_milestone = milestone
	
end	subroutine ds_log_milestone

function det_qindex(d)
integer det_qindex, d
det_qindex = 2*d
end function det_qindex

function det_uindex(d)
integer det_uindex, d
det_uindex = 2*d+1
end function det_uindex

function ds_get_time() result(seconds)
	integer count, rate, cmax
	real(4) seconds
	call system_clock(count, rate,cmax)
	seconds = 1.0*count
	seconds = seconds/rate
end function ds_get_time

function ds_get_max_time() result(cmax)
	integer count, rate, cmax
	call system_clock(count, rate,cmax)
end function ds_get_max_time


subroutine ds_logTime(message)
	character(*) message	
	write(*,*) message," time = ", ds_get_time()
end subroutine ds_logTime

subroutine ds_log(message,threshold)
character(*) :: message
integer :: threshold
if (ds_global_feedback .ge. threshold) write(*,'(A)') trim(adjustl(message))
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

ELEMENTAL function ds_isfinite_R8(x) result(n)
	real(8), intent(in) :: x
	logical :: n
#ifdef USE_IEEE
	n = ieee_is_finite(x)
#else
	n=.true.
	if (x .ne. x .or. .not. abs(x)<=huge(x) .or. 1/x==0.0) n=.false.
#endif
end function ds_isfinite_R8

ELEMENTAL function ds_isfinite_R4(x) result(n)
	real(4), intent(in) :: x
	logical :: n
#ifdef USE_IEEE
	n = ieee_is_finite(x)
#else
	n=.true.
	if (x .ne. x .or. .not. abs(x)<=huge(x)  .or. 1/x==0.0) n=.false.
#endif
end function ds_isfinite_R4

function standard_deviation(data)
	real(8), dimension(:) :: data
	real(8) :: standard_deviation
	real(8) :: xsum, x2sum
	integer i, n
	xsum = 0
	x2sum = 0
	n = 0
	do i=lbound(data,1),ubound(data,1)
		xsum = xsum + data(i)
		x2sum = x2sum + data(i)**2
	enddo
	
	n = ubound(data,1) - lbound(data,1) + 1
	
	standard_deviation = sqrt(x2sum/n + (xsum/n)**2 )
	
end function standard_deviation

function string_ends_with(str,ending) result(r)
	logical :: r
	character(*) :: str, ending
	integer n,m,i
	
	n=len_trim(str)
	m=len_trim(ending)
	
	r=.true.
	do i=0,m-1
		if (str(n-i:n-i) .ne. ending(m-i:m-i)) then
			r=.false.
			exit
		endif
	enddo
end function string_ends_with
end module ds_utils
