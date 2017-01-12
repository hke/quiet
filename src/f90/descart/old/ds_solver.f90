module ds_solver
use ds_precision
use ds_types
use ds_utils
use ds_multidetector
use ds_simple_prior
implicit none

contains

subroutine pcg(correlator,pointing,whiteNoise,npix,b,x)
! solve the system K^{-1} A x = K^{-1} b
! b= F^T Z y   where y= signal
! the arrays are distributed, 1 detector per process.
! Prepare FFTQ must have been run already.

!dummy arguments
integer nd
type(ds_correlator),intent(in) :: correlator
type(ds_detpointing),intent(in), dimension(0:correlator%my_ndet/2-1) :: pointing
type(ds_offsets),intent(inout), dimension(0:correlator%my_ndet-1) :: b
type(ds_offsets),intent(out), dimension(0:correlator%my_ndet-1) :: x
real(dp),intent(in), dimension(0:correlator%ndet-1) :: whiteNoise						!required for apply_A
integer,intent(in) :: npix							!required for dummyMap and dummyTime respectively
integer nt
integer :: i,j
real(dp) :: rho, rho_min1, rho_min2, beta, beta_min1, alpha, err
type(ds_offsets), dimension(0:correlator%ndet-1) :: r, p, q
type(ds_offsets), dimension(0:correlator%ndet-1) :: dummyOff
type(ds_map)  :: dummyQMap, dummyUMap
type(ds_timestream), dimension(0:correlator%ndet-1) :: dummyTime

integer,parameter :: itmax= 1000					!Add error message?
real(dp),parameter :: tol= 1.E-6_dp 
real(dp) :: target_err
integer d
real(dp),parameter :: units= 1.E9
character(len=125) :: message


nd = correlator%my_ndet
do d=0,nd-1

!b(d)%values= b(d)%values*units

nt = correlator%ntod(correlator%my_det(d))

call prepareOffsets(x(d),b(d))
call prepareOffsets(r(d),b(d))
call prepareOffsets(p(d),b(d))
call prepareOffsets(q(d),b(d))

!initialise dummyspaces for apply_A
call prepareOffsets(dummyOff(d),b(d))

!precondition rhs
r(d)%values=b(d)%values
enddo


call apply_prec(r,correlator)						!matrix ops - apply_prec

!initialise scalars
rho_min1= 0.0_dp
rho_min2= 0.0_dp
beta_min1= 0.0_dp

!set convergence critereon
target_err= tol * sqrt(inner3(correlator,r,r))

do i=1,itmax
	call assert_it(i,itmax,err)						!terminate if itmax reached
	rho_min2= rho_min1
	rho_min1 = inner3(correlator,r,r)	!inner product


	if(i==1) then
		do d=0,nd-1
			p(d)%values= r(d)%values
		enddo
	else	
		beta_min1= rho_min1 / rho_min2
		do d=0,nd-1
			p(d)%values= r(d)%values+ beta_min1 * p(d)%values	!update p
		enddo
!		print*,'beta_min1 = ',beta_min1
	endif
	do d=0,nd-1
		q(d)%values= p(d)%values
	enddo
	
!	if(correlator%proc==0) print*,"calling apply_A"
	
	call apply_A(q,pointing,correlator,dummyOff,whiteNoise,npix)			!matrix ops - apply_A
	call apply_prec(q,correlator)					!matrix ops - apply_prec
	
!	if(correlator%proc==0) print*,"apply_A and prec called successfully"
	
	if(rho_min1 .gt. 1.E-16_8) then
		alpha= rho_min1/ inner3(correlator,p,q)		!inner product
	else
		alpha = 0.0_8
	endif

	do d=0,nd-1
		x(d)%values= x(d)%values + alpha * p(d)%values			!update x	
		r(d)%values= r(d)%values - alpha * q(d)%values			!update r
	enddo
	
!	if(correlator%proc==0) print*,"calling inner3"
	
	!check convergence			
	err= sqrt(inner3(correlator,r,r))	!sqrt of inner product - should we normalise err to length of vector
	if(correlator%proc==0) then						!only the master node
		write(message,*),i,err,target_err							!inform user of iterative progress
		call ds_log(message,ds_feedback_quiet)
	endif	
	if(err.le. target_err) then								!exit iterations if convergence reached
		exit
	endif
	
end do	

do d=0,nd-1
call destroyOffsets(r(d))
call destroyOffsets(p(d))
call destroyOffsets(q(d))
call destroyOffsets(dummyOff(d))
enddo



end subroutine

subroutine assert_it(i,itmax,err)
	integer,intent(in)::i,itmax
	real(dp),intent(in)::err
	if(i==itmax) then
		print*,'PCG: max iterations received without convergence'
		print*,'PCG: itmax=',itmax,'error=',err
		print*,'PCG: terminating...'
		stop
	endif	
end subroutine

subroutine apply_A(offsets,pointing,correlator,dummyOffsets,whiteNoise,npix)
type(ds_correlator) :: correlator
type(ds_offsets), dimension(0:correlator%my_ndet-1) :: offsets, dummyOffsets, dummyOffsets2
type(ds_detpointing), dimension(0:correlator%my_ndet/2-1) :: pointing
!type(ds_timestream), dimension(0:correlator%my_ndet-1) :: dummyTimestream
!type(ds_map) :: dummyQMap, dummyUMap
real(dp),dimension(0:correlator%ndet-1) :: whiteNoise
integer,intent(in) :: npix
integer i,d,nd, ierr
real(dp),parameter :: units= 1.E9

!do d=0,correlator%my_ndet-1
!    offsets(d)%values= offsets(d)%values * units
!enddo

nd = correlator%my_ndet

!Apply the A matrix, A = C_w^{-1} F^T Z F + C_a^{-1}, in-place to the offset-space vector in offsets.
!The dummy spaces are used to avoid continually allocating and deallocating the memory - at the
!start of the function they should be allocated with the right sizes.
!The data in correlator should have been initialized with the prepareFFTQ routine.

!Copy the offsets into the space.  This might be as simple as dummyOFfsets=offsets
!but that might cause a memory leak (do not know).
!We apply the first term in the A matrix to the original and add in the 2nd term with the dummy

!The first term in the matrix:  C_w^{-1} F^T Z F.   We apply the operators in turn.

!call mpi_barrier(correlator%comm,ierr)
!if(correlator%proc==0) print*,"calling F operator"

!do d=0,nd-1
!    dummyTimestream(d)%timestream= 0.0_8  !DWPS: because the below now projects onto the input
!	call projectOffsetOntoTimestream(offsets(d),dummyTimestream(d),'add')   !F operator
!enddo

!!call mpi_barrier(correlator%comm,ierr)
!!if(correlator%proc==0) print*,"F operator called"

!call removeSignalNaive(dummyTimestream,pointing,dummyQMap,dummyUmap,correlator)   !Z operator

!!if(correlator%proc==0) print*,"signal naive removed"

!do d=0,nd-1
!	call deprojectTimestreamOntoOffset(dummyTimestream(d),offsets(d)) !F^T operator
!	offsets(d)%values = offsets(d)%values / whiteNoise(correlator%my_det(d)) 					!C_w^{-1} operator
!enddo

!!call mpi_barrier(correlator%comm,ierr)
!!if(correlator%proc==0) print*,"deprojetion complete"

do d=0,nd-1
    call copyOffsets(offsets(d),dummyOffsets(d))
enddo
call apply_FtZF(correlator,pointing,whiteNoise,npix,offsets,dummyOffsets)


if (.not. correlator%traditionalMode) then

!The second term in the matrix, C_a^{-1}.  This is the term where the different detectors interact,
!so it calls MPI routines.  The dummyOffsets at this point contain the input offsets.
!The final false argument indicates that we apply the main correlator, not the preconditioner.
!do d=0,nd-1
!	dummyOffsets(d)%values = offsets(d)%values
!enddo

do d=0,nd-1
    call copyOffsets(offsets(d),dummyOffsets2(d))
enddo

call d_applyCorrelator(dummyOffsets2,correlator,.false.) 

endif

!Sum the two terms together in the offsets vector.
!This is the useful output.
do d=0,nd-1
	offsets(d)%values = dummyOffsets(d)%values  !data temr
	if(.not.correlator%traditionalMode) then
	    offsets(d)%values= offsets(d)%values + dummyOffsets2(d)%values   !optional prior term
    endif
enddo

end subroutine



subroutine apply_prec(offsets,correlator)
type(ds_correlator) :: correlator
type(ds_offsets), dimension(0:correlator%my_ndet-1) :: offsets
!Apply the preconditioner matrix K^{-1}, where K= ( C_w^{-1} F^T F + C_a^{-1} )

!The final argument indicates that this is the preconditioner and not the main correlator 
if (correlator%traditionalMode) return
call d_applyCorrelator(offsets,correlator,.true.)

end subroutine

function inner1(correlator,vec1,vec2)
type(ds_correlator) :: correlator
real(8),dimension(:) :: vec1,vec2
real(8) :: inner1

inner1= sum(vec1 * vec2)

end function inner1

end module ds_solver
