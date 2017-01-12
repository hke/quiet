module ds_simple_prior
use ds_types
use ds_utils
use ds_multidetector
implicit none

contains

!new prepare_fftq2
!to operate on offsets, using noise info and correlator info
subroutine d_prepare_fftq(correlator,lc,nd,comm,sigma,fknee,alpha,nyquist,offsets,no_correlation_in)
type(ds_correlator) :: correlator
integer nd, lc
logical, optional :: no_correlation_in
logical  :: no_correlation
type(ds_offsets),dimension(0:correlator%my_ndet-1) :: offsets
real(dp),dimension(0:nd-1) :: whiteNoise, sigma, alpha
real(dp),dimension(0:nd-1) :: fknee
integer,dimension(0:correlator%ndet-1) ::  na
!real(dp),allocatable,dimension(:) :: timeCorrelator
integer :: comm
real(dp) :: nyquist
integer :: d, iglobal, fftlength, ierror

complex(dpc),allocatable,dimension(:),target :: Qinv, Pinv

if (present(no_correlation_in)) then
	no_correlation = no_correlation_in
else
	no_correlation = .false.
endif

ierror= 0 
correlator%comm = MPI_COMM_WORLD
!figure out the mpi context - size and rank.
call MPI_COMM_SIZE(comm,correlator%nproc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_size')
call MPI_COMM_RANK(comm,correlator%proc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_rank')

if (correlator%traditionalMode) then
	if (correlator%proc==0) call ds_log("Traditional Mode: No correlators.",ds_feedback_quiet)
	return
endif


!call ds_assert(.false., "correlated mode is disabled until you rewrite the - UNDERSTOOD!")

allocate(correlator%beta(0:correlator%ndet-1))

!nb: these are vectors
na=correlator%ntod/lc
correlator%ndet = nd
whiteNoise = sigma**2


!go through the offsets and build the required Pinv and Qinv
do d= 0, correlator%my_ndet-1
    iglobal= correlator%my_det(d)  !use this for noise info
    fftlength= na(iglobal)/2+1
        
    offsets(d)%ownscorrelators= .true.
    !build Qinv using pe-existing function
    allocate(offsets(d)%Qinv(fftlength))
    call buildOneCorrelator(correlator%ntod(iglobal), lc, whiteNoise(iglobal),1.0_8, &
        fknee(iglobal),alpha(iglobal),nyquist,offsets(d)%Qinv, .true.)
     
    offsets(d)%Qinv= 1.0_8 / offsets(d)%Qinv

    correlator%beta(iglobal) = real(lc,kind=8) / (regularisationParameter * whiteNoise(iglobal) )

    !build preconditioner Pinv
    allocate(offsets(d)%Pinv(fftlength))
    offsets(d)%Pinv= offsets(d)%Qinv + correlator%beta(iglobal)
    offsets(d)%Pinv= 1.0_8 / offsets(d)%Pinv
    
enddo

!if(correlator%proc==0) then
!open(unit=1001,file="Qinvfile.dat",status="replace")
!write(1001,*) real(offsets(0)%Qinv,kind=8)
!close(1001)
!endif
!stop

end subroutine d_prepare_fftq





!new apply_correlator
!to operate on offsets. goes through list of timestreams on this process, ffts the input offset
!values, divides by the input offset's Qinv (or Pinv), then ffts back.

subroutine d_applyCorrelator(offsets,correlator,isPreconditioner)
    type(ds_correlator),intent(in) :: correlator
    type(ds_offsets),dimension(0:correlator%my_ndet-1),intent(inout) :: offsets
    logical,intent(in) :: isPreconditioner

    complex(dpc),allocatable,dimension(:) :: dummySpace
    integer :: d, fftlength

    fftlength= 0
    do d= 0,correlator%my_ndet-1
        !do we need a different sized dummySpace?
        if(fftlength .ne. offsets(d)%na/2+1) then
            if(allocated(dummySpace)) deallocate(dummySpace)
            !update fftlength to new value
            fftlength= offsets(d)%na / 2 + 1
            allocate(dummySpace(fftlength))
        endif
    
        !fft the offsets
        call ds_fft(offsets(d)%values,dummySpace)

        call ds_assert(size(offsets(d)%Qinv) == fftlength,"length error in d_applyCorrelator")
        
        !covolve with the Qinv or Pinv
        if(isPreconditioner) then
            dummySpace= dummySpace * offsets(d)%Pinv
        else
            dummySpace= dummySpace * offsets(d)%Qinv
        endif
        
        !fft the offsets back
        call ds_ifft(dummySpace,offsets(d)%values)
    enddo

    deallocate(dummySpace)

end subroutine d_applyCorrelator


end module ds_simple_prior
