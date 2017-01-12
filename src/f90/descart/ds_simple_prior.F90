module ds_simple_prior

use ds_types
use ds_utils
use ds_multidetector
implicit none


contains

subroutine apply_prior(moduleScan,correlator,isPreconditioner)

    type(ds_correlator),intent(in) :: correlator
	type(ds_moduleScan), dimension(0:correlator%my_nmodules-1) :: moduleScan
    type(ds_offsets),dimension(:),pointer :: offsets
    logical,intent(in) :: isPreconditioner

    complex(dpc),allocatable,dimension(:) :: dummySpace
    integer :: d, fftlength, m

	do m=0,size(moduleScan)-1
		offsets=>moduleScan(m)%offsets_2
	    fftlength= 0
	    do d= 1,modulescan(m)%ndiodes
	        if(moduleScan(m)%flags(d)==0) cycle !ignore if bad
	    
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
	enddo

    deallocate(dummySpace)

end subroutine apply_prior



!operates on main data modulescan array (before any data is put in).
subroutine prepare_prior(correlator,lc,nyquist,modulescan)

type(ds_correlator) :: correlator
integer :: lc
real(dp) :: nyquist
type(ds_modulescan),dimension(0:correlator%my_nmodules-1),intent(inout) :: modulescan

integer :: d, m, fftlength, ierror, na
real(dp) :: beta
type(ds_offsets),pointer,dimension(:) :: offsets  !dummy for ease of reading

!#####################################################!
!##                                                 ##!
!## This chunk of code prepares the correlator type ##!
!##                                                 ##!
!#####################################################!

ierror= 0 
!figure out the mpi context - size and rank.
call MPI_COMM_SIZE(correlator%comm,correlator%nproc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_size')
call MPI_COMM_RANK(correlator%comm,correlator%proc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_rank')

if (correlator%traditionalMode) then
	if (correlator%proc==0) call ds_log("Traditional Mode: No correlators.",ds_feedback_quiet)
	return
endif


!##################################################################!
!##                                                              ##!
!## This chunk of code prepares fourier prior and preconditioner ##!
!##                                                              ##!
!##################################################################!


do m= 0, correlator%my_nmodules-1
    na= modulescan(m)%ntod/lc
    fftlength= na/2+1
    
    !check to make sure offsets not already used
    call ds_assert(associated(modulescan(m)%offsets_2), &
        "offsets in modulescan not associated in prepare_fftq2")
    
    !point offsets to the modulescan offsets
    offsets => modulescan(m)%offsets_2
    
    do d=1,modulescan(m)%ndiodes
        if(modulescan(m)%flags(d)==0) cycle !ignore if bad
        
        offsets(d)%ownscorrelators= .true.
        
        !build Qinv using pe-existing function
        allocate(offsets(d)%Qinv(fftlength))
        call buildOneCorrelator(modulescan(m)%ntod, lc, modulescan(m)%noise%sigma(d)**2 ,1.0_8, &
        modulescan(m)%noise%fknee(d),modulescan(m)%noise%alpha(d),nyquist,offsets(d)%Qinv, .true.)
     
        offsets(d)%Qinv= 1.0_8 / offsets(d)%Qinv

        beta = real(lc,kind=8) / ( modulescan(m)%noise%sigma(d)**2 )

        !build preconditioner Pinv
        allocate(offsets(d)%Pinv(fftlength))
        offsets(d)%Pinv= offsets(d)%Qinv + beta
        offsets(d)%Pinv= 1.0_8 / offsets(d)%Pinv
    enddo
enddo


end subroutine prepare_prior



subroutine prepare_one_data_prior(correlator,lc,modulescan)

type(ds_correlator) :: correlator
integer :: lc
type(ds_modulescan),intent(inout) :: modulescan
real(dp):: beta
integer :: d, m, fftlength, ierror, na
type(ds_offsets),pointer,dimension(:) :: offsets  !dummy for ease of reading

!#####################################################!
!##                                                 ##!
!## This chunk of code prepares the correlator type ##!
!##                                                 ##!
!#####################################################!

ierror= 0 
!figure out the mpi context - size and rank.
call MPI_COMM_SIZE(correlator%comm,correlator%nproc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_size')
call MPI_COMM_RANK(correlator%comm,correlator%proc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_rank')

if (correlator%traditionalMode) then
	if (correlator%proc==0) call ds_log("Traditional Mode: No correlators.",ds_feedback_quiet)
	return
endif


!##################################################################!
!##                                                              ##!
!## This chunk of code prepares fourier prior and preconditioner ##!
!##                                                              ##!
!##################################################################!


    na = modulescan%ntod/lc
    fftlength= na/2+1
    
    !check to make sure offsets not already used
    call ds_assert(associated(modulescan%offsets_2), &
        "offsets in modulescan not associated in prepare_fftq2")
    
    !point offsets to the modulescan offsets
    offsets => modulescan%offsets_2
    
    do d= 1,modulescan%ndiodes
        if(modulescan%flags(d)==0) cycle !ignore if bad
        
        offsets(d)%ownscorrelators= .true.
        
        !build Qinv using pe-existing function
        allocate(offsets(d)%Qinv(fftlength))
        call buildOneCorrelatorFromData(modulescan%ntod, lc , modulescan%timestreams(d)%timestream,  offsets(d)%Qinv,correlator%proc)
     
        offsets(d)%Qinv= 1.0_8 / offsets(d)%Qinv

        beta = real(lc,kind=8) / (modulescan%noise%sigma(d)**2)

        !build preconditioner Pinv
        allocate(offsets(d)%Pinv(fftlength))
        offsets(d)%Pinv= offsets(d)%Qinv + beta
        offsets(d)%Pinv= 1.0_8 / offsets(d)%Pinv
    enddo


end subroutine prepare_one_data_prior




subroutine buildOneCorrelator(nt,lc,whiteNoiseVariance,rho,fknee,alpha,nyquist,offsetPower,includeWhite)
real(dp) :: whiteNoiseVariance,fknee,alpha
integer lc
integer nt
integer na
logical includeWhite
real(dp), allocatable, dimension(:) :: timeCorrelator, offsetCorrelator
character(128) :: message
complex(dpc), dimension(:) :: offsetPower
complex(dpc), dimension(:), allocatable :: timePower
integer k,a,i,t1,t2, j
real(dp) :: dk, sigma_min1, sigma, rho, nyquist

	na = nt/lc

	call ds_assert(size(offsetPower)==na/2+1,"Wrong sized offsetPower")
	write(message,*) "Wrong sign convention for alpha; should be less than zero (",alpha,")"
	call ds_assert(alpha<0,message)
	if (rho==0) then
!		allocate(offsetPower(1:na/2+1) )
		offsetPower = dcmplx(0,0)
	   if (includeWhite) offsetPower = whiteNoiseVariance/lc
		return
	endif


   allocate(timePower(0:nt/2) )

   !Build the power spectrum
   timePower(0)=cmplx(0.0_dp,0.0_dp,kind=dpc)
   dk = 2.0_dp * nyquist / nt
   do k=1,nt/2
      timePower(k) = dcmplx( (1.0_dp / k / dk)**(-alpha), 0.0_dp)
   enddo

   !Convert power to time correlator
   allocate(timeCorrelator(0:nt-1) )
   call ds_ifft(timePower,timeCorrelator)
   deallocate(timePower)

   !convert time correlator to offset correlator, for off-diagonal components
   allocate(offsetCorrelator(1:na) )
   offsetCorrelator=0
   do a= 1,na					!a indexes the offset
      do i= 0,lc-1				!i indexes the time element in the offset
         t1= (a-1)*lc+i 		!t1 is the time index of the i'th element in offset a
         if(i == 0) then							!only need to do sum once per a
         
         	!replace with a forall with mask?
            do t2= 0,lc-1
               if(abs(t1 - t2).le.nt) then
                  offsetCorrelator(a)= offsetCorrelator(a) + timeCorrelator(abs(t1 - t2))
               endif
            enddo
            sigma_min1= offsetCorrelator(a)				!first row's sum
         else									
         	!following rows sum is sigma_min1 + new end element - old start element
            sigma = sigma_min1 - timeCorrelator(abs((-(a-1)+1)*lc-i)) + timeCorrelator(abs(-((a-1)*lc)-(i+1)+1)) !ichunk is 1 based,i and timeCorrelator are zero based			 
            sigma_min1 = sigma
            offsetCorrelator(a) = offsetCorrelator(a) + sigma
         endif
      enddo
   enddo

	deallocate(timeCorrelator)
   
   !Do (F^t F)^-2
   offsetCorrelator = offsetCorrelator/(lc**2)
!   allocate(offsetPower(1:na/2+1) )
   call ds_fft(offsetCorrelator,offsetPower)
   deallocate(offsetCorrelator)

   !Add this lengths's correlators to Qinv
   offsetPower = offsetPower * whiteNoiseVariance * rho * (fknee**(-alpha))
	if (includeWhite) offsetPower =  offsetPower + whiteNoiseVariance/lc

	call ds_assert(all(ds_isfinite(real(offsetPower))),"NaN or infinity detected in noise power spectrum.  Bad parameters?")


end subroutine buildOneCorrelator


subroutine buildOneCorrelatorFromData(nt,lc,timestream,offsetPower,rank)
real(dp) :: whiteNoiseVariance,fknee,alpha
integer lc
integer nt
integer na
logical includeWhite
real(dp), allocatable, dimension(:) :: timeCorrelator, offsetCorrelator
character(128) :: message
real(dp), dimension(:) :: timestream
complex(dpc), dimension(:) :: offsetPower
complex(dpc), dimension(:), allocatable :: timePower
integer k,a,i,t1,t2, j
real(dp) :: dk, sigma_min1, sigma, rho
integer, save :: ncalls = 1
integer rank
real(dp) :: minpow

	rho=1.0
	na = nt/lc

	call ds_assert(size(offsetPower)==na/2+1,"Wrong sized offsetPower in buildOneCorrelatorFromData")
	call ds_assert(size(timestream)==nt,"Wrong sized timestream in buildOneCorrelatorFromData")


   allocate(timePower(0:nt/2) )

   !Build the power spectrum
	call ds_fft(timestream,timePower)
	timePower=abs(timePower)**2
!	write(*,*) "min timepower = ", minval(real(timePower)), "lc = ", lc
   !Convert power to time correlator
   allocate(timeCorrelator(0:nt-1) )
   call ds_ifft(timePower,timeCorrelator)
   deallocate(timePower)

   !convert time correlator to offset correlator, for off-diagonal components
   allocate(offsetCorrelator(1:na) )
   offsetCorrelator=0
   do a= 1,na					!a indexes the offset
      do i= 0,lc-1				!i indexes the time element in the offset
         t1= (a-1)*lc+i 		!t1 is the time index of the i'th element in offset a
         if(i == 0) then							!only need to do sum once per a
         
         	!replace with a forall with mask?
            do t2= 0,lc-1
               if(abs(t1 - t2).le.nt) then
                  offsetCorrelator(a)= offsetCorrelator(a) + timeCorrelator(abs(t1 - t2))
               endif
            enddo
            sigma_min1= offsetCorrelator(a)				!first row's sum
         else									
         	!following rows sum is sigma_min1 + new end element - old start element
            sigma = sigma_min1 - timeCorrelator(abs((-(a-1)+1)*lc-i)) + timeCorrelator(abs(-((a-1)*lc)-(i+1)+1)) !ichunk is 1 based,i and timeCorrelator are zero based			 
            sigma_min1 = sigma
            offsetCorrelator(a) = offsetCorrelator(a) + sigma
         endif
      enddo
   enddo

	deallocate(timeCorrelator)
   
   !Do (F^t F)^-2
   offsetCorrelator = offsetCorrelator/(lc**2)
	if(rank==0) then
!		write(*,*) ncalls
		ncalls=ncalls+1
		open(file="correlator.dat",unit=12)
		do a=1,na
			write(12,*) offsetCorrelator(a)
		enddo
	endif
!   allocate(offsetPower(1:na/2+1) )
   call ds_fft(offsetCorrelator,offsetPower)
   deallocate(offsetCorrelator)
!	minpow = minval(abs(offsetPower))
!	maxpow = maxval(real(offsetPower))
!write(*,*) minval(abs(offsetPower)), maxval(abs(offsetPower))
offsetPower=real(offsetPower)

	
	
!	write(*,*) "minpow = ", minpow, "maxpow = ", maxval(real(offsetPower))
	do a=1,na/2+1
		if (real(offsetPower(a))<=0) offsetPower(a) = 1.0e-12
	enddo
   !Add this lengths's correlators to Qinv
!write(*,*) minval(real(offsetPower)), maxval(real(offsetPower))
!write(*,*) minval(imag(offsetPower)), maxval(imag(offsetPower))
	call ds_assert(all(ds_isfinite(real(offsetPower))),"NaN or infinity detected in noise power spectrum.  Bad data?")
	if (rank==0) call ds_assert(minval(real(offsetPower))>0.0,"Zero detected in noise power spectrum.  Hmm.")


end subroutine buildOneCorrelatorFromData



!#########################################!
!##										##!
!##	FFT routines calling FFTW library	##!
!##										##!
!#########################################!

subroutine ds_fft(inArr,outArr)
!Perform an FFT of a real input into a complex output (of half the length).
!The input array must be real(dp) and can have any length n.
!The output array must be complex(dpc) and must have length n/2+1
!NOTE WELL:
!	The transform is **NORMALIZED**

complex(dpc), dimension(:) :: outArr
real(8), dimension(:) :: inArr
integer(8) plan
integer arrSize


arrSize = size(inArr)
call ds_assert(size(outArr) == arrSize/2+1,'Wrong sized arrays in ds_fft')

!FFTW stuff - create a plan, execute it, then destroy it. 
!Subsequent plans should be generated much faster.
#ifndef NO_FFTW
call dfftw_plan_dft_r2c_1d(plan,arrSize,inArr,outArr,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
#else
stop 'FFTW DISABLED - NO FFTs POSSIBLE. Recompile without -DNO_FFTW.'
#endif
end subroutine ds_fft


subroutine ds_ifft(inArr,outArr)
!Perform an inverse Fourier transform from a complex array into real data.
!The input array must be complex(dpc) and must have length n/2+1
!The input array must be real(dp) and has length n.
!NOTE WELL:
!	The transform is **NORMALIZED**

complex(dpc), dimension(:) :: inArr
real(8), dimension(:) :: outArr
integer(8) plan
integer arrSize

arrSize = size(outArr)
call ds_assert(size(inArr) == arrSize/2+1,'Wrong sized arrays in ds_ifft')

!FFTW stuff - create a plan, execute it, then destroy it. 
!Subsequent plans should be generated much faster
#ifndef NO_FFTW
call dfftw_plan_dft_c2r_1d(plan,arrSize,inArr,outArr,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
#else
stop 'FFTW DISABLED - NO FFTs POSSIBLE. Recompile without -DNO_FFTW.'
#endif
!normalise
outArr= outArr / real(arrSize,kind=8)

end subroutine ds_ifft



end module ds_simple_prior
