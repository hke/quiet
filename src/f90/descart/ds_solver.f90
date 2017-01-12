module ds_solver
	use ds_precision
	use ds_types
	use ds_utils
	use ds_multidetector
	use ds_simple_prior
	implicit none

	contains

	subroutine pcg(correlator,b,x,map_covariance, tol, sigma_cut,cut_when)
		! solve the system K^{-1} A x = K^{-1} b
		! b= F^T Z y   where y= signal
		! the arrays are distributed, 1 detector per process.
		! Prepare FFTQ must have been run already.
		!inner products are over distribited arrays using inner3
		!b is destroyed (r=>b then r is destroyed)

		!dummy arguments
		type(ds_correlator) :: correlator
		type(ds_modulescan), target, dimension(0:) :: b
		type(ds_modulescan), pointer,dimension(:) :: x
		real(dp), optional :: sigma_cut
		integer, optional :: cut_when
		type(ds_modulescan), pointer, dimension(:) :: r, p, q
		type(ds_covariance) :: map_covariance

		real(dp) :: rho, rho_min1, rho_min2, beta, beta_min1, alpha, err, target_err,initial_err
		integer :: nd, d, idiode, i, j, checkunit
		integer :: cut_iteration
		character(len=125) :: message, checkname
		real(4) startTime,currentTime,timeTaken,timePerIteration
		integer remainingIterationsEstimate

		integer,parameter :: itmax= 1000					
		real(dp) :: tol

		integer :: ierror

		!do d=0,correlator%my_nmodules-1
		!	do idiode= 1,ndiodes_max
		!		if(b(d)%diode_flag(idiode)==0) cycle
		!		print*,correlator%proc,d,idiode,b(d)%offsets(idiode)%az_offset%hits(:)
		!	enddo
		!enddo

		
!		write(*,*) "size(b)=", size(b)
!		call ds_assert(size(b)==correlator%my_nmodules,"wrong sized b in ds_solver")

		nd = size(b)
		correlator%my_nmodules=nd

		allocate(x(0:correlator%my_nmodules-1))
		allocate(p(0:correlator%my_nmodules-1))
		allocate(q(0:correlator%my_nmodules-1))

		do d=0,nd-1
			call setup_moduleScan_ndiodes(x(d),b(d)%ndiodes_T,b(d)%ndiodes_P)
			call setup_moduleScan_ndiodes(p(d),b(d)%ndiodes_T,b(d)%ndiodes_P)
			call setup_moduleScan_ndiodes(q(d),b(d)%ndiodes_T,b(d)%ndiodes_P)
			call copy_modulescan(b(d),q(d),duplicateoffsets=.true.)
			call copy_modulescan(b(d),x(d),duplicateoffsets=.true.)

			!initialise the test solution
			do idiode= 1,x(d)%ndiodes
				if(x(d)%flags(idiode)==0) cycle
				x(d)%offsets_2(idiode)%values= 0.0_8

				if(x(d)%offsets_2(idiode)%azimuth_flag) x(d)%offsets_2(idiode)%az_offset%amplitudes= 0.0_8
			enddo
		enddo


		!print*,'b',b(1)%offsets(1)%az_offset%amplitudes
		!precondition rhs
		r => b
		
		call apply_prec(r,correlator)													!matrix ops - apply_prec

		!initialise scalars
		rho_min1= 0.0_dp
		rho_min2= 0.0_dp
		beta_min1= 0.0_dp
		
!		write(*,*) "my nmodules = ", correlator%my_nmodules
		
!		if (correlator%proc==0) then
!			do i=1,r(0)%na
!				
!			enddo
!		endif

		initial_err = sqrt(inner3(correlator,r,r))
		!set convergence critereon
		target_err= tol * initial_err
		if(correlator%proc==0) then							
			write(message,*),"Target error = ",target_err, ", tol = ", tol										!inform user of iterative progress
			call ds_log(message,ds_feedback_quiet)
		endif	
		
		cut_iteration = 30
		if (present(cut_when)) cut_iteration = cut_when

		if (correlator%proc==0 .and. present(sigma_cut)) then
			write(message,'(A,F5.2,A,I3)') "Will cut out ",sigma_cut,"-sigma residuals at iteration ",cut_iteration
			call ds_log(message,ds_feedback_quiet)
		endif
		
		startTime = ds_get_time()
		do i=1,itmax
			
			
			call assert_it(i,itmax,err)													!terminate if itmax reached
			
			if (present(sigma_cut)) then
				if (i==cut_iteration) call residual_cut(correlator,sigma_cut,r,p,q,x)
			endif
			
			rho_min2= rho_min1
			rho_min1 = inner3(correlator,r,r)											!inner product

			if(i==1) then
				do d=0,nd-1
					call copy_modulescan(r(d),p(d),duplicateoffsets=.true.)			
				enddo
			else	
				beta_min1= rho_min1 / rho_min2
				call update(beta_min1,p,'+',1.0_dp,r)
			endif
			call update(0.0_dp,q,'+',1.0_dp,p)
			
			call apply_A(q,map_covariance,correlator)													!matrix ops - apply_A
			call apply_prec(q,correlator)												!matrix ops - apply_prec

			if(rho_min1 .gt. 1.E-16_8) then
				alpha= rho_min1/ inner3(correlator,p,q)									!inner product
			else
				alpha = 0.0_8
			endif


			call update(1.0_dp,x,'+',alpha,p)
			call update(1.0_dp,r,'-',alpha,q)

			!check convergence			
			err= sqrt(inner3(correlator,r,r))											!sqrt of inner product


			!inform user of iterative progress	
			if(correlator%proc==0) then
				if (i>=3) then							
					currentTime = ds_get_time()
					timeTaken = currentTime - startTime
					!Deal with time overflow
					if (timeTaken<0) timeTaken = timeTaken + ds_get_max_time()
					timePerIteration = timeTaken/i
					remainingIterationsEstimate = ceiling( log(target_err/initial_err)/log(err/initial_err)*i-i)
					write(message,*),i,err/target_err!,"(about",timePerIteration*remainingIterationsEstimate," seconds remain)"
					call ds_log(message,ds_feedback_quiet)
				else
					write(message,*),i,err/target_err
					call ds_log(message,ds_feedback_quiet)
				endif
			endif	

			!print*,'stopping for debug',alpha !alpha should be 
			!stop

			if(err.le. target_err) then													!exit iterations if convergence reached
				exit
			endif

		end do	

		!clear up memory
		call destroy_modulescan(p,deallocateoffsets=.true.)
		call destroy_modulescan(q,deallocateoffsets=.true.)


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


	


	subroutine apply_A(modulescan,covariance,correlator)
		!Applies the LHS of the destriping equation to the offsets in modulescan.
		type(ds_covariance) :: covariance
		type(ds_correlator) :: correlator
		type(ds_modulescan), dimension(0:correlator%my_nmodules-1),intent(inout) :: modulescan

		type(ds_modulescan),pointer,dimension(:) :: dummy
		integer i,d,nd, ierr


		nd= correlator%my_nmodules

		!creates dummy modulescan if we are using Ca^-1 prior 

		if(.not. correlator%traditionalMode) then
			allocate(dummy(0:nd-1))
			do d=0,nd-1
				call setup_moduleScan_ndiodes(dummy(d),modulescan(d)%ndiodes_T,modulescan(d)%ndiodes_P)
				call copy_modulescan(modulescan(d),dummy(d),duplicateoffsets=.true.)
			enddo
		endif

		!Here we apply the destriping matrix, F^t Cw^-1 Z F, where Z is the naive signal removal
		!operator. This saves memory by allocating module timestreams only when needed - the input and
		!outputs are offsets.  In traditional mode, this is the only operation in apply_A

		!print*,'not calling destriping matrix'
		call apply_destriping_matrix(correlator,covariance,modulescan)


		!do d= 0,correlator%my_nmodules-1
		!	do i=1,ndiodes_max
		!		if(dummy(d)%diode_flag(i)==0) cycle
		!!!    	modulescan(d)%offsets(i)%values= modulescan(d)%offsets(i)%values / modulescan(d)%offsets(i)%length
		!		dummy(d)%offsets(i)%az_offset%amplitudes= dummy(d)%offsets(i)%az_offset%amplitudes * dummy(d)%offsets(i)%az_offset%hits
		!		
		!		
		!		print*,modulescan(d)%offsets(i)%az_offset%amplitudes-dummy(d)%offsets(i)%az_offset%amplitudes
		!	enddo
		!enddo

		if (.not. correlator%traditionalMode) then

			!Here we apply the optional prior matrix, Ca^{-1}, if traditionalMode is false.  The matrix becomes: 
			!    A= F^t Cw^-1 Z F +  Ca^{-1}
			!In this case, we add the offsets in the dummy to the offsets in the modulescan.




!!DWPS: 5/3/10 added this to weight azimuth templates vs offsets through priors.
!!It uses the un-priored dummy, without changing it.  It only changes the azmimuth offsets
!!of modulescan array
!!The math is this:  Prior_matrix= ( Ca^-1     0
!                                     0   var(az)^-1 )
!			do d= 0,nd-1
!				do i= 1,ndiodes_max
!					if(modulescan(d)%diode_flag(i)==0) cycle !ignore if bad
!				
!					if(modulescan(d)%offsets(i)%azimuth_flag)  modulescan(d)%offsets(i)%az_offset%amplitudes= &
!						modulescan(d)%offsets(i)%az_offset%amplitudes + &
!						modulescan(d)%inv_Cw(i,i)/4.0_8 * dummy(d)%offsets(i)%az_offset%amplitudes !assumes std_az = 2 std_noise
!				enddo
!			enddo




			call apply_prior(dummy,correlator,.false.) 

			do d= 0,nd-1
				do i= 1,modulescan(d)%ndiodes
					if(modulescan(d)%flags(i)==0) cycle !ignore if bad

					!add the dummy to the modulescan
					modulescan(d)%offsets_2(i)%values = modulescan(d)%offsets_2(i)%values + &
					dummy(d)%offsets_2(i)%values
				enddo
			enddo

			call destroy_modulescan(dummy,deallocateoffsets=.true.)
		endif


	end subroutine





	subroutine apply_prec(modulescan,correlator)
		type(ds_correlator) :: correlator
		type(ds_modulescan), dimension(0:correlator%my_nmodules-1) :: modulescan
		integer :: d, idiode

		!Apply diagonal preconditioner K= hits * sigma_w^-2 for optional azimuth offsets
		do d= 0,correlator%my_nmodules-1
			do idiode= 1, modulescan(d)%ndiodes
				if(modulescan(d)%flags(idiode)==0) cycle !cycle if bad

				if(modulescan(d)%offsets_2(idiode)%azimuth_flag) then
					modulescan(d)%offsets_2(idiode)%az_offset%amplitudes= modulescan(d)%offsets_2(idiode)%az_offset%amplitudes &
					/ (modulescan(d)%inv_Cw(idiode,idiode) * modulescan(d)%offsets_2(idiode)%az_offset%hits)
				endif
			enddo
		enddo

		!Apply the preconditioner matrix K^{-1}, where K= ( sigma^{-2} F^T F + C_a^{-1} )

		if (correlator%traditionalMode) then
			do d= 0,correlator%my_nmodules-1
				do idiode= 1, modulescan(d)%ndiodes
					if(modulescan(d)%flags(idiode)==0) cycle !cycle if bad
					
					modulescan(d)%offsets_2(idiode)%values= modulescan(d)%offsets_2(idiode)%values &
						/ ( modulescan(d)%inv_Cw(idiode,idiode) * modulescan(d)%offsets_2(idiode)%length )
				enddo
			enddo
		else
			!The final argument indicates that this is the preconditioner and not the main correlator 
			call apply_prior(modulescan,correlator,.true.)
		endif

	end subroutine




	subroutine update(coeff_sub,subject,op,coeff_dir,direction)
		type(ds_modulescan),dimension(0:),intent(inout) :: subject
		character(len=1) :: op
		real(dp) :: coeff_sub, coeff_dir
		type(ds_modulescan),dimension(0:),intent(in) :: direction

		integer :: nd, d, idiode

		nd= size(subject)
		do d= 0, nd-1
			do idiode=1,subject(d)%ndiodes
				if(subject(d)%flags(idiode)==0) cycle !ignore if bad

				!update timestream offsets
				if(op=='+') then
					subject(d)%offsets_2(idiode)%values= coeff_sub * subject(d)%offsets_2(idiode)%values + &
					coeff_dir * direction(d)%offsets_2(idiode)%values	

					elseif(op=='-') then
					subject(d)%offsets_2(idiode)%values= coeff_sub * subject(d)%offsets_2(idiode)%values - &
					coeff_dir * direction(d)%offsets_2(idiode)%values	

				else
					print*,"Illegal operator in update (+ or -)"
				endif


				if(subject(d)%offsets_2(idiode)%azimuth_flag) then

					!update azimuth offsets			
					if(op=='+') then
						subject(d)%offsets_2(idiode)%az_offset%amplitudes= coeff_sub * subject(d)%offsets_2(idiode)%az_offset%amplitudes + &
						coeff_dir * direction(d)%offsets_2(idiode)%az_offset%amplitudes

					else
						subject(d)%offsets_2(idiode)%az_offset%amplitudes= coeff_sub * subject(d)%offsets_2(idiode)%az_offset%amplitudes - &
						coeff_dir * direction(d)%offsets_2(idiode)%az_offset%amplitudes				
					endif
				endif
			enddo
		enddo

	end subroutine update

	subroutine residual_cut(correlator,sigma,residuals,a,b,c,d,e,f,g,h,i,j)
		type(ds_correlator),intent(in) :: correlator
		real(dp) :: sigma
		real(dp), dimension(:,:), allocatable :: residual_contribution
		type(ds_modulescan), pointer, dimension(:) :: residuals
		type(ds_modulescan), pointer, dimension(:), optional :: a,b,c,d,e,f,g,h,i,j
		real(dp) :: residual_mean, residual_stdev, residual_limit
		integer m,diode,N
		character(256) :: message

		allocate(residual_contribution(0:correlator%my_nmodules-1,residuals(0)%ndiodes))
		residual_contribution = 0.0
		call inner_product_map(correlator,residuals,residuals,residual_contribution)

		N=0
		do m = 0,correlator%my_nmodules-1
		    do diode= 1,residuals(m)%ndiodes
				if (residuals(m)%flags(diode) .ne. 0) N=N+1
			enddo
		enddo
		
		residual_mean = sum(residual_contribution)/N
		residual_stdev = sqrt( sum((residual_contribution-residual_mean)**2) / N )
		residual_limit = residual_mean + sigma*residual_stdev
		write(message,'(A,F12.4)') "Residual Mean = ",residual_mean
		call ds_log(message,ds_feedback_quiet)
		write(message,'(A,F12.4)') "Residual Std. dev. = ",residual_stdev
		call ds_log(message,ds_feedback_quiet)
		write(message,'(A,F12.4)') "Residual limit = ",residual_limit
		call ds_log(message,ds_feedback_quiet)
		
		do m = 0,correlator%my_nmodules-1
		    do diode= 1,residuals(m)%ndiodes
				if (residual_contribution(m,diode) > residual_limit) then
					write(*,*) "RESIDUAL_CUT",residuals(m)%id,diode,"[",residual_contribution(m,diode),"]"
!					residuals(m)%diode_flag(diode)=0
!					if (present(a)) residuals(m)%diode_flag(diode)=0
!					if (present(b)) residuals(m)%diode_flag(diode)=0
!					if (present(c)) residuals(m)%diode_flag(diode)=0
!					if (present(d)) residuals(m)%diode_flag(diode)=0
!					if (present(e)) residuals(m)%diode_flag(diode)=0
!					if (present(f)) residuals(m)%diode_flag(diode)=0
!					if (present(g)) residuals(m)%diode_flag(diode)=0
!					if (present(h)) residuals(m)%diode_flag(diode)=0
!					if (present(i)) residuals(m)%diode_flag(diode)=0
!					if (present(j)) residuals(m)%diode_flag(diode)=0
				endif
			enddo
		enddo
		deallocate(residual_contribution)
	end subroutine residual_cut


end module ds_solver
