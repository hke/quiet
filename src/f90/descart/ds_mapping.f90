!Once pointing information has been worked out and loaded into the modulescan type:
!  call make_naivecov to produce the naive covariance matrix
!When a naive map is required, initialise q and u map arrays and
!  call add2rhs in a loop to add each modulescan's data to the sum
!  call cov_mult on the maps to finish the map
!When you are finished making naive maps
!  call deallocate_cov

!NB: the inv_Cw matrix will have to be used in place of whiteNoise in simple_prior and driver too.
!Use invCW_mult to apply it, only after call to make_naivecov

!For data read-in, require pointing info (and angles) to be read first, to run make_naivecov.  Can then read-in
!the data and add2rhs as we go along, producing the first map. We must also convert each timestream into offset
!space at read-in. After these 2 ops are done on it, we can forget it. Z op is factorised thus:
! Zd = ( F^t - F^t P naive) d
!where naive includes Cw^-1 ops where appropriate.

!Support for bad diodes: set the diode_flag array of modulescan structure to zero AND set the sigma 
!of that diode to zero.  This will cause it to (i) be rejected from the mapping ops and (ii) have no
!effect on the accepted timestreams when inv_Cw is convolved with them. There might be a numerical 
!problem with using zero though.

module ds_mapping
use ds_types
implicit none

!real(dp),allocatable,dimension(:,:) :: naive_cov !make private later
!integer :: p_npix
integer,parameter :: bad_pixel= -1

contains


!NOISE TOOLS
subroutine make_inv_Cw(modulescan,cross_correlated) !called per modulescan

    type(ds_modulescan),intent(inout) :: modulescan
    logical :: cross_correlated
    
    integer :: i, j, ierror
    logical :: corr
    real(dp),dimension(:,:),pointer :: rho  !only care about the upper off-diagonals
	real(dp),dimension(:),pointer :: sigma
	

   	modulescan%inv_Cw= 0.0_8
    
    rho => modulescan%noise%corr
    sigma => modulescan%noise%sigma

    !make diagonals
    do i= 1,modulescan%ndiodes
        if(modulescan%flags(i)==1) then
            modulescan%inv_Cw(i,i) = sigma(i)**2
		else
			modulescan%inv_Cw(i,i) = 0.0d0 !JZ ???  Will this work with the inversion?
       endif
    enddo
    
    if(cross_correlated) then
        !make off diagonals
        do i= 1,modulescan%ndiodes-1
            if(modulescan%flags(i)==0) cycle
            do j= i+1,modulescan%ndiodes
                if(modulescan%flags(j)==0) cycle    !ignore if bad
                modulescan%inv_Cw(i,j) = sigma(i) * sigma(j) * rho(i,j)
                modulescan%inv_Cw(j,i) = modulescan%inv_Cw(i,j)
            enddo
        enddo
        !inversion
        ierror= 0 
		stop 'JZor DWPS sort this out'
        call dpotrf('U',modulescan%ndiodes,modulescan%inv_Cw,modulescan%ndiodes,ierror)
        call dpotri('U',modulescan%ndiodes,modulescan%inv_Cw,modulescan%ndiodes,ierror)
        do i=1,modulescan%ndiodes-1
            do j=i+1,modulescan%ndiodes
                modulescan%inv_Cw(j,i)= modulescan%inv_Cw(i,j)
            enddo
        enddo
    else
        !simple inversion
        do i= 1,modulescan%ndiodes
            if(modulescan%flags(i)==0) cycle    !ignore if bad
            modulescan%inv_Cw(i,i)= 1.0_8 / modulescan%inv_Cw(i,i)
        enddo
    endif


modulescan%inv_Cw_TT => modulescan%inv_Cw(1:modulescan%ndiodes_T, 1:modulescan%ndiodes_T)
modulescan%inv_Cw_PP => modulescan%inv_Cw(modulescan%ndiodes_T+1:modulescan%ndiodes, modulescan%ndiodes_T+1:modulescan%ndiodes)
modulescan%inv_Cw_TP => modulescan%inv_Cw(1:modulescan%ndiodes_T, modulescan%ndiodes_T+1:modulescan%ndiodes)
end subroutine make_inv_Cw


subroutine invCw_mult_tod(modulescan)
    type(ds_modulescan),intent(inout) :: modulescan

    integer :: i, k
    real(dp),dimension(modulescan%ndiodes) :: vec
    type(ds_timestream),pointer,dimension(:) :: timestream
	integer length_check
    timestream => modulescan%timestreams
    !check timestream lengths
	do k=1,modulescan%ndiodes
		if (modulescan%flags(k)==1) then
			length_check=timestream(k)%nt
			exit
		endif
	enddo
	
    do k=1,modulescan%ndiodes
        if(modulescan%flags(k)==1) call ds_assert(timestream(k)%nt == length_check, "Error in timestream lengths in invCw_mult_tod")
    enddo

    !This is not very optimised
    do i = 1,timestream(1)%nt
        do k = 1,modulescan%ndiodes
            if(modulescan%flags(k)==0) then
                vec(k) = 0.0_8 !set to zero if bad
                cycle
            endif
            vec(k) = timestream(k)%timestream(i)
        enddo
        vec = matmul(moduleScan%inv_Cw,vec)
        do k= 1,modulescan%ndiodes
            if(modulescan%flags(k)==0) cycle    !ignore diode if bad
            timestream(k)%timestream(i) = vec(k)
        enddo
    enddo
    
end subroutine invCw_mult_tod


subroutine invCw_mult_offsets(matinv,modulescan)
    type(ds_modulescan),intent(inout) :: modulescan
    real(dp),dimension(modulescan%ndiodes,modulescan%ndiodes),intent(in) :: matinv

    integer :: i, k
    real(dp),dimension(modulescan%ndiodes) :: vec
    type(ds_offsets),pointer,dimension(:) :: offsets

    offsets => modulescan%offsets_2

    !check timestream lengths
    do k=1,modulescan%ndiodes-1
        call ds_assert(offsets(k)%na == offsets(k+1)%na, "Error in offset stream lengths in invCw_mult_tod")
    enddo

    !This is not very optimised
    do i= 1, offsets(1)%na
        do k= 1,modulescan%ndiodes
            if(modulescan%flags(k)==0) then
                vec(k)= 0.0_8 !set to zero if bad
                cycle
            endif
            
            vec(k) = offsets(k)%values(i)
        enddo
    
        vec = matmul(matinv,vec)
        
        do k = 1,modulescan%ndiodes
            if(modulescan%flags(k)==0) cycle
            offsets(k)%values(i) = vec(k)
        enddo
    enddo
   
end subroutine invCw_mult_offsets



!MAPPING TOOLS

subroutine add2rhs(modulescan,rhs)   !called once per modulescan  !need an allreduce op after end of use
    !NB: do the C_w^-1 op separately
    type(ds_modulescan),intent(in) :: modulescan !one modulescan per module per scan
	type(ds_trimap) :: rhs

    integer :: t, pixel, idiode, i

	if (modulescan%has_p) call ds_assert(associated(modulescan%timestreams_P),"polarization timestream unassociated in call to add2rhs")
	if (modulescan%has_t) call ds_assert(associated(modulescan%timestreams_T), "temperauture timestream unassociated in call to add2rhs")
	
	if (modulescan%has_p) call ds_assert(rhs%has_p, "no polarization map available in call to add2rhs")
	if (modulescan%has_t) call ds_assert(rhs%has_t, "no temperature map available in call to add2rhs")
	
	if (rhs%has_t) call ds_assert(allocated(rhs%T%map),"imap not allocated in add2rhs!")
	if (rhs%has_p) call ds_assert(allocated(rhs%Q%map),"qmap not allocated in add2rhs!")
	if (rhs%has_p) call ds_assert(allocated(rhs%U%map),"umap not allocated in add2rhs!")


	!Temperature Modules
	if (modulescan%has_t) then
		do i = 1,modulescan%ndiodes_T
	        if(modulescan%flags_T(i)==0) cycle  !if diode bad, ignore
			do t = 1, modulescan%ntod
	            pixel = modulescan%pointing(t)
				if (pixel==bad_pixel) cycle
				rhs%T%map(pixel) = rhs%T%map(pixel) + modulescan%timestreams_T(i)%timestream(t)
	        enddo
		enddo
	endif

	!Polarization Modules
	if (modulescan%has_p) then
		do i = 1,modulescan%ndiodes_P
	        if(modulescan%flags_P(i)==0) cycle  !if diode bad, ignore
			do t = 1, modulescan%ntod
	            pixel = modulescan%pointing(t)
				if (pixel==bad_pixel) cycle
	            rhs%Q%map(pixel) = rhs%Q%map(pixel) + cos(2.0_8 * (modulescan%theta(t) + &
	                modulescan%dpsi(i))) * modulescan%timestreams_P(i)%timestream(t)

	            rhs%U%map(pixel)= rhs%U%map(pixel) + sin(2.0_8 * (modulescan%theta(t) + &
	                modulescan%dpsi(i))) * modulescan%timestreams_P(i)%timestream(t)        

	        enddo
		enddo
	endif



end subroutine add2rhs



subroutine add2cov(modulescan, cov) !called once per module (per code)
    type(ds_modulescan),intent(in) :: modulescan
	type(ds_covariance) :: cov
    integer :: t, p, i, j
	real(dp) :: theta_i, theta_j

	if (modulescan%has_t) call ds_assert(cov%has_t, "Naive cov has no temperature in add2cov")
	if (modulescan%has_p) call ds_assert(cov%has_p, "Naive cov has no temperature in add2cov")
	if (modulescan%has_t) then
		do i=1,modulescan%ndiodes_T
			if (modulescan%flags_T(i)==0) cycle
			do j=1,modulescan%ndiodes_T
				if (modulescan%flags_T(j)==0) cycle
				do t=1,modulescan%ntod
					p=modulescan%pointing(t)
					if (p==bad_pixel) cycle
					cov%TT(p) = cov%TT(p) + modulescan%inv_Cw_TT(i,j)
				enddo
			enddo
		enddo
	endif


	if (modulescan%has_p) then
		do i=1,modulescan%ndiodes_P
			if (modulescan%flags_P(i)==0) cycle
			do j=1,modulescan%ndiodes_P
				if (modulescan%flags_P(j)==0) cycle
				do t=1,modulescan%ntod
					
					p=modulescan%pointing(t)
					if (p==bad_pixel) cycle
					theta_i = 2*(modulescan%theta(t) + modulescan%dpsi(i))
					theta_j = 2*(modulescan%theta(t) + modulescan%dpsi(j))
					
					cov%QQ(p) = cov%QQ(p) + 					&
	                    cos(theta_i) * 							&
	                    modulescan%inv_Cw_PP(i,j) * 			&
	                    cos(theta_j)
	
					cov%QU(p) = cov%QU(p) + 					&
	                    cos(theta_i) * 							&
	                    modulescan%inv_Cw_PP(i,j) * 			&
	                    sin(theta_j)

					cov%UU(p) = cov%UU(p) + 					&
					    sin(theta_i) * 							&
					    modulescan%inv_Cw_PP(i,j) * 			&
					    sin(theta_j)
				enddo
			enddo
		enddo
	endif

	if (modulescan%has_p .and. modulescan%has_t) then
		do i=1,modulescan%ndiodes_T
			if (modulescan%flags_T(i)==0) cycle
			do j=1,modulescan%ndiodes_P
				if (modulescan%flags_P(j)==0) cycle
				do t=1,modulescan%ntod
					p=modulescan%pointing(t)
					if (p==bad_pixel) cycle
					theta_j = 2*(modulescan%theta(t) + modulescan%dpsi(j))
					
					cov%TQ(p) = cov%TQ(p) + 					&
	                    modulescan%inv_Cw_TP(i,j) * 			&
	                    cos(theta_j)
	
					cov%TU(p) = cov%TU(p) + 					&
	                    modulescan%inv_Cw_TP(i,j) * 			&
	                    sin(theta_j)
				enddo
			enddo
		enddo
	endif


        

end subroutine add2cov




subroutine map2tod(modulescan,maps) !called once per module per projection
    type(ds_modulescan),intent(inout) :: modulescan
	type(ds_trimap) :: maps
    integer :: t, pixel, i
    
    if (modulescan%has_p) call ds_assert(associated(modulescan%timestreams_P),"polarization timestream unassociated in call to map2tod")
    if (modulescan%has_t) call ds_assert(associated(modulescan%timestreams_T), "temperauture timestream unassociated in call to map2tod")

	if (modulescan%has_p) call ds_assert(maps%has_p, "no polarization map available in call to map2tod")
	if (modulescan%has_t) call ds_assert(maps%has_t, "no temperature map available in call to map2tod")

    if (maps%has_t) call ds_assert(allocated(maps%T%map),"imap not allocated in map2modulescan!")
    if (maps%has_p) call ds_assert(allocated(maps%Q%map),"qmap not allocated in map2modulescan!")
    if (maps%has_p) call ds_assert(allocated(maps%U%map),"umap not allocated in map2modulescan!")



	if (modulescan%has_t) then
	    do t=1, modulescan%ntod
	        pixel = modulescan%pointing(t)
	        if(pixel==bad_pixel) cycle	
	        do i=1,modulescan%ndiodes_T
	            if(modulescan%flags_T(i)==0) cycle  !if diode bad, ignore
	            modulescan%timestreams_T(i)%timestream(t) = maps%T%map(pixel)
	        enddo
		enddo
	endif
	
	if (modulescan%has_p) then
	    do t= 1, modulescan%ntod
	        pixel = modulescan%pointing(t)
	        if(pixel==bad_pixel) cycle	
	        do i=1,modulescan%ndiodes_P
	            if(modulescan%flags_P(i)==0) cycle  !if diode bad, ignore
	            modulescan%timestreams_P(i)%timestream(t) = &
	                cos(2.0_8 * (modulescan%theta(t) + modulescan%dpsi(i))) * maps%Q%map(pixel) + &
	                sin(2.0_8 * (modulescan%theta(t) + modulescan%dpsi(i))) * maps%U%map(pixel)
	        enddo
		enddo
	endif
end subroutine map2tod



subroutine cov_mult(maps,cov)
    !cov is only the lower triangle!
    !map has indices and npix
    type(ds_trimap),intent(inout) :: maps
	type(ds_covariance) :: cov
    real(dp) :: t,q,u

    integer :: p, npix
    
	if (maps%has_p) call ds_assert(cov%has_p,"Inconsistent cov and maps in cov_mult (pol)!")
	if (maps%has_t) call ds_assert(cov%has_t,"Inconsistent cov and maps in cov_mult (temp)!")
                                   
	if (maps%has_t) call ds_assert(cov%npix==maps%T%npix,"Inconsistent cov and maps in cov_mult (I size)!")
	if (maps%has_p) call ds_assert(cov%npix==maps%Q%npix,"Inconsistent cov and maps in cov_mult (I size)!")
	if (maps%has_p) call ds_assert(cov%npix==maps%U%npix,"Inconsistent cov and maps in cov_mult (I size)!")
	
	npix = maps%T%npix
    
	if (cov%has_t .and. cov%has_p) then
		do p=1,npix
			t = maps%T%map(p)
			q = maps%Q%map(p)
			u = maps%U%map(p)
			maps%T%map(p) = cov%TT(p)*t + cov%TQ(p)*q + cov%TU(p)*u
			maps%Q%map(p) = cov%TQ(p)*t + cov%QQ(p)*q + cov%QU(p)*u
			maps%U%map(p) = cov%TU(p)*t + cov%QU(p)*q + cov%UU(p)*u
		enddo
	else if (cov%has_p) then
		do p=1,npix
			q = maps%Q%map(p)
			u = maps%U%map(p)
			maps%Q%map(p) = cov%QQ(p)*q + cov%QU(p)*u
			maps%U%map(p) = cov%QU(p)*q + cov%UU(p)*u
		enddo		
	else if (cov%has_t) then
		do p=1,npix
			t = maps%T%map(p)
			maps%T%map(p) = cov%TT(p)*t
		enddo
	else
		call ds_assert(.false.,"In cov_mult: neither temperature nor polarization found")
	endif

end subroutine cov_mult



subroutine invert_weight(cov)  !to be used once per code. cov= W^-1
    integer :: p
	type(ds_covariance) :: cov
    real(dp) :: det
	real(dp) :: a11,  a12,  a13
	real(dp) :: a21,  a22,  a23
	real(dp) :: a31,  a32,  a33
	real(dp) :: x1,x2,x3
	

	if (cov%has_t) call ds_assert(allocated(cov%TT),"TT not allocated in invert_cov!")
	if (cov%has_p) call ds_assert(allocated(cov%QQ),"QQ not allocated in invert_cov!")
	if (cov%has_p) call ds_assert(allocated(cov%QU),"QU not allocated in invert_cov!")
	if (cov%has_p) call ds_assert(allocated(cov%UU),"UU not allocated in invert_cov!")
	if (cov%has_t .and. cov%has_p) call ds_assert(allocated(cov%TQ),"TQ not allocated in invert_cov!")
	if (cov%has_t .and. cov%has_p) call ds_assert(allocated(cov%TU),"TU not allocated in invert_cov!")

	if (cov%has_t .and. cov%has_p) then  !The covariance is a 3x3 matrix.
	    do p= 1, cov%npix
			a11 = cov%TT(p)
			a12 = cov%TQ(p)
			a13 = cov%TU(p)
			a21 = a12
			a22 = cov%QQ(p)
			a23 = cov%QU(p)
			a31 = a13
			a32 = a23
			a33 = cov%UU(p)
			x1 = a22*a33 - a23*a32
			x2 = a23*a31 - a21*a33
			x3 = a21*a32 - a22*a31
			det = a11*x1 + a22*x2 + a33*x3
			det = 1.0/det
			
			cov%TT(p) = x1*det
			cov%TQ(p) = x2*det
			cov%TU(p) = x3*det
			cov%QQ(p) = (a11*a33-a13*a31)*det
			cov%QU(p) = (a13*a21-a11*a23)*det
			cov%UU(p) = (a11*a22-a12*a21)*det
	    enddo
	else if (cov%has_t) then  !The covariance is a 1x1 matrix
		do p=1,cov%npix
			cov%TT(p) = 1/cov%TT(p)
		enddo
    else if (cov%has_p) then  !The covariance is a 2x2 matrix
		do p=1,cov%npix
			a11 = cov%QQ(p)
			a12 = cov%QU(p)
			a21 = a12
			a22 = cov%UU(p)
			det = a11*a22-a12*a21
			det = 1.0/det
			cov%QQ(p) =  a22*det
			cov%QU(p) = -a12*det
			cov%UU(p) =  a11*det
		enddo
	else
		call ds_assert(.false.,"In invert_weight: neither temperature nor polarization found")	
	endif


end subroutine invert_weight




!I/O routines

subroutine writePixelsFile(cov,filename,nside,maps)
!writes input file for the PPCL program, containing pixel list and rmsnoise per pixel
	character(len=*) :: filename
	integer :: nside
	type(ds_trimap),intent(in) :: maps
	type(ds_covariance) :: cov
	
	character(len=50) :: mess
	integer :: i, npix, err

	!open file
	err=0
	open(unit=5004,file=filename,status="replace",iostat=err)
	if(err.ne.0) then
		write(*,*) "Could not open file:",filename
		return
	endif

	!write file header
	write(mess,*) nside
	mess= "nside="//trim(adjustl(mess))
	write(5004,"(a)",iostat=err) trim(adjustl(mess))

	npix= maps%T%npix
	write(mess,*) npix
	mess= "npix="//trim(adjustl(mess))
	write(5004,"(a)") trim(adjustl(mess))
	if(err.ne.0) then
		write(*,*) "Error writing header in file:",filename
	endif


	!write file
	do i= 1,npix
		!Is the rms noise correct - like the magnitude of the complex intensity of polarised noise
		! rms_noise = sqrt( sigma_qq^2 + sigma_uu^2 )
		write(5004,"('  ',I6,A)") maps%T%indices(i), "NOISE OUTPUT BROKEN" !sqrt(naive_cov(i,1) + naive_cov(i,3)) !naive_cov is a variance
	enddo
	
	!close file
	close(5004)

end subroutine writePixelsFile


subroutine saveCovariance(cov,filename)!,map)
!Saves naive covariance matrix, pixel by pixel
	type(ds_covariance) :: cov
	character(len=*) :: filename
!	type(ds_map),intent(in) :: map
	integer,parameter :: unitno= 12987
	integer :: i

	print*,"Saving covariance to	",trim(adjustl(filename))
	open(unit= unitno, file= trim(adjustl(filename)), form= "formatted", status= "replace")
	if (cov%has_T .and. cov%has_P) then
		do i= 1,cov%npix
			!format: tt,tq,tu,qq, qu, uu
			write(unitno,'(I8,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8)') i,cov%TT(i),cov%TQ(i),cov%TU(i),cov%QQ(i),cov%QU(i),cov%UU(i)
		enddo
	else if (cov%has_T) then
		do i= 1,cov%npix
			!format: tt,tq,tu,qq, qu, uu
			write(unitno,'(I8,E16.8)') i,cov%TT(i)
		enddo

	else if (cov%has_P) then
		do i= 1,cov%npix
			!format: tt,tq,tu,qq, qu, uu
			write(unitno,'(I8,E16.8,E16.8,E16.8)') i,cov%QQ(i),cov%QU(i),cov%UU(i)
		enddo
		
	endif
	
	close(unitno)

end subroutine saveCovariance


!UTILS


end module ds_mapping

