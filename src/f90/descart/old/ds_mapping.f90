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
!space at read-in. After these 2 ops are done on it, we can forget it. - what about Z op?

!Support for bad diodes: set the diode_flag array of modulescan structure to zero AND set the sigma 
!of that diode to zero.  This will cause it to (i) be rejected from the mapping ops and (ii) have no
!effect on the accepted timestreams when inv_Cw is convolved with them. There might be a numerical 
!problem with using zero though.

module ds_mapping
use ds_types
!use ds_multidetector
implicit none

real(dp),allocatable,dimension(:,:) :: naive_cov !make private later
integer :: p_npix

contains

!subroutine to make the Npix 2*2 covariance matrices !maybe not here - need mpi calls and correlator type
!subroutine make_naivecov(npix,sigma,rho,correlator,modulescan)
!    type(ds_correlator) :: correlator
!    type(ds_modulescan),dimension(0:correlator%my_nmodules-1) :: modulescan
!    real(dp),intent(in) :: sigma(0:correlator%my_nmodules-1,4), rho(0:correlator%my_nmodules,6) !only half of the off-diagonals
!    integer :: npix
!    integer :: m, ierror
!
!    !initialise naive cov
!    call allocate_cov(npix)
!    
!    !go through list of modulescans
!    do m= 0,correlator%my_nmodules-1
!        !build Cw^-1
!        call make_inv_Cw(sigma(m,:),rho(m,:),modulescan(m)%inv_Cw)
!        
!        !add this modulescan's info to the covariance
!        call add2cov(modulescan(m))
!    enddo
!    
!    !call mpiallreduce for naivecov
!    call mpi_allreduce(MPI_IN_PLACE,naive_cov,3*p_npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)	!send to master proc, sum, then redistribute final result to every proc
!    call ds_assert(ierror==0, 'Error in mpi_allreduce in make_naivecov')	
!
!    !invert the 2*2 weight matrix for each pixel to get the cov
!    call invert_weight()
!    
!end subroutine


!NOISE TOOLS
subroutine make_inv_Cw(sigma,rho,mat) !called per modulescan
    real(dp),dimension(4) :: sigma
    real(dp),dimension(6) :: rho  !rho is only the upper off diagonals
    real(dp),dimension(4,4) :: mat
    
    integer :: i, j, k, ierror
    logical :: corr
    
   	mat= 0.0_8
    
    !make diagonals
    do i= 1,4
        mat(i,i)= sigma(i)**2
    enddo
    
#warning !must check that the lower/upper conventions here are consistent and work
    if(maxval(rho).ne.0.0_8) then
        !make off diagonals
        do i= 1,3
            do j= i+1,4
                if(i==1) k= j-i
                if(i==2) k= j-i + 3
                if(i==3) k= 6
                mat(i,j)= sigma(i) * sigma(j) * rho(k)
                mat(j,i)= mat(i,j)
            enddo
        enddo
        !inversion
        ierror= 0 
        call dpotrf('U',4,mat,4,ierror)
        call dpotri('U',4,mat,4,ierror)
        do i=1,3
            do j=i+1,4
                mat(j,i)= mat(i,j)
            enddo
        enddo
    else
        !simple inversion
        do i= 1,4
            mat(i,i)= 1.0_8 / mat(i,i)
        enddo
    endif
    
end subroutine make_inv_Cw


subroutine invCw_mult(matinv,time1,time2,time3,time4)
    real(dp),dimension(4,4),intent(in) :: matinv
    type(ds_timestream),intent(inout) :: time1, time2, time3, time4

    integer :: i
    real(dp),dimension(4) :: vec

    !check timestream lengths
    call ds_assert(time1%nt == time2%nt,"Error in timestream lengths in invCw_mult")
    call ds_assert(time2%nt == time3%nt,"Error in timestream lengths in invCw_mult")
    call ds_assert(time3%nt == time4%nt,"Error in timestream lengths in invCw_mult")

    !This is not very optimised
    do i= 1, time1%nt
        vec= (/time1%timestream(i),time2%timestream(i),time3%timestream(i),time4%timestream(i)/)
        vec= matmul(matinv,vec)
        time1%timestream(i)= vec(1)
        time2%timestream(i)= vec(2)
        time3%timestream(i)= vec(3)
        time4%timestream(i)= vec(4)
    enddo
   
end subroutine invCw_mult


!MAPPING TOOLS

subroutine add2rhs(idiode,modulescan,tod,rhs_q,rhs_u)   !called once per diode  !need an allreduce op after end of use
    !NB: do the C_w^-1 op separately
    integer,intent(in) :: idiode                 !diode number
    type(ds_modulescan),intent(in) :: modulescan !one modulescan per module per scan
    type(ds_timestream),intent(in) :: tod
    type(ds_map),intent(inout) :: rhs_q, rhs_u

    integer :: t, tpix

    if(modulescan%diode_flag(idiode)==0) return !if diode is bad, remove from sum

    call ds_assert(allocated(rhs_q%map),"rhs_q not allocated in add2rhs!")
    call ds_assert(allocated(rhs_u%map),"rhs_u not allocated in add2rhs!")

    do t= 1, modulescan%ntod
        tpix= modulescan%pointing(t)
        
        rhs_q%map(tpix)= rhs_q%map(tpix) + cos(2.0_8 * (modulescan%theta(t) + modulescan%dpsi(idiode))) * tod%timestream(t)
        rhs_u%map(tpix)= rhs_u%map(tpix) + sin(2.0_8 * (modulescan%theta(t) + modulescan%dpsi(idiode))) * tod%timestream(t)        
    enddo

end subroutine add2rhs



subroutine add2cov(modulescan) !called once per module (per code)
    type(ds_modulescan),intent(in) :: modulescan
    integer :: t, tpix, i, j

    call ds_assert(allocated(naive_cov),"naive_cov not allocated in add2cov!")

    !loop through the TOD. For the pixel, add to sum the permutations
    do t= 1, modulescan%ntod
        tpix= modulescan%pointing(t)
        
        do i= 1,4
            if(modulescan%diode_flag(i)==1) then   !only include if diode is good
                do j= 1,4
                    if(modulescan%diode_flag(j)==1) then   !only include if diode is good
                        !QQ element
                        naive_cov(tpix,1)= naive_cov(tpix,1) + &
                                (cos(2.0_8*(modulescan%theta(t) + modulescan%dpsi(i))) * &
                                 modulescan%inv_Cw(i,j) * &
                                 cos(2.0_8*(modulescan%theta(t) + modulescan%dpsi(j))) )
                        
                        !QU element
                        naive_cov(tpix,2)= naive_cov(tpix,2) + &
                                (cos(2.0_8*(modulescan%theta(t) + modulescan%dpsi(i))) * &
                                modulescan%inv_Cw(i,j) * &
                                sin(2.0_8*(modulescan%theta(t) + modulescan%dpsi(j))) )
                         
                        !UU element
                        naive_cov(tpix,3)= naive_cov(tpix,3) + &
                                (sin(2.0_8*(modulescan%theta(t) + modulescan%dpsi(i))) * &
                                 modulescan%inv_Cw(i,j) * &
                                sin(2.0_8*(modulescan%theta(t) + modulescan%dpsi(j))) )                         
                    endif
                enddo
            endif
        enddo
        
    enddo

end subroutine add2cov




subroutine map2tod(idiode,modulescan,qmap,umap,tod) !called once per diode per projection
    integer,intent(in) :: idiode   !diode number within the module
    type(ds_modulescan),intent(in) :: modulescan
    type(ds_map),intent(in) :: qmap, umap
    type(ds_timestream),intent(inout) :: tod

    integer :: t, tpix

    if(modulescan%diode_flag(idiode)==0) return !if diode is bad, remove from projection

    call ds_assert(allocated(qmap%map),"qmap not allocated in map2modulescan!")
    call ds_assert(allocated(umap%map),"umap not allocated in map2modulescan!")
    
    do t= 1, modulescan%ntod
        tpix= modulescan%pointing(t)
        
        !project the map onto this tod
        tod%timestream(t) = &
            cos(2.0_8 * (modulescan%theta(t) + modulescan%dpsi(idiode))) * qmap%map(tpix) + &
            sin(2.0_8 * (modulescan%theta(t) + modulescan%dpsi(idiode))) * umap%map(tpix)
    enddo

end subroutine map2tod



subroutine cov_mult(q,u)
    !cov is only the lower triangle!
    !map has indices and npix
    type(ds_map),intent(inout) :: q, u
    real(dp) :: qbuf, ubuf
    integer :: i
    
    call ds_assert(allocated(naive_cov),"naive_cov not allocated in cov_mult!")
    call ds_assert(q%npix==p_npix .and. u%npix==p_npix,"Wrong sized maps or weight matrix in cov_mult")
    
    do i=1,p_npix
        qbuf= q%map(i)
        ubuf= u%map(i)
        q%map(i)= naive_cov(i,1) * qbuf + naive_cov(i,2) * ubuf
        u%map(i)= naive_cov(i,2) * qbuf + naive_cov(i,3) * ubuf
    enddo

end subroutine cov_mult



subroutine invert_weight()  !to be used once per code. cov= W^-1
    integer :: i
    real(dp) :: det, a, b, d

    call ds_assert(allocated(naive_cov),"naive_cov not allocated in invert_cov!")
    do i= 1, p_npix
      !analytic solution faster than cholesky!
       a= naive_cov(i,1)
       b= naive_cov(i,2)
       d= naive_cov(i,3)
       det= a*d - b**2
       naive_cov(i,:) = (/ d,-b,a /) / det
    enddo
    
end subroutine invert_weight




!UTILS

subroutine allocate_cov(npix)
    integer :: npix
    call deallocate_cov()
    allocate(naive_cov(npix,3))
    naive_cov= 0.0_8
    p_npix= npix
end subroutine allocate_cov

subroutine deallocate_cov()
    if(allocated(naive_cov)) deallocate(naive_cov)
end subroutine deallocate_cov

end module ds_mapping

