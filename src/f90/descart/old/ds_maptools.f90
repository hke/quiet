module ds_maptools
use ds_types
use ds_precision
use ds_utils
!use ds_azimuth
implicit none
contains


subroutine addMap(source, dest)
type(ds_map) :: source, dest
call ds_assert(source%npix==dest%npix, "Incompatible maps added")
dest%map = dest%map + source%map
end subroutine addMap

subroutine projectMapOntoTimestream(map,timestream,pointing)
type(ds_map) :: map
type(ds_timestream) :: timestream
type(ds_detpointing) :: pointing
integer t,p
	call ds_assert(pointing%npix == map%npix, 'Pointing npix is not the same as map npix')
	call ds_assert(pointing%nt == timestream%nt, 'Pointing nt is not the same as timestream nt')

	do t = 1,pointing%nt
		p = pointing%pixel(t)
		timestream%timestream(t) = map%map(p)
	enddo
	
end subroutine projectMapOntoTimestream



subroutine addOffsetToTimestream(offsets,timestream)
type(ds_timestream), intent(inout) :: timestream
type(ds_offsets), intent(in) :: offsets
integer a,basisStart, basisEnd

	do a = 1, offsets%na
		basisStart = (a-1)*offsets%length+1
		basisEnd = a*offsets%length
!		timestream%timestream(basisStart:basisEnd) = offsets%values(a)
    		timestream%timestream(basisStart:basisEnd) = timestream%timestream(basisStart:basisEnd) + offsets%values(a)
	enddo
end subroutine addOffsetToTimestream

subroutine subtractOffsetFromTimestream(offsets,timestream)
type(ds_timestream), intent(inout) :: timestream
type(ds_offsets), intent(in) :: offsets
integer a,basisStart, basisEnd

	do a = 1, offsets%na
		basisStart = (a-1)*offsets%length+1
		basisEnd = a*offsets%length
!		timestream%timestream(basisStart:basisEnd) = offsets%values(a)
    		timestream%timestream(basisStart:basisEnd) = timestream%timestream(basisStart:basisEnd) - offsets%values(a)
	enddo
end subroutine subtractOffsetFromTimestream




subroutine deprojectTimestreamOntoOffset(timestream,offsets)
!This is the operator F^T in the PCG.
	type(ds_timestream), intent(in) :: timestream
	type(ds_offsets), intent(inout) :: offsets
	integer(dp) a,basisStart,basisEnd
	do a = 1,offsets%na
		basisStart = 1+offsets%length*(a-1)
		basisEnd = a*offsets%length
		offsets%values(a) = sum( timestream%timestream(basisStart:basisEnd) )
	enddo
	
!	if(offsets%azimuth_flag) call az_deproject_tod_onto_offset(timestream,offsets%az)
	
end subroutine 	deprojectTimestreamOntoOffset
	
	
	
!Alternate naive binning approach	(wrong place for them)
	
!subroutine sum_naive_rhs(correlator,pointing,timestream,rhs_q,rhs_u)
!!pointing should include cos2phi and a sin2phi array
!type(ds_correlator) :: correlator
!type(ds_map) :: rhs_q, rhs_u
!type(ds_detpointing), dimension(0:correlator%my_ndet/2-1) :: pointing
!type(ds_timestream),dimension(0:correlator%my_ndet-1) :: timestream
!integer :: d, p, i, ierr
!
!rhs_q%map= 0.0_8
!rhs_u%map= 0.0_8
!
!!Sum naive RHS for all timestreams on this process
!do d= 0, correlator%my_ndet-1
!    p= d/2  !which pointing array
!    do i= 1,pointing(p)%nt
!        rhs_q%map(pointing(p)%pix(i))= rhs_q%map(pointing(p)%pix(i)) + pointing(p)%cos2phi(i) * timestream(d)%timestream(i)
!        rhs_q%map(pointing(p)%pix(i))= rhs_q%map(pointing(p)%pix(i)) + pointing(p)%sin2phi(i) * timestream(d)%timestream(i)
!    enddo
!enddo
!
!!mpi_reduce the maps onto the 
!ierr=0
!call mpi_reduce(MPI_IN_PLACE,rhs_q%map,rhs_q%npix,MPI_DOUBLE_PRECISION,MPI_SUM,0,correlator%comm,ierr)
!call mpi_reduce(MPI_IN_PLACE,rhs_u%map,rhs_u%npix,MPI_DOUBLE_PRECISION,MPI_SUM,0,correlator%comm,ierr)
!
!end subroutine sum_naive_rhs
	
	
!subroutine sum_naive_weight_matrix(correlator,pointing,mat_qq,mat_qu,mat_uq,mat_uu)
!!similar to the above except summing a 2*2 matrix for each pixel
!type(ds_correlator) :: correlator
!type(ds_detpointing), dimension(0:correlator%my_ndet/2-1) :: pointing
!type(ds_timestream),dimension(0:correlator%my_ndet-1) :: timestream()
!real(dp),allocatable,dimension(:),intent(out) :: mat_qq,mat_qu,mat_uq,mat_uu
!
!!allocate the naive weight matrices
!allocate(mat_qq(1:pointing(0)%npix))	
!allocate(mat_qu(1:pointing(0)%npix))
!allocate(mat_uq(1:pointing(0)%npix))
!allocate(mat_uu(1:pointing(0)%npix))
!
!!initialise the matrices
!mat_qq= 0.0_8
!mat_qu= 0.0_8
!mat_uq= 0.0_8
!mat_uu= 0.0_8
!
!!Sum naive naive weight matrix for all timestreams on this process
!do d= 0, correlator%my_ndet-1
!    p= d/2  !which pointing array
!    do i= 1,pointing(p)%nt
!        mat_qq(pointing(p)%pix(i))= mat_qq(pointing(p)%pix(i)) + pointing(p)%cos2phi(i) * timestream(d)%timestream(i)
!        mat_uu(pointing(p)%pix(i))= mat_uu(pointing(p)%pix(i)) + pointing(p)%sin2phi(i) * timestream(d)%timestream(i)
!    enddo
!enddo
!
!end subroutine sum_naive_weight_matrix

end module ds_maptools





