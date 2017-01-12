module ds_maptools
use ds_types
use ds_precision
use ds_utils
implicit none
contains


subroutine addMap(source, dest)
type(ds_map) :: source, dest
integer i
call ds_assert(source%npix==dest%npix, "Incompatible maps added")
do i=lbound(source%map,1), ubound(source%map,1)
	dest%map(i) = dest%map(i) + source%map(i)
enddo
end subroutine addMap



subroutine subtractMap(subtractee, subtractor)
	!subtractee= subtractee- subtractor
	type(ds_map),intent(inout) :: subtractee
	type(ds_map),intent(in) :: subtractor
	integer i
	
	call ds_assert(subtractee%npix==subtractor%npix, "Incompatible maps subtracted")
	
	do i=lbound(subtractee%map,1), ubound(subtractee%map,1)
		subtractee%map(i) = subtractee%map(i) - subtractor%map(i)
	enddo
end subroutine subtractMap

subroutine subtractTriMap(a, b)
	!  a -= b
	type(ds_trimap),intent(inout) :: a
	type(ds_trimap),intent(in) :: b
	
	if (a%has_T) then
		call ds_assert(b%has_T, "Tried to subtracted incompatible trimaps (no Temperature in one map).")
		call subtractMap(a%T,b%T)
	endif
	
	if (a%has_P) then
		call ds_assert(b%has_P, "Tried to subtracted incompatible trimaps (no Temperature in one map).")
		call subtractMap(a%Q,b%Q)
		call subtractMap(a%U,b%U)
	endif

end subroutine subtractTriMap



!subroutine projectMapOntoTimestream(map,timestream,pointing)
!type(ds_map) :: map
!type(ds_timestream) :: timestream
!type(ds_detpointing) :: pointing
!integer t,p
!	call ds_assert(pointing%npix == map%npix, 'Pointing npix is not the same as map npix')
!	call ds_assert(pointing%nt == timestream%nt, 'Pointing nt is not the same as timestream nt')
!
!	do t = 1,pointing%nt
!		p = pointing%pixel(t)
!		timestream%timestream(t) = map%map(p)
!	enddo
!	
!end subroutine projectMapOntoTimestream
!
!

subroutine addOffsetToTimestream(offsets,timestream)
type(ds_timestream), intent(inout) :: timestream
type(ds_offsets), intent(inout) :: offsets
integer a,basisStart, basisEnd

	do a = 1, offsets%na
		basisStart = (a-1)*offsets%length+1
		basisEnd = a*offsets%length
!!		timestream%timestream(basisStart:basisEnd) = offsets%values(a)
    		timestream%timestream(basisStart:basisEnd) = timestream%timestream(basisStart:basisEnd) + offsets%values(a)
	enddo
	
!	if(offsets%azimuth_flag) call az_add_offset_to_tod(offsets%az_offset,timestream)

	
end subroutine addOffsetToTimestream

!subroutine subtractOffsetFromTimestream(offsets,timestream)
!type(ds_timestream), intent(inout) :: timestream
!type(ds_offsets), intent(in) :: offsets
!integer a,basisStart, basisEnd
!
!	do a = 1, offsets%na
!		basisStart = (a-1)*offsets%length+1
!		basisEnd = a*offsets%length
!!		timestream%timestream(basisStart:basisEnd) = offsets%values(a)
!    		timestream%timestream(basisStart:basisEnd) = timestream%timestream(basisStart:basisEnd) - offsets%values(a)
!	enddo
!end subroutine subtractOffsetFromTimestream


subroutine subtractTimestreamFromOffset(timestream,offsets)
	type(ds_timestream), intent(in) :: timestream
	type(ds_offsets), intent(inout) :: offsets
	integer(dp) a,basisStart,basisEnd
!	type(ds_az_offset) :: buffer
	do a = 1,offsets%na
		basisStart = 1+offsets%length*(a-1)
		basisEnd = a*offsets%length
		offsets%values(a) = offsets%values(a) - sum( timestream%timestream(basisStart:basisEnd) )
	enddo
	
!    if(offsets%azimuth_flag) then
!        call az_copy(offsets%az_offset,buffer)
!        buffer%amplitudes= 0.0_dp
!        call az_deproject_tod_onto_offset(timestream,buffer)
!        offsets%az_offset%amplitudes = offsets%az_offset%amplitudes - buffer%amplitudes
!        call az_destroy(buffer)
!    endif
	
end subroutine



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
	
!	if(offsets%azimuth_flag) call az_deproject_tod_onto_offset(timestream,offsets%az_offset)
	
end subroutine 	deprojectTimestreamOntoOffset
	
	


end module ds_maptools





