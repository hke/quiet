module ds_multidetector
use ds_precision
use ds_types
use ds_utils
use ds_mapping !remove this when not needed any more
use ds_maptools
!use ds_azimuth
!use mpi
implicit none

include 'mpif.h'
include 'fftw3.f'


!test


!#error This module is not finished


type ds_correlator
!integer fftlength
integer :: comm, stokes_comm
complex(dpc), allocatable, dimension(:) :: space
!complex(dpc), allocatable, dimension(:,:) :: masterSpace  !, Pinv, Qinv
real(dp), allocatable, dimension(:) :: beta
real(dp), allocatable, dimension(:) :: realspace
integer, allocatable, dimension(:) :: owner, my_det
integer, allocatable, dimension(:) :: ntod  !JZ The TOD length for all of the detectors
integer :: my_ndet, ndistinct, my_nmodules

!The idea is that pair i will only have a correlator if it needs to (ie isCorrelated is true). fftlength gives its length
logical :: traditionalMode
!For each detector pair (i,j) this Qinv variable is the inverse of  the fourier space (cross-)power spectrum of
!the two detector offset-space noise.
!Note that since the correlator is symmetric we do note use a 2x2 array of Qinvs.
!Instead Q_ij(w) should be in Qinv( triangleIndex(i,j,n),w )
!The variable should only be allocated for the master node, since the other nodes send their data to
!the master for the summation.  In addition it is a zero-based array.
!realspace is an npix sized space useful in the mapmaking phase.
integer proc, det
integer nproc, ndet  !The number of detectors and the number of processors. For the moment, these
!are the same, but that may change.
end type ds_correlator

type ds_timeCorrelator
integer :: p, i, j
real(dp),allocatable,dimension(:) :: vals
end type

type ds_fdCa
integer :: p
complex(dpc),allocatable,dimension(:) :: vals
end type
 

contains


subroutine destroy_correlator(C)
type(ds_correlator) :: C
integer i

!write this function
i=1
call ds_assert(i==0,"Called unwritten function destroy_correlator")

end subroutine destroy_correlator

subroutine inner_product_map(correlator,a,b,map)
	type(ds_correlator) :: correlator
	type(ds_modulescan), dimension(0:correlator%my_nmodules-1) :: a, b
	real(dp), dimension(0:,1:) :: map
	integer i,d

	call ds_assert(size(map,1)==correlator%my_nmodules,"Wrong size map 1 passed to inner_product_map")
	call ds_assert(size(map,2)==a(0)%ndiodes,"Wrong size map 2 passed to inner_product_map")

	map = 0.0_dp
	

	do d = 0,correlator%my_nmodules-1
	    do i= 1, a(d)%ndiodes
	    	if(a(d)%flags(i)==0) cycle !ignore if bad
	    	map(d,i) = sum(a(d)%offsets_2(i)%values * b(d)%offsets_2(i)%values)
!	    	if(a(d)%offsets(i)%azimuth_flag) map(d,i) = map(d,i) + az_inner_product(a(d)%offsets(i)%az_offset , b(d)%offsets(i)%az_offset)
if(a(d)%offsets_2(i)%azimuth_flag) stop 'azimuths broken in inner_product_map'
		enddo
	enddo
   
end subroutine inner_product_map

function inner3(correlator,a,b)
type(ds_correlator) :: correlator
type(ds_modulescan), dimension(0:correlator%my_nmodules-1) :: a, b
real(dp) :: inner3
integer :: ierror, d, i

inner3 = 0.0_dp
do d = 0,correlator%my_nmodules-1
    do i= 1, a(d)%ndiodes
    	if(a(d)%flags(i)==0) cycle !ignore if bad
 
    	inner3 = inner3 + sum(a(d)%offsets_2(i)%values * b(d)%offsets_2(i)%values)
!    	if(a(d)%offsets(i)%azimuth_flag) inner3= inner3 + &
!    	az_inner_product(a(d)%offsets(i)%az_offset , b(d)%offsets(i)%az_offset)
	enddo
	
enddo

call mpi_allreduce(MPI_IN_PLACE,inner3,1,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)	!send to master proc, sum, then redistribute final result to every proc
call ds_assert(ierror==0, 'Error in mpi_allreduce in inner')	

end function inner3



 
 subroutine ds_mpiTest(retcode, mesg)
 character(*) mesg
 character(len=1024) mpi_mesg
 integer retcode, mpi_mesg_len, ierror
 
 if (retcode .ne. 0) then
 call Mpi_error_string(retcode, mpi_mesg, mpi_mesg_len,ierror)
 write(*,*) "Error"
 write(*,*) mesg
 write(*,*) "MPI ERROR:"
 write(*,*) mpi_mesg(1:mpi_mesg_len+1)

 stop
endif 

 end subroutine ds_mpiTest
 

subroutine apply_destriping_matrix(correlator,cov,out)  !applies destriping matrix to modulescan 'out'.
!Add this to output of apply_prior if traditionalMode is false
!The routine does not require a dummyTimestream input and does not create a full dummy timestream 
!internally - only one modules timestreams are allocated per processsor at any given time.
type(ds_correlator),intent(in) :: correlator
type(ds_modulescan),intent(inout),dimension(0:) :: out
type(ds_covariance) :: cov
!these dummy maps held internally
type(ds_trimap) :: buffer
integer :: itod, m, ierror, i, npix
npix = cov%npix
!npix= p_npix !this is dirty because p_npix is public variable in ds_mapping

!initialise map sums
!call prepareMap(qbuf,npix)  !sets it to zero
!call prepareMap(ubuf,npix)
call prepareTriMap(buffer,npix, cov%has_t, cov%has_p)

do m= 0,correlator%my_nmodules-1

    !project this offset onto a timestream of appropriate size
    do i=1,out(m)%ndiodes
        if(out(m)%flags(i)==0) cycle
		call check_modulescan_initialized(out(m), "Scan uninitialized in apply_destriping matrix")
        call prepareTimestream(out(m)%timestreams(i),out(m)%ntod)
        call addOffsetToTimestream(out(m)%offsets_2(i),out(m)%timestreams(i))
    enddo
    
    !apply Cw^-1
    call invCw_mult_tod(out(m))
    
    !deproject the timestream into the offset
    !the offsets now contain F^t Cw^-1 F a
    do i=1,out(m)%ndiodes
        if(out(m)%flags(i)==0) cycle !ignore if bad
        call deprojectTimestreamOntoOffset(out(m)%timestreams(i),out(m)%offsets_2(i)) !F^T operator
    enddo
    
    !sum into map space
    call add2rhs(out(m),buffer)

    !deallocate timestream memory
    do i= 1,out(m)%ndiodes
        if(out(m)%flags(i)==0) cycle
        call destroyTimestream(out(m)%timestreams(i))
    enddo
enddo

!mpi_allreduce the map sums then average
	call share_maps(buffer, correlator%comm)

!average the T/Q/U map sums
call cov_mult(buffer,cov)

do m=0,correlator%my_nmodules-1
    !project map onto this timestream of appropriate size
    
    do i= 1,out(m)%ndiodes
        if(out(m)%flags(i)==0) cycle
        call prepareTimestream(out(m)%timestreams(i),out(m)%ntod)
    enddo
    
    !project map onto timestream
    call map2tod(out(m),buffer)

    call invCw_mult_tod(out(m))

    !subtract deprojected timestream from the offsets
    !The timestream is Cw^-1 P NAIVE(F a)
    !the offsets are  F^t Cw^-1 F a
    !so the subtraction is F^t Cw^-1 Z F a
    do i= 1,out(m)%ndiodes
        if(out(m)%flags(i)==0) cycle !ignore if bad
        call subtractTimestreamFromOffset( out(m)%timestreams(i),out(m)%offsets_2(i) )
        call destroyTimestream(out(m)%timestreams(i))
    enddo
    

enddo

!clear the map buffer memory
call destroyTriMap(buffer)

end subroutine apply_destriping_matrix



subroutine makeNaiveMap(correlator,modulescan,maps,cov,have_T,have_P)
!Makes naive map by projecting the offsets in the modulescan object. This means the whole timestream
!is not allocated - only one modulescan instance at a time.  All processes know the global maps at
!the end.


type(ds_correlator),intent(in) :: correlator
type(ds_modulescan),dimension(0:correlator%my_nmodules-1),intent(inout) :: modulescan
type(ds_trimap),intent(out) :: maps
type(ds_covariance) :: cov
logical :: have_P, have_T

integer :: m,i, ierror

call prepareTriMap(maps,cov%npix,have_T,have_P)

	!project this offset onto a timestream of appropriate size
do m=0,correlator%my_nmodules-1
    do i=1,modulescan(m)%ndiodes
        if(modulescan(m)%flags(i)==0) cycle
        call prepareTimestream(modulescan(m)%timestreams(i),modulescan(m)%ntod)
        call addOffsetToTimestream(modulescan(m)%offsets_2(i),modulescan(m)%timestreams(i))
    enddo
    
    !apply Cw^-1
    call invCw_mult_tod(modulescan(m))

	!add to the accumulated sum map
    call add2rhs(modulescan(m),maps)
    
    !deallocate timestream memory
    do i= 1,modulescan(m)%ndiodes
        if(modulescan(m)%flags(i)==0) cycle
        call destroyTimestream(modulescan(m)%timestreams(i))
    enddo
enddo

!reduce the maps
if (have_P) call mpi_allreduce(MPI_IN_PLACE,maps%Q%map,maps%Q%npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)
if (have_P) call mpi_allreduce(MPI_IN_PLACE,maps%U%map,maps%U%npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)
if (have_T) call mpi_allreduce(MPI_IN_PLACE,maps%T%map,maps%T%npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)

!apply the single pixel covariance matrices
call cov_mult(maps,cov)


end subroutine makeNaiveMap





!subroutine to make the Npix 3*3 covariance matrices 
subroutine make_naivecov(npix, cov, comm, modulescan, have_T, have_P)
	type(ds_covariance) :: cov
	integer nscan
	integer comm
    type(ds_modulescan),dimension(0:) :: modulescan
	logical :: have_T, have_P

    integer :: npix
    integer :: m, i, ierror


	nscan = size(modulescan)

    !initialise naive cov
    call prepare_covariance(cov, npix, have_T, have_P)
    
    !go through list of modulescans
	
    do m= 0,nscan-1
!		write(*,*) modulescan(m)%flags
!		write(*,*) modulescan(m)%inv_cw
    	!assert that unflagged diagonals of inv_Cw~=0. Zero is their default and is a unphysical 
    	!value, so if zero they must have not been pre-set
    	do i=1,modulescan(m)%ndiodes
			call ds_assert(modulescan(m)%flags(i)==0 .or. modulescan(m)%inv_Cw(i,i) .gt. 0, "Cw^-1 not valid in make_naivecov")
    	end do
        !add this modulescan's info to the covariance
        call add2cov(modulescan(m), cov)
    enddo
    
!    call copy_covariance(cov,modulescan%covariance)  !Copy the covariance into the scan so that we retain the unshared, uninverted covariance.

    !call mpiallreduce for naivecov
	call share_covariance(cov, comm)

    !invert the 2*2 weight matrix for each pixel to get the cov
    call invert_weight(cov)
    
end subroutine make_naivecov


subroutine share_covariance(cov, comm)
	type(ds_covariance) :: cov
	integer comm
	integer ierror

	if (cov%has_T) then
		call mpi_allreduce(MPI_IN_PLACE,cov%TT,cov%npix,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)
    	call ds_assert(ierror==0, 'Error in TT in share_covariance')	
	endif
	
	if (cov%has_P) then
		call mpi_allreduce(MPI_IN_PLACE,cov%QQ,cov%npix,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)
    	call ds_assert(ierror==0, 'Error in QQ in share_covariance')	
		call mpi_allreduce(MPI_IN_PLACE,cov%QU,cov%npix,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)
    	call ds_assert(ierror==0, 'Error in QU in share_covariance')	
		call mpi_allreduce(MPI_IN_PLACE,cov%UU,cov%npix,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)
    	call ds_assert(ierror==0, 'Error in UU in share_covariance')	
	endif
	if (cov%has_T .and. cov%has_P) then
		call mpi_allreduce(MPI_IN_PLACE,cov%TQ,cov%npix,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)
    	call ds_assert(ierror==0, 'Error in TQ in share_covariance')	
		call mpi_allreduce(MPI_IN_PLACE,cov%TU,cov%npix,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)
    	call ds_assert(ierror==0, 'Error in TU in share_covariance')	
	endif
end subroutine share_covariance

subroutine share_maps(maps, comm)
	type(ds_trimap) :: maps
	integer :: comm
	integer ierror
	
	if (maps%has_T) then
		call mpi_allreduce(MPI_IN_PLACE, maps%T%map, maps%T%npix, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
		call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
	endif
	
	if (maps%has_P) then
		call mpi_allreduce(MPI_IN_PLACE, maps%Q%map, maps%Q%npix, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
		call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
		call mpi_allreduce(MPI_IN_PLACE, maps%U%map, maps%U%npix, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierror)
		call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
	endif
end subroutine share_maps

	
end module
 
