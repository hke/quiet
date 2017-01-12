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

type ds_correlation_data
#warning SHTOP ... This type is not finished yet
integer n
end type ds_correlation_data

type ds_masterSpace
   complex(dpc),allocatable,dimension(:) :: vec
end type ds_masterSpace

type ds_timestream_list
   integer :: ntimestreams, ntod
   integer,allocatable,dimension(:) :: list
end type ds_timestream_list

type ds_Qinv
   logical :: isCorrelated
   integer :: fftlength
   complex(dpc),allocatable,dimension(:) :: Qinv, Pinv
end type ds_Qinv

type ds_correlator
!integer fftlength
integer :: comm, stokes_comm
complex(dpc), allocatable, dimension(:) :: space
!complex(dpc), allocatable, dimension(:,:) :: masterSpace  !, Pinv, Qinv
real(dp), allocatable, dimension(:) :: beta
real(dp), allocatable, dimension(:) :: realspace
integer, allocatable, dimension(:) :: owner, my_det
integer, allocatable, dimension(:) :: ntod  !JZ The TOD length for all of the detectors
integer :: my_ndet, ndistinct

!The idea is that pair i will only have a correlator if it needs to (ie isCorrelated is true). fftlength gives its length
type(ds_Qinv),allocatable,dimension(:) :: pair
type(ds_masterSpace),allocatable,dimension(:) :: masterSpace
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

#warning !deal with this parameter
!regularisation parameter used to weight Ca prior matrix.  The M-L derivation suggests this be unity,
!but in Sutton et al 2009, we used 2 as it reduced error in maps. What is the optimal value?
real(dp),parameter :: regularisationParameter= 1.0_dp  

contains


subroutine destroy_correlator(C)
type(ds_correlator) :: C
integer i
if (allocated(C%space)) deallocate(C%space)
!if (allocated(C%Qinv)) deallocate(C%Qinv)
if (allocated(C%masterSpace)) deallocate(C%masterSpace)
if (allocated(C%beta)) deallocate(C%beta)
!if (allocated(C%Pinv)) deallocate(C%Pinv)
if (allocated(C%realspace)) deallocate(C%realspace)
if (allocated(C%pair)) then
   do i= 0,C%ndet-1
      if(allocated(C%pair(i)%Qinv)) deallocate(C%pair(i)%Qinv)
      if(allocated(C%pair(i)%Pinv)) deallocate(C%pair(i)%Pinv)
   enddo
   deallocate(C%pair)
endif
C%proc = -1
C%det = -1
C%ndet = -1
C%nproc = -1
C%ndistinct= -1
C%comm = -1
C%stokes_comm = -1
end subroutine destroy_correlator

function triangleIndex(i,j,nd)
!Since the a-space covariance matrix is symmetric positive definite (in its blocks) we only need to store 
!the upper triangle of blocks. This function converts a pair of detector indices into a 1D index
!that moves along the upper triangle.  It expects, and returns, 0-based indices.  The rationale for this
!is that mpi ranks are zero-based.
!It only accepts i<=j.
integer :: triangleIndex,i,j,nd
call ds_assert((i <= j .and. j<nd .and. i>=0),'Bad indices in triangleIndex')
triangleIndex = i*nd-((i+1)*i)/2+j
end function triangleIndex

subroutine squareIndex(ndet,p,i,j)
!Function to return square matrix indices (i,j) from triangle index p. Needed by prepareFFTQ to decide which 
!timeCorrelator pair file to upload in loadCorrelator. 
!NOTES:
!(p,i,j) are zero-based
!are there precision issues with using floor?
integer, intent(in) :: ndet
integer,intent(in) :: p
integer,intent(out) :: i,j
integer :: N, k
real(dp) :: l

!return the row number
N= (ndet**2 + ndet )/2				!total number of triangle elements
l= real(N-1-p,kind=8)										!number of triagle elements greater than p
k= floor((sqrt(8.0_8*l+1.0_8)-1.0_8)/2.0_8)					!find largest integer i where (i^2+i)/2 <=l (number of rows of higher index than i)
i= ndet-k-1										!convert to row number (-1 makes the rows zero-based)
!return the column number
j= ndet - (k+1) + ((k+1)*(k+2)/2 - anint(l)) -1	!k+1 is also the number of triangle elements on row k

!print*,'p=',p,'k=',k,'h=',(k+1)*(k+2)/2,'l=',anint(l),'i=',i,'j=',j
end subroutine squareIndex

function numberOfPairs(n)
!return the number of detector pairs - the size of the upper triangular part of a matrix.
integer n,numberOfPairs
numberOfPairs = n*(n+1)/2
end function numberOfPairs



subroutine loadCorrelator(correlatorFilenameRoot,i,j,source,realSize)
!JZ Read in a real-space timestream (auto-)?correlation function and transform it to the first row 
!of the fourier-space.
!The two detectors being correlated have indices i and j, the length of the correlation is assumed
!to be the full length of the data, realSize, and the result is put into the destination space.
!The filenames are assumed to be correaltorFilenameRoot_i_j

integer i,j,size,realSize,t
character(*) :: correlatorFilenameRoot
character(512) :: filename
character(20) :: ipart, jpart
!real(dp), dimension(realSize) :: destination
real(dp), dimension(0:realSize-1) :: source
integer lun
!Build the filename for the correlator.
write(ipart,*) i
write(jpart,*) j
filename = correlatorFilenameRoot // "_" // trim(adjustl(ipart)) // "_" // trim(adjustl(jpart))


!read in the file into our temporary space
lun = ds_get_lun()
!open( unit=lun, file=filename, form='binary', status='old')
!do t=0,realSize-1
!read(lun) source(t)
!enddo
open(unit=lun, file=filename,form='BINARY', status='OLD')
read(lun) source
close(lun)

!Fourier transform into the target and correct FFTW's normalization - DWPS: this is now done in offsets space by transform_timestream_correlator_to_offset
!so we now just output source directly
!call ds_fft(source, destination)
!destination = destination / sqrt(real(realSize))
!Is there something else you have to do here? like square it?
!#warning Is there something else to do here?
end subroutine loadCorrelator












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
call dfftw_plan_dft_r2c_1d(plan,arrSize,inArr,outArr,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

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
call dfftw_plan_dft_c2r_1d(plan,arrSize,inArr,outArr,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

!normalise
outArr= outArr / real(arrSize,kind=8)

end subroutine ds_ifft







function skyToDetectorQ(theta,Q,U) result(detQ)
real(dp) :: theta, Q, U, detQ
detQ = cos(2*theta)*Q + sin(2*theta)*U
end function skyToDetectorQ

function skyToDetectorU(theta,Q,U) result(detU)
real(dp) :: theta, Q, U, detU
detU = -sin(2*theta)*Q + cos(2*theta)*U
end function skyToDetectorU


function detectorToSkyQ(theta,Q,U) result(skyQ)
real(dp) :: theta, Q, U, skyQ
skyQ = cos(2*theta)*Q - sin(2*theta)*U
end function detectorToSkyQ

function detectorToSkyU(theta,Q,U) result(skyU)
real(dp) :: theta, Q, U, skyU
skyU = sin(2*theta)*Q + cos(2*theta)*U
end function detectorToSkyU



	


function inner(correlator,vec1,vec2)
!perform inner product on distributed array

type(ds_correlator) :: correlator
real(8),dimension(:) :: vec1, vec2
real(8) :: inner																	!function result
integer :: ierror


inner= sum(vec1 * vec2)														!sum the element-wise products for each processor

call mpi_allreduce(MPI_IN_PLACE,inner,1,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)	!send to master proc, sum, then redistribute final result to every proc
call ds_assert(ierror==0, 'Eror in mpi_allreduce in inner')	

end function



function inner3(correlator,a,b)
type(ds_correlator) :: correlator
type(ds_offsets), dimension(0:correlator%my_ndet-1) :: a, b
real(dp) :: inner3
integer :: ierror, d

inner3 = 0.0_dp
do d = 0,correlator%my_ndet-1
	inner3 = inner3 + sum(a(d)%values * b(d)%values)
	
!	if(a(d)%azimuth_flag) inner3= inner3 + az_inner_product(a(d)%az)
enddo

call mpi_allreduce(MPI_IN_PLACE,inner3,1,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)	!send to master proc, sum, then redistribute final result to every proc
call ds_assert(ierror==0, 'Error in mpi_allreduce in inner')	

end function inner3




subroutine invert_fourier_long_matrix(matrix,fftlength,npair)
integer :: npair, fftlength
complex(dpc),dimension(fftlength,0:npair-1) :: matrix

integer :: a, ndet, ierror
complex(dpc) :: determ, q00
complex(dpc),allocatable,dimension(:) :: inversionSpace

ndet= floor(sqrt(2.*real(npair)))


!The master node now has all the Q and must invert them to Qinv
!This is done in Fourier space, inverting each mode separately
!calls to hardware optimised lapack routines:
!zpptrf - computes double complex cholesky factorization of matrix Q
!zpptri - computes double complex inverse of Q from output of zpptrf
!these routines run on packed upper triangle of Q
allocate(inversionSpace(0:npair-1))
ierror=0
do a=2,fftlength
	inversionSpace = matrix(a,:)
	call zpptrf('L',ndet,inversionSpace,ierror)
	if(ierror.ne.0) print*,ierror,a,inversionSpace
	call ds_assert(ierror==0, "Cholesky factorization failure 1")
	call zpptri('L',ndet,inversionSpace,ierror)
	call ds_assert(ierror==0, "Cholesky factorization failure 2 ")
	matrix(a,:)= inversionSpace
enddo
!We set the zeroth mode to zero by hand
matrix(1,:)= 0.0_8
deallocate(inversionSpace)

end subroutine invert_fourier_long_matrix





subroutine create_mpi_fdCa_type(thetype,bsize,mesg_mpi_t,ierror)

type(ds_fdCa)   thetype
integer		 bsize
integer      mesg_mpi_t , ierror
 
integer      ierr 
integer      block_lengths(2)
integer      displacements(2)
integer(8)      address(3)
integer      typelist(2)

!  First specify the types.
typelist(1) = MPI_INTEGER
typelist(2) = MPI_DOUBLE_COMPLEX

 
!  Specify the number of elements of each type.
block_lengths(1) = 1
block_lengths(2) = bsize

!  Calculate the displacements of the members relative to indata.
call MPI_Address(thetype,   address(1), ierr)
call ds_assert(ierr==0, 'Failed to get address 1')
call MPI_Address(thetype%p, address(2), ierr)
call ds_assert(ierr==0, 'Failed to get address 2')
call MPI_Address(thetype%vals, address(3), ierr)
call ds_assert(ierr==0, 'Failed to get address 3')
displacements(1) = address(2) - address(1)
displacements(2) = address(3) - address(1)

!  Build the derived datatype  
call MPI_TYPE_STRUCT(2, block_lengths, displacements,typelist, mesg_mpi_t, ierr)
 call ds_assert(ierr==0,"Build of type struct failed in create_mpi_fdca_type")
!  Commit it -- tell system we'll be using it for communication.
call MPI_TYPE_COMMIT(mesg_mpi_t, ierr)
 call ds_assert(ierr==0,"Commit failed in create_mpi_fdca_type")

ierror=ierr

end subroutine create_mpi_fdCa_type

 
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
 

 
!This might be a good place to build the global hit count map.


subroutine addHitPixels(correlator, pointing, inputMapHit,currentPixel,nt)
  type(ds_detpointing) :: pointing
  type(ds_correlator) :: correlator
integer ::  currentPixel, t
  integer, dimension(*) :: inputMapHit
integer :: nt
  do t = 1,nt
    if (inputMapHit(pointing%pixel(t)) == 0) then
      inputMapHit(pointing%pixel(t)) = currentPixel
      currentPixel = currentPixel +1
    endif
  enddo
end subroutine addHitPixels




subroutine applyCorrelator(offsets,correlator,isPreconditioner)
type(ds_correlator) :: correlator
type(ds_offsets), dimension(0:correlator%my_ndet-1) :: offsets
integer status(MPI_STATUS_SIZE)
logical isPreconditioner
real(dp) :: masterSpace
integer i,j,p, imode, d, my_j
integer ierror



!Perform the Fourier space summation: w_i = SUM_j Q_ij v_j
!The sum is performed on the master node.  Each of the w_i is a vector.
!One by one the other nodes send their segment v_j to the master node, into the array
!correlator%space.The master node adds that contribution to the result in correlator%masterSpace.
!to w_i.  Once the summation is complete, the master node returns the segments to the other nodes.

if (correlator%proc == 0) then								!the master

!allocate and initialise masterSpace instances for summation
do i= 0,correlator%ndet-1
   p= triangleIndex(i,i,correlator%ndet) !use this to get correct fftlength
   allocate(correlator%masterSpace(i)%vec( correlator%pair(p)%fftlength ))
   correlator%masterSpace(i)%vec = 0
enddo

my_j=0
do j=0,correlator%ndet-1    								!loop through all the detectors
        !allocate correlator%space correctly for this timestream
        p= triangleIndex(j,j,correlator%ndet)
        allocate(correlator%space(correlator%pair(p)%fftlength))

	if (correlator%owner(j)==0) then  						!if this detector is owned by the master node 
		call ds_fft(offsets(my_j)%values,correlator%space)	                !then do the FFT myself 
		my_j=my_j+1  								!while keeping track of how many of my detectors I have done
		do i=0,correlator%ndet-1    						!and perform the matrix multiplication
			if(i.ge.j) then
				p = triangleIndex(j,i,correlator%ndet)
			else
				p=  triangleIndex(i,j,correlator%ndet)
			endif
                        if(correlator%pair(p)%isCorrelated) then
                           if (isPreconditioner) then
                              correlator%masterSpace(i)%vec = correlator%masterSpace(i)%vec + correlator%pair(p)%Pinv * correlator%space
                           else !forward matrix application
                              correlator%masterSpace(i)%vec = correlator%masterSpace(i)%vec + correlator%pair(p)%Qinv*correlator%space !correlator%space is v_j
                           endif
                        endif
                 enddo
              else 						
           !otherwise this detector is not owned by the master.  get it from whoever does own it.
                p= triangleIndex(j,j,correlator%ndet)  !get the p for the fftlength
		call MPI_RECV(correlator%space, correlator%pair(p)%fftlength, MPI_DOUBLE_COMPLEX, correlator%owner(j),j, correlator%comm,status,ierror)
		call ds_assert(ierror==0, 'Error in MPI_RECV')
		do i=0,correlator%ndet-1							!and then do the same matrix multiplication 
			if(i.ge.j) then
				p = triangleIndex(j,i,correlator%ndet)
			else
				p=  triangleIndex(i,j,correlator%ndet)
			endif
                        if(correlator%pair(p)%isCorrelated) then
                           if (isPreconditioner) then
                              correlator%masterSpace(i)%vec = correlator%masterSpace(i)%vec + correlator%pair(p)%Pinv * correlator%space
                           else !forward matrix application
                              correlator%masterSpace(i)%vec = correlator%masterSpace(i)%vec + correlator%pair(p)%Qinv * correlator%space !correlator%space is v_j
                           endif
                        endif
		enddo !inner matrix loop
     endif
    deallocate(correlator%space)
enddo ! outer matrix loop
else  !not the master processor - they just need to send the data.
	do d=0,correlator%my_ndet-1
		j = correlator%my_det(d)
                p= triangleIndex(j,j,correlator%ndet)
                allocate(correlator%space(correlator%pair(p)%fftlength))
		call ds_fft(offsets(d)%values,correlator%space)	
		call MPI_SEND(correlator%space,correlator%pair(p)%fftlength,MPI_DOUBLE_COMPLEX,0,j,correlator%comm,ierror)
		call ds_assert(ierror==0, 'Error in MPI_SEND')
                deallocate(correlator%space)
	enddo
endif

!We have now finished the matrix multiplication on the master node and send out the results (where appropriate)
call MPI_Barrier(correlator%comm,ierror)
call ds_assert(ierror==0, 'Error in MPI_barrier')

!master sends out summed array.
if (correlator%proc == 0) then

my_j=0
do j=0,correlator%ndet-1
	if (correlator%owner(j)==0) then
		call ds_ifft(correlator%masterSpace(j)%vec, offsets(my_j)%values)
		my_j = my_j + 1
	else
                p= triangleIndex(j,j,correlator%ndet)
		call MPI_SEND(correlator%masterSpace(j)%vec,correlator%pair(p)%fftlength,MPI_DOUBLE_COMPLEX,correlator%owner(j),j,correlator%comm,ierror)
		call ds_assert(ierror==0, 'Error in MPI_SEND - 2')
	endif
        deallocate(correlator%masterSpace(j)%vec)
enddo
else !slave node
	do j=0,correlator%my_ndet-1
!this needs to be changed
                p= triangleIndex( correlator%my_det(j), correlator%my_det(j), correlator%ndet)
                allocate(correlator%space(correlator%pair(p)%fftlength))
		call MPI_RECV(correlator%space, correlator%pair(p)%fftlength, MPI_DOUBLE_COMPLEX, 0,correlator%my_det(j), correlator%comm,status,ierror)
		call ds_assert(ierror==0, 'Error in MPI_RECV - 2')
		call ds_ifft(correlator%space,offsets(j)%values)
                deallocate(correlator%space)
	enddo
endif

end subroutine applyCorrelator





subroutine makeNaiveMap(timestreams,qmap,umap,pointings,correlator)
!Put a naive map from the timestream and pointing arguments into the map argument.
!The map should already have been allocated.
!A naive map is simply an average of all the measurements of a pixel
!with no regard for noise.
	type(ds_correlator), intent(in) :: correlator
	type(ds_timestream), intent(in), dimension(0:correlator%my_ndet-1) :: timestreams
	type(ds_detpointing), intent(in), dimension(0:correlator%my_ndet/2-1) :: pointings
	type(ds_map), intent(inout) :: qmap, umap
!	real(dp),dimension(0:correlator%my_ndet-1),intent(in) :: var
	integer t,p, ierror,status(MPI_STATUS_SIZE), n, i, my_pairs, d, dq, du
	real(dp) :: c,s

	my_pairs = correlator%my_ndet / 2

	call ds_assert(mod(correlator%my_ndet,2)==0, 'Need even number of detectors - q and u')
	qmap%map = 0
	umap%map = 0

	do d = 0,my_pairs-1
		dq = det_qindex(d)
		du = det_uindex(d)
        do t = 1,pointings(d)%nt
                p = pointings(d)%pixel(t)
                qmap%map(p) = qmap%map(p) + detectorToSkyQ(pointings(d)%theta(t), timestreams(dq)%timestream(t), timestreams(du)%timestream(t) )
                umap%map(p) = umap%map(p) + detectorToSkyU(pointings(d)%theta(t), timestreams(dq)%timestream(t), timestreams(du)%timestream(t) )
        enddo
    enddo


	!mpi_reduce the Q and U maps from all the modules separately
	if(correlator%ndet>2) then
		call mpi_allreduce(MPI_IN_PLACE, qmap%map, qmap%npix, MPI_DOUBLE_PRECISION, MPI_SUM, correlator%comm, ierror)
		call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
		call mpi_allreduce(MPI_IN_PLACE, umap%map, umap%npix, MPI_DOUBLE_PRECISION, MPI_SUM, correlator%comm, ierror)
		call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
	endif

	!average the Q/U map sums
	qmap%map = qmap%map / pointings(0)%GlobalHitMap%map
	umap%map = umap%map / pointings(0)%GlobalHitMap%map

end subroutine makeNaiveMap

!function det_qindex(d)
!integer det_qindex, d
!det_qindex = 2*d
!end function det_qindex

!function det_uindex(d)
!integer det_uindex, d
!det_uindex = 2*d+1
!end function det_uindex




subroutine removeSignalNaive(timestreams,pointings,qmap,umap,correlator)
!The operator Z in the PCG.  This turns a timestream into a naive map, re-projects
!the map into a timestream, and then deducts the projected timestream from the original one.
!The map argument is an input to this subroutine so we do now have to keep reallocating it.
!It will be over-written.
!	type(ds_timestream), intent(inout), dimension(correlator%my_ndet) :: timestreams
!	type(ds_detpointing), intent(in), dimension(correlator%my_ndet/2) :: pointings
	type(ds_map), intent(inout) :: qmap,umap
	type(ds_correlator), intent(in) :: correlator
	integer t,p, ierror, status(MPI_STATUS_SIZE), d, dq, du, my_pairs
	real(dp) :: c,s
	type(ds_timestream), dimension(0:correlator%my_ndet-1) :: timestreams
	type(ds_detpointing), intent(in), dimension(0:correlator%my_ndet/2-1) :: pointings


	my_pairs = correlator%my_ndet / 2

!	call ds_assert(size(timestreams)==correlator%my_ndet, 'Correlator%my_ndet and size of timestreams vector mismatch')
!	call ds_assert(size(pointings)==my_pairs, 'Correlator%my_ndet and size of pointings vector mismatch')
!	call ds_assert(lbound(timestreams,1) == 0, 'Expecting a zero-based array for timestreams in removeSignalNaiveMulti')
!	call ds_assert(lbound(pointings,1) == 0, 'Expecting a zero-based array for timestreams in removeSignalNaiveMulti')		
	call ds_assert(mod(correlator%my_ndet,2)==0, 'Need even number of detectors - q and u')

  call makeNaiveMap(timestreams,qmap,umap,pointings,correlator)


	do d = 0,my_pairs-1
		dq = det_qindex(d)
		du = det_uindex(d)
        do t = 1,pointings(d)%nt
                p = pointings(d)%pixel(t)
				timestreams(dq)%timestream(t) = timestreams(dq)%timestream(t) - skyToDetectorQ(pointings(d)%theta(t), qmap%map(p), umap%map(p) )
				timestreams(du)%timestream(t) = timestreams(du)%timestream(t) - skyToDetectorU(pointings(d)%theta(t), qmap%map(p), umap%map(p) )
        enddo
    enddo

end subroutine removeSignalNaive






subroutine makeSharedPointingContinuous(correlator,pointings,maxIndex,filename,originalIndices)
!This subroutine resets the pointing of a ds_pointing instance so that the 
!pixel indices it uses run from 1..npix (and can therefore be used as array indices)
!The maxIndex argument specifies the number of pixels in the input class,
!for example if you are converting a set of healpix indices then maxIndex should be 12*Nside^2
  type(ds_correlator) :: correlator
  type(ds_detpointing), dimension(0:correlator%my_ndet/2-1) :: pointings
type(ds_detpointing) :: tempPointing  
integer :: maxIndex, currentPixel,p,t
integer ierror, outfile, full_npix, d, my_d, q, new_p, owner
integer max_ntod, u,p2
integer , dimension(MPI_STATUS_SIZE) :: status
integer, allocatable, dimension(:) :: originalIndices
  integer,allocatable, dimension(:) :: inputMapHit
  character(*) :: filename
  character(128) :: filename2
character(128) :: message

  integer next_proc
  
  integer det_pairs, my_pairs
  integer nSelf,nRec
  max_ntod = maxval(correlator%ntod)
  det_pairs = correlator%ndet/2 
  my_pairs = correlator%my_ndet/2

  
  !Send all the pointings from each node to the master, which then enumerates the hit pixels in them.  
  currentPixel=1
  
nSelf=0
nRec=0





if (correlator%proc == 0) then
	allocate(inputMapHit(maxIndex)) !JZ - do healpix indices start at zero or one?
    inputMapHit = 0
	call preparePointing(tempPointing, max_ntod, pointings(0)%npix)
	my_d=0
	do d=0,det_pairs-1
		call ds_log(message,ds_feedback_debug)
		q = det_qindex(d)
		owner = correlator%owner(q)
!		write(message,*) "Working on pointing for pair", d, "of ", det_pairs, " owner = ",owner
		if (owner == correlator%proc) then
           call ds_assert(pointings(my_d)%nt == correlator%ntod(q), "pointing ntod mismatch")
           call addHitPixels(correlator,pointings(my_d),inputMapHit,currentPixel,correlator%ntod(q))
           my_d = my_d + 1
           nSelf=nSelf+1
		else
!			write(*,*) correlator%proc," receiving ",q," from ",owner," size = ",correlator%ntod(q)
           call MPI_Recv(tempPointing%pixel,correlator%ntod(q),MPI_INTEGER, owner,q,correlator%comm,status,ierror)
!		write(*,*) "received"
           call ds_mpiTest(ierror,"Recv error 1 in makeSharedPointingContinuous")
           call addHitPixels(correlator,tempPointing,inputMapHit,currentPixel,correlator%ntod(q))	
!		write(*,*) "added"
           nRec=nRec+1
		endif
	enddo
else	
	my_d=0
	do d=0,det_pairs-1
        q = det_qindex(d)
		owner = correlator%owner(q)
		if (owner == correlator%proc) then
           call ds_assert(pointings(my_d)%nt == correlator%ntod(q), "pointing ntod mismatch 2")
!			write(*,*) correlator%proc," sending ",q
			call MPI_Send(pointings(my_d)%pixel, pointings(my_d)%nt, MPI_INTEGER, 0, q, correlator%comm, status, ierror)
           call ds_mpiTest(ierror,"Send error 1 in makeSharedPointingContinuous")
			my_d = my_d + 1
		endif
	enddo
endif


full_npix = currentPixel-1



!The master has now calculated the total number of pixels.  It shares that info with the other nods
call MPI_Bcast(full_npix, 1, MPI_INTEGER, 0, correlator%comm,ierror)
call ds_mpiTest(ierror, "BCast error in makeSharedPointingContinuous")

do d=0,my_pairs-1
	pointings(d)%npix = full_npix
enddo
!JZ Each processor needs only one copy of the globalHitMap.  Maybe it should really be an element of ds_correlator
call prepareMap(pointings(0)%globalHitMap, pointings(0)%npix)


!The inputMapHit array, which has size of a full map (ie 12*nside**2 in the case of healpix)
!now contains a global indexing scheme.  We now do the above loop again, but this time after the
!root gets each pointing it re-indexes it and sends it back.
!We do not want to pass around the inputMapHit because it can be very large - 3 Gig in the case of nside=8192



if (correlator%proc == 0) then
	call preparePointing(tempPointing, max_ntod, pointings(0)%npix)
	my_d=0
	do d=0,det_pairs-1
!		write(message,*) "Returning pointing for pair", d, " of ", det_pairs
!		call ds_log(message,ds_feedback_debug)
		q = det_qindex(d)
		owner = correlator%owner(q)
		if (owner == correlator%proc) then
			do t=1,pointings(my_d)%nt
				new_p = inputMapHit(pointings(my_d)%pixel(t))
				pointings(my_d)%pixel(t) = new_p
				pointings(0)%globalHitMap%map(new_p) = pointings(0)%globalHitMap%map(new_p) + 1
			enddo
			my_d = my_d + 1
		else
			call MPI_Recv(tempPointing%pixel,correlator%ntod(q),MPI_INTEGER, owner,q,correlator%comm,status,ierror)
			call ds_mpiTest(ierror,"Recv error 3 in makeSharedPointingContinuous")
			do t=1,correlator%ntod(q)
				new_p = inputMapHit(tempPointing%pixel(t))
				tempPointing%pixel(t) = new_p
				pointings(0)%globalHitMap%map(new_p) = pointings(0)%globalHitMap%map(new_p) + 1
			enddo
			call MPI_Send(tempPointing%pixel,correlator%ntod(q),MPI_INTEGER, owner,q,correlator%comm,status,ierror)
			call ds_mpiTest(ierror,"Send error 3 in makeSharedPointingContinuous")
		endif
	enddo
	call destroyPointing(tempPointing)
else	 !This is not the master node.
	my_d=0
	do d=0,det_pairs-1
		q = det_qindex(d)
		owner = correlator%owner(q)
		if (owner == correlator%proc) then
			call MPI_Send(pointings(my_d)%pixel, pointings(my_d)%nt, MPI_INTEGER, 0, q, correlator%comm, status, ierror)
			call ds_mpiTest(ierror,"send error 4 in makeSharedPointingContinuous")
			call MPI_Recv(pointings(my_d)%pixel, pointings(my_d)%nt, MPI_INTEGER, 0, q, correlator%comm, status, ierror)			
			call ds_mpiTest(ierror,"Recv error 4 in makeSharedPointingContinuous")
			my_d = my_d + 1
		endif
	enddo
endif

!Prepare the reverse indexing so we can convert back at the end. -surely output this to the driver for the master node? -DWPS  Yes - JZ
if (correlator%proc == 0) then
	allocate(originalIndices(currentPixel-1))!We add +1 when we finish the final pixel-DWPS
	do p=1,maxIndex
		if(inputMapHit(p) .gt. 0 ) originalIndices(inputMapHit(p)) = p
	enddo
!	outfile = ds_get_lun()
!	open(unit=outfile, file="pixels.old.dat")
!	do p=1,pointings(0)%npix
!		write(outfile,*) originalIndices(p)
!	enddo
	deallocate(inputMapHit)
!	close(outfile)
!We no longer divide by two as we have only one pointing per QU pair
endif
  
!#warning !added code to output hit map
!if(correlator%proc==0) then
!print*,"writing hit map"
!open(unit=1016,file="global_hits.old.dat",status="replace")
!do p= 1,pointings(0)%npix
!   write(1016,*) originalIndices(p),pointings(0)%globalHitMap%map(p)
!enddo
!close(1016)
!stop
!endif


call MPI_Bcast(pointings(0)%globalHitMap%map, pointings(0)%npix, MPI_DOUBLE_PRECISION, 0, correlator%comm,ierror)
call ds_mpiTest(ierror,"Bcast in makeSharedPointingContinuous")		
  
  
end subroutine


subroutine assign_detectors(correlator,nd)
type(ds_correlator) :: correlator
integer d,nd,np,c,i,ierror

call MPI_COMM_SIZE(correlator%comm,correlator%nproc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_size')
call MPI_COMM_RANK(correlator%comm,correlator%proc,ierror)
call ds_assert(ierror==0,'Error in mpi_comm_rank')

allocate(correlator%owner(0:nd-1))
correlator%ndet = nd
!Loop through the Q/U pairs assigning them to processors
c=0
correlator%my_ndet = 0
do i=0,nd/2-1  !loop through the Q/U pairs
correlator%owner(det_qindex(i)) = c
correlator%owner(det_uindex(i)) = c
if (c==correlator%proc) correlator%my_ndet = correlator%my_ndet+2
c = c+1
if (c==correlator%nproc) c = 0
enddo

allocate(correlator%my_det(0:correlator%my_ndet-1))

c=0
do i=0,nd/2-1
if (correlator%owner(det_qindex(i))==correlator%proc)  then
	call ds_assert(correlator%owner(det_uindex(i))==correlator%proc, "Processor seems to own Q but not U")
	correlator%my_det(det_qindex(c)) = det_qindex(i)
	correlator%my_det(det_uindex(c)) = det_uindex(i)
	c = c+1
	endif
enddo

end subroutine assign_detectors




subroutine prepareFFTQ2(correlator,lc,nd,comm,sigma,fknee,alpha,correlationData,nyquist,no_correlation_in)
!Populate the structure correlator, a ds_correlator instance, to be run on the correlator comm.
!The correlatorFilenameRoot variable describes where the timestream correlators are stored. (see above)
!The offsets are need primarily for the size and number of the offsets.
!The master node stores a lot more data than the other nodes, since it performs the summation 
!when the correlator is applied.  This is to avoid a many-to-many MPI operation.
!DWPS: na and nt must now be arrays
!DWPS: lc must be given as an argument
!DWPS: the TOD and nts must have been cut to integer multiples of lc ALREADY
type(ds_correlation_data) :: correlationData
type(ds_correlator) :: correlator
integer nd, c
logical, optional :: no_correlation_in
logical  :: no_correlation
real(dp),dimension(0:nd-1) :: whiteNoise, sigma, alpha
real(dp),dimension(0:nd-1) :: fknee
!real(dp),dimension(0:nd-1,0:nd-1) :: rho
integer nf
integer,dimension(0:correlator%ndet-1) ::  na
!real(dp), allocatable, dimension(:) :: timeCorrelator
type(ds_timeCorrelator),allocatable,dimension(:) :: timeCorrelator
integer :: comm
integer n,i,j,p, lc, npair
integer ierror, newcomm_size
real(dp) :: nyquist
if (present(no_correlation_in)) then
	no_correlation = no_correlation_in
else
	no_correlation = .false.
endif

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


call ds_assert(.false., "correlated mode is disabled until you rewrite the ")

!nb: these are vectors
na=correlator%ntod/lc
correlator%ndet = nd
whiteNoise = sigma**2


npair = numberOfPairs(correlator%ndet)
!The master node performs the summation over the detectors, so it needs the Q^{-1}_{ij} fourier-space
! inverse correlators


call prepare_Qinv_array(correlator, correlationData, na, no_correlation) !DWPS added this subroutine


if (correlator%proc .ne. 0) return
!this is for the MASTER only!

allocate(correlator%beta(0:correlator%ndet-1))
do i=0,correlator%ndet-1
	correlator%beta(i) = real(lc,kind=8) / (regularisationParameter * whiteNoise(i) )
enddo


allocate(correlator%masterSpace(0:correlator%ndet-1))


call ds_log("Building Covariance",ds_feedback_quiet)
call buildCorrelators2(correlator,na,lc,whiteNoise,fknee,alpha,correlationData,nyquist,no_correlation)
call ds_log("Building Preconditioner",ds_feedback_quiet)

	
!build preconditioner in correlator%Pinv
do p=0,npair-1
   if(correlator%pair(p)%isCorrelated) then
      correlator%pair(p)%Pinv= correlator%pair(p)%Qinv
   endif
enddo

do i= 0, correlator%ndet-1
	p= triangleIndex(i,i,correlator%ndet)
	if(correlator%pair(p)%isCorrelated) correlator%pair(p)%Pinv= correlator%pair(p)%Qinv + correlator%beta(i)
enddo

call invert_correlator_long_matrix(correlator,.true.)
!call invert_fourier_long_matrix(correlator%Pinv,correlator%fftlength,numberOfPairs(correlator%ndet))

end subroutine prepareFFTQ2




subroutine read_noisefile(noisefile,nd,nyquist,variance,fknee,rho)
use inifile
!DWPS: modified to read white noise variance instead of std
character(*) :: noisefile
integer nd
real(dp) nyquist, corr_value
real(dp), dimension(0:nd-1) :: variance, fknee
real(dp), dimension(0:nd-1,0:nd-1) :: rho
integer f,i,j
logical :: param_status, set_corr, bad_corr
type(TIniFile) :: param_file

f = ds_get_lun()
call Ini_Open_File(param_file, noisefile, f, param_status, .false.)
call ds_assert(.not. param_status, 'Could not open parameter file '// noisefile)
	
nyquist = ini_read_double_file(param_file, 'nyquist')
variance = ini_read_double_file(param_file, 'variance')
fknee = ini_read_double_file(param_file, 'fknee')
	
! Now set the individual detector values.  The default is their current value (ie the value set for)
! all of them.
do i=0,nd-1
	variance(i) = ini_read_double_array_file(param_file, 'variance', i, variance(i) )
	fknee(i) = ini_read_double_array_file(param_file, 'fknee', i, fknee(i) )
enddo


!Set the detector correlations.
!First, set the default value for all the pairs
!If not found, this is zero.
rho = ini_read_double_file(param_file,'correlation',0.0_8)

!Then, set the remaining pairs
do i=0,nd-1
	do j=i,nd-1
		corr_value = ini_read_double_array_2d_file(param_file,'correlation',i,j,222.0D0)
		set_corr = (.not. (corr_value > 221.99) .and. (corr_value < 222.01) )
		bad_corr = (abs(corr_value) > 1.0)
		if (set_corr .and. (.not. bad_corr)) then
			rho(i,j) = corr_value
			rho(j,i) = corr_value
		endif
		if (set_corr .and. bad_corr) write(*,*) "Warning: Ignored funny correlation value ", corr_value, " for ", i,",",j
		if (set_corr .and. (i==j) .and. (rho(i,i) .ne. 0.0)) write (*,*) "Warning: 1/f autocorrelation ignored for parameter ", i
	enddo
enddo

! Finally, set the auto-correlations to 1.0
do i=0,nd-1
	rho(i,i) = 1.0
enddo


call Ini_Close_File(param_file)



end subroutine read_noisefile


subroutine buildOneCorrelator(nt,lc,whiteNoiseVariance,rho,fknee,alpha,nyquist,offsetPower,includeWhite)
real(dp) :: whiteNoiseVariance,fknee,alpha
integer lc
integer nt
#warning check start of qinv zero or one
integer na
logical includeWhite
real(dp), allocatable, dimension(:) :: timeCorrelator, offsetCorrelator
complex(dpc), dimension(:) :: offsetPower
complex(dpc), dimension(:), allocatable :: timePower
integer k,a,i,t1,t2, j
real(dp) :: dk, sigma_min1, sigma, rho, nyquist

	na = nt/lc

	call ds_assert(size(offsetPower)==na/2+1,"Wrong sized offsetPower")

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
            do t2= 0,lc-1
               if(abs(t1 - t2).le.nt) then
                  offsetCorrelator(a)= offsetCorrelator(a) + timeCorrelator(abs(t1 - t2))
               endif
            enddo
            sigma_min1= offsetCorrelator(a)				!first row's sum
         else										!following rows sum is sigma_min1 + new end element - old start element
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
	if (includeWhite) offsetPower = offsetPower + whiteNoiseVariance/lc



end subroutine buildOneCorrelator

subroutine buildCorrelators2(correlator,na,lc,sigma2,fknee,alpha,correlationData,nyquist,no_correlation)
	type(ds_correlator) :: correlator
	real(dp), dimension(0:correlator%ndet-1) :: sigma2, fknee,alpha
!	real(dp), dimension(0:correlator%ndet-1,0:correlator%ndet-1) :: rho

	type(ds_correlation_data) :: correlationData
	integer,dimension(0:correlator%ndet-1),intent(in) :: na
	integer,intent(in) :: lc
	logical, optional :: no_correlation
	real(dp) :: nyquist, s2, fk,a
	integer :: nt


	integer :: i_chunk, i, i_tod, j_tod, bsize, ierror , npair_proc, p, k,j, itstr, jtstr
	real(dp) :: sigma, sigma_min1, denom, dk
	character(128) :: msg

	type(ds_timestream_list),allocatable,dimension(:) :: timestream_list

	integer status(MPI_STATUS_SIZE)
	complex(dpc) :: q00, determ
	integer :: t, i_length!number of timestreams with this observation length

	!DWPS: We only need to build the correlator once using ffts for a given chunk length. This means we only need do it once for all pairs in a given observation length
	!get list of detectors with a given observation length
	
	!Doesn't that assume the same alpha for all of the timestreams?
	call get_pairs(correlator,timestream_list)

	call ds_assert(correlator%proc==0, "Only the master should call buildCorrelators")

	
	call ds_log("Warning: Assuming noise 1/f**alpha has alpha = 1.0",ds_feedback_silent)

	do i_length= 1, correlator%ndistinct
		write(msg,'("Building correlators for scan ",I4," of ",I4)') i_length,correlator%ndistinct
		call ds_log(trim(msg),ds_feedback_noisy)


	   !Add this lengths's correlators to Qinv
	   do i= 1,timestream_list(i_length)%ntimestreams
	      do j= i,timestream_list(i_length)%ntimestreams
	         itstr= timestream_list(i_length)%list(i)
	         jtstr= timestream_list(i_length)%list(j)
	         p = triangleIndex(itstr,jtstr,correlator%ndet)
	         if( correlator%pair(p)%isCorrelated) then
	            !the pair indexed by p should have the same ntod as this length instance
	            call ds_assert( correlator%ntod(itstr) == timestream_list(i_length)%ntod ,'Error: correlator lengths do not agree in buidCorrelators')
	            call ds_assert( correlator%ntod(jtstr) == timestream_list(i_length)%ntod ,'Error: correlator lengths do not agree in buidCorrelators')
				s2 = sqrt(sigma2(itstr)*sigma2(jtstr))
				fk = sqrt(fknee(itstr)*fknee(jtstr))
				a = alpha(itstr)
				nt = correlator%ntod(itstr)
!				call buildOneCorrelator(nt, lc, s2,rho(itstr,jtstr),fk,a,nyquist,correlator%pair(p)%Qinv, .true.)
#warning buildOneCorrelator call removed
	         endif
	      enddo
	   enddo

!	subroutine buildOneCorrelator(nt,lc,whiteNoiseVariance,rho,fknee,alpha,nyquist,offsetPower,includeWhite)
!correlator%pair(p)%Qinv= fd_Ca * sqrt(white(itstr) * white(jtstr)) * rho(itstr,jtstr) * sqrt(fknee(itstr)*fknee(jtstr))
!if(itstr==jtstr) correlator%pair(p)%Qinv=  correlator%pair(p)%Qinv + white(itstr) / lc

	enddo

	print*,'calling inversion'

	call invert_correlator_long_matrix(correlator,.false.)

	print*,'inversion complete'



end subroutine buildCorrelators2

subroutine buildCorrelators(correlator,na,lc,white,fknee,rho,nyquist,no_correlation)
!Transforms real-space circulant time correlator to circulant offset correlator, then fourier transforms and inverts.
!The Fourier transforms are distributed with nipairs for each process (nipairs can be different for each process where 
!nipairs <= npairs/nproc).  The final result correlator%Qinv is held only by the master process. The transformation 
!from time to offset uses C_a= (F^T F)^{-2} F^T C_n F, where the first term (F^T F)^{-2} = (lc)^{-2} I, is proportional 
!to the identity matrix.  The second term sums all terms in an lc*lc block into a single term. The algorithm uses the 
!circulance property of C_n to avoid explicitly summing more than once per block - this will have to be modified when 
!we switch to real data with non-circulant toeplitz C_n.
!
!Notes:
!target is correlator%Qinv(:,p)
!lc is chunk length - offset%length
!timeCorrelator is real and zero-based
!timeCorrelator's vals are destroyed within the subroutine
!nipairs  is the number of pairs the calling process is responsible for

type(ds_correlator) :: correlator
real(dp), dimension(0:correlator%ndet-1) :: white, fknee
real(dp), dimension(0:correlator%ndet-1,0:correlator%ndet-1) :: rho
integer,dimension(0:correlator%ndet-1),intent(in) :: na
integer,intent(in) :: lc
logical, optional :: no_correlation
real(dp) :: nyquist


integer :: i_chunk, i, i_tod, j_tod, bsize, ierror , npair_proc, p, a, k,j, itstr, jtstr
real(dp) :: sigma, sigma_min1, denom, dk

type(ds_timestream_list),allocatable,dimension(:) :: timestream_list

integer status(MPI_STATUS_SIZE)
complex(dpc) :: q00, determ
real(dp), allocatable,dimension(:) :: Ca !(na)
complex(dpc),allocatable,dimension(:) :: fd_Ca !(na/2+1)
real(dp), allocatable, dimension(:) :: timeCorrelator
complex(dpc), allocatable, dimension(:) :: power
integer :: t, i_length!number of timestreams with this observation length

!DWPS: We only need to build the correlator once using ffts for a given chunk length. This means we only need do it once for all pairs in a given observation length
!get list of detectors with a given observation length
call get_pairs(correlator,timestream_list)

call ds_assert(correlator%proc==0, "Only the master should call buildCorrelators")


call ds_log("Warning: Assuming noise 1/f**alpha has alpha = 1.0",ds_feedback_silent)

do i_length= 1, correlator%ndistinct

   allocate(timeCorrelator(0:timestream_list(i_length)%ntod-1))
   allocate(power(0:timestream_list(i_length)%ntod/2))

   !Build the power spectrum
   power(0)=cmplx(0.0_dp,0.0_dp,kind=dpc)
   dk = 2.0_dp * nyquist / timestream_list(i_length)%ntod
   do k=1,timestream_list(i_length)%ntod/2
      power(k) = dcmplx( (1.0_dp / k / dk)**1.0, 0.0_dp)
   enddo

   !Convert power to time correlator
   call ds_ifft(power,timeCorrelator)
   deallocate(power)

   !convert time correlator to offset correlator, for off-diagonal components
   allocate(Ca(1:timestream_list(i_length)%ntod/lc)) !na
   allocate(fd_Ca(1:(timestream_list(i_length)%ntod/lc)/2+1))  !na/2+1
   Ca=0
   do i_chunk= 1,timestream_list(i_length)%ntod/lc !ie: na
      do i= 0,lc-1
         i_tod= (i_chunk-1)*lc+i
         if(i == 0) then							!only need to do sum once per i_chunk
            do j_tod= 0,lc-1
               if(abs(i_tod - j_tod).le.timestream_list(i_length)%ntod) then
                  Ca(i_chunk)= Ca(i_chunk) + timeCorrelator(abs(i_tod - j_tod))
               endif
            enddo
            sigma_min1= Ca(i_chunk)				!first row's sum
         else										!following rows sum is sigma_min1 + new end element - old start element
            sigma= sigma_min1 - timeCorrelator(abs((-(i_chunk-1)+1)*lc-i)) + timeCorrelator(abs(-((i_chunk-1)*lc)-(i+1)+1)) !ichunk is 1 based,i and timeCorrelator are zero based			 
            sigma_min1= sigma
            Ca(i_chunk)= Ca(i_chunk) + sigma
         endif
      enddo
   enddo
   
   !Do (F^t F)^-2
   denom= real(lc,kind=8)**2
   Ca= Ca/denom
   call ds_fft(Ca,fd_Ca)
   deallocate(Ca)

   !Add this lengths's correlators to Qinv
   do i= 1,timestream_list(i_length)%ntimestreams
      do j= i,timestream_list(i_length)%ntimestreams
         itstr= timestream_list(i_length)%list(i)
         jtstr= timestream_list(i_length)%list(j)
         p = triangleIndex(itstr,jtstr,correlator%ndet)
         if( correlator%pair(p)%isCorrelated) then
            !the pair indexed by p should have the same ntod as this length instance
            call ds_assert( correlator%ntod(itstr) == timestream_list(i_length)%ntod ,'Error: correlator lengths do not agree in buidCorrelators')
            call ds_assert( correlator%ntod(jtstr) == timestream_list(i_length)%ntod ,'Error: correlator lengths do not agree in buidCorrelators')
            correlator%pair(p)%Qinv= fd_Ca * sqrt(white(itstr) * white(jtstr)) * rho(itstr,jtstr) * sqrt(fknee(itstr)*fknee(jtstr))
            if(itstr==jtstr) correlator%pair(p)%Qinv=  correlator%pair(p)%Qinv + white(itstr) / lc
         endif
      enddo
   enddo

!   correlator%Qinv = 0.0
!   do i=0,correlator%ndet-1
!      do j=i,correlator%ndet-1
!         p = triangleIndex(i,j,correlator%ndet)
!         if (rho(i,j) .ne. 0 .and. (i==j .or. .not. no_correlation) ) correlator%Qinv(:,p)= fd_Ca * sqrt(white(i) * white(j)) * rho(i,j) * sqrt(fknee(i)*fknee(j))
!         if (i==j) correlator%Qinv(:,p) = correlator%Qinv(:,p) + white(i) / lc
!      enddo
!   enddo

 deallocate(fd_Ca)
 deallocate(timeCorrelator)
enddo

print*,'calling inversion'

call invert_correlator_long_matrix(correlator,.false.)

print*,'inversion complete'
!call invert_fourier_long_matrix(correlator%Qinv,correlator%fftlength,numberOfPairs(correlator%ndet))

!!$!dump to file
!!$#warning !dumping Q_inv to file in prepareFFTQ2
!!$print*,'dumping Q_inv to file'
!!$open(unit=1100,file="Q_inv_new.bin",form="binary",status="replace")
!!$write(1100) correlator%Qinv
!!$close(1100)
!!$stop

end subroutine buildCorrelators

subroutine writeNoiseRmsFile(filename,pointings,originalIndices,whitenoise,nside)
!whitenoise is variance. Assume it is the same for all streams
   character(len=100),intent(in) :: filename
   type(ds_detpointing),intent(in) :: pointings(0:)
   integer,intent(in) :: originalIndices(pointings(0)%npix), nside
   real(dp),intent(in) :: whitenoise(0:)
   integer :: p

   print*,"writing hit map",filename
   open(unit=1016,file=trim(adjustl(filename)),status="replace")
   write(1016,*) "nside=", nside
   write(1016,*) "npix=", pointings(0)%npix
   do p= 1,pointings(0)%npix
      write(1016,*) originalIndices(p), sqrt( whitenoise(0) / real(pointings(0)%globalHitMap%map(p),kind=8) )
   enddo
   close(1016)

end subroutine writeNoiseRmsFile


!Creates array of the Qinv stucture, allocating space only when there is non-zero correlation
!For instances of Qinv with memory allocated, isCorrelated= .true.
subroutine prepare_Qinv_array(C,correlation_data,na,no_correlation)
	type(ds_correlation_data) :: correlation_data
  type(ds_correlator),intent(inout) :: C
!  real(dp),dimension(0:C%ndet-1,0:C%ndet-1) :: rho
  integer,dimension(0:C%ndet-1) :: na
  logical,intent(in) :: no_correlation
  integer :: p, i, j, npairs
  npairs= (C%ndet**2+C%ndet)/2
  allocate(C%pair(0:npairs-1))
  do p= 0,npairs-1
     C%pair(p)%isCorrelated= .false.     !default is that there is no correlation
     call squareIndex(C%ndet,p,i,j)
     if(i==j) then                       !diagonals by definition are correlated
        C%pair(p)%isCorrelated=.true.    !pair is correlated
        C%pair(p)%fftlength= na(i)/2 +1  !NB: must make na and ntod arrays over detectors
        if(C%proc==0) then
           allocate(C%pair(p)%Qinv(C%pair(p)%fftlength))
           allocate(C%pair(p)%Pinv(C%pair(p)%fftlength))
        endif
!     elseif((.not. no_correlation) .and. (rho(i,j) .ne. 0.0_8)) then
#warning line replaced until make correlation_data type
     elseif(.not. no_correlation) then
        !off diagonals with non-zero correlation
        C%pair(p)%isCorrelated= .true.
        C%pair(p)%fftlength= na(i)/2 +1  
        if(C%proc==0) then
           allocate(C%pair(p)%Qinv(C%pair(p)%fftlength))
           allocate(C%pair(p)%Pinv(C%pair(p)%fftlength))
        endif
     endif
  enddo
end subroutine prepare_Qinv_array


!DWPS: rename this
subroutine get_pairs(correlator,timestream_list)
!get list of all timestreams, in groups of ntod
  type(ds_correlator),intent(inout) :: correlator
  !integer,dimension(0:correlator%ndet-1),intent(in) :: ntod
  type(ds_timestream_list),allocatable,dimension(:),intent(out) :: timestream_list

  integer,allocatable,dimension(:,:) :: buf
  integer :: p, i, j, nt,d, ndistinct

  !count number of distinct ntods
  allocate(buf(2,0:correlator%ndet-1))
  buf= -1
  correlator%ndistinct=0
  do d= 0, correlator%ndet-1
     nt= correlator%ntod(d)
     !cycle through list of ntods
     i=0
     count: do
        if(buf(2,i).eq.-1) then
           !this is the first instance of a new ntod
           correlator%ndistinct= correlator%ndistinct+1
           buf(1,i)= nt
           buf(2,i)= 1
           exit count
        elseif(nt==buf(1,i)) then
           !this is another instance of an existing ntod
           buf(2,i)= buf(2,i)+1
           exit count
        endif
        i=i+1
     enddo count
  enddo
  
  !put this info into a nice ds_timestream_list structure
  allocate(timestream_list(1:correlator%ndistinct))
  do i=1,correlator%ndistinct
     timestream_list(i)%ntod = buf(1,i-1)
     timestream_list(i)%ntimestreams= buf(2,i-1)
     allocate(timestream_list(i)%list(1:timestream_list(i)%ntimestreams))
     timestream_list(i)%list= -1
  enddo

  deallocate(buf)

  !now put detector numbers into the list
  do d=0,correlator%ndet-1
     nt= correlator%ntod(d)
     list: do i= 1,correlator%ndistinct
        if(nt .eq. timestream_list(i)%ntod) then
           element: do j= 1,timestream_list(i)%ntimestreams
              if(timestream_list(i)%list(j) .eq. -1) then
                 timestream_list(i)%list(j)= d
                 exit element
              endif
           enddo element
           exit list
        endif
     enddo list
  enddo

end subroutine get_pairs


subroutine invert_correlator_long_matrix(C,prec) !matrix,fftlength,npair)

!DWPS: This routine replaces invert_fourier_long_matrix. The Cholesky decomp
!is done mode-by-mode only for each set of correlated pairs, which by 
!definition have the same length. We don't consider correlated pairs with
!different lengths. This also splits the potentially very large (ntimestream)^2
!Cholesky decomp into a more managable set by scan (or by "length").
!prec is a flag replacing Qinv with Pinv to invert the preconditioner instead

type(ds_correlator) :: C
!integer,dimension(0:C%ndet-1) :: ntod
logical :: prec
integer::iscan

integer :: i_scan, a, npairs, ndet, ierror
integer :: i,j,p,p_thischunk,itstr,jtstr
complex(dpc) :: determ, q00
complex(dpc),allocatable,dimension(:) :: inversionSpace
type(ds_timestream_list),allocatable,dimension(:) :: timestream_list

!get list of timestreams in a scan ("of a length")
call get_pairs(C,timestream_list)

do i_scan= 1, C%ndistinct
  
   !get the number of timestreams in this inversion
   ndet= timestream_list(i_scan)%ntimestreams !need special case if ==1
   npairs= (ndet**2+ndet)/2

   !The master node now has all the Q and must invert them to Qinv
   !This is done in Fourier space, inverting each mode separately
   !calls to hardware optimised lapack routines:
   !zpptrf - computes double complex cholesky factorization of matrix Q
   !zpptri - computes double complex inverse of Q from output of zpptrf
   !these routines run on packed upper triangle of Q

   allocate(inversionSpace(0:npairs-1))
   ierror=0

   do a=2, C%pair(timestream_list(i_scan)%list(1))%fftlength
      !build the matrix for inversion
      do i= 0,ndet-1
         do j= i,ndet-1
            p_thischunk= triangleIndex(i,j,ndet) !triangle index for THIS inversion

            itstr= timestream_list(i_scan)%list(i+1) !list is 1, not 0, based.
            jtstr= timestream_list(i_scan)%list(j+1)
            p= triangleIndex(itstr,jtstr,C%ndet) !triangle index for global data
            if(C%pair(p)%isCorrelated) then
               if(prec) then
                  inversionSpace(p_thischunk)= C%pair(p)%Pinv(a)
               else
                  inversionSpace(p_thischunk)= C%pair(p)%Qinv(a)
               endif
            else
               inversionSpace(p_thischunk)= 0.0_dpc
            endif
         enddo
      enddo
      
      !apply Cholesky LAPACK routines
      call zpptrf('L',ndet,inversionSpace,ierror)
      if(ierror.ne.0) print*,ierror,a,inversionSpace
      call ds_assert(ierror==0, "Cholesky factorization failure 1")
      call zpptri('L',ndet,inversionSpace,ierror)
      call ds_assert(ierror==0, "Cholesky factorization failure 2 ")

      !put inverted matrix back
      do i= 0,ndet-1
         do j= i,ndet-1
            p_thischunk= triangleIndex(i,j,ndet) !triangle index for THIS inversion

            itstr= timestream_list(i_scan)%list(i+1) !list is 1 based, not zero based
            jtstr= timestream_list(i_scan)%list(j+1)
            p= triangleIndex(itstr,jtstr,C%ndet) !triangle index for global data

            if(C%pair(p)%isCorrelated) then
               if(prec) then
                  C%pair(p)%Pinv(a)= inversionSpace(p_thischunk)
               else
                  C%pair(p)%Qinv(a)= inversionSpace(p_thischunk)
               endif
            endif
         enddo
      enddo
   enddo
   !We set the zeroth mode to zero by hand
   do i= 0,timestream_list(i_scan)%ntimestreams-1
      do j= i,timestream_list(i_scan)%ntimestreams-1
         p= triangleIndex( timestream_list(i_scan)%list(i+1),timestream_list(i_scan)%list(j+1),C%ndet) 
         if(C%pair(p)%isCorrelated) then
            if(prec) then
               C%pair(p)%Pinv(1)= 0.0_8
            else   
               C%pair(p)%Qinv(1)= 0.0_8
            endif
         endif
      enddo
   enddo

   deallocate(inversionSpace)
enddo

end subroutine invert_correlator_long_matrix



subroutine apply_FtZF(correlator,pointing,whiteNoise,npix,offset_in,offset_out)
!This takes and returns an offset vector, applying 'F^t ZF' to the trial offset solution 'a'
!The routine does not require a dummyTimestream input and does not create a full dummy timestream 
!internally - only one 'detectors' timestream' is allocated per processsor at any given time.
type(ds_correlator),intent(in) :: correlator
type(ds_detpointing),intent(in),dimension(0:correlator%my_ndet/2-1) :: pointing
real(dp),dimension(0:correlator%ndet-1) :: whiteNoise
integer,intent(in) :: npix
type(ds_offsets),dimension(0:correlator%my_ndet-1) :: offset_in, offset_out

!these dummy maps and timestreams held internally
type(ds_map) :: qbuf, ubuf
type(ds_timestream) :: dummyTime1, dummyTime2 !NOTE: only 1 instance of each of these!
integer :: pair, my_pairs, dq, du, itod, p, ierror

my_pairs = correlator%my_ndet / 2
call ds_assert(mod(correlator%my_ndet,2)==0, 'Need even number of detectors - q and u')

!initialise map sums
call prepareMap(qbuf,npix)
call prepareMap(ubuf,npix)

do pair= 0,my_pairs-1
    dq= det_qindex(pair)
    du= det_uindex(pair)
    !project this offset onto a timestream of appropriate size
    call prepareTimestream(dummyTime1,pointing(pair)%nt)
    call prepareTimestream(dummyTime2,pointing(pair)%nt)
!    call projectOffsetOntoTimestream(offset_in(dq),dummyTime1,'add')   !F operator    
!    call projectOffsetOntoTimestream(offset_in(du),dummyTime2,'add')   !F operator 
    call addOffsetToTimestream(offset_in(dq),dummyTime1)
    call addOffsetToTimestream(offset_in(du),dummyTime2)

    !sum into map space
    do itod = 1,pointing(pair)%nt
            p = pointing(pair)%pixel(itod)
            qbuf%map(p) = qbuf%map(p) + detectorToSkyQ(pointing(pair)%theta(itod), dummyTime1%timestream(itod), dummyTime2%timestream(itod) )
            ubuf%map(p) = ubuf%map(p) + detectorToSkyU(pointing(pair)%theta(itod), dummyTime1%timestream(itod), dummyTime2%timestream(itod) )
    enddo    
    
    !clear the timestream memory
    call destroyTimestream(dummyTime1)
    call destroyTimestream(dummyTime2)
enddo

!mpi_allreduce the map sums then average
if(correlator%ndet>2) then
	call mpi_allreduce(MPI_IN_PLACE, qbuf%map, qbuf%npix, MPI_DOUBLE_PRECISION, MPI_SUM, correlator%comm, ierror)
	call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
	call mpi_allreduce(MPI_IN_PLACE, ubuf%map, ubuf%npix, MPI_DOUBLE_PRECISION, MPI_SUM, correlator%comm, ierror)
	call ds_assert(ierror==0, 'Error in mpi_allreduce in makeNaiveMap')
endif

!average the Q/U map sums
qbuf%map = qbuf%map / pointing(0)%GlobalHitMap%map
ubuf%map = ubuf%map / pointing(0)%GlobalHitMap%map

do pair=0,my_pairs-1
    !project map onto this timestream of appropriate size
    call prepareTimestream(dummyTime1,pointing(pair)%nt)
    call prepareTimestream(dummyTime2,pointing(pair)%nt)
    
    !deproject map onto this timestream and then onto an offset vector
    dq= det_qindex(pair)
    du= det_uindex(pair)
    do itod= 1, pointing(pair)%nt
        p= pointing(pair)%pixel(itod)
        dummyTime1%timestream(itod) = skyToDetectorQ(pointing(pair)%theta(itod), qbuf%map(p), ubuf%map(p) )
		dummyTime2%timestream(itod) = skyToDetectorU(pointing(pair)%theta(itod), qbuf%map(p), ubuf%map(p) )
    enddo
    call deprojectTimestreamOntoOffset(dummyTime1,offset_out(dq)) !F^T operator
    call deprojectTimestreamOntoOffset(dummyTime2,offset_out(du)) !F^T operator
    
    !Now have F^t P naive() F a in offset_out. Subtract this from F^t F a because
    ! F^t Z F a = (F^t F - F^t P naive() F) a
    !where naive()= (P^tP)^-1 P^t  and F^t F = offset length
    offset_out(dq)%values= (offset_in(dq)%length * offset_in(dq)%values) - offset_out(dq)%values
    offset_out(du)%values= (offset_in(du)%length * offset_in(du)%values) - offset_out(du)%values    

	offset_out(dq)%values = offset_out(dq)%values / whiteNoise(correlator%my_det(dq)) 					!C_w^{-1} operator
	offset_out(du)%values = offset_out(du)%values / whiteNoise(correlator%my_det(du)) 					!C_w^{-1} operator	
	
    !clear the timestream memory
    call destroyTimestream(dummyTime1)
    call destroyTimestream(dummyTime2)
enddo

!clear the map buffer memory
call destroyMap(qbuf)
call destroyMap(ubuf)

end subroutine apply_FtZF


subroutine d_test_naive(timestreams,qmap,umap,pointings,correlator)
!subroutine to test new naive mapping engine.  The result should be the same as make_naive_map
type(ds_correlator) :: correlator
type(ds_timestream),dimension(0:correlator%my_ndet-1) :: timestreams
type(ds_detPointing),dimension(0:correlator%my_ndet/2-1) :: pointings
type(ds_map) :: qmap, umap
integer :: d,i,j, ierror, nmodules, imodule, idiode
type(ds_modulescan),allocatable,dimension(:) :: modulescan

!create the modulescan structure
nmodules= correlator%my_ndet/2 !we are pretending the modules have only 2 diodes
allocate(modulescan(0:nmodules-1))
do d= 0,nmodules-1
	modulescan(d)%ntod= pointings(d)%nt
    modulescan(d)%diode_flag=0
    modulescan(d)%diode_flag(1:2)= 1
    modulescan(d)%dpsi(1)= 0.0_8
    modulescan(d)%dpsi(2)= PI / 4.0_8
    allocate(modulescan(d)%theta(pointings(d)%nt))
    allocate(modulescan(d)%pointing(pointings(d)%nt))
    modulescan(d)%theta= pointings(d)%theta
    modulescan(d)%pointing= pointings(d)%pixel
    !set C_w^-1 = I
    modulescan(d)%inv_Cw=0.0_8
    modulescan(d)%inv_Cw(1,1)=1.0_8
    modulescan(d)%inv_Cw(2,2)=1.0_8
enddo

!make the cov matrix
call allocate_cov(qmap%npix)

do d= 0, nmodules-1
    call add2cov(modulescan(d))
enddo

call mpi_allreduce(MPI_IN_PLACE,naive_cov,3*p_npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)

call invert_weight()

do d= 0,correlator%my_ndet-1
    imodule= d/2
    idiode= modulo(d,2)+1
    !assumes that timestreams i and i+1 belong to module i/2
    call add2rhs(idiode,modulescan(imodule),timestreams(d),qmap,umap)
enddo

!reduce the maps
call mpi_allreduce(MPI_IN_PLACE,qmap%map,p_npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)
call mpi_allreduce(MPI_IN_PLACE,umap%map,p_npix,MPI_DOUBLE_PRECISION,MPI_SUM,correlator%comm,ierror)

!apply the single pixel covariance matrices
call cov_mult(qmap,umap)

!destroy the cov matrix
call deallocate_cov()

end subroutine d_test_naive

end module
 
