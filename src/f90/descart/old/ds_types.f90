module ds_types
  use ds_precision
  use ds_utils
  implicit none
  

!!                                                  !!
!!  CONSTANT DEFINITIONS (MOST IN HEALPIX TYPES)    !!
!!                                                  !!

real(kind=DP), parameter, public :: HOUR2RAD = TWOPI / 24.0_DP

  
  
!!!!  DWPS
!!!!  These are the azimuth offset structures.  One holds offset info and other is a pointer based
!!!!  memory cache so that we only store the Nt length azimuth pointings once
!!!!  The ds_offset type should contain an instance of ds_az_offset (and an azimuth_flag)
  
type ds_az_retriever
    !!holds pointing memory cache OR points to the cache
    integer,allocatable,dimension(:) :: pointing
end type ds_az_retriever

  
  !the azimuth offsets object (1 per CES scan per diode stream) - MOVED TO DS_TYPES
type ds_az_offset
    integer :: n_az,ntod
    type(ds_az_retriever),pointer :: cache_retrieve !size ntod vector mapping tod to offset
    real(dp),allocatable,dimension(:) :: amplitudes  !size n_az
    integer,allocatable,dimension(:) :: hits     !size n_az
end type ds_az_offset

  
type ds_map
  real(dp), allocatable, dimension(:) :: map
  integer, allocatable, dimension(:) :: indices
  integer :: npix
end type ds_map
  

type ds_detpointing
  integer, allocatable, dimension(:) :: pixel
  real(dp), allocatable, dimension(:) :: theta
  integer :: npix
  integer :: nt
  type(ds_map) :: globalHitMap
end type ds_detpointing
  

type ds_timestream
  integer :: nt, my_module
  real(dp), allocatable, dimension(:) :: timestream
end type ds_timestream
  
  
type ds_offsets
  integer :: na, my_module
  real(dp), allocatable, dimension(:) :: values
  integer :: length
  logical :: azimuth_flag
  type(ds_az_offset) :: az 
  complex(dpc), pointer, dimension(:) :: Qinv
  complex(dpc), pointer, dimension(:) :: Pinv
  logical ownsCorrelators
end type ds_offsets
  
  
type ds_converter
  integer, allocatable, dimension(:) :: pixelTable
  integer length
  integer currentIndex_
  logical final_
end type ds_converter


type ds_modulescan
    integer :: ntod
    real(dp),dimension(4) :: dpsi !radians! !offset of diode from theta
    real(dp) :: inv_Cw(1:4,1:4)
    integer,dimension(4) :: diode_flag        !1 to use, 0 to reject whole diode
    integer,allocatable,dimension(:) :: pointing
    real(dp),allocatable,dimension(:) :: theta !radians!
end type ds_modulescan

interface prepareOffsets
    module procedure prepareOffsets_stan, prepareOffsets_withaz
end interface
 
  
  
contains
  

	!JZ Notes
	!Using the offset-owned correlators.
	!In your code to build/load the correlators, allocate
	!complex(dpc) vectors with the "target" attribute
	!for Pinv and Qinv.
	!Call the setOffsetCorrelators(offset,Pinv,Qinv) routine
	!Do not deallocate your Pinv and Qinv.  Just let them be.
	!They will be deallocated when the destroyOffsets is called.
	

  
  !!                                            !!
  !!    SUBROUTINES TO DO WITH THE DATA TYPES   !!
  !!                                            !!
  
  
  subroutine prepareOffsets_stan(offsets, na,offsetLength)
    type(ds_offsets) :: offsets
    integer, intent(in) :: na, offsetLength
    offsets%length = offsetLength
    offsets%na = na
	call ds_checkAllocate(offsets%values, offsets%na)
	
    offsets%azimuth_flag= .false. !default is false
	nullify(offsets%Qinv)
	nullify(offsets%Pinv)
	offsets%ownsCorrelators=.false.
	
  end subroutine prepareOffsets_stan

  subroutine prepareOffsets_withaz(offsets, source)
    type(ds_offsets) :: offsets, source
    offsets%length = source%length
    offsets%na = source%na
	call ds_checkAllocate(offsets%values, offsets%na)
	
    offsets%azimuth_flag=source%azimuth_flag
    if (offsets%azimuth_flag) call az_copy(source%az,offsets%az)
    if (associated(source%Qinv)) offsets%Qinv=>source%Qinv
	if (associated(source%Pinv)) offsets%Pinv=>source%Pinv
	offsets%ownsCorrelators=.false.
	
  end subroutine prepareOffsets_withaz

  
  subroutine copyOffsets(source, dest)
      type(ds_offsets) :: source, dest
	  integer i
      call prepareOffsets(dest, source%na, source%length)
	  do i= 1, source%na
			dest%values(i) = source%values(i)
	  enddo
	  
	  dest%azimuth_flag= source%azimuth_flag
      if(dest%azimuth_flag) call az_copy(source%az, dest%az)
    if (associated(source%Qinv)) dest%Qinv=>source%Qinv
	if (associated(source%Pinv)) dest%Pinv=>source%Pinv
	dest%ownsCorrelators=.false.

  end subroutine copyOffsets

  subroutine setOffsetCorrelators(offset,pinv,qinv)
	complex(dpc), target, dimension(:) :: pinv,qinv
	type(ds_offsets) :: offset
	integer n
	
	n=size(pinv)
	call ds_assert(2*n + 1 == offset%na, "Wrong size correlators in setOffsetCorrelators")
	call ds_assert(size(qinv) == n, "Different size correlators in setOffsetCorrelators")
	
	offset%pinv => pinv
	offset%qinv => qinv
	offset%ownsCorrelators = .true.
	
  end subroutine setOffsetCorrelators
  
  
  subroutine destroyOffsets(offsets)
    type(ds_offsets),intent(inout) :: offsets
    if(allocated(offsets%values)) deallocate(offsets%values)
    offsets%na= -1
    offsets%length= -1
	if (offsets%ownsCorrelators) then
		deallocate(offsets%qinv)
		deallocate(offsets%pinv)
	endif
    if (associated(offsets%Qinv)) nullify(offsets%Qinv)
    if (associated(offsets%Pinv)) nullify(offsets%Pinv)

   if(offsets%azimuth_flag) call az_destroy(offsets%az)
  end subroutine destroyOffsets
  
  subroutine prepareTimestream(timestream,nt)
    type(ds_timestream) :: timestream
    integer :: nt
    
    call ds_checkAllocate(timestream%timestream, nt)
    timestream%nt = nt
  end subroutine prepareTimestream
  
  subroutine copyTimestream(source, dest)
    !Copy a timestream into another.
    !If the destination is not allocated with the right size, (re-)allocate it.
    type(ds_timestream) :: source, dest
	integer i
    call prepareTimestream(dest,source%nt)
	do i= 1, source%nt
		dest%timestream(i) = source%timestream(i)
	enddo
  end subroutine copyTimestream
  
  subroutine destroyTimestream(timestream)
    type(ds_timestream) :: timestream
    if (allocated(timestream%timestream)) deallocate(timestream%timestream)
    timestream%nt = -1
  end subroutine destroyTimestream
  
  
  subroutine prepareMap(map, npix)
    type(ds_map) :: map
    integer npix
    
    call ds_checkAllocate(map%map,npix)
    call ds_checkAllocate(map%indices,npix)
        
    map%npix = npix
  end subroutine prepareMap
  
  subroutine destroyMap(map)
    type(ds_map) :: map
    if (allocated(map%map)) deallocate(map%map)
    if (allocated(map%indices)) deallocate(map%indices)
    map%npix=-1
  end subroutine destroyMap
  
  
  subroutine copyMap(source, dest)
    type(ds_map) :: source, dest
	integer i
    call prepareMap(dest, source%npix)
	do i= 1, source%npix
		dest%indices(i) = source%indices(i)
		dest%map(i) = source%map(i)
	enddo
    
  end subroutine copyMap
  
  subroutine preparePointing(pointing, nt, np,do_hitmap)
    type(ds_detpointing) :: pointing
	integer nt
	integer, optional :: np
	logical, optional :: do_hitmap
	integer npix
	
	if (present(np)) then
		npix = np
	else
		npix = 0
	endif
	
	call ds_checkAllocate(pointing%pixel, nt)
	call ds_checkAllocate(pointing%theta, nt)
	pointing%nt = nt
	pointing%npix = npix
    if (present(do_hitmap)) then 
    	if (do_hitmap) then
		    call prepareMap(pointing%globalHitMap, npix)
		    pointing%globalHitMap%map = 0
	    endif
	endif 
	
  end subroutine preparePointing
  
  subroutine destroyPointing(pointing)
      type(ds_detpointing) :: pointing
    if (allocated(pointing%pixel)) deallocate(pointing%pixel)
    if (allocated(pointing%theta)) deallocate(pointing%theta)
    pointing%npix=-1
    pointing%nt = -1
    call destroyMap(pointing%globalHitMap)
end subroutine destroyPointing


subroutine copyPointing(source, dest, do_hitmap)
      type(ds_detpointing) :: source,dest
      logical,optional :: do_hitmap
	  integer i

	call ds_checkAllocate(dest%pixel,source%nt)
	call ds_checkAllocate(dest%theta,source%nt)
    
	do i= 1,source%nt
		dest%pixel(i) = source%pixel(i)
		dest%theta(i) = source%theta(i)
	enddo
    dest%npix = source%npix
    dest%nt = source%nt
    
!    if (present(do_hitmap)) then 
!    	if (do_hitmap) then
!		    call prepareMap(dest%globalHitMap, dest%npix)
!			do i= 1,dest%npix
!				dest%globalHitMap = source%globalHitMap
!			enddo
!	    endif
!	endif 
	
end subroutine copyPointing
  

!                                                                   !
!   utilities related to scan-synchronous offset type ds_az_offset  !
!                                                                   !

!subroutine to initialise a ds_az_offset instance
subroutine az_initialise(ntod,n_az,az)
integer,intent(in) :: ntod, n_az
type(ds_az_offset),intent(inout) :: az

az%n_az= n_az
az%ntod= ntod
!if(allocated(az%cache_retrieve)) deallocate(az%cache_retrieve)
if(allocated(az%amplitudes)) deallocate(az%amplitudes)
if(allocated(az%hits)) deallocate(az%hits)
allocate(az%amplitudes(n_az))
allocate(az%hits(n_az))

az%amplitudes= 0.0_8
az%hits= 0
nullify(az%cache_retrieve)

end subroutine az_initialise




!subroutine to destroy a ds_az_offset instance
subroutine az_destroy(az)
type(ds_az_offset),intent(inout) :: az

az%n_az= -1
az%ntod= -1
if(allocated(az%amplitudes)) deallocate(az%amplitudes)
if(associated(az%cache_retrieve)) nullify(az%cache_retrieve)
if(allocated(az%hits)) deallocate(az%hits)

end subroutine az_destroy



!subroutine to copy an azimuth pointing (except the amplitudes field which is allocated and blank)
subroutine az_copy(az_in,az_out)
!Initialise az_out using array sizes from az_in AND copies the cache_retrieve and hits from az_in
type(ds_az_offset),intent(in) :: az_in
type(ds_az_offset),intent(inout) :: az_out

call az_initialise(az_in%ntod,az_in%n_az,az_out)
if(associated(az_in%cache_retrieve)) then
    az_out%cache_retrieve => az_in%cache_retrieve
else
    nullify(az_out%cache_retrieve)
endif
az_out%hits = az_in%hits

end subroutine az_copy
  
end module
