module ds_types
  use ds_precision
  use ds_utils
  implicit none
  

!!                                                  !!
!!  CONSTANT DEFINITIONS (MOST IN HEALPIX TYPES)    !!
!!                                                  !!

real(kind=DP), parameter, public :: HOUR2RAD = TWOPI / 24.0_DP
!integer, parameter, public :: ndiodes_T = 2
!integer, parameter, public :: ndiodes_P = 4
!integer, parameter, public :: ndiodes_max = ndiodes_T + ndiodes_P
integer, parameter :: MAGIC_VALUE = 250582

  
!!!!  DWPS
!!!!  These are the azimuth offset structures.  One holds offset info and other is a pointer based
!!!!  memory cache so that we only store the Nt length azimuth pointings once
!!!!  The ds_offset type should contain an instance of ds_az_offset (and an azimuth_flag)
  

type ds_noiseinfo
  real(dp),pointer,dimension(:) :: sigma, fknee, alpha
  real(dp),pointer,dimension(:,:) :: corr
end type ds_noiseinfo


type ds_az_retriever
    !!holds pointing memory cache OR points to the cache
    integer,allocatable,dimension(:) :: pointing
end type ds_az_retriever

  
  !the azimuth offsets object (1 per CES scan per diode stream) - MOVED TO DS_TYPES
type ds_az_offset
    integer :: n_az,ntod
!    type(ds_az_retriever),pointer :: cache_retrieve !size ntod vector mapping tod to offset
	integer,dimension(:),pointer :: pointing
    real(dp),allocatable,dimension(:) :: amplitudes  !size n_az
    integer,allocatable,dimension(:) :: hits     !size n_az
end type ds_az_offset

  
type ds_map
  real(dp), allocatable, dimension(:) :: map
  integer, allocatable, dimension(:) :: indices
  integer :: npix
end type ds_map
  

type ds_trimap
	type(ds_map) :: T,Q,U
	logical :: has_p, has_t
end type ds_trimap



type ds_detpointing
  integer, allocatable, dimension(:) :: pixel
  real(dp), allocatable, dimension(:) :: theta
  integer :: npix
  integer :: nt
  type(ds_map) :: globalHitMap
end type ds_detpointing
  


type ds_covariance
	real(dp), allocatable, dimension(:) :: TT
	real(dp), allocatable, dimension(:) :: QQ
	real(dp), allocatable, dimension(:) :: UU
	real(dp), allocatable, dimension(:) :: TQ
	real(dp), allocatable, dimension(:) :: TU
	real(dp), allocatable, dimension(:) :: QU
	logical :: has_p, has_t
	integer npix
end type ds_covariance


type ds_modulescan
    integer :: ntod
	integer :: ndiodes_T, ndiodes_P, ndiodes
	integer magic_check
    real(dp),dimension(:), pointer :: dpsi !radians! !offset of diode from theta
    real(dp), dimension(:,:), pointer  :: inv_Cw !? Do we need this?
    real(dp), pointer :: inv_Cw_TT(:,:)
    real(dp), pointer :: inv_Cw_PP(:,:)
    real(dp), pointer :: inv_Cw_TP(:,:)
    integer,dimension(:), pointer :: flags_T        !1 to use, 0 to reject whole diode
    integer,dimension(:), pointer :: flags        !1 to use, 0 to reject whole diode
    integer,dimension(:), pointer :: flags_P        !1 to use, 0 to reject whole diode
    integer,pointer,dimension(:) :: pointing
	integer, pointer, dimension(:) :: az_pointing, az_hits
	logical :: owns_az
    real(dp),pointer,dimension(:) :: theta !radians!
	logical :: has_p, has_t  !has
    type(ds_timestream),pointer,dimension(:) :: timestreams !will have dimension ndiodes_max
    type(ds_timestream),pointer,dimension(:) :: timestreams_P !will have dimension ndiodes_P  !pointer to subarray of timestreams
    type(ds_timestream),pointer,dimension(:) :: timestreams_T !will have dimension ndiodes_T  !pointer to subarray of timestreams
    type(ds_offsets),pointer,dimension(:) :: offsets_2  !Will rename later.
    type(ds_offsets),pointer,dimension(:) :: offsets_P
    type(ds_offsets),pointer,dimension(:) :: offsets_T
	type(ds_noiseinfo) :: noise
	integer, dimension(3) :: id !JAZ - run, scan, seg
	real(dp) :: scanfreq
	
	!Map and covariance for this scan only - we are going to perform multiple null tests so need to keep these.
	type(ds_trimap) :: maps
	type(ds_covariance) :: covariance
	
end type ds_modulescan


type ds_timestream
  integer :: nt, my_module
  real(dp), allocatable, dimension(:) :: timestream
end type ds_timestream
  
  
type ds_offsets
  integer :: na, my_module
  real(dp), allocatable, dimension(:) :: values
  integer :: length
  logical :: azimuth_flag
  type(ds_az_offset) :: az_offset
  complex(dpc), pointer, dimension(:) :: Qinv
  complex(dpc), pointer, dimension(:) :: Pinv
  logical ownsCorrelators
!  type(ds_modulescan),pointer :: my_modulescan    !points to the modulescan that owns this timestream
!  integer :: idiode                               !the diode number of this timestream in the above scan
end type ds_offsets
  
  
type ds_converter
  integer, allocatable, dimension(:) :: pixelTable
  integer length
  integer currentIndex_
  logical final_
end type ds_converter


interface prepareOffsets
    module procedure prepareOffsets_stan, prepareOffsets_withaz
end interface
 
interface saveOffsetsToFile
  module procedure saveOffsetsToFileName, saveOffsetsToFileNum
end interface
  
contains
	subroutine setup_moduleScan_ndiodes(ms,nt,np)
		type(ds_modulescan) :: ms
		integer nt, np
		integer nd 
		nd = nt+np
		
		ms%ndiodes_T = nt
		ms%ndiodes_P = np
		ms%ndiodes = nt+np
		
		
		ms%magic_check = MAGIC_VALUE

		allocate(ms%inv_Cw(nd,nd))
		ms%inv_Cw_TT => ms%inv_Cw(1:nt,1:nt)
		ms%inv_Cw_PP => ms%inv_Cw(nt+1:nd,nt+1:nd)
		ms%inv_Cw_TP => ms%inv_Cw(1:nt,nt+1:nd)

		allocate(ms%timestreams(nd))
		ms%timestreams_T => ms%timestreams(1:nt)
		ms%timestreams_P => ms%timestreams(nt+1:nd)

		allocate(ms%offsets_2(nd))
		ms%offsets_T => ms%offsets_2(1:nt)
		ms%offsets_P => ms%offsets_2(nt+1:nd)
		
		allocate(ms%flags(nd))
		ms%flags_T => ms%flags(1:nt)
		ms%flags_P => ms%flags(nt+1:nd)
		
		allocate(ms%dpsi(1:np))
		
		nullify(ms%theta)
		nullify(ms%pointing) !DWPS: always nullify pointers when you make them
		
		call prepare_noise(ms%noise, nd)
		
	end subroutine setup_moduleScan_ndiodes
  

	subroutine prepare_noise(noise,ndiodes)
		type(ds_noiseinfo) :: noise
		integer ndiodes
		
		allocate(noise%sigma(ndiodes))
		allocate(noise%alpha(ndiodes))
		allocate(noise%fknee(ndiodes))
		allocate(noise%corr(ndiodes,ndiodes))
	end subroutine prepare_noise
	
	subroutine destroy_noise(noise)
		type(ds_noiseinfo) :: noise
		if (associated(noise%sigma)) deallocate(noise%sigma)
		if (associated(noise%sigma)) deallocate(noise%alpha)
		if (associated(noise%sigma)) deallocate(noise%fknee)
		if (associated(noise%sigma)) deallocate(noise%corr)
	
	end subroutine destroy_noise

	subroutine prepare_covariance(cov, npix, has_t, has_p)
		type(ds_covariance) :: cov
		logical has_t, has_p
		integer npix
		
		cov%has_t = has_t
		cov%has_p = has_p
		cov%npix = npix
		

		if (has_t) then
			allocate(cov%TT(npix))
			cov%TT = 0
		endif
		
		if (has_p) then
			allocate(cov%QQ(npix))
			allocate(cov%UU(npix))
			allocate(cov%QU(npix))
			cov%QQ = 0
			cov%UU = 0
			cov%QU = 0
		endif

		if (has_p .and. has_t) then
			allocate(cov%TQ(npix))
			allocate(cov%TU(npix))
			cov%TQ = 0
			cov%TU = 0
		endif

	end subroutine prepare_covariance

	subroutine destroy_covariance(cov)
		type(ds_covariance) :: cov
		cov%npix = -1
		if (allocated(cov%TT)) deallocate(cov%TT)
		if (allocated(cov%TQ)) deallocate(cov%TQ)
		if (allocated(cov%TU)) deallocate(cov%TU)
		if (allocated(cov%QQ)) deallocate(cov%QQ)
		if (allocated(cov%QU)) deallocate(cov%QU)
		if (allocated(cov%UU)) deallocate(cov%UU)
		
		cov%has_t = .false.
		cov%has_p = .false.
	end subroutine destroy_covariance


	subroutine copy_covariance(source, dest)
		type(ds_covariance) :: source, dest
		
		call prepare_covariance(dest, source%npix, source%has_t, source%has_p)
		if (source%has_t) dest%TT=source%TT
		if (source%has_p) then
			dest%QQ=source%QQ
			dest%UU=source%UU
			dest%QU=source%QU
			if (source%has_t) then
				dest%TQ=source%TQ
				dest%TU=source%TU
			endif
		endif
	
	end subroutine copy_covariance




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

  subroutine prepareOffsets_withaz(offsets, na,offsetLength,azimuth_flag)
    type(ds_offsets) :: offsets
    integer, intent(in) :: na, offsetLength
    logical :: azimuth_flag
    
    offsets%length = offsetLength
    offsets%na = na
	call ds_checkAllocate(offsets%values, offsets%na)
	
    offsets%azimuth_flag= azimuth_flag
	nullify(offsets%Qinv)
	nullify(offsets%Pinv)
	offsets%ownsCorrelators=.false.
	
  end subroutine prepareOffsets_withaz

subroutine saveOffsetsToFileName(offsets,filename)
	type(ds_offsets) :: offsets
    character(*) :: filename
	integer unit
	unit=ds_get_lun()
	open(unit=unit,file=filename,form='unformatted',action='write')
	call saveOffsetsToFileNum(offsets,unit)
	close(unit)
end subroutine saveOffsetsToFileName

subroutine saveOffsetsToFileNum(offsets,unit)
	type(ds_offsets) :: offsets
	integer unit
	write(unit) offsets%length
	write(unit) offsets%na
	write(unit) offsets%values
	write(unit) offsets%azimuth_flag
	if (offsets%azimuth_flag) then
		write(unit) offsets%az_offset%n_az
		write(unit) offsets%az_offset%amplitudes
	endif
end subroutine saveOffsetsToFileNum



subroutine loadOffsetsFromFile(offsets,filename)
	type(ds_offsets) :: offsets
    character(*) :: filename
	integer unit
	integer na,lc
	unit=ds_get_lun()
	open(unit=unit,file=filename,form='unformatted',action='read')
	read(unit) lc
	read(unit) na
	call prepareOffsets(offsets,na,lc)
	read(unit) offsets%values
	close(unit)
end subroutine loadOffsetsFromFile

  
  subroutine copyOffsets(source, dest)
      type(ds_offsets) :: source, dest
	  integer i
      call prepareOffsets(dest, source%na, source%length)
	  do i= 1, source%na
			dest%values(i) = source%values(i)
	  enddo
	  
	  dest%azimuth_flag= source%azimuth_flag
      if(dest%azimuth_flag) call az_copy(source%az_offset, dest%az_offset)
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

   if(offsets%azimuth_flag) call az_destroy(offsets%az_offset)
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
  
  
  subroutine prepareMap(map, npix, zero_based)
    type(ds_map) :: map
    integer npix
	logical, optional :: zero_based
	integer :: start, finish
	
	start = 1
	if (present(zero_based)) then
		if (zero_based) start=0
	endif
	finish=start+npix-1
    
	allocate(map%map(start:finish))
	allocate(map%indices(start:finish))
	map%map = 0.0D0
	map%indices = 0
    map%npix = npix

  end subroutine prepareMap

	subroutine prepareTriMap(maps,npix,temperature, polarization, zero_based)
		type(ds_trimap) :: maps
		integer npix
		logical temperature, polarization
		logical, optional :: zero_based
		logical :: zero_based_
		maps%has_t = temperature
		maps%has_p = polarization
		
		zero_based_ = .false.
		if (present(zero_based)) zero_based_ = zero_based
		
		if (temperature)  call prepareMap(maps%T,npix,zero_based_)
		if (polarization) call prepareMap(maps%Q,npix,zero_based_)
		if (polarization) call prepareMap(maps%U,npix,zero_based_)
		
	end subroutine prepareTriMap



  subroutine destroyMap(map)
    type(ds_map) :: map
    if (allocated(map%map)) deallocate(map%map)
    if (allocated(map%indices)) deallocate(map%indices)
    map%npix=-1
  end subroutine destroyMap
  
	subroutine destroyTriMap(maps)
  	type(ds_trimap) :: maps
	call destroyMap(maps%T)
	call destroyMap(maps%Q)
	call destroyMap(maps%U)
	
	maps%has_t = .false.
	maps%has_p = .false.

	end subroutine destroyTriMap
	
	subroutine copyTriMap(source,dest)
		type(ds_trimap) :: source, dest
		dest%has_t = source%has_t
		dest%has_p = source%has_p
		if (source%has_t) call copyMap(source%T,dest%T)
		if (source%has_p) call copyMap(source%Q,dest%Q)
		if (source%has_p) call copyMap(source%U,dest%U)
	end subroutine copyTriMap

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
nullify(az%pointing)

end subroutine az_initialise




!subroutine to destroy a ds_az_offset instance
subroutine az_destroy(az)
type(ds_az_offset),intent(inout) :: az

az%n_az= -1
az%ntod= -1
if(allocated(az%amplitudes)) deallocate(az%amplitudes)
!if(allocated(az%pointing)) deallocate(az%pointing)
if(allocated(az%hits)) deallocate(az%hits)

end subroutine az_destroy



!subroutine to copy an azimuth pointing (except the amplitudes field which is allocated and blank)
subroutine az_copy(az_in,az_out)
!Initialise az_out using array sizes from az_in AND copies the hits from az_in
!Will point to the pointing of cache_in (not duplicate)
type(ds_az_offset),intent(in) :: az_in
type(ds_az_offset),intent(inout) :: az_out

call az_initialise(az_in%ntod,az_in%n_az,az_out)
if(associated(az_in%pointing)) then
    az_out%pointing => az_in%pointing
else
    nullify(az_out%pointing)
endif
az_out%hits = az_in%hits
az_out%amplitudes= az_in%amplitudes

end subroutine az_copy

subroutine copy_modulescan(source,dest,duplicateoffsets,duplicatetimestreams,shareoffsets,sharetimestreams, duplicatemaps, duplicatepointing)
!copies modulescan object.  The pointing and theta pointers are shared.  By default, the offsets and
!timestream pointers are nullified. To share either, set shareoffsets and/or sharetimestreams to 1.
!To duplicate either (make own hard copy with allocated memory), set duplicateoffsets and/or 
!duplicatetimestreams to 1.

    type(ds_modulescan),intent(in) :: source
    type(ds_modulescan),intent(inout) :: dest
    logical,optional :: duplicateoffsets
    logical,optional :: duplicatetimestreams
    logical,optional :: duplicatepointing
    logical,optional :: shareoffsets
    logical,optional :: sharetimestreams
    logical,optional :: duplicatemaps
	integer npix
	integer i

	call check_modulescan_initialized(dest, "Copy modulescan")
    !eventually add an assertion that sharing and copying cannot both be true.
	if (present(sharetimestreams).and. present(duplicatetimestreams))  then
		call ds_assert(.not.(sharetimestreams .and. duplicatetimestreams),"Cannot both copy and share timestreams")
	endif
	if (present(shareoffsets).and. present(duplicateoffsets))  then
		call ds_assert(.not.(shareoffsets .and. duplicateoffsets),"Cannot both copy and share timestreams")
	endif

    dest%ntod    = source%ntod
	dest%ndiodes_T = source%ndiodes_T
	dest%ndiodes_P = source%ndiodes_P
	dest%ndiodes = source%ndiodes
	dest%id      = source%id
    dest%dpsi    = source%dpsi
    dest%inv_Cw  = source%inv_Cw
    dest%flags   = source%flags
    dest%flags_T = source%flags_T
    dest%flags_P = source%flags_P
    dest%theta => source%theta
	dest%has_T = source%has_T
	dest%has_P = source%has_P
	dest%noise%sigma = source%noise%sigma
	dest%noise%alpha = source%noise%alpha
	dest%noise%fknee = source%noise%fknee
	dest%noise%corr = source%noise%corr
	
	
	if (present(duplicatepointing)) then
		if (duplicatepointing) then
			allocate(dest%pointing(size(source%pointing) ))
			dest%pointing = source%pointing
		else
			dest%pointing => source%pointing
		endif
	else
		dest%pointing => source%pointing
    endif
	
	
	
!    nullify(dest%timestreams)
    if(present(sharetimestreams)) then
        if(sharetimestreams) then
 			call ds_assert(.false., "Coding Error: sharing of timestreams is broken in copy_modulescan")
           dest%timestreams => source%timestreams
        endif
    endif
    if(present(duplicatetimestreams)) then
        if(duplicatetimestreams) then
!            nullify(dest%timestreams)
!            allocate(dest%timestreams(1:ndiodes_max))
            do i=1,source%ndiodes
                if(source%flags(i)==0) cycle !ignore if bad
                call copytimestream(source%timestreams(i),dest%timestreams(i))
            enddo
        endif
    endif
    
!    nullify(dest%offsets_2)
    if(present(shareoffsets)) then
		call ds_assert(.false., "Coding Error: sharing of offsets is broken in copy_modulescan")
        if(shareoffsets) dest%offsets_2 => source%offsets_2
    endif

    if(present(duplicateoffsets)) then
        if(duplicateoffsets) then
!            nullify(dest%offsets_2)
!            allocate(dest%offsets_2(1:ndiodes_max))
            do i=1,source%ndiodes
                if(source%flags(i)==0) cycle !ignore if bad
                call copyOffsets(source%offsets_2(i),dest%offsets_2(i))
            enddo
        endif
    endif

    if(present(duplicatemaps)) then
        if(duplicatemaps) then
			call destroytrimap(dest%maps)
			npix = 0
			if (source%maps%has_t) npix = source%maps%T%npix
			if (source%maps%has_p) npix = source%maps%Q%npix
			
!			call prepareTriMap(dest%maps, npix, source%has_t, source%has_p)
			call copyTriMap(source%maps, dest%maps)
        endif
    endif


!	dest%offsets_T => dest%offsets_2(1:ndiodes_T)
!	dest%offsets_P => dest%offsets_2(ndiodes_T+1:ndiodes_max)
!	dest%timestreams_T => dest%timestreams(1:ndiodes_T)
!	dest%timestreams_P => dest%timestreams(ndiodes_T+1:ndiodes_max)
    

end subroutine copy_modulescan


subroutine destroy_modulescan(arg,deallocateoffsets,deallocatepointing)
!destroys pointer array of modulescans. Should work on an instance too.
!Default is to nullify offset pointer. Set deallocateoffsets=.true. to deallocate
    type(ds_modulescan),pointer,dimension(:) :: arg
    logical,optional :: deallocateoffsets, deallocatepointing
    integer :: m, nmodules
    logical :: deall, depoint
	type(ds_modulescan), pointer :: A

    deall= .false.
	depoint=.false.
    if(present(deallocateoffsets)) deall= deallocateoffsets
    if(present(deallocatepointing)) depoint = deallocatepointing

    if(.not.associated(arg)) return

    nmodules= size(arg)
    do m= 0,nmodules-1
		A => arg(m)
		A%id=-1
		A%magic_check = 0
        !destroy instance of arg
        nullify(A%theta)
		if (depoint) then
			deallocate(A%pointing)
		endif
        nullify(A%pointing)
        if(deall) then
            if(associated(A%offsets_2)) deallocate(A%offsets_2)
			nullify(A%offsets_T)
			nullify(A%offsets_P)
        else
            nullify(A%offsets_2)
            nullify(A%offsets_T)
            nullify(A%offsets_P)
        endif
		if (associated(A%flags)) deallocate(A%flags)
		nullify(A%flags_T)
		nullify(A%flags_P)

		if (associated(A%inv_Cw)) deallocate(A%inv_Cw)
		nullify(A%inv_Cw)
		nullify(A%inv_Cw_TT)
		nullify(A%inv_Cw_PP)
		nullify(A%inv_Cw_TP)

		if (associated(A%timestreams)) deallocate(A%timestreams)
		nullify(A%timestreams)
		nullify(A%timestreams_T)
		nullify(A%timestreams_P)
		call destroyTriMap(A%maps)
		call destroy_noise(A%noise)
    enddo
    
    deallocate(arg)

end subroutine destroy_modulescan

subroutine check_modulescan_initialized(ms, message)
	type(ds_modulescan) :: ms
	character(*) :: message
	
	call ds_assert(ms%magic_check==MAGIC_VALUE, message)
end subroutine check_modulescan_initialized


subroutine saveAllOffsetsToFiles(moduleScans,dir)
	type(ds_modulescan), pointer, dimension(:) :: moduleScans
	character(*) dir
    type(ds_modulescan), pointer :: moduleScan
    integer i,d
    character(512) :: filename
    do i=0,size(moduleScans)-1
       moduleScan=>moduleScans(i)
		do d=1,moduleScan%ndiodes
       		filename=filenameForSavedOffset(dir,moduleScan%id(1),moduleScan%id(2),d)
			if (modulescan%flags(d) .ne. 0) call saveOffsetsToFile(moduleScan%offsets_2(d),filename)
       enddo
    enddo
end subroutine saveAllOffsetsToFiles



function filenameForSavedOffset(directory,run,scan,diode) result(f)
	character(512) :: f
	character(*) :: directory
	character(*), parameter :: fmt = '( A, "/", I5.5, "_",I5.5,"_", I5.5, ".off"  )'
	integer :: run, scan, diode
	write(f,fmt) trim(directory),run,scan,diode
end function filenameForSavedOffset



end module
