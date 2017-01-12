module ds_azimuth
use ds_types
use ds_utils
implicit none


!keep all of the azimuth pointing targets in this array
!type(ds_az_retriever),pointer,dimension(:) :: cache

contains


!subroutine to include azimuth offsets in the inner products
function az_inner_product(az1,az2) result(inner)
type(ds_az_offset) :: az1, az2
real(dp) :: inner

!#warning !big nbug in az_inner_product
!print*,'Big bug here: this should take inner product of two vectors, not just one input'

inner= 0.0_8
inner= sum( az1%amplitudes * az2%amplitudes )
!to use this, add to the sum in inner before the square root

end function az_inner_product



!subroutine to project az offset onto tod
subroutine az_add_offset_to_tod(az,tod)
type(ds_az_offset) :: az
type(ds_timestream) :: tod
integer :: itod

!print*,"PROJECT"

call ds_assert(az%ntod == tod%nt,"Error in azimuth offset projection: arrays don't match")
    
!Add appropriate offset value to TOD
do itod=1,tod%nt
    tod%timestreams(itod) = tod%timestreams(itod) + az%amplitudes( az%pointing(itod) )
enddo

end subroutine az_add_offset_to_tod



!subroutine to deproject tod into az offset
subroutine az_deproject_tod_onto_offset(tod,az)
type(ds_timestream) :: tod
type(ds_az_offset) :: az
integer :: itod

!print*,"BIN"

call ds_assert(az%ntod == tod%nt,"Error in azimuth offset deprojection: arrays don't match")

!sum the TOD into the offsets
az%amplitudes= 0.0_8
do itod = 1,tod%nt
    az%amplitudes( az%pointing(itod) )= &
    	az%amplitudes( az%pointing(itod) ) + tod%timestreams(itod)
enddo

!complete the average - no, this is a summation
!az%amplitudes = az%amplitudes / az%hits

end subroutine az_deproject_tod_onto_offset





!DWPS: this will now have to "eat" a modulescan object. Point each diode's az_offset to the correct
!az pointing cache for the horn.

subroutine az_create_pointing(ms,memory_diode,corrlength,n_bins,azimuth_scan)

type(ds_modulescan),intent(inout) :: ms
integer,intent(in) :: corrlength,n_bins,memory_diode
real(dp),dimension(ms%ntod),intent(inout) :: azimuth_scan

!!subroutine to create an instance of the azimuth pointing from azimuth data
!subroutine az_create_pointing(iscan,ntod,corrlength,n_bins,azimuth_scan,az_offset)
!iscan is the global index of this timestream
!integer :: iscan,ntod, corrlength, n_bins
!real(dp),dimension(ntod),intent(in) :: azimuth_scan
!type(ds_az_offset),intent(inout) :: az_offset


real(dp) :: minv,maxv,binsize
integer :: idiode, itod, ioffset, iset, ibin, n_sets, n_az, leftover, ntod, t
logical :: absorb
character(len=120) :: mssg

!corrlength is the number of TOD to bin into a particular azimuth bin set
ntod= ms%ntod
print*,ntod,corrlength
n_sets= ntod / corrlength  !DWPS: corrlength MUST be longer than two full sweeps

if(ntod < corrlength) then
	n_sets=1
endif

!leftover= mod(ntod,corrlength)
!if(leftover .le. corrlength/2) then
!    !absorb
!    absorb= .true.
!else
!    !make new ?and reset penultimate boundary?
!    absorb= .false.
!    n_sets= n_sets+1
! endif

absorb= .true.
n_az= n_bins * n_sets      !total number of bins for this timestream

!initialise the azimuth offsets for each diode
do idiode = 1,ms%ndiodes
	if(ms%flags(idiode)==0) cycle
	call az_initialise(ntod,n_az,ms%offsets_2(idiode)%az_offset)
enddo


minv = minval(azimuth_scan)
maxv = maxval(azimuth_scan)

if(minv .le. 0.1*pi .and. maxv .ge. 1.9 * pi) then
    !scan straddles azimuth discontinuity - redefine coords
	do t=1,ntod
		if (azimuth_scan(t) .le. pi) azimuth_scan(t) = 2*pi + azimuth_scan(t)
	enddo
!    azimuth_scan= azimuth_scan - 180.0_dp

endif

!Set up the bin grid using min, max and nbins
binsize= (maxv - minv) / real(n_bins,kind=dp)

!create an azimuth pointing instance in the cache for this modulescan
!call ds_assert(associated(cache),'cache not allocated yet in az_create_pointing')
!if(allocated(cache(ims)%pointing)) deallocate(cache(ims)%pointing)
!allocate(cache(ims)%pointing(ntod))

allocate(ms%offsets_2(memory_diode)%az_offset%pointing(ntod))


!Loop through the azimuth scan and assign each datum to a bin, counting #hits as we go
do itod= 1,ntod
    if(azimuth_scan(itod) .le. minv) then
        ibin=1
    elseif(azimuth_scan(itod) .ge. maxv) then
    	ibin= n_bins
    else
        ibin = ceiling((azimuth_scan(itod) - minv) / binsize)   !which azimuth bin
    endif
    call ds_assert(ibin .le. n_bins,"ibin > n_bins")
    
    iset= (itod-1) / corrlength           !which set of azimuth bins (ZERO BASED)
    
    if(absorb .and. iset==n_sets) iset= iset-1 !absorb into last bin
    
    ioffset= ibin + iset * n_bins         !this offset index
    mssg= 'Error in building azimuth scan: ioffset out of range: '//trim(adjustl(inttostr(ioffset)))//trim(adjustl(inttostr(n_az)))//trim(adjustl(inttostr(itod)))//trim(adjustl(inttostr(ntod)))
    call ds_assert((ioffset .le. n_az) .and. (ioffset .gt. 0),mssg)
    
    ms%offsets_2(memory_diode)%az_offset%pointing(itod)= ioffset
    ms%offsets_2(memory_diode)%az_offset%hits(ioffset)= ms%offsets_2(1)%az_offset%hits(ioffset) + 1
enddo

!do itod=1,n_az
!	print*,ms%id(1),ms%id(2),ms%id(3),itod,ms%offsets(memory_diode)%az_offset%hits(itod)
!enddo

call ds_assert(minval(ms%offsets_2(memory_diode)%az_offset%hits) .gt. 0,"Zero hits in an azimuth bin")

!direct each diode's az_offsets pointing to correct modules scan pointing cache
do idiode=1,ms%ndiodes
	if(ms%diode_flag(idiode)==0) cycle  !ignore if bad
    if(idiode == memory_diode) cycle

	ms%offsets_2(idiode)%az_offset%hits= ms%offsets_2(memory_diode)%az_offset%hits
	ms%offsets_2(idiode)%az_offset%pointing => ms%offsets_2(memory_diode)%az_offset%pointing
enddo

end subroutine az_create_pointing



!Utilities: 

!create azimuth cache structure so it can start holding pointings
!subroutine create_azimuth_cache(n)
!n is the total number of timestreams globally
!    integer :: n, i
!    allocate(cache(0:n-1)) !timestreams zero indexed
!    do i=0,n-1
!        nullify(cache(i)%pointing)
!    enddo
!end subroutine create_azimuth_cache

!subroutine to destroy the azimuth cache structure
!subroutine destroy_azimuth_cache()
!    integer :: i,n
!    n= size(cache)
!    do i=0,n-1
!        if(allocated(cache(i)%pointing)) deallocate(cache(i)%pointing)    
!    enddo
!    deallocate(cache)
!end subroutine destroy_azimuth_cache

subroutine az_build_pointing_hits(pointing,hits,n_bins,corrlength,azimuth_scan)
	real(dp), dimension(1:)  :: azimuth_scan
	real(dp) :: az
	integer, dimension(:), pointer :: pointing, hits
	integer :: n_bins, n_sets,corrlength,n_az
	real(dp) :: binsize
	logical, parameter :: absorb= .true.
	integer ntod
	character(256) :: mssg
	real(dp) :: minv,maxv
	
	integer :: t,bin,a,set

	ntod=size(azimuth_scan)

	n_sets= ntod / corrlength  !DWPS: corrlength MUST be longer than two full sweeps
	n_az= n_bins * n_sets      !total number of bins for this timestream

	allocate(hits(n_az))
	allocate(pointing(ntod))

	
	minv=minval(azimuth_scan)
	maxv=maxval(azimuth_scan)

	binsize= (maxv - minv) / real(n_bins,kind=dp)
	
	do t= 1,ntod
		az=azimuth_scan(t)
		!Determine azimuth bin
	    if(az .le. minv) then
	        bin=1
	    elseif(az .ge. maxv) then
	    	bin= n_bins
	    else
	        bin = ceiling((az - minv) / binsize)
	    endif
	    call ds_assert(bin .le. n_bins,"ibin > n_bins")
    

		!which set of azimuth bins (ZERO BASED)
	    set= (t-1) / corrlength           
    
	    if(absorb .and. set==n_sets) set= set-1 !absorb into last bin
    
	    a= bin + set * n_bins         !this offset index
		write(mssg,*) "Error in building azimuth scan: ioffset out of range:",a,n_az,t,ntod
	    call ds_assert((a .le. n_az) .and. (a .gt. 0),mssg)
    
	    pointing(t)=a
	    hits(a)= hits(a) + 1
	enddo
end subroutine az_build_pointing_hits

end module ds_azimuth