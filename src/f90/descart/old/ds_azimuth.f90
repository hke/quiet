module ds_azimuth
use ds_types
use ds_utils
implicit none


!keep all of the azimuth pointing targets in this array
type(ds_az_retriever),pointer,dimension(:) :: cache

contains


!subroutine to include azimuth offsets in the inner products
function az_inner_product(az) result(inner)
type(ds_az_offset) :: az
real(dp) :: inner

inner= 0.0_8
inner= sum( az%amplitudes**2 )
!to use this, add to the sum in inner before the square root

end function az_inner_product



!subroutine to project az offset onto tod
subroutine az_project_offset_onto_tod(az,tod)
type(ds_az_offset) :: az
type(ds_timestream) :: tod
integer :: itod

call ds_assert(az%ntod == tod%nt,"Error in azimuth offset projection: arrays don't match")
    
!Set TOD to appropriate offset value
do itod=1,tod%nt
    tod%timestream(itod)= az%amplitudes( az%cache_retrieve%pointing(itod) )
enddo

end subroutine az_project_offset_onto_tod



!subroutine to deproject tod into az offset
subroutine az_deproject_tod_onto_offset(tod,az)
type(ds_timestream) :: tod
type(ds_az_offset) :: az
integer :: itod

call ds_assert(az%ntod == tod%nt,"Error in azimuth offset deprojection: arrays don't match")

!sum the TOD into the offsets
az%amplitudes= 0.0_8
do itod= 1,tod%nt
    az%amplitudes( az%cache_retrieve%pointing(itod) )= az%amplitudes( az%cache_retrieve%pointing(itod) ) + tod%timestream(itod)
enddo

!complete the average - no, this is a summation
!az%amplitudes = az%amplitudes / az%hits

end subroutine az_deproject_tod_onto_offset





!subroutine to create an instance of the azimuth pointing from azimuth data
subroutine az_create_pointing(iscan,ntod,corrlength,n_bins,azimuth_scan,az_offset)
!iscan is the global index of this timestream
integer :: iscan,ntod, corrlength, n_bins
real(dp),dimension(ntod),intent(in) :: azimuth_scan
type(ds_az_offset),intent(inout) :: az_offset
real(dp) :: minv,maxv,binsize
integer :: itod, ioffset, iset, ibin, n_sets, n_az, leftover
logical :: absorb

!corrlength is the number of TOD to bin into a particular azimuth bin set
n_sets= ntod / corrlength  !DWPS: corrlength MUST be a factor of ntod and longer than two full sweeps

leftover= mod(ntod,corrlength)
if(leftover .le. corrlength/2) then
    !absorb
    absorb= .true.
else
    !make new ?and reset penultimate boundary?
    absorb= .false.
    n_sets= n_sets+1
endif

n_az= n_bins * n_sets      !total number of bins for this timestream

!#warning !uncomment this call eventually
call az_initialise(ntod,n_az,az_offset)
!print*,"warning: missed commented call to az_initalise"
!stop

!Set up the bin grid using min, max and nbins
minv= minval(azimuth_scan) !DWPS: what about discontinuities - like the 0 or 360 deg meridian?
maxv= maxval(azimuth_scan)
binsize= (maxv - minv) / real(n_bins,kind=dp)

!create an azimuth pointing instance in the cache
call ds_assert(associated(cache),'cache not allocated yet in az_create_pointing')
if(allocated(cache(iscan)%pointing)) deallocate(cache(iscan)%pointing)
allocate(cache(iscan)%pointing(ntod))

!Loop through the tod and assign each one to a bin, counting #hits as we go
do itod= 1,az_offset%ntod

    ibin = ceiling((azimuth_scan(itod) - minv) / binsize)   !which azimuth bin
    call ds_assert(ibin .le. n_bins,"ibin > n_bins")
    
    iset= (itod-1) / corrlength           !which set of azimuth bins
    
    if(absorb .and. iset==n_sets) iset= iset-1 !absorb into last bin
    
    ioffset= ibin + iset * n_bins         !this offset index
    call ds_assert((ioffset .le. n_az) .and. (ioffset .gt. 0), 'Error in building azimuth scan: ioffset out of range')
    
    cache(iscan)%pointing(itod)= ioffset
    az_offset%hits(ioffset)= az_offset%hits(ioffset) + 1
enddo

!direct this az_offsets pointing to correct pointing array
az_offset%cache_retrieve => cache(iscan)!%pointing

end subroutine az_create_pointing



!Utilities: 

!create azimuth cache structure so it can start holding pointings
subroutine create_azimuth_cache(n)
!n is the total number of timestreams globally
    integer :: n, i
    allocate(cache(0:n-1)) !timestreams zero indexed
!    do i=0,n-1
!        nullify(cache(i)%pointing)
!    enddo
end subroutine create_azimuth_cache

!subroutine to destroy the azimuth cache structure
subroutine destroy_azimuth_cache()
    integer :: i,n
    n= size(cache)
    do i=0,n-1
        if(allocated(cache(i)%pointing)) deallocate(cache(i)%pointing)    
    enddo
    deallocate(cache)
end subroutine destroy_azimuth_cache


end module ds_azimuth