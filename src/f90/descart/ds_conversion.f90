module ds_conversion
  
  use ds_types
  use ds_utils
  implicit none

real, parameter ::  converterLengthFactor = 2

  
  contains


subroutine makePointingContinuous(pointing,maxIndex,reverser)
!This subroutine resets the pointing of a ds_pointing instance so that the 
!pixel indices it uses run from 1..npix (and can therefore be used as array indices)
!The maxIndex argument specifies the number of pixels in the input class,
!for example if you are converting a set of healpix indices then maxIndex should be 12*Nside^2
  type(ds_detpointing) :: pointing
  type(ds_converter), intent(inout) :: reverser
  integer(dp), allocatable, dimension(:) :: inputMapHit
  integer maxIndex
  integer t
  integer uniquePix,currentPixel
  allocate(inputMapHit(maxIndex))
  
  inputMapHit = 0
  uniquePix = 0
  currentPixel = 1
  
  call ds_converter_clear(reverser)
  call ds_converter_init(reverser,1000)
  
  do t = 1,pointing%nt
    if (inputMapHit(pointing%pixel(t)) == 0) then
      inputMapHit(pointing%pixel(t)) = currentPixel
      call ds_converter_append(reverser,pointing%pixel(t))
      pointing%pixel(t) = currentPixel !inputMapHit(pointing%pixel(t))
      currentPixel = currentPixel +1
    endif
  enddo
  
  pointing%npix = currentPixel-1
  deallocate(inputMapHit)
  call ds_converter_finalize(reverser)
end subroutine



 subroutine resetOriginalMapIndices(map, reverser)
 !Not sure if we want this routine.
 !It resets the map indices.
  type(ds_map), intent(inout) :: map
  type(ds_converter), intent(in) :: reverser !I think this should be intent in
  integer p
  
  call ds_assert(map%npix == reverser%length,'Pointing and index-reverser have different npix')
  do p=1,map%npix
  	map%indices(p) = reverser%pixelTable( map%indices(p) )
  enddo 
end subroutine resetOriginalMapIndices


 subroutine resetOriginalPointingIndices(pointing, reverser)
  type(ds_detpointing) :: pointing
  type(ds_converter), intent(in) :: reverser
  call ds_assert(pointing%npix == reverser%length,'Pointing and index-reverser have different npix')
  call ds_converter_apply(reverser,pointing)
   
end subroutine resetOriginalPointingIndices




!----------------------------------------------------!
!Methods to initialize and use the ds_converter type.
!The only one 

  subroutine ds_converter_apply(converter,pointing)
  !Apply a pointing converter.
  !Typical usage is to use with the output reverser from the makePointingContinuous subroutine
  !To return to the original pointing or map.
    type(ds_converter) :: converter
    type(ds_detpointing) :: pointing
    integer t
    do t = 1,pointing%nt
      pointing%pixel(t) = converter%pixelTable(pointing%pixel(t))
    end do
  end subroutine ds_converter_apply
  


subroutine ds_converter_init(converter,initialLength)
  type(ds_converter) :: converter
  integer :: initialLength
  converter%length = initialLength
  converter%currentIndex_ = 1
  converter%final_ = .false.
  allocate(converter%pixelTable(initialLength))
end subroutine


subroutine ds_converter_finalize(converter)
  type(ds_converter) :: converter
  if(converter%final_) return
  call ds_converter_realloc(converter,converter%currentIndex_ - 1)
  converter%final_ = .true.
end subroutine


subroutine ds_converter_realloc(converter,newLength)
!This subroutine lengthens a ds_converter's internal pixel table.
!It should only be used internally by the converter methods - 
!You should not use it yourself.
  type(ds_converter) :: converter
  integer, allocatable, dimension(:) :: tempArray
  integer oldLength, newLength, copyLength
  oldLength = converter%length-1
  copyLength = min(oldLength,newLength)
  allocate(tempArray(copyLength))
  tempArray = converter%pixelTable(1:copyLength)
  deallocate(converter%pixelTable)
  allocate(converter%pixelTable(newLength))
  converter%pixelTable(1:copyLength) = tempArray
  converter%length = newLength
  deallocate(tempArray)
end subroutine


subroutine ds_converter_append(converter,pixel)
  type(ds_converter) :: converter
  integer :: pixel
  call ds_assert(.not. converter%final_,'Converter is finalized')
  if (converter%currentIndex_ == converter%length+1) then
    call ds_converter_realloc(converter,nint(converter%length*converterLengthFactor))
  endif
  converter%pixelTable(converter%currentIndex_) = pixel
  converter%currentIndex_ = converter%currentIndex_ + 1
end subroutine



subroutine ds_converter_clear(converter)
  type(ds_converter) :: converter
  if (allocated(converter%pixelTable)) deallocate(converter%pixelTable)
  converter%length = 0
  converter%currentIndex_ = 1
  converter%final_ = .false.

end subroutine





end module ds_conversion



