module rms_mod
  use quiet_utils
  use quiet_fileutils
  use alm_tools
  implicit none

contains
  
  subroutine make_rmsmap(file_in, file_mask, file_rms, rms_default)

    character(len=512)                          :: file_in, file_mask, file_rms
    logical(lgt)                                :: rms_default
    integer(i4b)                                :: order, nside, lmin, lmax, n_rms, pix_rms, i, nlist, n_masked
    integer(i4b), dimension(:), allocatable     :: listpix
    real(dp)                                    :: radius, mu
    real(dp), dimension(3)                      :: vector
    real(dp), dimension(:,:), allocatable       :: fgmap, rmsmap, mask
    complex(dpc), dimension(:,:,:), allocatable :: alm
    
    write(*,*) 'RMS-computing in progress...'

    ! Reading input map
    call read_map(fgmap, order, file_in)
    nside = npix2nside(size(fgmap,1))
    if (order == nest) call convert_nest2ring(nside, fgmap)
    
    ! Converting to harmonics 
    lmin = 25
    lmax = 300
    allocate(alm(1:3, 0:lmax, 0:lmax))
     
    call map2alm(nside, lmax, lmax, fgmap(:,1:3), alm)

    ! Including only l's between 25 and 300
    alm(:,0:25,:) = 0.d0
    
    ! Converting back
    call alm2map(nside, lmax, lmax, alm, fgmap(:,1:3))
    deallocate(alm)

    ! Storing P = sqrt(U^2 + Q^2) in 1st map layer (since this is redundant)
    fgmap(:,1) = sqrt(fgmap(:,2)**2 + fgmap(:,3)**2)

    ! Allocating the output map with nside=64
    n_rms = 64
    pix_rms = nside2npix(n_rms)
    radius = 10.d0 * DEG2RAD 
    
    allocate(rmsmap(0:pix_rms-1,1))
    allocate(listpix(0:size(fgmap,1)-1))

    ! Masking the foreground
    call read_map(mask, order, file_mask)
    if (order==2) call convert_nest2ring(nside, mask)
    fgmap(:,1) = fgmap(:,1)*mask(:,1)
    
    ! Storing RMS value of P on a disk of 10 degrees radius for each pixel in rmsmap
    do i=0,pix_rms-1
       vector = pix2vec(n_rms, ring, i)
       call query_disc(nside, vector, radius, listpix, nlist)
       n_masked = sum(mask(listpix(0:nlist-1),1))
       if (n_masked.gt. 0.0) then
          mu = sum(fgmap(listpix(0:nlist-1),1))/n_masked
          rmsmap(i,1) = sqrt(sum((fgmap(listpix(0:nlist-1),1)-mu)**2) / n_masked)
       else
          rmsmap(i,1) = -1.6375d30
       endif
       
    end do
    write(*,*) "Done - writing map"

    if (rms_default) file_rms = 'rms.fits'
    call write_map(rmsmap, ring, file_rms)
    
    deallocate(fgmap)
    deallocate(rmsmap)
    deallocate(listpix)
    deallocate(mask)
    
  end subroutine make_rmsmap
end module rms_mod

