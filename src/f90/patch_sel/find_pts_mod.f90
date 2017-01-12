module find_pts_mod
  use quiet_fileutils
  use quiet_utils
  use quiet_pointing_mod
  implicit none

contains
  
  subroutine find_low_fg_pts(file_in, npix_in, elevation_limit)
    character(len=512), intent(in)      :: file_in
    character(len=512)                  :: file_list, file_map
    integer(i4b)                        :: order, npix_in, nside, smallind, i, unit
    real(dp)                            :: elevation_limit, theta, phi, theta_cel, phi_cel, phi_rot, theta_rot, psi_rot, smallval
    real(dp), dimension(3)              :: vector, vector_cel
    real(dp), dimension(3,3)            :: euler_matrix
    real(dp), dimension(:), allocatable :: rmsmap, coldpts
    logical(lgt), dimension(:), allocatable :: mask
    
    ! input is file in galactic coordinates - insert an option for rotating if it is not?

    ! Read precomputed foreground RMS    
    call read_map(rmsmap, order, file_in)
    nside = npix2nside(size(rmsmap,1))

    ! Computing euler matrix for conversion from galactic to equatorial coordinates
    call compute_euler_matrix_zyz(-j2000_euler_gamma*deg2rad, -j2000_euler_beta*deg2rad, -j2000_euler_alpha*deg2rad, euler_matrix)
    
    ! Mask to use during sorting
    allocate(mask(size(rmsmap,1)))
    mask = rmsmap .gt. 0.d0
    
    allocate(coldpts(0:size(rmsmap,1)-1))
    coldpts = -1.6375d30
    
    unit = getlun()
    file_list = 'coords'//trim(itoa(npix_in))//'.dat'
    open(unit, file = file_list)
    write(unit,"(7G15.5)") 'No.', 'Pixel', 'Gal.long', 'Gal.lat', 'RA', 'Dec', 'RMS'
    
    i=1
    do while (i<=npix_in)
       ! Find smallest, store, update mask
       smallval = minval(rmsmap,1,mask)
       smallind = minloc(rmsmap,1,mask) ! Minloc counts from 1 even though rmsmap doesn't
       mask(smallind) = .false.
       vector = pix2vec(nside, order, smallind-1)
       vector_cel = matmul(euler_matrix, vector)
       
       ! Converting to galactic and equatorial angles
       call vec2ang(vector_cel, theta_cel, phi_cel)
       call vec2ang(vector, theta, phi)
       
       ! Points with declination above a certain limit are not visible from QUIET site and thus not stored
       if (theta_cel*RAD2DEG .gt. elevation_limit+23.d0) then
          coldpts(smallind-1) = rmsmap(smallind-1)
          write(unit,"(7G15.5)") i, smallind-1, phi*RAD2DEG, 90.d0-theta*RAD2DEG, phi_cel*RAD2DEG, 90.d0-theta_cel*RAD2DEG, smallval
          i=i+1
       end if
    end do
    
    close(unit)

    file_map = 'coldpts'//trim(itoa(npix_in))//'.fits'
    call write_map(coldpts, ring, file_map)
    
    deallocate(rmsmap)
    deallocate(mask)
    deallocate(coldpts)
    
  end subroutine find_low_fg_pts
end module find_pts_mod

