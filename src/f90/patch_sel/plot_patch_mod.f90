module plot_patch_mod
  use quiet_fileutils
  use quiet_utils
  implicit none

  type patch
     real(dp)           :: pos(2)
     character(len=512) :: name
  end type patch

contains

  subroutine plot_patch_and_elevation(file_in, file_out, file_rms, ndays, rms_default, make_map)

    ! input file with list of points, plot to .fits map if make_map==T
    ! must be of format galactic long, lat in first two columns
    
    character(len=512)            :: file_in, file_out, file_rms, line, root
    logical(lgt)                  :: make_map, rms_default
    integer(i4b)                  :: ndays, unit, length, i, order, ind, nside, pixel
    real(dp), dimension(:), allocatable    :: rmsmap
    type(patch), dimension(:), allocatable :: patches

    if(rms_default) file_rms = 'rms.fits'
    call read_map(rmsmap, order, file_rms)
    nside = npix2nside(size(rmsmap,1))

    ! Read through file to find number of lines
    unit = getlun()
    open(unit, file=file_in, action="read", status="old")
    length = 0
    do
       read(unit, '(a)', end=1) line
       length = length+1
    end do
1   rewind(unit)
    write(*,*) length
    allocate(patches(length))
    do i = 1, length
       read(unit,*) patches(i)%pos, patches(i)%name ! Long, lat, name
       patches(i)%pos(1) = patches(i)%pos(1)*DEG2RAD
       patches(i)%pos(2) = (90.d0-patches(i)%pos(2))*DEG2RAD
    end do
    close(unit)
    
        
    ! Making a root for output filenames
    ind = index(file_out, '.dat', .false.)
    write(root, *) file_out(1:ind-1)
    
    unit = getlun()
    open(unit, file=file_out, action="write", status="replace")
    write(unit,'(2(A10), A14, 3X, A)') "long", "lat", "rms", "Name"

    do i=1,length
       call write_elevation(ndays, patches(i), root)
       pixel = ang2pix(nside, order, patches(i)%pos(2), patches(i)%pos(1))
       write(unit,'(2(F10.2), ES14.5,3X, A)') patches(i)%pos(1)*RAD2DEG, 90.d0 - patches(i)%pos(2)*RAD2DEG, rmsmap(pixel), trim(patches(i)%name)

       if (make_map) then
          rmsmap(pixel) = 0.d0
          rmsmap(pixel-1) = -1.6375d30
          rmsmap(pixel+1) = -1.6375d30
       endif
    end do
    
    call sunmoon(ndays)
    
    ! Writing .fits file with patches marked
    if (make_map) then
       ! Splicing together a reasonable filename
       file_in = trim(root)//'_patchmap.fits'
       call write_map(rmsmap, order, file_in)
    endif
    
    deallocate(rmsmap)
    deallocate(patches)    
  end subroutine plot_patch_and_elevation
  
  
  ! This sub writes el/az coordinates for a given number of days for a given point in galactic long/lat, for the observation
  ! site in Chile.
  subroutine write_elevation(ndays, p, root)
    use quiet_fileutils
    use quiet_pointing_mod
    implicit none
    
    character(len=512)           :: root, filename
    integer(i4b)                 :: ndays, i, n, unit
    real(dp)                     :: theta, phi, mjd0, hour, mjd
    real(dp), dimension(3)       :: position, hor_pos
    real(dp), dimension(3,3)     :: euler_matrix
    type(patch)                  :: p

    ! Converting to celestial coordinates
    call ang2vec(p%pos(2), p%pos(1), position)
    call compute_euler_matrix_zyz(-j2000_euler_gamma*deg2rad, -j2000_euler_beta*deg2rad, -j2000_euler_alpha*deg2rad, euler_matrix)
    position = matmul(euler_matrix, position)
    
    ! Number of hours/time samples
    n = 24*ndays
    
    filename = trim(root)//'_el_'//trim(p%name)//'.dat'
    unit = getlun()
    open(unit, file=filename, action="write", status="replace")
    write(unit,'(12A, 10A, 10A)') 'MJD', 'El', 'Az'

    mjd0 = 55927d0 ! 1st jan 2012 
    hour = 1.d0/24.d0
    do i=0,n-1
       mjd = mjd0+i*hour
       euler_matrix = rot_equ2hor(mjd)
       hor_pos = matmul(euler_matrix, position)
       call vec2ang(hor_pos, theta, phi)
       write(unit,'(F12.2, 2(F10.2))') mjd, 90.d0-theta*RAD2DEG, phi*RAD2DEG
    end do
    
    close(unit)
    
  end subroutine write_elevation

  subroutine sunmoon(ndays)
    use quiet_ephem_mod
    use quiet_pointing_mod
    implicit none

    character(len=512)       :: file_sun, file_moon
    integer(i4b)             :: ndays, unit1, unit2, n, i
    real(dp)                 :: mjd, hour, ang_sun(2), xyz_sun(3), ang_moon(2), xyz_moon(3)
    real(dp), dimension(3,3) :: euler_matrix

    n = 24*ndays

    file_sun = 'sun_el_'//trim(adjustl(itoa(ndays)))//'.dat'
    file_moon = 'moon_el_'//trim(adjustl(itoa(ndays)))//'.dat'

    unit1 = getlun()
    open(unit1, file=file_sun)
    unit2 = getlun()
    open(unit2, file=file_moon)
    write(unit1,'(12A, 10A, 10A)') 'MJD', 'El', 'Az'
    write(unit2,'(12A, 10A, 10A)') 'MJD', 'El', 'Az'

    mjd = 55927d0 ! 1st jan 2012 00:00 
    hour = 1.d0/24.d0
    do i=0,n-1
       mjd = mjd+i*hour
       ang_sun = ephem(EPHEM_SUN, mjd)
       ang_moon = ephem(EPHEM_MOON, mjd)
       call ang2vec(ang_sun(2), ang_sun(1), xyz_sun) ! ang_sun = (phi, theta) 
       call ang2vec(ang_moon(2), ang_moon(1), xyz_moon)
       euler_matrix = rot_equ2hor(mjd)
       xyz_sun = matmul(euler_matrix, xyz_sun)
       xyz_moon = matmul(euler_matrix, xyz_moon)
       call vec2ang(xyz_sun, ang_sun(2), ang_sun(1)) 
       call vec2ang(xyz_moon, ang_moon(2), ang_moon(1)) 

       write(unit1,'(F12.2, 2(F10.2))') mjd, 90.d0 - ang_sun(2)*RAD2DEG, ang_sun(1)*RAD2DEG
       write(unit2,'(F12.2, 2(F10.2))') mjd, 90.d0 - ang_moon(2)*RAD2DEG, ang_moon(1)*RAD2DEG
    end do

    close(unit1)
    close(unit2)

  end subroutine sunmoon

end module plot_patch_mod

