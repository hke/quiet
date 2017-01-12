! small programs that are useful

program convert_Planck_noise
!  use quiet_pixspace_mod
  use quiet_fileutils
  use rngmod
  use quiet_utils
  implicit none
  character(len=512) :: command


  character(len=512), dimension(5)            :: infiles, outmaps, outinvQU, outinvTT
  real(dp),     dimension(:,:), allocatable   :: map, map2
  real(dp),     dimension(:,:,:), allocatable :: invQUmap, invTTmap
  real(dp)    , dimension(5)                  :: sigma0QU, sigma0I  
  integer(i4b), dimension(:), allocatable        :: pixels
  integer(i4b)       :: nside, ordering, npix, i, n, j
  real (dp)           :: nullval, scale, matrix(2,2)
  logical (lgt)       :: anynull

  call getarg(1, command)
!  if (command=="beam") call command_beam

  infiles(1) = "LFI_SkyMap_044_1024_DX11D_full_nomp.fits"

  outmaps(1) = "Planck_44_map_n1024.fits"

  !these hdf-files are variances  (sigma**2)
  outinvQU(1) = "Planck_invQU_n1024_44.hdf"

  !these hdf-files are variances (sigma**2)
  outinvTT(1) = "Planck_invTT_n1024_44.hdf"


  write(*,*), "reading ", trim(infiles(1))

  do i=1,size(infiles)

     !read the maps and put into fitsfiles
     call read_map(map, pixels, nside, ordering, trim(infiles(i)))
     npix = size(map,1)

!     allocate(map(0:npix-1,10))
!     call read_bintab( trim(infiles(i)), map, npix, 10, nullval, anynull , extno=0)
     ! we have maps in microKelvin, Planck maps are in K
     map(:,1:3) = map(:,1:3)*1000*1000
     
     write(*,*) map(34,1)
     call write_map(map(:,1:3), pixels,nside, ordering, outmaps(i))

     !read the QQ, QU and UU in extension 1
     allocate(invQUmap(npix,2,2))
!     call read_bintab( trim(infiles(i)), map2, npix, 4, nullval, anynull , extno=1)
!     scale = (sigma0QU(i)*1e3)**2

     map(:,8:10)=map(:,8:10)*1e12 !from K to muK
     write(*,*) "QQ ", map(35,8)
     do j=1,npix
        matrix(1,1)=map(j,8)
        matrix(1,2)=map(j,9)
        matrix(2,1)=matrix(1,2)
        matrix(2,2)=map(j,10)
        call invert_matrix(matrix)

        invQUmap(j,1,1) = matrix(1,1)
        invQUmap(j,1,2) = matrix(1,2)
        invQUmap(j,2,1) = invQUmap(j,1,2)
        invQUmap(j,2,2) = matrix(2,2)
     end do

     write(*,*) "invQQ ", invQUmap(35,1,1)
     call write_map(invQUmap, ordering, trim(outinvQU(i)))


     !read the TT 
!!$     allocate(TTmap(0:npix-1,1,1))
!!$     scale = (sigma0I(i)*1e3)**2
!!$     TTmap(:,1,1) = scale/map(:,4)
!!$     call write_map(TTmap, order, outTT(i))

     deallocate(map, invQUmap, invTTmap)
  end do


contains


  subroutine command_readwritemap
    implicit none

    character(len=512) :: infile, outfile
    integer(i4b)       :: nside, order, npix
    real(dp),     dimension(:,:), allocatable :: map
    integer(i4b), dimension(:),   allocatable :: pixels
    real (dp)           :: nullval
    logical (lgt)       :: anynull

!    npix = getsize_fits("sim_K23_30.fits", nside=nside, ordering=ordering)
!    allocate(map(0:npix-1,1:nmap))
!  call read_bintab( trim(wmapfile), mapwmap, npix, nmap, nullval, anynull )

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)

    npix = getsize_fits(trim(infile), nside=nside, ordering=order)
    allocate(map(0:npix-1,1:3))
    call read_bintab( trim(infile), map, npix, 3, nullval, anynull )
    pixels = irange(npix)-1
!    call read_map(map, order, infile)
    call write_map(map, order, outfile)
    deallocate(map)


  end subroutine command_readwritemap


  subroutine command_fitsmake1map
    implicit none

    character(len=512) :: infile, outfile
    integer(i4b)       :: nside, order
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: map

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)

    call read_map(map, pixels, nside, order, infile)
    call write_map(map(:,1), pixels, nside, order, outfile)
    deallocate(map, pixels)
  end subroutine


end program convert_Planck_noise
