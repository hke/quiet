! small programs that are useful

program convert_WMAP9_noise
!  use quiet_pixspace_mod
  use quiet_fileutils
  use rngmod
  use quiet_utils
  implicit none
  character(len=512) :: command


  character(len=512), dimension(5)            :: infiles, outmaps, outQU, outTT
  real(dp),     dimension(:,:), allocatable   :: map, map2
  real(dp),     dimension(:,:,:), allocatable :: QUmap, TTmap
  real(dp)    , dimension(5)                  :: sigma0QU, sigma0I  
  integer(i4b)       :: nside, order, npix, i, n, j
  real (dp)           :: nullval, scale, matrix(2,2)
  logical (lgt)       :: anynull

  call getarg(1, command)
!  if (command=="beam") call command_beam

  infiles(1) = "wmap_band_iqumap_r9_9yr_K_v5.fits"
  infiles(2) = "wmap_band_iqumap_r9_9yr_Ka_v5.fits"
  infiles(3) = "wmap_band_iqumap_r9_9yr_Q_v5.fits"
  infiles(4) = "wmap_band_iqumap_r9_9yr_V_v5.fits"
  infiles(5) = "wmap_band_iqumap_r9_9yr_W_v5.fits"

  outmaps(1) = "wmap9_K_map_n512.fits"
  outmaps(2) = "wmap9_Ka_map_n512.fits"
  outmaps(3) = "wmap9_Q_map_n512.fits"
  outmaps(4) = "wmap9_V_map_n512.fits"
  outmaps(5) = "wmap9_W_map_n512.fits"

  !these hdf-files are variances  (sigma**2)
  outQU(1) = "wmap9_QU_n512_K.hdf"
  outQU(2) = "wmap9_QU_n512_Ka.hdf"
  outQU(3) = "wmap9_QU_n512_Q.hdf"
  outQU(4) = "wmap9_QU_n512_V.hdf"
  outQU(5) = "wmap9_QU_n512_W.hdf"

  !these hdf-files are variances (sigma**2)
  outTT(1) = "wmap9_TT_n512_K.hdf"
  outTT(2) = "wmap9_TT_n512_Ka.hdf"
  outTT(3) = "wmap9_TT_n512_Q.hdf"
  outTT(4) = "wmap9_TT_n512_V.hdf"
  outTT(5) = "wmap9_TT_n512_W.hdf"

  sigma0QU(1) = 1.435
  sigma0QU(2) = 1.472 	
  sigma0QU(3) = 2.197
  sigma0QU(4) = 3.141
  sigma0QU(5) = 6.560

  sigma0I(1) = 1.429
  sigma0I(2) = 1.466
  sigma0I(3) = 2.188 	
  sigma0I(4) = 3.131
  sigma0I(5) = 6.544

  !in the wmap_band_iqumap_r9_9yr_K_v5.fits files, there are two N_obs
  ! map(:,4)  has min and max 3.8e5 - 1.8e7
  ! map2(:,1) has min and max 3.8e2 - 1.8e4
  !use map(:,4) as N-obs for pol, the other for temp   ???

  write(*,*), trim(infiles(3))


  do i=1,size(infiles)
     npix = getsize_fits(trim(infiles(i)), nside=nside, ordering=order)

     !read the maps and put into fitsfiles
     allocate(map(0:npix-1,4))
     call read_bintab( trim(infiles(i)), map, npix, 4, nullval, anynull , extno=0)
     map(:,1:3) = map(:,1:3)*1000  ! we have maps in microKelvin, lambda maps are in mK
     write(*,*) map(34,1)
     call write_map(map(:,1:3), order, outmaps(i))

     !read the QQ, QU and UU in extension 1
     allocate(QUmap(0:npix-1,2,2),map2(0:npix-1,4))
     call read_bintab( trim(infiles(i)), map2, npix, 4, nullval, anynull , extno=1)
     scale = (sigma0QU(i)*1e3)**2
     
     do j=1,npix
        matrix(1,1)=map2(j,2)
        matrix(1,2)=map2(j,3)
        matrix(2,1)=matrix(1,2)
        matrix(2,2)=map2(j,4)
        call invert_matrix(matrix)

        QUmap(j,1,1) = scale*matrix(1,1)
        QUmap(j,1,2) = scale*matrix(1,2)
        QUmap(j,2,1) = QUmap(j,1,2)
        QUmap(j,2,2) = scale*matrix(2,2)
     end do

!     QUmap(:,1,1) = scale/map2(:,2)
!     QUmap(:,1,2) = scale/map2(:,3)
!     QUmap(:,2,1) = QUmap(:,1,2)
!     QUmap(:,2,2) = scale/map2(:,4)

     write(*,*) QUmap(35,2,1)
     call write_map(QUmap, order, outQU(i))


     !read the TT 
     allocate(TTmap(0:npix-1,1,1))
     scale = (sigma0I(i)*1e3)**2
     TTmap(:,1,1) = scale/map(:,4)
     call write_map(TTmap, order, outTT(i))

     deallocate(map, map2, QUmap, TTmap)
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


end program convert_WMAP9_noise
