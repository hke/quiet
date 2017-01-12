MODULE DS_FITSTOOLS

USE DS_PRECISION
USE DS_TYPES
USE DS_UTILS

IMPLICIT NONE

PRIVATE
PUBLIC ::  ds_write_fits_map, ds_write_planar_map

TYPE fits_header
	INTEGER :: nrecords
	INTEGER,ALLOCATABLE,DIMENSION(:) :: record_int
	REAL(DP),ALLOCATABLE,DIMENSION(:) :: record_dp
	CHARACTER(len=80),ALLOCATABLE,DIMENSION(:) :: name, record_str 
	CHARACTER(len=1),ALLOCATABLE,DIMENSION(:) :: htype !see FITSIO guide (p.17) for conventions: J=integer, D=real(8), S=string
END TYPE fits_header


CONTAINS


subroutine write_hpmap_clover(filename, obs_npix, pixel, tmap, qmap, umap, &
     &                     nside, order, units, coord, extra_cards)
!=======================================================================
! writes a fits file for cut sky data set with 4 columns:
!  PIXEL, T, Q, U
! 
! write_hpmap(filename, obs_npix, pixel, tmap, qmap, umap
!                 nside, order, coord)
!=======================================================================
 use pix_tools

  character(len=*), intent(in) :: filename
  integer(4) :: obs_npix
  integer(4), dimension(obs_npix), intent(in) :: pixel
  real(8), dimension(obs_npix), intent(in) :: tmap, qmap, umap
  integer(4), intent(in), optional :: nside, order
  character(len=*), intent(in), optional :: coord, units
  character(len=80), dimension(:), optional :: extra_cards
  character(len=*), parameter :: routine = 'write_hpmap'

  ! --- healpix/fits related variables ----
  integer(4)     :: ncol, grain
  integer(4)     :: npix_hd, nside_hd, nside_final, npix_final, i
  integer(4)     :: maxpix, minpix
  character(len=1) :: char1, coord_usr

  character(len=20) :: units_usr
  logical :: done_nside, done_order, done_coord

  ! --- cfitsio related variables ----
  integer(4) ::  status, unit, blocksize, bitpix, naxis, naxes(1)
  logical ::  simple, extend
  character(LEN=80) :: svalue, comment

  integer(4), parameter :: MAXDIM = 20 !number of columns in the extension
  integer(4) :: nrows, tfields, varidat
  integer(4) :: frow,  felem, repeat, repeatg
  character(len=20) :: ttype(MAXDIM), tform(MAXDIM), tunit(MAXDIM), extname
  character(len=80) ::  card
  character(len=4) :: srepeat, srepeatg
  !=======================================================================

  ncol = 4
  grain = 1
  units_usr = ' '
  if (present(units)) units_usr = units

  status=0
  unit = 100

  !     create the new empty FITS file
  blocksize=1
  call ftinit(unit,filename,blocksize,status)

  !     -----------------------------------------------------
  !     initialize parameters about the FITS image
  simple=.true.
  bitpix=32     ! integer*4
  naxis=0       ! no image
  naxes(1)=0
  extend=.true. ! there is an extension

  !     ----------------------
  !     primary header
  !     ----------------------
  !     write the required header keywords
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !     write the current date
  call ftpdat(unit,status) ! format ccyy-mm-dd

  !     ----------------------
  !     image : none
  !     ----------------------
  
  !     ----------------------
  !     extension
  !     ----------------------

  !     creates an extension
  call ftcrhd(unit, status)

  !     writes required keywords
  repeat = 1
  nrows    = (obs_npix + repeat - 1)/ repeat ! naxis1
  if (obs_npix < repeat) repeat = 1
  write(srepeat,'(i4)') repeat
  srepeat = adjustl(srepeat)

  repeatg = repeat * grain
  write(srepeatg,'(i4)') repeatg
  srepeatg = adjustl(srepeatg)

  tfields  = ncol
  ttype(1:4) = (/ 'PIXEL ','T_MAP ','Q_MAP ','U_MAP ' /)
  tform(1) = trim(srepeat)//'J'
  tform(2) = trim(srepeatg)//'E'
  tform(3) = trim(srepeatg)//'E'
  tform(4) = trim(srepeatg)//'E'

  tunit =  ' '      ! optional, will not appear
  tunit(2) = units_usr
  tunit(3) = units_usr
  tunit(4) = units_usr
  extname  = 'DESCART_MAP'  ! default, will be overide by user provided one if any
  varidat  = 0
  

  call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
       &     extname, varidat, status)

  call ftpcom(unit,'------------------------------------------',status)
  call ftpcom(unit,'          Pixelisation Specific Keywords    ',status)
  call ftpcom(unit,'------------------------------------------',status)
  call ftpkys(unit,'PIXTYPE','HEALPIX ',' HEALPIX Pixelisation',status)      
  call ftpkyu(unit,'NSIDE',   ' ',status) ! place holder, will be updated later on
  call ftpkyu(unit,'ORDERING',' NESTED',status)
  call ftpkys(unit,'COORDSYS',' ',' ',status)
  call ftpcom(unit,'  G = Galactic, E = ecliptic, C = celestial = equatorial', status)   
  call ftpcom(unit,'------------------------------------------',status)
  call ftpcom(unit,'          Data Specific Keywords    ',status)
  call ftpcom(unit,'------------------------------------------',status)
  call ftpkys(unit,'INDXSCHM','EXPLICIT',' Indexing : IMPLICIT or EXPLICIT', status)
  call ftpkyj(unit,'GRAIN',  grain,     ' Grain of pixel indexing',status)
  call ftpcom(unit,'GRAIN=0 : no indexing of pixel data (IMPLICIT) ',status)
  call ftpcom(unit,'GRAIN=1 : 1 pixel index -> 1 pixel data (EXPLICIT)',status)
  call ftpcom(unit,'GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)',status)
  call ftpkys(unit,'OBJECT','PARTIAL ',' Sky coverage represented by data',status)     
  call ftpkyj(unit,'OBS_NPIX',obs_npix, ' Number of pixels observed and recorded',status)
	if (present(extra_cards)) then
		do i=1,size(extra_cards)
			call ftprec(unit,extra_cards(i),status) 
		enddo
	endif

  ! add required Healpix keywords (NSIDE, ORDER) if provided by user
  done_order = .false.
  if (present(order)) then

     if (order == 1) then
        call ftukys(unit, 'ORDERING','RING',' Pixel ordering scheme, either RING or NESTED',status)
        done_order = .true.
     elseif (order == 2) then
        call ftukys(unit, 'ORDERING','NESTED',' Pixel ordering scheme, either RING or NESTED ',status)
        done_order = .true.
     else
        print*,'Invalid ORDER given : ',order, ' instead of 1 (RING) or 2 (NESTED)'
     endif
  endif

  done_nside = .false.
  if (present(nside)) then
     if (nside2npix(nside) > 0) then ! valid nside
        call ftukyj(unit,'NSIDE',nside,' Resolution parameter for HEALPIX', status)
        done_nside = .true.
        nside_final = nside
     else
        print*,'Invalid NSIDE given : ',nside
     endif
  endif

  ! add non required Healpix keyword (COORD)
  done_coord = .false.
  if (present(coord)) then
     coord_usr = adjustl(coord)
     char1 = coord_usr(1:1)
     call ftupch(char1) ! uppercase
     if (char1 == 'C' .or. char1 == 'Q') then
        coord_usr = 'C'
        done_coord = .true.
     elseif (char1 == 'G') then
        coord_usr = 'G'
        done_coord = .true.
     elseif (char1 == 'E' ) then
        coord_usr='E'
        done_coord = .true.
     else
        print*,'Unrecognised COORD given : ',coord,' instead of C, G, or E'
        print*,'Proceed at you own risks '
        coord_usr = char1
        done_coord = .true.
     endif
     if (done_coord) then
        call ftukys(unit, 'COORDSYS',coord_usr,' Pixelisation coordinate system',status)
     endif
  endif

  !     write the extension one column by one column
  frow   = 1  ! starting position (row)
  felem  = 1  ! starting position (element)
  call ftpclj(unit, 1, frow, felem, obs_npix, pixel , status)
  call ftpcld(unit, 2, frow, felem, obs_npix, tmap, status)
  if (tfields >= 3) call ftpcld(unit, 3, frow, felem, obs_npix, qmap, status)
  if (tfields >= 4) call ftpcld(unit, 4, frow, felem, obs_npix, umap, status)

  !     ----------------------
  !     close and exit
  !     ----------------------

  !     close the file and free the unit number
  call ftclos(unit, status)

  !     check for any error, and if so print out error messages
  if (status > 0) write (6,*) "Problem with healpix map output:",filename, status

  return
end subroutine write_hpmap_clover


! Cookbook routine to delete file

subroutine deletefile(filename,status)

integer(4), intent(out) :: status
character(len=80), intent(in) :: filename
integer(4) :: unit,blocksize


! Simply return if status is greater than zero
if (status .gt. 0)return

! Get an unused Logical Unit Number to use to open the FITS file
call ftgiou(unit,status)

! Try to open the file, to see if it exists
call ftopen(unit,filename,1,blocksize,status)

if (status .eq. 0)then!
! file was opened;  so now delete it 
   call ftdelt(unit,status)
else if (status .eq. 103)then
! file doesn't exist, so just reset status to zero and clear errors
   status=0
   call ftcmsg
else
! there was some other error opening the file; delete the file anyway
   status=0
   call ftcmsg
   call ftdelt(unit,status)
end if

! Free the unit number for later reuse
call ftfiou(unit, status)

end subroutine


subroutine ds_write_planar_map(filename,maps,xpix,ypix,do_T, do_P)
	character(*) :: filename
	type(ds_trimap) :: maps
	logical :: do_T, do_P
	integer :: xpix,ypix
	integer :: unit
	integer, parameter :: blocksize = 1
	integer status
	integer pixel
	integer nelements
	status=0
	call deletefile(filename,status)
	
	
	!Get a logical unit number for the file.
	call FTGIOU(unit, status)
	
	!Create the FITS file
	call FTINIT(unit,filename,blocksize,status)
	
	if (do_T) then
		call ds_add_hdu_with_map(unit,xpix,ypix,maps%T,status)
	endif
	if (do_P) then
		call ds_add_hdu_with_map(unit,xpix,ypix,maps%Q,status)
		call ds_add_hdu_with_map(unit,xpix,ypix,maps%U,status)
	endif
	
	
	!Close the FITS file
	call FTCLOS(unit,status)
	
	!Free the logical unit number
	call FTFIOU(unit, status)


end subroutine ds_write_planar_map

subroutine ds_add_hdu_with_map(unit,xpix,ypix,map,status)
	integer :: unit, xpix, ypix
	integer status
	type(ds_map) :: map
	integer :: npix, index
	integer, parameter :: naxis = 2
	integer, parameter :: bitpix = -64
	integer, parameter :: first_pixel = 1
	integer, parameter :: group = 0
	integer, dimension(2) :: naxes
	real(dp), dimension(:), allocatable :: values
	integer p
	

	npix = xpix*ypix
	naxes(1) = xpix
	naxes(2) = ypix

	!Create the primary image HDU
	call FTIIMG(unit,bitpix,naxis,naxes,status)

	!Create the array of values from the map
	allocate(values(npix))
	values=0.0
	do p=1,map%npix
		index = map%indices(p)
		values(index) = map%map(p)
	enddo

	!Write the data into the image
	call FTPPRD(unit,group,first_pixel,npix,values,status)
	deallocate(values)
	
end subroutine ds_add_hdu_with_map

subroutine ds_write_fits_map(filename,maps,nside,partial_out,units,isRing,extra_headers)
use fitstools, only : output_map
use head_fits, only : write_minimal_header
type(ds_trimap) :: maps
character(len=80), dimension(1:), optional :: extra_headers

integer nside
logical :: partial_out
real(dp), dimension(:), pointer :: dummy_map
character(*) :: filename
integer card_pos
logical write_card
integer status, i
logical, optional :: isRing
integer(4) order_code
integer(4), dimension(:), pointer :: indices
character(8), optional :: units
character(8) :: unitsCode
integer :: nheal
!integer,parameter :: map_type=2
real(sp),allocatable,dimension(:,:) :: map
character(len=80),dimension(1:300) :: header
real(dp), parameter :: hpx_dbadval = -1.6375e30_dp
character(len=20) :: dtype, coordsys
character(len=8) :: ordering
integer npix
integer t
call ds_assert(maps%has_T .or. maps%has_P, "Maps with neither T nor P passed to write_fits_maps")

if (maps%has_t) then
	npix = maps%T%npix
else
	npix = maps%Q%npix	
endif

if (.not. (maps%has_T .and. maps%has_P)) then
	allocate(dummy_map(npix))
	dummy_map = 0
endif



status=0
call deletefile(filename,status)

if (.not. present(isRing)) then
   order_code = 2
	ordering='NESTED'
elseif (isRing) then
   order_code = 1
ordering='RING'
else
	ordering='NESTED'
   order_code=2
endif

if (.not. present(units)) then
   unitsCode = 'unknown '
else
   unitsCode= units
endif

if(partial_out) then
	if (present(extra_headers)) then
		if (maps%has_T .and. maps%has_P) then
			call write_hpmap_clover(filename, npix, maps%T%indices, maps%T%map, maps%Q%map, maps%U%map, nside,order=order_code,units=unitsCode,extra_cards=extra_headers)
		else if (maps%has_T) then
			call write_hpmap_clover(filename, npix, maps%T%indices, maps%T%map, dummy_map,dummy_map, nside,order=order_code,units=unitsCode,extra_cards=extra_headers)
		else  !Polarization only
			call write_hpmap_clover(filename, npix, maps%Q%indices, dummy_map, maps%Q%map, maps%U%map, nside,order=order_code,units=unitsCode,extra_cards=extra_headers)
		endif
	else
		if (maps%has_T .and. maps%has_P) then
			call write_hpmap_clover(filename, npix, maps%T%indices, maps%T%map, maps%Q%map, maps%U%map, nside,order=order_code,units=unitsCode)
		else if (maps%has_T) then
			call write_hpmap_clover(filename, npix, maps%T%indices, maps%T%map, dummy_map,dummy_map, nside,order=order_code,units=unitsCode)
		else  !Polarization only
			call write_hpmap_clover(filename, npix, maps%Q%indices, dummy_map, maps%Q%map, maps%U%map, nside,order=order_code,units=unitsCode)
		endif
	endif
else
	
	nheal= nside**2 *12
	allocate(map(0:nheal-1,1:3))
	map= hpx_dbadval
	if (maps%has_T) map(maps%T%indices,1)= maps%T%map
	if (maps%has_P) map(maps%Q%indices,2)= maps%Q%map
	if (maps%has_P) map(maps%U%indices,3)= maps%U%map
	dtype= 'MAP'
	
	call write_minimal_header( header, dtype, nside=nside, &
		ordering=ordering, coordsys='E',polar=.true.,units=unitsCode)

	if (present(extra_headers)) then
		
		card_pos = size(header)+1
		write_card=.true.
		do
			card_pos=card_pos-1
			if (header(card_pos) .ne. '') exit
			if (card_pos==0) then
				write(*,*) "FAILED TO WRITE ANY EXTRA HEADERS TO OUTPUT"
				write_card=.false.
			endif
		enddo
		header(card_pos) = "COMMENT  "
		card_pos = card_pos + 1
		header(card_pos) = "COMMENT  List of cuts applied follows"
		card_pos = card_pos + 1
		header(card_pos) = "COMMENT  "
		card_pos = card_pos + 1
		if (write_card) then
			do i=1,size(extra_headers)
				if (card_pos>size(header)) then
					write(*,*) "FAILED TO WRITE SOME EXTRA HEADERS TO OUTPUT"
					exit
				endif
				header(card_pos) = extra_headers(i)
				card_pos = card_pos + 1
			enddo
		endif
	endif
		
	call output_map(map,header,filename)

	deallocate(map)
	if (.not. (maps%has_T .and. maps%has_P)) deallocate(dummy_map)	
endif

end subroutine ds_write_fits_map






END MODULE DS_FITSTOOLS

