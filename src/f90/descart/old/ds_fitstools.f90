MODULE DS_FITSTOOLS

USE DS_PRECISION
USE DS_TYPES
USE DS_UTILS

IMPLICIT NONE

PRIVATE
PUBLIC :: READ_POINTING_HEADER, FAST_FITS_READ, fits_header, savePolMaps, read_fits_header, readDetPointing, read_simple_tod, read_ntod, write_simple_tod, read_simple_tod_1col

TYPE fits_header
	INTEGER :: nrecords
	INTEGER,ALLOCATABLE,DIMENSION(:) :: record_int
	REAL(DP),ALLOCATABLE,DIMENSION(:) :: record_dp
	CHARACTER(len=80),ALLOCATABLE,DIMENSION(:) :: name, record_str 
	CHARACTER(len=1),ALLOCATABLE,DIMENSION(:) :: htype !see FITSIO guide (p.17) for conventions: J=integer, D=real(8), S=string
END TYPE fits_header



CONTAINS



SUBROUTINE READ_POINTING_HEADER(filename,naxis1,naxis2,ob_num,tot_obs,datarate,obs_el,hwp_vel,bore_step,bore_ang,longitude,latitude)

  TYPE(fits_header) :: pointhead
  CHARACTER(len=*),INTENT(IN) :: filename
  INTEGER,INTENT(OUT) :: ob_num, tot_obs,datarate,bore_step,naxis1,naxis2
  REAL(DP),INTENT(OUT) :: obs_el, hwp_vel, bore_ang,longitude,latitude
  

  pointhead%nrecords= 11!34 total number in header
  ALLOCATE(pointhead%name(pointhead%nrecords))
  ALLOCATE(pointhead%htype(pointhead%nrecords))
  ALLOCATE(pointhead%record_dp(pointhead%nrecords))
  ALLOCATE(pointhead%record_int(pointhead%nrecords))
  
  pointhead%name(1)= 'OB_NUM'
  pointhead%htype(1)= 'J'
  pointhead%name(2)= 'TOT_OBS'
  pointhead%htype(2)= 'J'
  pointhead%name(3)= 'DATARATE'
  pointhead%htype(3)= 'J'
  pointhead%name(4)= 'OBS_EL'
  pointhead%htype(4)= 'D'
  pointhead%name(5)= 'HWP_VEL'
  pointhead%htype(5)= 'D'
  pointhead%name(6)= 'BORE_STP'!ask CEN about this
  pointhead%htype(6)= 'J'
  pointhead%name(7)= 'BORE_ANG'
  pointhead%htype(7)= 'D'
  pointhead%name(8)= 'NAXIS1'
  pointhead%htype(8)= 'J'
  pointhead%name(9)= 'NAXIS2'
  pointhead%htype(9)= 'J'
  pointhead%name(10)= 'SITE_LON'
  pointhead%htype(10)= 'D'
  pointhead%name(11)= 'SITE_LAT'
  pointhead%htype(11)= 'D'
  
  
  CALL READ_FITS_HEADER(filename,pointhead)

  ob_num= pointhead%record_int(1)
  tot_obs= pointhead%record_int(2)
  datarate= pointhead%record_int(3)
  obs_el= pointhead%record_dp(4)
  hwp_vel= pointhead%record_dp(5)
  bore_step= pointhead%record_int(6)
  bore_ang= pointhead%record_dp(7)
  naxis1=pointhead%record_int(8)
  naxis2=pointhead%record_int(9)
  longitude=pointhead%record_dp(10)
  latitude=pointhead%record_dp(11)

  DEALLOCATE(pointhead%name)
  DEALLOCATE(pointhead%htype)
  DEALLOCATE(pointhead%record_dp)
  DEALLOCATE(pointhead%record_int)

END SUBROUTINE READ_POINTING_HEADER




SUBROUTINE FAST_FITS_READ(filename, naxis1, naxis2, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10)
 
  CHARACTER(len=*),INTENT(IN) :: filename
  INTEGER,INTENT(IN) :: naxis1, naxis2
  REAL(DP),ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: col1
  REAL(DP),ALLOCATABLE,DIMENSION(:),INTENT(OUT),OPTIONAL :: col2, col3, col4, col5, col6, col7, col8, col9, col10
  !nb: naxis must be the number of columns in the table, regardless of the number of columns you read.

  INTEGER :: status, unit, nrows, rwmode, blocksize, icall, tot_calls, group, nullval, fpixel, nelements, first, last
  LOGICAL :: anyf
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: values

  status=0
  CALL FTGIOU(unit,status)
  CALL ASSERT_STATUS(status)

  rwmode= 0 !read-only
  CALL FTOPEN(unit,trim(filename),rwmode,blocksize,status)
  CALL ASSERT_STATUS(status)

  !determine the optimal number of rows
!  CALL FTGRSZ(unit,nrows,status)
!  CALL ASSERT_STATUS(status)
  nrows= naxis2

  tot_calls= ceiling(real(naxis2)/real(nrows))

  group=0
  nullval=0
  fpixel=1
  nelements= nrows*naxis1

  ALLOCATE(values(naxis1,nrows)) ! (ncolumns,nrows)

  ALLOCATE(col1(naxis2))
  IF(PRESENT(col2)) ALLOCATE(col2(naxis2))
  IF(PRESENT(col3)) ALLOCATE(col3(naxis2))
  IF(PRESENT(col4)) ALLOCATE(col4(naxis2))
  IF(PRESENT(col5)) ALLOCATE(col5(naxis2))
  IF(PRESENT(col6)) ALLOCATE(col6(naxis2))
  IF(PRESENT(col7)) ALLOCATE(col7(naxis2))
  IF(PRESENT(col8)) ALLOCATE(col8(naxis2))
  IF(PRESENT(col9)) ALLOCATE(col9(naxis2))
  IF(PRESENT(col10)) ALLOCATE(col10(naxis2))


  last=0
  DO icall= 1, tot_calls     
     if(icall==tot_calls) nelements = (naxis2-last) * naxis1
     CALL FTGPVD(unit, group,fpixel,nelements,nullval,values,anyf,status)
     CALL ASSERT_STATUS(status)
     fpixel= fpixel + nelements

     first= last + 1
     last= min(icall * nrows , naxis2)

     col1(first:last) = values(1,:)
     IF(PRESENT(col2)) col2(first:last) = values(2,:)
     IF(PRESENT(col3)) col3(first:last) = values(3,:)
     IF(PRESENT(col4)) col4(first:last) = values(4,:)
     IF(PRESENT(col5)) col5(first:last) = values(5,:)
     IF(PRESENT(col6)) col6(first:last) = values(6,:)
     IF(PRESENT(col7)) col7(first:last) = values(7,:)
     IF(PRESENT(col8)) col8(first:last) = values(8,:)
     IF(PRESENT(col9)) col9(first:last) = values(9,:)
     IF(PRESENT(col10)) col10(first:last) = values(10,:)

  ENDDO

  DEALLOCATE(values)

  CALL FTCLOS(unit,status)
  CALL ASSERT_STATUS(status)

  CALL FTFIOU(unit, status)
  CALL ASSERT_STATUS(status)

END SUBROUTINE FAST_FITS_READ


SUBROUTINE READ_FITS_HEADER(filename,fitshead)
   
  CHARACTER(len=*),INTENT(IN) :: filename
  TYPE(fits_header),INTENT(INOUT) :: fitshead
	
  INTEGER :: unit, status
  INTEGER :: blocksize, rwmode
  INTEGER :: irecord,vali
  CHARACTER(len=80) :: comment, valstr
  REAL(DP) :: vald

  status=0
  CALL FTGIOU(unit,status)
  CALL ASSERT_STATUS(status)

  rwmode= 0 !read-only
  CALL FTOPEN(unit,trim(filename),rwmode,blocksize,status)
  CALL ASSERT_STATUS(status)

  DO irecord=1, fitshead%nrecords
     SELECT CASE(fitshead%htype(irecord))
        CASE('J')
           CALL FTGKYJ(unit,trim(fitshead%name(irecord)),vali,comment,status)
           CALL ASSERT_STATUS(status)
           fitshead%record_int(irecord)= vali
        CASE('D')
           CALL FTGKYD(unit,trim(fitshead%name(irecord)),vald,comment,status)
           CALL ASSERT_STATUS(status)
           fitshead%record_dp(irecord)= vald
        CASE('S')	
        	CALL FTGKYS(unit,trim(fitshead%name(irecord)),valstr,comment,status)
        	CALL ASSERT_STATUS(status)
        	fitshead%record_str(irecord)= valstr
     END SELECT
  END DO

  CALL FTCLOS(unit,status)
  CALL ASSERT_STATUS(status)

  CALL FTFIOU(unit, status)
  CALL ASSERT_STATUS(status)
	

END SUBROUTINE READ_FITS_HEADER





SUBROUTINE ASSERT_STATUS(status)

  INTEGER,INTENT(IN):: status
  CHARACTER(len=30) :: errtext

  IF(status .gt. 0) THEN
     PRINT*,'FITS ERROR',status
     CALL FTGERR(status,errtext)
     PRINT*,errtext
     STOP
  ENDIF
  
END SUBROUTINE ASSERT_STATUS









subroutine write_hpmap_clover(filename, obs_npix, pixel, tmap, qmap, umap, &
     &                     nside, order, units, coord)
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
  ttype(1:4) = (/ 'PIXEL ','T_MAP ','Q_MAP ','U_MAP' /)
  tform(1) = trim(srepeat)//'J'
  tform(2) = trim(srepeatg)//'E'
  tform(3) = trim(srepeatg)//'E'
  tform(4) = trim(srepeatg)//'E'

  tunit =  ' '      ! optional, will not appear
  tunit(2) = units_usr
  tunit(3) = units_usr
  tunit(4) = units_usr
  extname  = 'CLOVER_SIMULATION'  ! default, will be overide by user provided one if any
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



subroutine savePolMaps(filename,qmap,umap,nside)
type(ds_map) :: qmap, umap
integer nside
real(dp), dimension(:), allocatable :: tmap
character(*) :: filename
integer status

call ds_assert(qmap%npix == umap%npix,'Error in save pol maps')
allocate(tmap(qmap%npix))
tmap = 0.0

status=0
call deletefile(filename,status)

call write_hpmap_clover(filename, qmap%npix, qmap%indices, tmap, qmap%map, umap%map, nside)


end subroutine savePolMaps



subroutine readDetPointing(filename,nside,pointing,datarate)
	use pix_tools
	integer :: status, unit, blocksize, ncol
	type(ds_detpointing) :: pointing
	real(dp), allocatable, dimension(:) :: ra, dec
	character(*), intent(in)  :: filename
	real(dp) :: theta_i
	real(dp) :: datarate
	integer :: nside, ntod,i,p
	logical :: bad
	integer, parameter :: read_only = 0
	character(80) :: comment
	character(8) :: keyword
	integer, parameter :: nfield = 3

	status=0

	call ftgiou(unit,status)
	call ds_assert(status==0, "Could not get a logical unit number to open FITS file")
	
	call ftdopn(unit,filename,read_only,status)
	call ds_assert(status==0, "Could not open FITS file called " // filename)
	
	call ftgkys(unit,"FIELD1  ", keyword,comment,status)
	call ds_assert((status==0) .and. (keyword == "RA"), "Keyword FIELD1 could not be loaded or was not 'RA'; wrong format. It was:" // keyword // " in file " // filename)

	call ftgkys(unit,"FIELD2  ", keyword,comment,status)
	call ds_assert((status==0) .and. (keyword == "DEC"), "Keyword FIELD2 could not be loaded or was not 'DEC'; wrong format. It was:" // keyword // " in file " // filename)

	call ftgkys(unit,"FIELD3  ", keyword,comment,status)
	call ds_assert((status==0) .and. (keyword == "THETA"), "Keyword FIELD3 could not be loaded or was not 'THETA'; wrong format. It was:" // keyword // " in file " // filename)

	call ftgkyj(unit,"NAXIS2", ntod,comment,status)
	call ds_assert(status==0, "Could not get the keyword NAXIS2 in FITS file " // filename // " - WTF?")

	call ftgkyd(unit,"DATARATE", datarate,comment,status)
	call ds_assert(status==0, "Could not get the keyword NAXIS2 in FITS file " // filename // " - WTF?")

	
	call ftgncl(unit, ncol, status)
	call ds_assert(status==0, "Could not get number of columns FITS file " // filename // " - WTF?")
	call ds_assert(ncol .ge. nfield, "Insufficient columns in FITS file "// filename)

	allocate(ra(ntod))
	allocate(dec(ntod))
	call preparePointing(pointing, ntod, 0, .false.)  !The last argument tells it not to build the globalhitmap - that happens after re-pixelization.	
	
	call ftgcvd(unit,1,1,1,ntod,0.0D0,ra,bad,status)
	call ds_assert(status==0, "Could not read RA column from file " // filename)

	call ftgcvd(unit,2,1,1,ntod,0.0D0,dec,bad,status)
	call ds_assert(status==0, "Could not read DEC column from file " // filename)
	
	call ftgcvd(unit,3,1,1,ntod,0.0D0,pointing%theta,bad,status)
	call ds_assert(status==0, "Could not read THETA column from file " // filename)

	call ftclos(unit, status)
	call ds_assert(status==0, "Could not close file " // filename)
	
	do i=1,ntod
		theta_i = halfpi - dec(i)
		call ang2pix_ring(nside, theta_i, ra(i), p)
		pointing%pixel(i) = p
!		write(*,*) trim(filename), i, p
	enddo
	
	deallocate(ra)
	deallocate(dec)

	call FTFIOU(unit, status)
	call ds_assert(status==0, "Could not free a logical unit number to open FITS file")


end subroutine readDetPointing

subroutine read_simple_tod(filename, timestream,noise_only)
	character(*) :: filename
	type(ds_timestream) :: timestream
	logical, optional :: noise_only
	integer :: status, unit, blocksize, column, ntod
	character(80) :: comment
	character(8) :: keyword
	logical :: bad
	real(dp),allocatable,dimension(:,:) :: table
	integer, parameter :: read_only = 0

	if (.not. present(noise_only)) then
		column = 1
	else
		if (noise_only) then
			column = 2
		else 
			column = 1
		endif
	endif
	
	status=0

	call ftgiou(unit,status)
	call ds_assert(status==0, "Could not get a logical unit number to open FITS file")
	call ftdopn(unit,filename,read_only,status)
	call ds_assert(status==0, "Could not open FITS file called " // filename)
	call ftgkys(unit,"FIELD1  ", keyword,comment,status)
	call ds_assert((status==0) .and. (keyword == "S+N"), "Keyword FIELD1 could not be loaded or was not 'S+N'; wrong format. It was:" // keyword // " in file " // filename)

	call ftgkys(unit,"FIELD2  ", keyword,comment,status)
	call ds_assert((status==0) .and. (keyword == "N"), "Keyword FIELD1 could not be loaded or was not 'N'; wrong format. It was:" // keyword // " in file " // filename)

	call ftgkyj(unit,"NAXIS2", ntod,comment,status)
	call ds_assert(status==0, "Could not get the keyword NAXIS2 in FITS file " // filename // " - WTF?")

	call prepareTimestream(timestream,ntod)

	allocate(table(2,ntod)) !more generally table(naxis1,naxis2)
	call ftgpvd(unit,0,1,2*ntod,0,table,bad,status)
	call ds_assert(status==0,"Could not read data from file"//filename)	
	timestream%timestream = table(column,:)
	deallocate(table)
	
!	call ftgcvd(unit,column,1,1,ntod,0.0D0,timestream%timestream,bad,status)
!	call ds_assert(status==0, "Could not read data column " // trim(inttostr(column)) // " from file " // filename)

	call ftclos(unit, status)
	call ds_assert(status==0, "Could not close file " // filename)

	call FTFIOU(unit, status)
	call ds_assert(status==0, "Could not free a logical unit number to open FITS file")


end subroutine read_simple_tod


subroutine read_simple_tod_1col(filename, timestream)
	character(*) :: filename
	type(ds_timestream) :: timestream
	integer :: status, unit, blocksize, column, ntod
	character(80) :: comment
	character(8) :: keyword
	logical :: bad
	integer, parameter :: read_only = 0
	
	column=1
	status=0

	call ftgiou(unit,status)
	call ds_assert(status==0, "Could not get a logical unit number to open FITS file")
	call ftdopn(unit,filename,read_only,status)
	call ds_assert(status==0, "Could not open FITS file called " // filename)

	call ftgkyj(unit,"NAXIS2", ntod,comment,status)
	call ds_assert(status==0, "Could not get the keyword NAXIS2 in FITS file " // filename // " - WTF?")

	call prepareTimestream(timestream,ntod)

	call ftgpvd(unit,0,1,ntod,0,timestream%timestream,bad,status)
	call ds_assert(status==0,"Could not read data from file"//filename)	

	call ftclos(unit, status)
	call ds_assert(status==0, "Could not close file " // filename)

	call FTFIOU(unit, status)
	call ds_assert(status==0, "Could not free a logical unit number to open FITS file")


end subroutine read_simple_tod_1col



subroutine read_ntod(filename,ntod)
	character(*),intent(in) :: filename
	integer,intent(inout) :: ntod
	integer :: status, unit
	character(80) :: comment
	integer,parameter :: read_only=0
	
	status= 0
	call ftgiou(unit,status)
	call ds_assert(status==0, "Could not get a logical unit number to open FITS file")
	
	call ftdopn(unit,filename,read_only,status)
	call ds_assert(status==0, "Could not open FITS file called " // filename)

	call ftgkyj(unit,"NAXIS2", ntod,comment,status)
	call ds_assert(status==0, "Could not get the keyword NAXIS2 in FITS file " // filename // " - WTF?")

	call ftclos(unit, status)
	call ds_assert(status==0, "Could not close file " // filename)

	call FTFIOU(unit, status)
	call ds_assert(status==0, "Could not free a logical unit number to open FITS file")

end subroutine read_ntod




subroutine write_simple_tod(filename,delta_t,nside,signalNoise,noise)
	character(*) :: filename
	integer :: nside
	real(dp) :: delta_t
	type(ds_timestream) :: signalNoise
	type(ds_timestream), optional :: noise
	integer :: status, unit, column, ntod
	integer, parameter :: blocksize = 1
	character(80) :: comment
	character(8) :: keyword
	logical :: bad
	logical, parameter :: simple=.TRUE., extend=.true.
	integer, parameter :: bitpix=-64 , naxis=2, decimals=6
	integer, dimension(2) :: naxes
	real(dp),allocatable,dimension(:,:) :: table
	integer, parameter :: read_only = 0
	
	naxes(1)=2
	naxes(2)=signalNoise%nt
  
	

	
	status=0

	call ftgiou(unit,status)
	call ds_assert(status==0, "Could not get a logical unit number to open FITS file")

	call ftinit(unit,trim(filename),blockSize,status)
	call ds_assert(status==0, "Could not open FITS file for writing:"//filename)
	
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	call ds_assert(status==0, "Could not write FITS basic header:"//filename)
	
	keyword = 'delta_t'
	comment="Sample time [s]"
	call ftpkyg(unit,keyword,delta_t,decimals,comment,status)

	keyword = 'datarate'
	comment="Sample rate [Hz]"
	call ftpkyg(unit,keyword,1.0_dp / delta_t,decimals,comment,status)

	keyword= 'nside'
	comment="Healpix nside"
	CALL FTPKYJ(unit,keyword, nside,comment,status)
	
	keyword= 'field1'
	comment="signal + noise [uK]"
	CALL FTPKYS(unit,keyword, "S+N",comment,status)
	
	keyword= 'field2'
	comment="noise [uK]"
	CALL FTPKYS(unit,keyword, "N",comment,status)

	call ds_assert(status==0, "Could not write FITS keywords:"//filename)


	allocate(table(2,signalNoise%nt)) !more generally table(naxis1,naxis2)

	table(1,:) = signalNoise%timestream
	if (present(noise)) then
		table(2,:) = noise%timestream
	else
		table(2,:) = 0.0_dp
	endif

  	CALL FTPPRD(unit,1,1,2*signalNoise%nt,table,status)
	call ds_assert(status==0, "Could not write FITS data:"//filename)	

	call ftclos(unit, status)
	call ds_assert(status==0, "Could not close file " // filename)

	call ftfiou(unit,status)
	call ds_assert(status==0, "Could not free unit")
	
	deallocate(table)


end subroutine write_simple_tod




END MODULE DS_FITSTOOLS

