module ds_focalplane
  USE ds_types

  IMPLICIT NONE


  TYPE DETECTOR
     REAL(DP) :: f, x, y, orient! [GHz,mm,mm,deg]
     CHARACTER(len=100) :: name
  END TYPE DETECTOR

TYPE DS_PLANE  !JZ I made this into a type to avoid the global variable.
  TYPE(DETECTOR),ALLOCATABLE,DIMENSION(:) :: det
  INTEGER(4) :: ndet 
  real(dp) :: focal_length
 END TYPE DS_PLANE

!  PUBLIC :: LOAD_FOCAL_PLANE

CONTAINS

  SUBROUTINE LOAD_FOCAL_PLANE(fp_file,plane,focal_length)
    real(8) focal_length
    CHARACTER(len=*),INTENT(IN) :: fp_file
     type(ds_plane) :: plane
    INTEGER(4) :: ncol, nrow
    CHARACTER(100),ALLOCATABLE,DIMENSION(:) :: columns, units, formats
    CHARACTER(100),ALLOCATABLE,DIMENSION(:,:) :: values
    INTEGER(4) :: iname, itype, istat, ifreq, ixpos, iypos, iorient, icol, irow, idet
    
    !read the focal plane info
    CALL READ_FOCAL_PLANE_PARAM(trim(fp_file), ncol, columns, units, formats, nrow, values)

    ! find out which column is type, status, names, x, y, orient etc..
    do icol = 1, ncol
       if (trim(columns(icol)).eq.'CHANNEL_NAME') iname = icol
       if (trim(columns(icol)).eq.'CHANNEL_TYPE') itype = icol
       if (trim(columns(icol)).eq.'CHANNEL_STATUS') istat = icol
       if (trim(columns(icol)).eq.'FREQUENCY') ifreq = icol
       if (trim(columns(icol)).eq.'X_POSITION') ixpos = icol
       if (trim(columns(icol)).eq.'Y_POSITION') iypos = icol
       if (trim(columns(icol)).eq.'DESIGN_ORIENTATION') iorient = icol
    enddo
    !  CHANNEL_TYPE - reject all detectors where this is not 'LIGHT'
    !  CHANNEL_STATUS - reject all detector where this is not 'OK'

    ! find out how many good detectors
    plane%ndet = 0 
    do irow = 1, nrow
       if (trim(values(itype, irow)).eq.'LIGHT' .and. trim(values(istat, irow)).eq.'OK') plane%ndet = plane%ndet + 1 
    enddo

    ! allocate arrays 
    ALLOCATE(plane%det(plane%ndet))
	plane%focal_length = focal_length
    ! fill arrays with working detector info
    idet = 0
    do irow = 1, nrow
       if (trim(values(itype, irow)).eq.'LIGHT' .and. trim(values(istat,irow)).eq.'OK') then
          idet = idet + 1
          plane%det(idet)%f = a2d(values(ifreq, irow))
          plane%det(idet)%x = a2d(values(ixpos, irow))
          plane%det(idet)%y = a2d(values(iypos, irow))
          plane%det(idet)%orient = a2d(values(iorient, irow)) * deg2rad
          plane%det(idet)%name = trim(values(iname,irow))
       endif
    enddo

    ! deallocate temporary arrays
    deallocate(columns)
    deallocate(units)
    deallocate(formats)
    deallocate(values)

  END SUBROUTINE LOAD_FOCAL_PLANE










  !###########################################!
  !# subroutine for reading a .csv file      #!
  !#   from Clover simulation code 03/09/08  #!
  !###########################################!

  subroutine READ_FOCAL_PLANE_PARAM(infile, ncol, columns, units, formats, nbody, values)

    integer(4), parameter :: line_max = 1024
    character(len=*), intent(in) :: infile
    integer(4), intent(out) :: ncol, nbody
    character(len=*), intent(inout), allocatable :: columns(:), formats(:)
    character(len=*), intent(inout), allocatable :: units(:), values(:,:)

    character(200), allocatable :: values_tmp(:,:)
    integer(4) :: hash_addr, comma_addr, old_addr
    integer(4) :: ibody, icol, ihead, nhead, nlines
    character ::  body
    character(line_max) :: line

    ! first, get the total number of lines in the file 
    open(10, file=infile, status='old')
    nlines = 0
    do 
       read (10,'(a)',end=10) line
       nlines = nlines + 1
    enddo
10  close(10)

    ! open the file
    body = 'h'
    ihead = 0
    ibody = 0
    open(10, file=infile, status='old')

    readlines: do ! read each line in turn
       read (10,'(a)',end=20) line
       hash_addr = scan(line,'#')

       ! get rid of leading # if present
       if (hash_addr .gt. 0) then
          line = trim( adjustl(line(hash_addr+1:)) )
       endif

       ! bail if line is empty
       if (len(trim(line)) .le. 2) cycle readlines

       ! see if we are now in the body
       if (line(1:4).eq."BODY") then 
          nhead = ihead
          body = 'c'
          cycle readlines
       endif

       select case(body)

       case('h') ! we're in the header
          ihead = ihead + 1

       case('c') ! parse the column headings
          ! first get the number of columns
          comma_addr = scan(line, ',')
          old_addr = comma_addr
          ncol = 1
          do while (comma_addr .ne. 0)
             comma_addr = scan(line(old_addr+1:), ',')
             if (comma_addr .eq. 0) cycle
             comma_addr = comma_addr + old_addr 
             old_addr = comma_addr
             ncol = ncol + 1
          enddo
          ncol = ncol + 1
          ! now fill the columns array
          allocate(columns(ncol))
          comma_addr = 0
          do icol=1, ncol
             old_addr = comma_addr
             comma_addr = scan(line(old_addr+1:), ',')
             comma_addr = comma_addr + old_addr
             columns(icol) = trim(line(old_addr+1:comma_addr-1))
             if (icol .eq. ncol) columns(icol) = trim(line(old_addr+1:))
          enddo
          body = 'u'

       case('u') ! parse the units
          allocate(units(ncol))
          comma_addr = 0
          do icol=1, ncol
             old_addr = comma_addr
             comma_addr = scan(line(old_addr+1:), ',')
             comma_addr = comma_addr + old_addr
             units(icol) = trim(line(old_addr+1:comma_addr-1))
             if (icol .eq. ncol) units(icol) = trim(line(old_addr+1:))
          enddo
          body = 'f'

       case('f') ! parse the formats
          allocate(formats(ncol))
          comma_addr = 0
          do icol=1, ncol
             old_addr = comma_addr
             comma_addr = scan(line(old_addr+1:), ',')
             comma_addr = comma_addr + old_addr
             formats(icol) = trim(line(old_addr+1:comma_addr-1))
             if (icol .eq. ncol) formats(icol) = trim(line(old_addr+1:))
          enddo
          body = 'b'
          allocate(values_tmp(ncol,nlines))

       case('b') ! parse the main body 
          ibody = ibody + 1
          comma_addr = 0
          do icol=1, ncol
             old_addr = comma_addr
             comma_addr = scan(line(old_addr+1:), ',')
             comma_addr = comma_addr + old_addr
             values_tmp(icol, ibody) = trim(line(old_addr+1:comma_addr-1))
             if (icol .eq. ncol) values_tmp(icol, ibody) = trim(line(old_addr+1:))
          enddo

       end select

    end do readlines

20  close(10)

    nbody = ibody
    allocate(values(icol, nbody))
    do icol=1,ncol
       do ibody=1,nbody
          values(icol,ibody) = values_tmp(icol,ibody)
       enddo
    enddo
    deallocate(values_tmp)

  end subroutine READ_FOCAL_PLANE_PARAM


  FUNCTION a2d(buf) RESULT(res)
    IMPLICIT NONE
    !// Arguments
    CHARACTER(LEN=*)               :: buf
    DOUBLE PRECISION               :: res
    !// Local variables
    INTEGER                        :: i,k,itmp, foffs, lng
    INTEGER                        :: p, q
    DOUBLE PRECISION               :: dtmp
!    INTEGER, EXTERNAL              :: front_trim
!    INTEGER, EXTERNAL              :: a2i
    LOGICAL                        :: neg
    LOGICAL                        :: exponential
    neg = .FALSE.
    exponential = .FALSE.
    lng = LEN_TRIM(buf)
    foffs = front_trim(buf)
    k = INDEX(buf,'.')
    p = INDEX(buf,'e')
    q = INDEX(buf,'E')
    IF(p /= -1 .AND. p /= 0) THEN
       exponential = .TRUE.
    END IF
    IF(q /= -1 .AND. q /= 0) THEN
       exponential = .TRUE.
       p = q
    END IF
    IF (k /= -1 .AND. k /= 0) THEN
       !// We have a floating point number
       itmp = 0
       !// Get the integer part of the number
       DO i = foffs, k - 1
          IF(buf(i:i) .EQ. '-') THEN
             neg = .TRUE.
             CONTINUE
          END IF
          IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND. &
               buf(i:i) .NE. ' ') THEN
             itmp = itmp * 10
             itmp = itmp + (IACHAR(buf(i:i)) - 48)
          END IF
       END DO
       res = DBLE(itmp)
       dtmp = 0.
       q = 0
       IF(.NOT. exponential) THEN
          !// We do not have an exponential number
          DO i = LEN_TRIM(buf), k+1, -1
             itmp = (IACHAR(buf(i:i)) - 48)
             dtmp = dtmp + FLOAT(itmp)
             dtmp = dtmp / 10.
          END DO
       ELSE
          !// We have an exponential number
          DO i = p-1, k+1, -1
             itmp = (IACHAR(buf(i:i)) - 48)
             dtmp = dtmp + FLOAT(itmp)
             dtmp = dtmp / 10.
          END DO
          q = a2i(buf(p+1:LEN_TRIM(buf)))
       END IF
       res = res + dtmp
       IF(exponential) THEN
          res = res * 10.**q
       END IF
       IF(neg) THEN
          res = res * -1
       END IF
    ELSE
       !// We have an integer
       res = DBLE(a2i(buf))
    END IF
  END FUNCTION a2d

  FUNCTION front_trim(buf) RESULT(res)
    IMPLICIT NONE
    !// Arguments
    CHARACTER(LEN=*), INTENT(IN)   :: buf
    INTEGER                        :: res
    !// Local variables
    INTEGER                        :: i, lng
    lng = LEN_TRIM(buf)
    DO i = 1, lng
       IF (buf(i:i) .NE. ' ') THEN
          res = i
          RETURN
       END IF
    END DO
  END FUNCTION front_trim
  
  FUNCTION a2i(buf) RESULT(res)
    IMPLICIT NONE
    !// Arguments
    CHARACTER(LEN=*), INTENT(IN)   :: buf
    INTEGER                        :: res
    !// Local variables
    INTEGER                        :: i, foffs, lng
!    INTEGER, EXTERNAL              :: front_trim
    LOGICAL                        :: neg
    neg = .FALSE.
    lng = LEN_TRIM(buf)
    foffs = front_trim(buf)
    res = 0
    DO i = foffs, lng
       IF(buf(i:i) .EQ. '-') THEN
          neg = .TRUE.
          CONTINUE
       END IF
       IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND. &
            buf(i:i) .NE. ' ') THEN
          res = res * 10
          res = res + (IACHAR(buf(i:i)) - 48)
       END IF
    END DO
    IF(neg) THEN
       res = res * -1
    END IF
  END FUNCTION a2i
  
end module ds_focalplane
