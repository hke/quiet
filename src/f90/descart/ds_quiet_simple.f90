MODULE DS_QUIET_SIMPLE
  USE DS_TYPES
  USE DS_FITSTOOLS
  USE PIX_TOOLS
  USE DS_MULTIDETECTOR
  USE DS_UTILS !for ds_assert	

  IMPLICIT NONE

!  INCLUDE 'mpif.h'

  CONTAINS

    SUBROUTINE READ_QUIET_TOD(filename,detPointing,detData,correlator,noise_only,nside,offsetLength)
      CHARACTER(len=*),INTENT(IN) :: filename
      TYPE(ds_detpointing),INTENT(OUT) :: detPointing
	  TYPE(ds_timestream),INTENT(OUT) :: detData
	  TYPE(ds_correlator),INTENT(INOUT) :: correlator
	  LOGICAL,INTENT(IN) :: noise_only
	  integer, intent(in) :: offsetLength
	integer, intent(in) :: nside
      INTEGER :: naxis1,naxis2, itod, ipix, size, rank, ierror, dummy_nside
      REAL(8) :: delta_t, theta, phi
      REAL(8),ALLOCATABLE,DIMENSION(:) :: dummysig1, dummysig2, dummyq, dummyu, ra, dec, dummyNoise
		integer excessData	  

	  !Determine whether this process is q or u
	  ierror=0
	  CALL MPI_COMM_SIZE(correlator%comm,size,ierror)
	  call ds_assert(ierror==0,'Error in mpi_comm_size')
	  CAll MPI_COMM_RANK(correlator%comm,rank,ierror)
	  call ds_assert(ierror==0,'Error in mpi_comm_rank')

      !read npix and nt from fits header
      CALL READ_QUIET_HEADER(TRIM(filename),naxis1,naxis2,dummy_nside,delta_t)
      excessData = mod(naxis2,offsetLength)
      if (excessData .ne. 0 ) then
      	write(*,*) "We are throwing away ", excessData, "datapoints"
      	write(*,*) "This is a percentage ", 100*real(excessData,kind=dp) / real(naxis2,kind=dp), "of the total."
      endif 
      
      naxis2 = naxis2 - excessData
      
	  detPointing%npix= 12*nside**2
	  detPointing%nt= naxis2
	  detData%nt= naxis2

      !allocate TOD
	  !if noise_only is called, need two dummy arrays for qsn and usn
	  !if this is U, need a further dummy array
!	  ALLOCATE(detPointing%theta(1:detPointing%nt))
!	  ALLOCATE(detData%timestream(1:detPointing%nt))
!	  IF(noise_only) THEN
!			ALLOCATE(dummysig1(1:detPointing%nt))
!			ALLOCATE(dummysig2(1:detPointing%nt))
!	  ENDIF	
!	  IF(.NOT.correlator%amIQ) THEN
!			ALLOCATE(dummyq(1:detPointing%nt))
!	  ENDIF

      !read tod
!	  ALLOCATE(ra(1:detPointing%nt))
!	  ALLOCATE(dec(1:detPointing%nt))
	write(*,*) "Need to rewrite the read-in code"
	stop
!      IF(noise_only) THEN
!			IF(correlator%amIQ) THEN
!				CALL FAST_FITS_READ(TRIM(filename),naxis1, naxis2,ra,dec,detPointing%theta,dummyq,dummyu,detData%timestream)
!			ELSE	
!				CALL FAST_FITS_READ(TRIM(filename),naxis1, naxis2,ra,dec,detPointing%theta,dummyq,dummyu,dummynoise,detData%timestream)
!				DEALLOCATE(dummynoise)
!			ENDIF	
!			DEALLOCATE(dummyq)
!			DEALLOCATE(dummyu)
!	  ELSE
!			IF(correlator%amIQ) THEN
!				CALL FAST_FITS_READ(TRIM(filename),naxis1, naxis2,ra,dec,detPointing%theta,detData%timestream)
!			ELSE	
!				CALL FAST_FITS_READ(TRIM(filename),naxis1, naxis2,ra,dec,detPointing%theta,dummyq,detData%timestream)
!				DEALLOCATE(dummyq)
!			ENDIF
!	  ENDIF
	  
	  !convert ra and dec to pixel
	  ALLOCATE(detPointing%pixel(1:detPointing%nt))	
	  DO itod= 1,detPointing%nt
	    theta= dec(itod) * DEG2RAD
		phi= ra(itod) * DEG2RAD
		CALL ANG2PIX_RING(nside,theta,phi,ipix)
		detPointing%pixel(itod) = ipix
	  ENDDO	

	  DEALLOCATE(ra)
	  DEALLOCATE(dec)

    END SUBROUTINE READ_QUIET_TOD



    SUBROUTINE READ_QUIET_HEADER(filename,naxis1,naxis2,nside,delta_t)

      TYPE(fits_header) :: pointhead
      CHARACTER(len=*),INTENT(IN) :: filename
      INTEGER,INTENT(OUT) :: naxis1,naxis2, nside
      REAL(DP),INTENT(OUT) :: delta_t


      pointhead%nrecords= 4 !34 total number in header
      ALLOCATE(pointhead%name(pointhead%nrecords))
      ALLOCATE(pointhead%htype(pointhead%nrecords))
      ALLOCATE(pointhead%record_dp(pointhead%nrecords))
      ALLOCATE(pointhead%record_int(pointhead%nrecords))

      pointhead%name(1)= 'DELTA_T'
      pointhead%htype(1)= 'D'
      pointhead%name(2)= 'NAXIS2'
      pointhead%htype(2)= 'J'
      pointhead%name(3)= 'NSIDE'
      pointhead%htype(3)= 'J'
      pointhead%name(4)= 'NAXIS1'
      pointhead%htype(4)= 'J'

      CALL READ_FITS_HEADER(filename,pointhead)

      delta_t= pointhead%record_dp(1)
      naxis2=pointhead%record_int(2)
      nside=pointhead%record_int(3)
      naxis1=pointhead%record_int(4)
	  

      DEALLOCATE(pointhead%name)
      DEALLOCATE(pointhead%htype)
      DEALLOCATE(pointhead%record_dp)
      DEALLOCATE(pointhead%record_int)

    END SUBROUTINE READ_QUIET_HEADER


  END MODULE DS_QUIET_SIMPLE
