MODULE ds_globalpointing
  USE ds_types
  USE ds_fitstools
  IMPLICIT NONE
  external  sla_dbear
 real(dp) sla_dbear

  TYPE FP_POINTING
     INTEGER :: n_obs
     REAL(DP) :: latitude, longitude ![deg],[deg]
     REAL(DP),ALLOCATABLE,DIMENSION(:) :: lst, az, el, fp_angle, hwp_angle ![hr],[deg],[deg],[deg],[deg]
     INTEGER,ALLOCATABLE,DIMENSION(:) :: flag ![see convention below]
  END TYPE FP_POINTING

  !Flag conventions:
  !  flag  > 0 : inside field
  ! |flag| = 1 : constant vel
  ! |flag| = 2 : accelerating
  ! |flag| = 3 : deccelerating
  
!JZ I moved the fp_point array into the subroutines so we have no global variable so that 
!DWPS: added lat and long to read_pointing_header and then to the fp_pointing type


  CONTAINS

    SUBROUTINE LOAD_GLOBAL_POINTING(n_point, point_path,fp_point)
  TYPE(FP_POINTING),ALLOCATABLE,DIMENSION(:) :: fp_point
      
      INTEGER,INTENT(IN) :: n_point
      CHARACTER(len=*),INTENT(IN) :: point_path

      INTEGER :: ipoint, ob_num, tot_obs, bore_step
      REAL(DP) :: obs_el, hwp_vel, bore_ang,latitude,longitude
      integer :: datarate  !JZ moved from real
      CHARACTER(len=200),ALLOCATABLE,DIMENSION(:) :: filename
      INTEGER,ALLOCATABLE,DIMENSION(:) :: naxis1, naxis2
      REAL(DP),ALLOCATABLE,DIMENSION(:) :: ra,dec
      REAL(8),ALLOCATABLE,DIMENSION(:) :: crossing_angle
      REAL(DP), allocatable, dimension(:) :: temp_flags  !JZ The fast_fits_read routine only reads doubles
      													 !JZ so we need to convert to ints
      
      ALLOCATE(filename(n_point))
      ALLOCATE(naxis1(n_point), naxis2(n_point))

      !Build filenames
      DO ipoint= 1, n_point
         write(filename(ipoint),'(a,i1,a)') trim(point_path),ipoint,'.fits'
      ENDDO

      !read data size from each file
      DO ipoint= 1, n_point
         CALL READ_POINTING_HEADER(filename(ipoint),naxis1(ipoint),naxis2(ipoint),ob_num,tot_obs,datarate,obs_el,hwp_vel,bore_step,bore_ang,longitude,latitude)
      ENDDO

      ALLOCATE(fp_point(n_point))

      fp_point(:)%n_obs= naxis2  !JZ modified to avoid rank-assignment error.  Nice fortran feature here.
      fp_point(:)%latitude= latitude
      fp_point(:)%longitude= longitude

      DO ipoint= 1,n_point
         ALLOCATE(fp_point(ipoint)%lst(naxis2(ipoint)))
         ALLOCATE(fp_point(ipoint)%az(naxis2(ipoint)))
         ALLOCATE(fp_point(ipoint)%el(naxis2(ipoint)))
         ALLOCATE(fp_point(ipoint)%hwp_angle(naxis2(ipoint)))
         ALLOCATE(fp_point(ipoint)%fp_angle(naxis2(ipoint)))
         ALLOCATE(fp_point(ipoint)%flag(naxis2(ipoint)))

         ALLOCATE(ra(naxis2(ipoint)), dec(naxis2(ipoint)))
         allocate(temp_flags(naxis2(ipoint)))

         CALL FAST_FITS_READ(filename(ipoint),naxis1(ipoint),naxis2(ipoint),fp_point(ipoint)%lst,ra,dec,fp_point(ipoint)%az,  &  !JZ changed arguments to match ds_fitstools
              fp_point(ipoint)%el,fp_point(ipoint)%fp_angle,fp_point(ipoint)%hwp_angle,temp_flags)
              fp_point(ipoint)%flag = floor(temp_flags)  !JZ We convert real(8) -> integer
         
         ra  =  ra * deg2rad
         dec = dec * deg2rad
         CALL SkyOrientation(naxis2(ipoint),ra,dec,crossing_angle)
         fp_point(ipoint)%fp_angle = fp_point(ipoint)%fp_angle + ( crossing_angle * rad2deg )!DWPS: both are -90<phi<+90, measured east from zenith (looking out)

         !Added this code to make sure angles are -180<phi<180
         WHERE(fp_point(ipoint)%fp_angle.lt.-180.0_8)
            fp_point(ipoint)%fp_angle = 360.0_8 + fp_point(ipoint)%fp_angle
         ELSEWHERE(fp_point(ipoint)%fp_angle.gt.180.0_8)
            fp_point(ipoint)%fp_angle = -360.0_8 + fp_point(ipoint)%fp_angle
         ENDWHERE

         DEALLOCATE(ra,dec,temp_flags)
         
      ENDDO


    END SUBROUTINE LOAD_GLOBAL_POINTING

    
    !###############################################################
    !#  Call this routine when finished with global pointing data  #
    !###############################################################
        !JZ When compiling this on gfortran this subroutine could be replaced with a simple
        !JZ deallocate(fp_point), but I do not know if that is true on other compilers.
    SUBROUTINE CLEAR_GLOBAL_POINTING(fp_point)
      INTEGER :: npoint, ipoint
        TYPE(FP_POINTING),ALLOCATABLE,DIMENSION(:) :: fp_point

      npoint= size(fp_point)

      DO ipoint= 1, npoint
         IF(ALLOCATED(fp_point(ipoint)%lst)) DEALLOCATE(fp_point(ipoint)%lst)
         IF(ALLOCATED(fp_point(ipoint)%az)) DEALLOCATE(fp_point(ipoint)%az)
         IF(ALLOCATED(fp_point(ipoint)%el)) DEALLOCATE(fp_point(ipoint)%el)
         IF(ALLOCATED(fp_point(ipoint)%fp_angle)) DEALLOCATE(fp_point(ipoint)%fp_angle)
      ENDDO

      IF(ALLOCATED(fp_point)) DEALLOCATE(fp_point)

    END SUBROUTINE CLEAR_GLOBAL_POINTING


    !##########################################################
    !#  Determine orientation of telescope from its movement  #
    !##########################################################
    !JZ changed pio2_d to halfpi
    !DWPS changed output to degrees and removed coord change
    SUBROUTINE SkyOrientation(ndata,ra,dec,crossing_angle)
      
      INTEGER(4),INTENT(IN):: ndata
      REAL(8),ALLOCATABLE,INTENT(OUT),DIMENSION(:) :: crossing_angle
      REAL(8),ALLOCATABLE,INTENT(INOUT),DIMENSION(:) :: ra, dec
      INTEGER(4) :: idata

      ! construct array of sky crossing angles
      allocate(crossing_angle(ndata))
      do idata = 2, ndata
         crossing_angle(idata) = sla_dbear(ra(idata-1), dec(idata-1), ra(idata), dec(idata))!Need to be sure about conventions
         ! rotate to coordinate system used for beam convolution
         !if (crossing_angle(idata) .ge. -HALFPI) then
         !   crossing_angle(idata) = HALFPI - crossing_angle(idata)
         !else
         !   crossing_angle(idata) = -3.d0*HALFPI - crossing_angle(idata)
         !endif
      enddo
      crossing_angle(1) = crossing_angle(2)
      crossing_angle= crossing_angle * rad2deg

    END SUBROUTINE SkyOrientation
      


END MODULE ds_globalpointing
