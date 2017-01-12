module ds_pointing

use ds_utils
use ds_types
use ds_precision 
use ds_focalplane
use ds_globalpointing

implicit none


contains

!JZ we can probably replace a bunch of these functions with things from slalib

function rotationMatrixAbout(axis,theta) result(R)
implicit none
real(dp), dimension(3) :: axis
real(dp), dimension(3,3) :: R
real(dp) :: theta
real(dp) :: c,s,t,wx,wy,wz
c = cos(theta)
s = sin(theta)
t = 1-c
wx = axis(1)
wy = axis(2)
wz = axis(3)

R(1,1) = c+wx*wx*t
R(2,1) = wz*s + wx*wy*t
R(3,1) = -wy*s+wx*wz*t

R(1,2) = wx*wy*t-wz*s
R(2,2) = c+wy*wy*t
R(3,2) = wx*s+wy*wz*t

R(1,3) = wy*s+wx*wz*t
R(2,3) = -wx*s+wy*wz*t
R(3,3) = c+wz*wz*t



end function rotationMatrixAbout


function crossProduct(a,b) result(n)
implicit none
real(dp), dimension(3) :: a,b,n
n(1) = a(2)*b(3) - a(3)*b(2)
n(2) = a(3)*b(1) - a(1)*b(3)
n(3) = a(1)*b(2) - a(2)*b(1)
end function crossProduct


subroutine normalize(x)
implicit none
real(dp), dimension(3) :: x
x = x/vectorNorm(x)
end subroutine normalize

function dotProduct(x,y)
implicit none
real(dp), dimension(3) ::x,y
real(dp) dotProduct
dotProduct = sum(x*y)
end function dotProduct

function vectorNorm(x)
implicit none
real(dp), dimension(3) :: x
real(dp) vectorNorm
vectorNorm = sqrt(dotProduct(x,x))
end function vectorNorm


function getDetectorAngle(focalX, focalY, focalScale, az, el,boresightAngle)
!Convert a focal plane location into an azimuth and elevation after any rotation of the telescope.
!The rotation is specified by a rotation about north by the boresight angle followed by a 
implicit none

real(dp), dimension(2) :: getDetectorAngle
real(dp) :: focalX, focalY, az, el, boresightAngle, focalScale
real(dp), dimension(3) :: x, boresightDirection,north,d,axis
real(dp), dimension(3,3) :: R1,R2
real(dp) :: cel

!Get the normal cartesian vector corresponding to the unrotated detector position.
x(1) = focalX / focalScale
x(2) = focalY / focalScale
x(3) = 1 - x(1)*x(1) - x(2)*x(2)

!Convert the boresight direction to cartesian.
cel = cos(el)
boresightDirection(1) = cel * sin(az) 
boresightDirection(2) = cel * cos(az)
boresightDirection(3) = sin(el)

!The first rotation, about the boresight, is about north, and has magnitude boresightAngle
north(1)=0
north(2)=0
north(3)=1
R2 = rotationMatrixAbout(north,boresightAngle)


!The second rotation moves north to the boresight direction.  It is therefore about the direction orthogonal to both.
axis = crossProduct(north,boresightDirection)
call normalize(axis)
R1 = rotationMatrixAbout(axis,HALFPI - el)

!Apply both rotations to the detector position.
d = matmul(matmul(R2,R1),x)

!Convert back to phi and theta angles.
getDetectorAngle(1) = atan2(d(2),d(1))   !phi
getDetectorAngle(2) = asin(d(3)/vectorNorm(d))   !theta
end function getDetectorAngle



subroutine getHealpixIndicesForDetector(focalX,focalY,focalScale,azimuth, elevation,boresightAngle,localSidereal,latitude,nside,useRing,pixelIndices)

use pix_tools
implicit none
real(dp), intent(in) :: focalX, focalY,focalScale
real(dp), intent(in) :: latitude
real(dp), intent(in), dimension(1:) :: azimuth, elevation,boresightAngle,localSidereal
integer(sp), intent(out), dimension(1:) :: pixelIndices
logical, intent(in) :: useRing
integer(sp), intent(in) ::  nside

real(dp), dimension(2) :: detectorAngle
real(dp) :: theta,phi
integer(dp) :: t


call ds_assert(size(azimuth) == size(elevation),'Azimuths and elevations have different lengths')
call ds_assert(size(azimuth) == size(boresightAngle),'Azimuths and boresight angles have different lengths') 
call ds_assert(size(azimuth) == size(pixelIndices),'Azimuths and pixel indices have different lengths')


do t=1,size(azimuth)
	detectorAngle = getDetectorAngle(focalX,focalY,focalScale,azimuth(t),elevation(t),boresightAngle(t))
        call sla_dh2e(detectorAngle(1),detectorAngle(2),latitude,phi,theta)  !Convert az,el to hour angle and dec.
        phi = localSidereal(t) - phi
        theta = HALFPI-theta
	if (useRing) then
		call ang2pix_ring(nside,theta,phi,pixelIndices(t))
	else
		call ang2pix_nest(nside,theta,phi,pixelIndices(t))
	endif
enddo

end subroutine getHealpixIndicesForDetector


!JZ This is the main driver routine
subroutine generateDetectorPointingsFromFiles(pointing_file_base, focal_plane_file, focal_length,nside,detector_pointings)
	character(*) pointing_file_base, focal_plane_file
	real(dp) focal_length
	type(FP_POINTING), allocatable, dimension(:) :: pointing_array
	type(ds_detpointing), allocatable, dimension(:), intent(OUT) :: detector_pointings
	type(ds_plane) :: plane
	integer nside
	
	call LOAD_FOCAL_PLANE(focal_plane_file,plane,focal_length)
	call LOAD_GLOBAL_POINTING(1, pointing_file_base,pointing_array)
	call generateDetectorPointings(pointing_array(1), plane, nside, detector_pointings)
	
end subroutine generateDetectorPointingsFromFiles

function hours_to_radians(hr)
real(dp) :: hr,hours_to_radians
hours_to_radians = hr / 24.0_dp * TWOPI 
end function

subroutine generateDetectorPointings(focal_plane_pointing, focal_plane_info, nside, detector_pointings)
use pix_tools
	type(ds_plane), intent(IN) :: focal_plane_info
	type(fp_pointing), intent(IN) :: focal_plane_pointing
	type(ds_detpointing), allocatable, dimension(:), intent(OUT) :: detector_pointings
	integer, intent(IN) :: nside
	real(dp) :: ra,dec,az,el,xoff,yoff, det_orient, phi,theta,slat,clat,orient_offset, hour_angle
	integer :: det, t
	integer pix
	
	allocate(detector_pointings(focal_plane_info%ndet))
	
	
	slat = sin(focal_plane_pointing%latitude*deg2rad)
	clat = cos(focal_plane_pointing%latitude*deg2rad)

	do det=1,focal_plane_info%ndet

		call preparePointing(detector_pointings(det), focal_plane_pointing%n_obs)
        xoff = focal_plane_info%det(det)%x / focal_plane_info%focal_length
        yoff = focal_plane_info%det(det)%y / focal_plane_info%focal_length
        orient_offset = focal_plane_info%det(det)%orient
       write(*,*) "Error here - rotate focal plane" 
        

		do t=1,focal_plane_pointing%n_obs

               el = focal_plane_pointing%el(t)*deg2rad + yoff
               az = focal_plane_pointing%az(t)*deg2rad + xoff / cos(el)
               
               ! {az,el} -> {ra,dec}
!               call azel2hadec(az, el, slat,clat, ra, dec, -pi)         

				call sla_dh2e(az,el,focal_plane_pointing%latitude*deg2rad,hour_angle,dec)

               ra = focal_plane_pointing%lst(det)*HOUR2RAD - hour_angle
               ! Get pixel
               phi = mod(ra, twopi)
               theta = (halfpi - dec)
               call ang2pix_ring(nside, theta, phi, pix)

               detector_pointings(det)%pixel(t) = pix
				!Get the orientation of the detector at this time step.
				!The sum of the plane rotation and the orientation of the detector in the plane
               det_orient = focal_plane_pointing%fp_angle(t) + orient_offset

               ! ensure we're in [-pi, pi]
				det_orient = mod(det_orient, twopi)
               if (det_orient .gt. twopi) det_orient = det_orient - twopi
               if (det_orient .lt. -twopi) det_orient = det_orient + twopi

               detector_pointings(det)%theta(t) = det_orient
              
        enddo
	enddo
	
end subroutine generateDetectorPointings




end module ds_pointing




