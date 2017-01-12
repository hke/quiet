program bscan
  use quiet_utils
  use quiet_mapfile_mod
  use quiet_pointing_mod
  use quiet_ephem_mod
  use rngmod
  use spline_1D_mod
  use quiet_module_mod
  implicit none

  type flight
     character(len=128)       :: name
     logical(lgt)             :: enable, output_hitcount
     real(dp)                 :: altitude
     real(dp)                 :: phi0, phi1, theta0, theta1
     real(dp)                 :: lon0, lat0
     real(dp)                 :: duration, mjd_start, mjd
     real(dp), dimension(3,3) :: M
  end type flight

  type patch
     character(len=128)       :: name
     real(dp)                 :: phi, theta, phi_gal, theta_gal
     real(dp), dimension(3,3) :: M
     logical(lgt)             :: enable
  end type patch

  integer(i4b)       :: myid, ierr, numprocs, root, i, j, k, l, p, n, d, iargc, npatch, nflight
  integer(i4b)       :: nside, npix, ordering, unit, listpix(0:3000000), state, npt, naz, counter, nmod
  real(dp)           :: psi, dt, mjd, az, el, dk, pos(2), pos2(2), sundist, min_sun_dist, vec1(3), vec2(3), el_min, el_max, dec, ra, theta, phi, this_el_max, this_el_min, dec_min, dec_max, obstime, suntime, eltime, samprate
  real(dp)           :: az_tele, el_tele, dk_tele, psize, scanspeed, patchdist, theta_tele, phi_tele, psi_tele, vec_tele(3), time_since_turn, turn_time, az_speed, accel, az_limit, az_repoint_limit, min_ces_time, time_since_repointing, cestime, pdist
  real(dp)           :: theta_diode, psi_diode, phi_diode, el_diode, az_diode, dk_diode, vec_diode(3), mat(3,3), ra_tele, dec_tele, daz
  character(len=256) :: paramfile, filename, outfile
  character(len=2)   :: itext
  character(len=4)   :: ctext
  real(dp), allocatable, dimension(:)     :: az_profile
  real(dp), allocatable, dimension(:,:)   :: mask, hitcount, hitcount_tot
  real(dp), allocatable, dimension(:,:,:) :: centerpath
  logical(lgt)       :: leftscan, ok
  type(planck_rng)   :: handle
  
  type(flight), allocatable, dimension(:) :: flights
  type(patch),  allocatable, dimension(:) :: patches

!  call mpi_init(ierr)
!  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
!  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  
!  root = 0
  unit = getlun()
  
  if (iargc() /= 1) then
     write(*,*) ''
     write(*,*) 'Usage: bscan [parameter file]'
     write(*,*) ''
     
!     call mpi_finalize(ierr)
     stop
  end if
  
  call getarg(1,paramfile)
  call get_parameter(unit, paramfile, 'NSIDE',                par_int=nside)
  call get_parameter(unit, paramfile, 'NUMPATCH',             par_int=npatch)
  call get_parameter(unit, paramfile, 'NUMFLIGHT',            par_int=nflight)
  call get_parameter(unit, paramfile, 'TIME_STEP_MINUTES',    par_dp=dt)
  call get_parameter(unit, paramfile, 'MINIMUM_SUN_DISTANCE', par_dp=min_sun_dist)
  call get_parameter(unit, paramfile, 'MINIMUM_ELEVATION',    par_dp=el_min)
  call get_parameter(unit, paramfile, 'MAXIMUM_ELEVATION',    par_dp=el_max)
  call get_parameter(unit, paramfile, 'SAMPLING_RATE',        par_dp=samprate)
  call get_parameter(unit, paramfile, 'PATCH_SIZE',           par_dp=psize)
  call get_parameter(unit, paramfile, 'SCAN_SPEED',           par_dp=scanspeed)
  call get_parameter(unit, paramfile, 'SCAN_ACCELERATION',    par_dp=accel)
  call get_parameter(unit, paramfile, 'AZ_REPOINT_LIMIT',     par_dp=az_repoint_limit)
  call get_parameter(unit, paramfile, 'MIN_CES_TIME',         par_dp=min_ces_time)
  psize = psize * DEG2RAD
  dt = dt / (24.d0*60.d0)
  scanspeed = scanspeed * DEG2RAD
  turn_time = 1.d0
  accel     = accel     * DEG2RAD
  az_repoint_limit = az_repoint_limit * DEG2RAD
  call rand_init(handle, 18471)

  ! Initialize modules
  call initialize_quiet_pointing_mod(paramfile, apply_mount_model=.false., apparent_correction=.false.)
  call ephem_set_db('unix.405')
  nmod      = get_num_modules()

  ! Initialize flights
  allocate(flights(nflight))
  do i = 1, nflight
     call int2string(i, itext)
     call get_parameter(unit, paramfile, 'FLIGHT_NAME'//itext,               par_string=flights(i)%name)
     call get_parameter(unit, paramfile, 'FLIGHT_LAUNCH_LATITUDE'//itext,    par_dp=flights(i)%lat0)
     call get_parameter(unit, paramfile, 'FLIGHT_LAUNCH_LONGITUDE'//itext,   par_dp=flights(i)%lon0)
     call get_parameter(unit, paramfile, 'FLIGHT_LANDING_LATITUDE'//itext,   par_dp=flights(i)%theta1)
     call get_parameter(unit, paramfile, 'FLIGHT_LANDING_LONGITUDE'//itext,  par_dp=flights(i)%phi1)
     call get_parameter(unit, paramfile, 'FLIGHT_DURATION'//itext,           par_dp=flights(i)%duration)
     call get_parameter(unit, paramfile, 'FLIGHT_LAUNCH_MJD'//itext,         par_dp=flights(i)%mjd_start)
     call get_parameter(unit, paramfile, 'FLIGHT_ALTITUDE'//itext,           par_dp=flights(i)%altitude)
     call get_parameter(unit, paramfile, 'FLIGHT_ENABLE'//itext,             par_lgt=flights(i)%enable)
     call get_parameter(unit, paramfile, 'FLIGHT_OUTPUT_HITCOUNT'//itext,    par_lgt=flights(i)%output_hitcount)
     flights(i)%theta0 = flights(i)%lat0 * DEG2RAD
     flights(i)%theta1 = flights(i)%theta1 * DEG2RAD
     flights(i)%phi0   = flights(i)%lon0   * DEG2RAD
     flights(i)%phi1   = flights(i)%phi1   * DEG2RAD
     call swap_coordinate_convention(flights(i)%phi0, flights(i)%theta0, psi, COORD_CEL)
     call swap_coordinate_convention(flights(i)%phi1, flights(i)%theta1, psi, COORD_CEL)
     pos = ephem(EPHEM_SUN, flights(i)%mjd_start)
     call coord_convert(COORD_CEL, pos(1), pos(2), 0.d0, COORD_GAL, phi, theta, psi)
     write(*,*) i, ', gal sun = ', modulo(phi*RAD2DEG,360.d0), 90.d0-theta*RAD2DEG

     allocate(mask(0:12*nside**2-1,1))
     mask = 1.d0
     pos = ephem(EPHEM_SUN, flights(i)%mjd_start)
     call coord_convert(COORD_CEL, pos(1), pos(2), 0.d0, COORD_GAL, phi, theta, psi)
     call ang2vec(theta, phi, vec1)
     call query_disc(nside, vec1, min_sun_dist*DEG2RAD, listpix, j)
     mask(listpix(0:j-1),1) = 0.d0
     call write_map(mask, 1, 'mask_sun'//itext//'.fits')
     deallocate(mask)
  end do

  allocate(mask(0:12*nside**2-1,1))
  mask = 0.d0 !-1.6375d30
  do i = 1, 12
     mjd = 57023.d0 + (i-1)*30
     pos = ephem(EPHEM_SUN, mjd)
!     write(*,*) 'sun     = ', pos(1)*RAD2DEG, (0.5d0*pi-pos(2))*RAD2DEG
!     call ang2vec(pos(2), pos(1), vec1)
!     call ang2vec(pos2(2), pos2(1), vec2)
!     call angdist(vec1, vec2, mjd)
!     write(*,*) mjd*RAD2DEG
     
     call coord_convert(COORD_CEL, pos(1), pos(2), 0.d0, COORD_GAL, phi, theta, psi)
     call ang2vec(theta, phi, vec1)
     call query_disc(nside, vec1, 5.d0*DEG2RAD, listpix, j)
     mask(listpix(0:j-1),1) = mask(listpix(0:j-1),1) + i
  end do
  call write_map(mask, 1, 'sunpos.fits')
  deallocate(mask)

  ! Initialize patches
  allocate(patches(npatch))
  do i = 1, npatch
     call int2string(i, itext)
     call get_parameter(unit, paramfile, 'PATCH_NAME'//itext,            par_string=patches(i)%name)
     call get_parameter(unit, paramfile, 'PATCH_GAL_LATITUDE'//itext,    par_dp=patches(i)%theta_gal)
     call get_parameter(unit, paramfile, 'PATCH_GAL_LONGITUDE'//itext,   par_dp=patches(i)%phi_gal)
     call get_parameter(unit, paramfile, 'PATCH_ENABLE'//itext,          par_lgt=patches(i)%enable)
     patches(i)%theta_gal = patches(i)%theta_gal * DEG2RAD
     patches(i)%phi_gal   = patches(i)%phi_gal   * DEG2RAD
     call swap_coordinate_convention(patches(i)%phi_gal, patches(i)%theta_gal, psi, COORD_GAL)
     call coord_convert(COORD_GAL, patches(i)%phi_gal, patches(i)%theta_gal, 0.d0, COORD_CEL, patches(i)%phi, patches(i)%theta, psi)
!     write(*,*) 'Celestial coordinates for patch ', i, ': ', real(patches(i)%phi*RAD2DEG,sp), real(90.d0-patches(i)%theta*RAD2DEG,sp)
  end do


  ! Compute paths
  open(58,file='az.dat')
  open(59,file='el.dat')
  open(60,file='sundist.dat')
  write(*,*) 'Patch    Flight    Observing time   Total accept    Sun accept    El accept'
  do j = 1, npatch
     if (.not. patches(j)%enable) cycle
     do i = 1, nflight
        if (.not. flights(i)%enable) cycle

        write(58,*) '# Flight             = ', trim(flights(i)%name)
        write(58,*) '#   Launch latitude  = ', flights(i)%lat0
        write(58,*) '#   Launch longitude = ', flights(i)%lon0
        write(58,*) '#   Altitude         = ', flights(i)%altitude
        mjd     = flights(i)%mjd_start
        obstime = 0.d0
        suntime = 0.d0
        eltime  = 0.d0
        npt = floor(flights(i)%duration/dt)+1
        allocate(centerpath(npt,3,2))
        p = 1
        do while (p <= npt)

           ! Update position
           QUIET_GEODETIC_LATITUDE  = flights(i)%lat0
           QUIET_GEODETIC_LONGITUDE = flights(i)%lon0
           QUIET_ALTITUDE           = flights(i)%altitude

           pos2 = ephem(EPHEM_JUPITER, flights(i)%mjd_start)
           pos = ephem(EPHEM_SUN, mjd)
!           call ang2vec(pos(2), pos(1), vec1)
!           write(*,*) 'sun     = ', pos(1)*RAD2DEG, (0.5*pi-pos(2))*RAD2DEG
!           write(*,*) 'jupiter = ', pos2(1)*RAD2DEG, (0.5*pi-pos2(2))*RAD2DEG
!           call ang2vec(pos2(2), pos2(1), vec2)
!           call angdist(vec1, vec2, mjd)
!           write(*,*) mjd*RAD2DEG
!           stop

           call ang2vec(pos(2), pos(1), vec1)
           call coord_convert(COORD_GAL, patches(j)%phi_gal, patches(j)%theta_gal, 0.d0, COORD_HOR, az, el, dk, mjd=mjd)           
           call ang2vec(patches(j)%theta, patches(j)%phi, vec2)
           centerpath(p,1,:) = mjd
           centerpath(p,2,1) = az
           centerpath(p,2,2) = el
           sundist = acos(sum(vec1*vec2)) * RAD2DEG
           write(58,fmt='(2f12.4)') (mjd-flights(i)%mjd_start)*24., az*RAD2DEG
           write(59,fmt='(2f12.4)') (mjd-flights(i)%mjd_start)*24., 90.-el*RAD2DEG
           write(60,fmt='(2f12.4)') (mjd-flights(i)%mjd_start)*24., sundist
           call swap_coordinate_convention(az, el, dk, COORD_HOR)
           el = el*RAD2DEG
           if ((sundist > min_sun_dist .or. pos(2) < 0.d0) .and. el > el_min .and. el < el_max) obstime = obstime + dt
           if ((sundist > min_sun_dist .or. pos(2) < 0.d0)) suntime = suntime + dt
           if (el > el_min .and. el < el_max) eltime = eltime + dt

           ! Move to next time step
           mjd = mjd + dt
           p   = p+1
        end do
        write(58,*) 
        write(59,*) 
        write(60,*) 

        ! Make sure that azimuth is continous
        do p = 2, npt
           if (abs(centerpath(p,2,1)-centerpath(p-1,2,1)) > 0.5d0*pi) then
              if (centerpath(p-1,2,1) < 0.d0) then
                 centerpath(p,2,1) = centerpath(p,2,1) - 2.d0*pi
              else
                 centerpath(p,2,1) = centerpath(p,2,1) + 2.d0*pi
              end if
           end if
        end do
        
        write(*,fmt='(2i6,f14.4,3f14.1)') j, i, 24*obstime, obstime / flights(i)%duration*100, suntime / flights(i)%duration*100, eltime / flights(i)%duration*100

        call spline(centerpath(:,1,1), centerpath(:,2,1), 1.d30, 1.d30, centerpath(:,3,1))
        call spline(centerpath(:,1,2), centerpath(:,2,2), 1.d30, 1.d30, centerpath(:,3,2))

        if (flights(i)%output_hitcount) then
           allocate(hitcount(0:12*nside**2-1,1))

           mjd      = flights(i)%mjd_start
           hitcount = 0.d0
           QUIET_GEODETIC_LATITUDE  = flights(i)%lat0
           QUIET_GEODETIC_LONGITUDE = flights(i)%lon0
           QUIET_ALTITUDE           = flights(i)%altitude
           dt                       = 1.d0 / samprate / (3600.d0*24.d0)
           daz                      = scanspeed * (dt*24.d0*3600.d0)
           npt = floor(flights(i)%duration/dt)+1
           open(18,file='scan.dat', recl=1024)
           !Initialize position
           counter = 1
           az_tele = 0.d0
           !el_tele = (90.d0 - (el_min + (el_max - el_min) * rand_uni(handle))) * DEG2RAD
           !el_tele = (90.d0 - el_min) * DEG2RAD
           !el_tele = 45.d0 * DEG2RAD
           el_tele = (90.d0-min(max(abs(flights(i)%lat0),el_min),el_max)) * DEG2RAD
           dk_tele = 0.d0
           write(*,*) 'counter = ', counter, (0.5d0*pi-el_tele)*RAD2DEG
           !call coord_convert(COORD_GAL, patches(j)%phi_gal, patches(j)%theta_gal, 0.d0, COORD_HOR, az_tele, el_tele, dk_tele, mjd=mjd)     
           !az_speed = 0.d0
!           write(*,*) az_tele, el_tele
!           write(*,*) splint(centerpath(:,1,1), centerpath(:,2,1), centerpath(:,3,1), mjd), splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), mjd)
           !az_limit = find_az_limit(mjd, psize, az_tele, el_tele, centerpath)
!           write(*,*) az_limit-az_tele, az_repoint_limit
!           stop
           !time_since_repointing = 0.d0
!           if (abs(az_limit-az_tele) < az_repoint_limit) then
           !   call find_optimal_elevation(mjd, psize, centerpath, el_min, el_max, el_tele, cestime)
           !   az_tele = splint(centerpath(:,1,1), centerpath(:,2,1), centerpath(:,3,1), mjd)
!           end if
           !call compute_az_profile(samprate, az_tele, az_limit, scanspeed, accel, turn_time, az_profile)
           obstime = 0.d0
           call int2string(counter,ctext)
           open(19,file='scan'//ctext//'.dat', recl=1024)
           do while (mjd < flights(i)%mjd_start + flights(i)%duration)
              pos = ephem(EPHEM_SUN, mjd) 
              call ang2vec(pos(2), pos(1), vec1)
              call coord_convert(COORD_GAL, patches(j)%phi_gal, patches(j)%theta_gal, 0.d0, COORD_HOR, az, el, dk, mjd=mjd)     
              call ang2vec(patches(j)%theta_gal, patches(j)%phi_gal, vec2)
              do k = 1, min(1.d0, flights(i)%mjd_start+flights(i)%duration-mjd)/dt
                 az_tele = mod(az_tele + daz,2.d0*pi)
                 call coord_convert(COORD_HOR, az_tele, el_tele, 0.d0, COORD_GAL, phi_tele, theta_tele, psi_tele, mjd=mjd+k*dt)
                 call ang2vec(theta_tele, phi_tele, vec_tele)
                 call angdist(vec_tele, vec2, pdist)
!                 if (el_tele >= (90.d0-el_max)*DEG2RAD .and. el_tele <= (90.d0-el_min)*DEG2RAD .and. pdist <= psize) then
                    call coord_convert(COORD_HOR, az_tele, el_tele, 0.d0, COORD_CEL, ra_tele,  dec_tele, psi_tele, mjd=mjd+k*dt)
                    call swap_coordinate_convention(ra_tele, dec_tele, psi_tele, COORD_CEL)
                    call vec2pix_ring(nside, vec_tele, l)
                    hitcount(l,1) = hitcount(l,1) + dt * 3600.d0*24.d0
                    write(18,*) (mjd+k*dt-flights(i)%mjd_start)*3600*24, az_tele*RAD2DEG, (0.5d0*pi-el_tele)*RAD2DEG, acos(sum(vec2*vec_tele))*RAD2DEG, ra_tele*RAD2DEG, dec_tele*RAD2DEG, mjd2lst(mjd+k*dt-flights(i)%mjd_start, QUIET_GEODETIC_LONGITUDE), QUIET_GEODETIC_LATITUDE
!                    write(19,*) (mjd+k*dt-flights(i)%mjd_start)*3600*24, az_tele*RAD2DEG, (0.5d0*pi-el_tele)*RAD2DEG, acos(sum(vec2*vec_tele))*RAD2DEG
                    obstime = obstime + dt
!                 end if
              end do
              mjd = mjd + 1.d0

              if (mjd < flights(i)%mjd_start + flights(i)%duration) then
                 ! Find new path
                 counter = counter+1
                 !el_tele = (90.d0 - (el_min + (el_max - el_min) * (counter-1)/(floor(flights(i)%duration)-1.d0))) * DEG2RAD
                 el_tele = (90.d0 - (el_min + (el_max - el_min) * (counter-2)/(floor(flights(i)%duration)-2.d0))) * DEG2RAD
                 call int2string(counter,ctext)
                 open(19,file='scan'//ctext//'.dat', recl=1024)                 
                 write(*,*) 'counter = ', counter, (0.5d0*pi-el_tele)*RAD2DEG
              end if

           end do
           close(18)

           where (hitcount == 0.d0) 
              hitcount = -1.6375d30
           end where
           call write_map(hitcount, 1, 'hitcount.fits')

           if (.false.) then
              allocate(hitcount_tot(0:12*nside**2-1,1))
              hitcount_tot = 0.d0
              do p = 0, 12*nside**2-1
                 if (mod(p,1000)==0) write(*,*) p, 12*nside**2
                 if (hitcount(p,1) > 0.d0) then
                    call pix2ang_ring(nside, p, theta, phi)
                    do d = 0, nmod-1
                       call coord_convert(COORD_TELE, phi, theta, 0.d0, COORD_HOR, phi_tele, theta_tele, psi_tele, mjd=flights(i)%mjd_start, mod=d, diode=0)
                       call ang2pix_ring(nside, theta_tele, phi_tele, l)
                       hitcount_tot(l,1) = hitcount_tot(l,1) + hitcount(p,1)
                    end do
                 end if
              end do

              where (hitcount_tot == 0.d0) 
                 hitcount_tot = -1.6375d30
              end where
              call write_map(hitcount_tot, 1, 'hitcount_tot.fits')
           end if


           write(*,*) 'Acceptable observation time = ', obstime / flights(i)%duration
           
           deallocate(hitcount)
        end if


     end do
     write(*,*)
  end do
  close(58)
  close(59)


  ! Output elevation and sun distance masks
  allocate(mask(0:12*nside**2-1,1))
  do i = 1, nflight
     if (.not. flights(i)%enable) cycle

!     if (flights(i)%lat0 > el_max) then
!        dec_max = max(90 - (flights(i)%lat0-el_max), 90 - (flights(i)%lat0-el_min))
!     else
!        dec_max = 90.d0
!     end if

!        dec_min = max(90 - (flights(i)%lat0-el_max), 90 - (flights(i)%lat0-el_min))
     dec_min = flights(i)%lat0 - (90.d0-el_min)
     dec_max = flights(i)%lat0 + (90.d0-el_min)
     if (flights(i)%lat0 >  el_max) dec_max =  90.d0 - (flights(i)%lat0-el_max)
     if (flights(i)%lat0 < -el_min) dec_min = -90.d0 - (flights(i)%lat0+el_max)
!     write(*,*) i, dec_min, dec_max

     call int2string(i, itext)

     ! Compute elevation mask
     mask = 0.d0
     do j = 0, 12*nside**2-1
        call pix2ang_ring(nside, j, theta, phi)
        call coord_convert(COORD_GAL, phi, theta, 0.d0, COORD_CEL, ra, dec, psi)
        call swap_coordinate_convention(ra, dec, psi, COORD_CEL) 
        dec = dec * RAD2DEG
        if (dec > dec_min .and. dec < dec_max) mask(j,1) = 1.d0
     end do
     
     call write_map(mask, 1, 'mask_el'//itext//'.fits')

  end do



  deallocate(patches, flights, mask)

contains

  function find_az_limit(mjd, patchsize, az, el, centerpath)  
    implicit none

    real(dp),                       intent(in) :: mjd, az, el, patchsize
    real(dp),     dimension(:,:,:), intent(in) :: centerpath
    real(dp)                                   :: find_az_limit

    real(dp) :: daz, dist, vec1(3), vec2(3), az_patch, el_patch, az_out, grace
    logical(lgt) :: positive

    daz = 0.1d0 * DEG2RAD ! Accuracy
    grace = 3d0 * DEG2RAD

    az_patch = splint(centerpath(:,1,1), centerpath(:,2,1), centerpath(:,3,1), mjd)
    el_patch = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), mjd)
    positive = az_patch > az

    az_out = az
    call ang2vec(el_patch, az_patch, vec1)
    call ang2vec(el,       az_out,   vec2)
    call angdist(vec1, vec2, dist)
    !write(*,*) 'dist = ', dist, patchsize
    do while (dist < patchsize .or. ((az_patch-az)*(az_patch-az_out) > 0.d0))
       if (positive) then
          az_out = az_out + daz
       else
          az_out = az_out - daz
       end if
       call ang2vec(el,       az_out,  vec2)
       call angdist(vec1, vec2, dist)
    end do
!    write(*,*) az, az_out, az_patch
!    stop

    find_az_limit = az_out
!    if (positive) then
!       find_az_limit = az + 2*patchsize !az_out
!    else
!       find_az_limit = az - 2*patchsize
!    end if

  end function find_az_limit

  subroutine compute_az_profile(samprate, az_tele, az_limit, v, a, turn_time, az_profile)  
    implicit none

    real(dp), intent(in) :: samprate, az_tele, az_limit, v, a, turn_time
    real(dp), allocatable, dimension(:) :: az_profile

    integer(i4b) :: i, j, n
    real(dp)     :: s, t1, t2, dt, t, u, sign

    ! 1) s = a*t1^2 + v*t2 = abs(az_limit-az_tele)
    ! 2) a*t1 = v
    ! => t2 = (s - a * (v/a)^2) / v = (s - v^2/a) / v = s/v - v/a
    ! => t1 = sqrt((s - v*t2)/a)
    ! 
    ! if (t2 < 0) then
    !    s = a*t1^2 => sqrt(s/a)

    if (az_tele > az_limit) then
       sign = -1.d0
    else
       sign = 1.d0
    end if
    s  = abs(az_tele - az_limit)
    t2 = s/v - v/a
    if (t2 > 0) then
       t1 = sqrt((s - v*t2)/a)
    else
       t2 = 0.d0
       t1 = sqrt(s/a)
    end if

    dt = 1.d0 / samprate
    n  = (2*t1 + t2 + turn_time) / dt+1
    allocate(az_profile(0:n))
    t = 0.d0
    i = 1
    u = 0.d0
!    write(*,*) 
!    write(*,*) az_tele, az_limit
!    write(*,*) t1, t2
    az_profile(0) = az_tele
    do while (t < t1)
       u             = u + sign*a*dt
       az_profile(i) = az_profile(i-1) + u * dt
       i             = i+1
       t             = t+dt
    end do
    do while (t < t1+t2)
       az_profile(i) = az_profile(i-1) + u * dt
       i             = i+1
       t             = t+dt
    end do
    do while (t < t1+t2+t1)
       u             = u - sign*a*dt
       az_profile(i) = az_profile(i-1) + u * dt
       i             = i+1
       t             = t+dt
    end do
    do while (t < t1+t2+t1+turn_time)
       az_profile(i) = az_profile(i-1)
       i             = i+1
       t             = t+dt
    end do

!    open(68,file='az_profile.dat')
!    do i = 0, size(az_profile)-1
!       write(68,*), i/samprate, az_profile(i)
!    end do
!    close(68)
!    stop

  end subroutine compute_az_profile

  subroutine find_optimal_elevation(mjd, patch_size, centerpath, el_min, el_max, el_tele, cestime)
    implicit none

    real(dp),                   intent(inout) :: mjd
    real(dp),                   intent(in)    :: patch_size, el_min, el_max
    real(dp), dimension(:,:,:), intent(in)    :: centerpath
    real(dp),                   intent(out)   :: el_tele, cestime
    
    integer(i4b) :: i
    real(dp) :: el_low, el_high, t, dt, el, el0, dist, az, vec1(3), vec2(3), time, rising, el_old
    

!!$    dist = 1.d30
!!$
!!$    do i = 1, 10
!!$       el_high = (90.d0 - el_min) * DEG2RAD
!!$       el_low  = (90.d0 - el_max) * DEG2RAD
!!$       el_tele = el_low + rand_uni(handle) * (el_high-el_low)
!!$       el      = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), mjd)
!!$       az      = splint(centerpath(:,1,1), centerpath(:,2,1), centerpath(:,3,1), mjd)
!!$       call ang2vec(el,      az,  vec1)
!!$       call ang2vec(el_tele, az,  vec2)
!!$       call angdist(vec1, vec2, dist)
!!$       if (dist < patch_size) exit
!!$    end do
!!$    return

    el_high = (90.d0 - el_min) * DEG2RAD
    el_low  = (90.d0 - el_max) * DEG2RAD
    
    dt = 1.d0 / (24.d0*60.d0)
    t  = mjd

    ! Check whether patch is rising or setting
    el0 = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), t)
    el  = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), t+dt)
    if (el > el0) then
       el = min(el_high, el0+patch_size)
    else
       el = max(el_low, el0-patch_size)
    end if
    el_tele = el

    time = mjd
    el = el0
    el_old = el
    do while (abs(el-el_tele) < patch_size .and. (el_old-el0)*(el-el0) >= 0.d0)
       time = time + dt
       el_old = el
       el     = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), time)
    end do
    if (el_tele == el_high .or. el_tele == el_low) then
       cestime = 1d0 / 24.d0
    else
       Cestime = max(time-mjd,10.d0/60.d0/24.d0)
    end if

    write(*,*) el0, el_tele, cestime*24.d0
    return

    time = 0.d0
    el0 = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), t)
    el  = el0
    do while (abs(el-el0) < patch_size .and. el < el_high .and. el > el_low)
       t   = t + dt
       el = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), t)
    end do

    if (t - mjd < min_ces_time/(24.d0*60.d0)) then
       ! Fast forward to acceptable time
       do while (el < el_low .or. el > el_high)
          t = t + dt
          el = splint(centerpath(:,1,2), centerpath(:,2,2), centerpath(:,3,2), t)
       end do
       mjd = t
    end if

    ! Add a random fluctuation
    el_tele = el 

  end subroutine find_optimal_elevation
  
end program bscan




