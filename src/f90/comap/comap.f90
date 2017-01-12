program comap
  use comap_utils
  use rngmod
  use quiet_hdf_mod
  use quiet_fft_mod
  implicit none

  integer(i4b)     :: numfreq, seed
  real(dp)         :: t_tot, samprate
  type(planck_rng) :: rng_handle

  type tod_type
     real(dp)     :: samprate, Tsys
     integer(i4b) :: numsamp, numdet, numfreq
     real(dp)     :: fmin, fmax, df
     real(dp), allocatable, dimension(:)     :: t, f              ! (time or freq)
     real(dp), allocatable, dimension(:,:,:) :: d, d_raw, g, rms  ! (time,freq,det)
     real(dp), allocatable, dimension(:,:,:) :: point             ! (time,3,numdet)
  end type tod_type

  type map_type
     integer(i4b) :: n_x, n_y, numfreq, n_k
     real(dp)     :: dtheta
     real(dp), allocatable, dimension(:)     :: x, y, f, k                 ! (n_x or n_y or numfreq or n_k)
     real(dp), allocatable, dimension(:,:,:) :: m, rms, dsum, nhit, div    ! (n_x,n_y,numfreq)
     real(dp), allocatable, dimension(:,:)   :: Pk, Pk_rms, Pk_mu          ! (n_k,numfreq)
  end type map_type

  type corr_type
     integer(i4b) :: numfreq, nbin
     real(dp)     :: dtheta, theta_min, theta_max
     real(dp), allocatable, dimension(:)     :: theta
     real(dp), allocatable, dimension(:,:)   :: C, n, rms
  end type corr_type

  type(tod_type)  :: tod
  type(map_type)  :: map
  type(corr_type) :: corrfunc

  integer(i4b)       :: i, j, k
  real(dp)           :: RA_cent, dec_cent, patch_cut, mu, w, rms_scale
  character(len=512) :: pointfile, todfile, outprefix


  ! Define parameters and initialize basic data structures
  numfreq   = 256
  seed      = 82719
  rms_scale = 1.d0 !/ sqrt(365.d0) ! Scale to one year integration

  RA_cent   = 180.
  dec_cent  =  50.
  patch_cut =   5.

!!$  ! Test set 1
!!$  t_tot    = 7.d0/60.d0           ! in hours
!!$  samprate  = 10.d0          ! in Hz
!!$  pointfile = 'timestream0.txt'
!!$  todfile   = 'starburstPacketDump1044_0_mf.h5'
!!$  outprefix = 'test'

  ! Test set 1
  t_tot     = 23.d0          ! in hours
  samprate  = 50.d0          ! in Hz
  call getarg(1,pointfile)
  call getarg(2,todfile)
  call getarg(3,outprefix)
  !pointfile = 'Pointing.dat'
  !todfile   = 'starburstPacketDump1046_0_mf.h5'
  !outprefix = 'testbed1046'


  call rand_init(rng_handle, seed)


  ! Set up time-ordered data
  write(*,*) 'Setting up TOD'
  if (.true.) then
     call simulate_tod(samprate, t_tot, numfreq, rng_handle, pointfile, todfile, tod)
  end if
  call output_tod(trim(outprefix)//'_tod.dat', 1, tod)
  call output_allan_std(trim(outprefix)//'_allan.dat', 1, tod)


  ! Filter time-ordered data

  ! Produce maps
  write(*,*) 'Making maps'
  call compute_maps(tod, map)
  call output_map(trim(outprefix), map)

  ! Estimate and output mean
  mu = 0.d0
  w  = 0.d0

  do k = 1, map%numfreq
     do i = 1, map%n_x
        do j = 1, map%n_y
           if (map%rms(i,j,k) > 0.d0) then
              mu = mu + map%m(i,j,k)/map%rms(i,j,k)**2
              w  = w  + 1.d0/map%rms(i,j,k)**2
           end if
        end do
     end do
     mu = mu/w
     write(*,fmt='(a,i4,a,f8.3,a)') 'Monopole freq ', k, ' = ', mu, ' uK'
     !map%m(:,:,k) = map%m(:,:,k) - mu
  end do

  ! Compute two-point correlation function
  write(*,*) 'Computing two-point correlation function'
  call compute_2pt_corr(map, corrfunc)
  call output_2pt_corr(trim(outprefix)//'_corr.dat', corrfunc)

  ! Compute power spectrum
  call compute_powspec(rng_handle, map)
  call output_powspec(trim(outprefix)//'_Pk.dat', map)


  ! Compute excess variance
  call compute_excess_variance(outprefix, map)

contains


  subroutine simulate_tod(samprate, t_tot, numfreq, rng_handle, pointfile, todfile, tod)
    implicit none
    real(dp),         intent(in)    :: samprate, t_tot
    integer(i4b),     intent(in)    :: numfreq
    type(planck_rng), intent(inout) :: rng_handle
    character(len=*), intent(in)    :: pointfile, todfile
    type(tod_type),   intent(inout) :: tod

    integer(i4b) :: i, j, k, unit, ext(3), nu_cut, numsamp, nchan, dtau
    real(dp)     :: sigma0_n, fknee_n, alpha_n, t, tmp1, tmp2
    real(dp), allocatable, dimension(:,:)   :: C, corr
    real(dp), allocatable, dimension(:,:,:) :: buffer
    type(hdf_file)               :: file

    ! Allocate data structures
    unit         = getlun()
    tod%numsamp  = int(t_tot*3600.d0 * samprate)
    tod%numdet   = 2
    tod%numfreq  = numfreq
    tod%samprate = samprate
    tod%fmin     = 26.d9 ! Hz
    tod%fmax     = 34.d9 ! Hz
    tod%df       = (tod%fmax-tod%fmin)/tod%numfreq
    tod%Tsys     = 35.d0 / rms_scale ! K
    nu_cut       = 10000

    allocate(tod%t(tod%numsamp), tod%d(tod%numsamp,tod%numfreq, tod%numdet), &
         & tod%d_raw(tod%numsamp,tod%numfreq, tod%numdet), &
         & tod%f(tod%numfreq), tod%point(tod%numsamp,3,tod%numdet), &
         & tod%g(tod%numsamp,tod%numfreq, tod%numdet), &
         & tod%rms(tod%numsamp,tod%numfreq, tod%numdet))
    
    do i = 1, tod%numsamp
       tod%t(i) = i / real(samprate,dp)
    end do
    do i = 1, tod%numfreq
       tod%f(i) = tod%fmin + (i-0.5d0) * tod%df
    end do

    ! Read pointing from file
    open(unit,file=trim(pointfile))
    read(unit,*)
    do i = 1, tod%numsamp
       read(unit,*) t, tmp1, tmp2, tod%point(i,1,1), tod%point(i,2,1)
       !read(unit,*) t, tod%point(i,1,1), tod%point(i,2,1)
    end do
    close(unit)

    
    if (.false.) then
       ! Simulate noise internally
       sigma0_n = 1.d8   ! in uK_RJ
       do k = 1, tod%numdet
          do j = 1, tod%numfreq
             do i = 1, tod%numsamp
                tod%d_raw(i,j,k) = sigma0_n * rand_gauss(rng_handle)
             end do
             
             ! Estimate RMS in ADU
             tod%rms(:,j,k) = sqrt(variance(tod%d_raw(:,j,k)))

             ! Estimate gain in uK/ADU
             tod%g(:,j,k) = sqrt(variance(tod%d_raw(:,j,k))) / (tod%Tsys*1d6/sqrt(tod%df/tod%samprate))

          end do
       end do
       tod%d = tod%d_raw

    else
       ! Read time stream from HDF file
       call open_hdf_file(todfile, file, "r")
       call get_size_hdf(file, "auto_py", ext)
       allocate(buffer(ext(1),ext(2),ext(3)))
       call read_hdf(file, "auto_py", buffer)

       open(58,file='tod.dat')
       do i = 1, tod%numsamp
          write(58,*) tod%t(i), buffer(208,2,i)
       end do
       close(58)
       stop

       nchan = 1 !int(256/numfreq)
       do k = 1, tod%numdet
          do j = 1, tod%numfreq
             do i = 1, tod%numsamp
                tod%d_raw(i,j,k) = mean(buffer((j-1)*nchan+1:j*nchan,k,i))
             end do

             ! Apply Fourier filter
             tod%d(:,j,k) = tod%d_raw(:,j,k)
             call apply_filter_fft(nu_cut, tod%d(:,j,k))
             
             ! Estimate RMS in ADU
             tod%rms(:,j,k) = sqrt(variance(tod%d(:,j,k)))
             
             ! Estimate gain in uK/ADU
             tod%g(:,j,k) = sqrt(variance(tod%d(:,j,k))) / (tod%Tsys*1d6/sqrt(tod%df/tod%samprate))
          end do
       end do
       call close_hdf_file(file)
       deallocate(buffer)

    end if

    ! Compute correlation matrix
    allocate(C(2*tod%numfreq,2*tod%numfreq), corr(2*tod%numfreq,2*tod%numfreq))
    do i = 1, tod%numfreq
       write(*,*) i
       do j = 1, tod%numfreq
          C(i,j)         = mean(tod%d(1:10*36000,i,1)*tod%d(1:10*36000,j,1))
          C(256+i,j)     = mean(tod%d(1:10*36000,i,2)*tod%d(1:10*36000,j,1))
          C(i,256+j)     = mean(tod%d(1:10*36000,i,1)*tod%d(1:10*36000,j,2))
          C(256+i,256+j) = mean(tod%d(1:10*36000,i,2)*tod%d(1:10*36000,j,2))
          !C(j,i) = C(i,j)
       end do
    end do

    do i = 1, 2*tod%numfreq
       do j = 1, 2*tod%numfreq
          corr(i,j) = C(i,j) / sqrt(C(i,i)*C(j,j))
       end do
    end do

    open(58,file='corrmat.dat',recl=100000)
    do i = 1, 2*tod%numfreq
       write(58,*) real(corr(i,:),sp)
    end do
    close(58)
    stop
    

!!$    ! Estimate RMS as a function of time
!!$    dtau = 5000
!!$    open(58,file='rms_vs_t.dat', recl=1024)
!!$    do i = 1, int(tod%numsamp/dtau)
!!$       write(58,*) (i-0.5)*dtau/tod%samprate, sqrt(variance(tod%d((i-1)*dtau+1:i*dtau,1,1))), &
!!$            & sqrt(variance(tod%d((i-1)*dtau+1:i*dtau,2,1))), &
!!$            & sqrt(variance(tod%d((i-1)*dtau+1:i*dtau,3,1)))
!!$    end do
!!$    close(58)
!!$    stop

!!$    open(85,file='tod_raw.dat', recl=100000)
!!$    do i = 1, 180000
!!$       write(85,*) tod%t(i), real(tod%d_raw(i,:,1),sp)
!!$    end do
!!$    close(85)
!!$
!!$
!!$    open(85,file='tod_filt.dat', recl=100000)
!!$    do i = 1, 180000
!!$       write(85,*) tod%t(i), real(tod%d(i,:,1),sp)
!!$    end do
!!$    close(85)
!!$
!!$    stop

    ! Add signal


    ! Select relevant data segments
    numsamp = 0
    do i = 1, tod%numsamp
       if (tod%point(i,1,1) > RA_cent-patch_cut  .and. tod%point(i,1,1) < RA_cent+patch_cut .and. &
           tod%point(i,2,1) > dec_cent-patch_cut .and. tod%point(i,2,1) < dec_cent+patch_cut) then
          numsamp = numsamp + 1
          tod%t(numsamp)         = tod%t(i)
          tod%point(numsamp,:,:) = tod%point(i,:,:)
          tod%d(numsamp,:,:)     = tod%d(i,:,:)
          tod%d_raw(numsamp,:,:) = tod%d_raw(i,:,:)
          tod%g(numsamp,:,:)     = tod%g(i,:,:)
          tod%rms(numsamp,:,:)   = tod%rms(i,:,:)
       end if
    end do
    write(*,*) '   Total number of accepted samples = ', numsamp, ' of ', tod%numsamp
    tod%numsamp = numsamp

  end subroutine simulate_tod


  subroutine output_tod(filename, det, tod)
    implicit none
    character(len=*), intent(in) :: filename
    integer(i4b),     intent(in) :: det
    type(tod_type),   intent(in) :: tod

    integer(i4b) :: unit, i, j

    unit = getlun()

    open(unit, file=trim(filename), recl=4096)
    do i = 1, tod%numsamp
       write(unit,fmt='(f16.8)',advance='no') tod%t(i)
       do j = 1, tod%numfreq
          write(unit,fmt='(2f24.8)',advance='no') tod%d(i,j,det), tod%d_raw(i,j,det)
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine output_tod


  subroutine output_allan_std(filename, det, tod)
    implicit none
    character(len=*), intent(in) :: filename
    integer(i4b),     intent(in) :: det
    type(tod_type),   intent(in) :: tod

    real(dp)     :: tau, tau_min, tau_max
    integer(i4b) :: i, j, k, n_t, dt, unit, n
    real(dp), allocatable, dimension(:)   :: allan
    real(dp), allocatable, dimension(:,:) :: d

    unit    = getlun()
    tau_min = max(1d-3, 2.d0/tod%samprate)
    tau_max = min(1d3,  tod%numsamp/tod%samprate)
    n_t     = 20

    open(unit,file=trim(filename), recl=4096)
    do i = 1, n_t-1
       tau = tau_min * (tau_max/tau_min)**(real(i-1,dp)/real(n_t-1,dp))
       dt  = tau * tod%samprate
       n   = tod%numsamp / dt
       allocate(d(n,tod%numfreq), allan(tod%numfreq))
       do k = 1, tod%numfreq
          do j = 1, n
             d(j,k) = mean(tod%d((j-1)*dt+1:j*dt,k,det))
          end do
          allan(k) = sqrt(variance(d(:,k)))
       end do
       write(unit,*) tau, allan
       deallocate(d, allan)
    end do
    close(unit)

  end subroutine output_allan_std



  subroutine compute_maps(tod, map)
    implicit none
    type(tod_type),   intent(in)    :: tod
    type(map_type),   intent(inout) :: map

    integer(i4b) :: i, j, k, p, q
    real(dp)     :: x_min, x_max, y_min, y_max, pad, chisq, nu

    ! Set up map grid
    pad         = 0.3d0       ! degrees
    map%numfreq = tod%numfreq
    !map%dtheta = 1.d0 * pi / 180.d0 / 60.d0  ! Arcmin -> radians
    map%dtheta = 5.d0 / 60.d0  ! Arcmin
    x_min = minval(tod%point(1:tod%numsamp,1,:))-pad; x_max = maxval(tod%point(1:tod%numsamp,1,:))+pad; map%n_x = (x_max-x_min)/map%dtheta+1
    y_min = minval(tod%point(1:tod%numsamp,2,:))-pad; y_max = maxval(tod%point(1:tod%numsamp,2,:))+pad; map%n_y = (y_max-y_min)/map%dtheta+1
    allocate(map%x(map%n_x), map%y(map%n_y))
    do i = 1, map%n_x
       map%x(i) = x_min + (i-1)*map%dtheta
    end do
    do i = 1, map%n_y
       map%y(i) = y_min + (i-1)*map%dtheta
    end do

    ! Set up map structures
    if (.not. allocated(map%dsum)) then
       allocate(map%m(map%n_x,map%n_y,map%numfreq), map%dsum(map%n_x,map%n_y,map%numfreq), &
            & map%nhit(map%n_x,map%n_y,map%numfreq), map%div(map%n_x,map%n_y,map%numfreq), &
            & map%rms(map%n_x,map%n_y,map%numfreq))
       map%dsum = 0.d0
       map%nhit = 0.d0
       map%div  = 0.d0
    end if

    ! Co-add into maps
    do k = 1, tod%numdet
       do i = 1, tod%numsamp
          p = min(max(int((tod%point(i,1,k)-x_min)/map%dtheta),1),map%n_x)
          q = min(max(int((tod%point(i,2,k)-y_min)/map%dtheta),1),map%n_y)
          do j = 1, tod%numfreq
             map%dsum(p,q,j) = map%dsum(p,q,j) + tod%g(i,j,k)    / tod%rms(i,j,k)**2 * tod%d(i,j,k) 
             map%div(p,q,j)  = map%div(p,q,j)  + tod%g(i,j,k)**2 / tod%rms(i,j,k)**2
             map%nhit(p,q,j) = map%nhit(p,q,j) + 1.d0
          end do
       end do
    end do
    where (map%nhit > 0) 
       map%m   = map%dsum / map%div
       map%rms = 1.d0 / sqrt(map%div)
    elsewhere
       map%m   = 0.d0
       map%rms = 0.d0
    end where

!!$    ! Replace with Gaussian
!!$    do k = 1, map%numfreq
!!$       do i = 1, map%n_x
!!$          do j = 1, map%n_y
!!$             map%m(i,j,k) = map%rms(i,j,k) * rand_gauss(rng_handle)
!!$          end do
!!$       end do
!!$    end do

    ! Report reduced chisquares
    do i = 1, map%numfreq
       nu    = count(map%nhit(:,:,i)>0)
       chisq = sum((map%m(:,:,i)/map%rms(:,:,i))**2,map%nhit(:,:,i)>0)
       write(*,fmt='(a,i5,a,f8.3,a,f8.3)') 'Freq = ', i, ', red chisq = ', chisq/nu, ', sigma = ', (chisq-nu)/sqrt(2.d0*nu)
    end do

  end subroutine compute_maps

  subroutine output_map(prefix, map)
    implicit none
    character(len=*), intent(in) :: prefix
    type(map_type),   intent(in) :: map

    integer(i4b)       :: i, j, k, unit
    character(len=4)   :: itext
    character(len=512) :: filename

    unit = getlun()
    do i = 1, map%numfreq
       call int2string(i, itext)
       filename = trim(prefix) //'_freq'// itext // '_map.dat'
       open(unit,file=trim(filename), recl=100000)
       write(unit,*) '# n_x = ', map%n_x
       write(unit,*) '# n_y = ', map%n_y
       write(unit,*) '# x   = ', real(map%x,sp)
       write(unit,*) '# y   = ', real(map%y,sp)
       do j = 1, map%n_x
          do k = 1, map%n_y
             write(unit,fmt='(f16.8)',advance='no') map%m(j,k,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

    unit = getlun()
    do i = 1, map%numfreq
       call int2string(i, itext)
       filename = trim(prefix) //'_freq'// itext // '_rms.dat'
       open(unit,file=trim(filename), recl=100000)
       write(unit,*) '# n_x = ', map%n_x
       write(unit,*) '# n_y = ', map%n_y
       write(unit,*) '# x   = ', real(map%x,sp)
       write(unit,*) '# y   = ', real(map%y,sp)
       do j = 1, map%n_x
          do k = 1, map%n_y
             write(unit,fmt='(f16.8)',advance='no') map%rms(j,k,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

    unit = getlun()
    do i = 1, map%numfreq
       call int2string(i, itext)
       filename = trim(prefix) //'_freq'// itext // '_nhit.dat'
       open(unit,file=trim(filename), recl=100000)
       write(unit,*) '# n_x = ', map%n_x
       write(unit,*) '# n_y = ', map%n_y
       write(unit,*) '# x   = ', real(map%x,sp)
       write(unit,*) '# y   = ', real(map%y,sp)
       do j = 1, map%n_x
          do k = 1, map%n_y
             write(unit,fmt='(f16.8)',advance='no') map%nhit(j,k,i)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

  end subroutine output_map


  subroutine compute_2pt_corr(map, corrfunc)
    implicit none
    type(map_type),   intent(in)  :: map
    type(corr_type),  intent(out) :: corrfunc

    integer(i4b) :: i, j, k, x1, x2, y1, y2
    real(dp)     :: b

    corrfunc%numfreq = map%numfreq
    corrfunc%dtheta  = map%dtheta
    corrfunc%nbin    = int(sqrt(real(map%n_x**2+map%n_y**2,dp)))
    
    allocate(corrfunc%C(corrfunc%nbin,corrfunc%numfreq), corrfunc%theta(corrfunc%nbin), &
         & corrfunc%n(corrfunc%nbin,corrfunc%numfreq), corrfunc%rms(corrfunc%nbin,corrfunc%numfreq))
    do i = 1, corrfunc%nbin
       corrfunc%theta(i) = (i-1)*corrfunc%dtheta
    end do

    ! Compute correlation functions
    corrfunc%C = 0.d0
    corrfunc%n = 0.d0
    do x1 = 1, map%n_x
       do y1 = 1, map%n_y
          if (map%nhit(x1,y1,1) == 0) cycle

          do x2 = 1, map%n_x
             do y2 = 1, map%n_y
                if (map%nhit(x2,y2,1) == 0) cycle
                
                b = min(int(sqrt(real(x1-x2,dp)**2 + real(y1-y2,dp)**2))+1,corrfunc%nbin)

                do i = 1, map%numfreq
                   corrfunc%C(b,i) = corrfunc%C(b,i) + (map%m(x1,y1,i)/map%rms(x1,y1,i)**2)*&
                        & (map%m(x2,y2,i)/map%rms(x2,y2,i)**2)
                   corrfunc%n(b,i) = corrfunc%n(b,i) + 1.d0 / (map%rms(x1,y1,i)**2 * map%rms(x2,y2,i)**2)
                end do
             end do
          end do

       end do
    end do
    do b = 1, corrfunc%nbin
       do i = 1, corrfunc%numfreq
          if (corrfunc%n(b,i) > 0.d0) then
             corrfunc%C(b,i) = corrfunc%C(b,i) / corrfunc%n(b,i) 
          end if
       end do
    end do

  end subroutine compute_2pt_corr


  subroutine output_2pt_corr(filename, corrfunc)
    implicit none
    character(len=*), intent(in) :: filename
    type(corr_type),   intent(in) :: corrfunc

    integer(i4b)       :: i, j, k, unit

    unit = getlun()
    open(unit,file=trim(filename), recl=100000)
    do i = 1, corrfunc%nbin
       write(unit,fmt='(f8.3)',advance='no') corrfunc%theta(i)
       do j = 1, corrfunc%numfreq
          write(unit,fmt='(f16.8)',advance='no') corrfunc%C(i,j)
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine output_2pt_corr


  subroutine compute_powspec(handle, map)
    implicit none
    type(planck_rng), intent(inout)  :: handle
    type(map_type),   intent(inout)  :: map

    integer(i4b) :: i

    map%n_k = 20
    allocate(map%k(map%n_k), map%Pk(map%n_k,map%numfreq), map%Pk_rms(map%n_k,map%numfreq), &
         & map%Pk_mu(map%n_k,map%numfreq))

    do i = 1, map%numfreq
       call compute_2d_powspec(handle, map%dtheta, map%m(:,:,i), map%rms(:,:,i), &
            & map%k, map%Pk(:,i), map%Pk_rms(:,i), map%Pk_mu(:,i))
    end do

  end subroutine compute_powspec

   subroutine compute_2d_powspec(handle, dtheta, map, rms, ks, Pk, Pk_rms, Pk_mu)
     implicit none
     type(planck_rng),           intent(inout) :: handle
     real(dp),                   intent(in)    :: dtheta
     real(dp), dimension(1:,1:), intent(in)    :: map, rms
     real(dp), dimension(1:),    intent(out)   :: ks, Pk, Pk_rms, Pk_mu

     integer(i4b) :: i, j, b, m, n, nk, numsim, s
     real(dp)     :: dk, kx, ky, kmin, kmax, k, rms_min, area
     integer*8    :: plan
     real(dp),     allocatable, dimension(:)   :: w
     real(dp),     allocatable, dimension(:,:) :: invN_map, sims
     complex(dpc), allocatable, dimension(:,:) :: fk

     numsim = 100
     
     m = size(map,1)
     n = size(map,2)
     allocate(invN_map(m,n), fk(m/2+1,n))

     ! Set up grid
     nk   = size(ks)
     kmin = 0.d0
     kmax = 2.d0/dtheta
     dk   = (kmax-kmin)/nk
     do i = 1, nk
        ks(i) = (i-0.5d0) *dk
     end do

     ! Inverse noise variance weighting and masking
     rms_min = minval(rms, rms > 0.d0)
     where (rms /= 0.d0)
        invN_map = map * rms_min**2 / rms**2
     elsewhere
        invN_map = 0.d0
     end where

     ! Compute 2D FFT
     call dfftw_plan_dft_r2c_2d(plan, m, n, invN_map, fk, fftw_estimate + fftw_unaligned)
     call dfftw_execute_dft_r2c(plan, invN_map, fk)
     fk = fk/sqrt(real(size(invN_map),dp))

     ! Compute power spectrum
     allocate(w(nk))
     Pk = 0.d0
     w  = 0.d0
     do i = 1, m/2+1
        kx = 2.d0/dtheta * real(i-1,dp)/ real(m,dp)
        do j = 1, n
           if (j <= n/2) then
              ky = 2.d0/dtheta * real(j-1,dp)/ real(n,dp)
           else
              ky = 2.d0/dtheta * real(n-j+1,dp)/ real(n,dp)
           end if
           k     = sqrt(kx**2 + ky**2)
           b     = min(max(int(k/dk)+1,1),nk)
           Pk(b) = Pk(b) + abs(fk(i,j))**2
           w(b)  = w(b)  + 1.d0
        end do
     end do

     ! Normalize
     where (w > 0.d0)
        Pk = Pk / w
     elsewhere
        Pk = 0.d0
     end where


     if (numsim > 0) then

        allocate(sims(nk,numsim))
        sims = 0.d0

        do s = 1, numsim

           ! Draw white noise simulated map
           do i = 1, m
              do j = 1, n
                 invN_map(i,j) = rms(i,j) * rand_gauss(handle)
              end do
           end do
     
           ! Inverse noise variance weighting and masking
           where (rms /= 0.d0)
              invN_map = invN_map * rms_min**2 / rms**2
           elsewhere
              invN_map = 0.d0
           end where

           ! Compute 2D FFT
           call dfftw_execute_dft_r2c(plan, invN_map, fk)
           fk = fk/sqrt(real(size(invN_map),dp))

           ! Compute power spectrum
           do i = 1, m/2+1
              kx = 2.d0/dtheta * real(i-1,dp)/ real(m,dp)
              do j = 1, n
                 if (j <= n/2) then
                    ky = 2.d0/dtheta * real(j-1,dp)/ real(n,dp)
                 else
                    ky = 2.d0/dtheta * real(n-j+1,dp)/ real(n,dp)
                 end if
                 k  = sqrt(kx**2 + ky**2)
                 if (k > kmin .and. k < kmax) then
                    b         = min(max(int(k/dk)+1,1),nk)
                    sims(b,s) = sims(b,s) + abs(fk(i,j))**2
                 end if
              end do
           end do

           ! Normalize
           where (w > 0.d0)
              sims(:,s) = sims(:,s) / w
           elsewhere
              sims(:,s) = 0.d0
           end where

        end do

        ! Debias spectrum, and compute error bars
        do i = 1, nk
           Pk_mu(i)  = mean(sims(i,:))
           !Pk(i)     = Pk(i) - Pk_mu(i)
           Pk_rms(i) = sqrt(variance(sims(i,:)))
        end do

        deallocate(sims)
     end if

     call dfftw_destroy_plan(plan)

     deallocate(invN_map, fk, w)

   end subroutine compute_2d_powspec



  subroutine output_powspec(filename, map)
    implicit none
    character(len=*), intent(in) :: filename
    type(map_type),   intent(in) :: map

    integer(i4b)       :: i, j, k, unit

    unit = getlun()
    open(unit,file=trim(filename), recl=100000)
    do i = 1, map%n_k
       if (any(map%Pk(i,:) /= 0.d0)) then
          write(unit,fmt='(f8.3)',advance='no') map%k(i)
          do j = 1, map%numfreq
             write(unit,fmt='(3f20.12)',advance='no') map%Pk(i,j), map%Pk_rms(i,j), map%Pk_mu(i,j)
          end do
          write(unit,*)
       end if
    end do
    close(unit)

  end subroutine output_powspec


  subroutine compute_excess_variance(prefix, map)
    implicit none
    character(len=*),              intent(in)  :: prefix
    type(map_type),                intent(in)  :: map

    integer(i4b) :: i, j, na, unit1, unit2
    character(len=4) :: jtext
    real(dp)     :: amin, amax, da, mu, sigma
    real(dp), allocatable, dimension(:) :: a, Pa

    unit1 = getlun()
    na = 101
    allocate(a(na), Pa(na))

    open(unit1,file=trim(prefix)//'_a.dat')
    unit2 = getlun()
    do j = 1, map%numfreq

       amin = -0.99d0*minval(map%rms(:,:,j), map%rms(:,:,j) > 0.d0)
       amax = -amin
       da   = (amax-amin)/(na-1.d0)
       do i = 1, na
          a(i) = amin + (i-1)*da
       end do

       call compute_lnL_A(map%m(:,:,j), map%rms(:,:,j), a, Pa)

!!$       do i = 1, na
!!$          if (a(i) < 0.d0) then
!!$             a(i) = -sqrt(-a(i))
!!$          else
!!$             a(i) = sqrt(a(i))
!!$          end if
!!$       end do
       
       ! Report mean and standard deviation
       mu    = sum(a*Pa) * da
       sigma = sqrt(sum((a-mu)**2*Pa) * da)
       write(*,fmt='(a,f16.3,a,f16.3)') '   A = ', mu, ' +/- ', sigma
       write(unit1,*) j, mu, sigma
       
       call int2string(j, jtext)
       open(unit2,file=trim(outprefix)//'_freq'//jtext//'_Pa.dat')
       do i = 1, na
          write(unit2,*) a(i), Pa(i)
       end do
       close(unit2)
    end do
    close(unit1)

    deallocate(a, Pa)

  end subroutine compute_excess_variance

   subroutine compute_lnL_A(map, rms, a, Pa)
     implicit none
     real(dp), dimension(1:,1:), intent(in)    :: map, rms
     real(dp), dimension(1:),    intent(in)    :: a
     real(dp), dimension(1:),    intent(out)   :: Pa

     integer(i4b) :: i, j, k, m, n
     real(dp)     :: da, mu, sigma, var

     m = size(map,1)
     n = size(map,2)

     Pa = 0.d0
     do k = 1, size(a)
        if (a(k) >= 0.d0) then
           var = a(k)**2
        else
           var = -a(k)**2
        end if
        do i = 1, m
           do j = 1, n
              if (rms(i,j) > 0.d0) then
                 Pa(k) = Pa(k) - 0.5d0 * (map(i,j)**2 / (var+rms(i,j)**2) + log(var+rms(i,j)**2))
              else
                 Pa(k) = Pa(k) - 1000.d0
              end if
           end do
        end do
     end do

     ! Exponentiate and normalize
     da = a(2)-a(1)
     Pa = exp(Pa-maxval(Pa))
     Pa = Pa / (sum(Pa)*da)

   end subroutine compute_lnL_A

end program comap
