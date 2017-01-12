program noise_tester
  use quiet_glitchy_noise_mod
  use quiet_hdf_mod
  use rngmod
  use quiet_mpi_mod
  use quiet_task_mod
  implicit none

  integer(i4b)       :: i, numsamp, nfreq, seed, s, numsim, numit, numglitch
  integer(i4b)       :: unit, myid, numprocs, ierr, root, tall, glength
  real(dp)           :: sigma0, alpha, fknee, var(3), skew(3), av(4), std(4)
  real(dp)           :: sigma0_in, alpha_in, fknee_in
  real(dp)           :: sigma0_guess, alpha_guess, fknee_guess
  real(dp)           :: nu, samprate, mu, sigma, t1, t2
  character(len=512) :: outfile, outprefix
  character(len=5) :: txtnum
  logical(lgt)       :: dosample
  type(planck_rng)   :: rng_handle
  type(task_list)    :: tasks
  logical(lgt), dimension(:),   allocatable :: mask
  real(sp),     dimension(:),   allocatable :: tod, powspec, templates
  complex(spc), dimension(:),   allocatable :: ft
  real(dp),     dimension(:,:), allocatable :: params, samleparams

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

  dosample=.true.
  numsim = 100
  glength = 30000
  numglitch = 1
  if(iargc() > 1) then
     call get_numsamp(numsamp,samprate)
  else
     numsamp  = 100000
     samprate = 50.d0
  end if
  nfreq  = numsamp/2+1

  seed      =  4767714  !104875 ! 4767714!374685 !4767714!5784598
  sigma0_in =  1.d-5
  alpha_in  = -1.8d0
  fknee_in  =  0.1d0
  sigma0_guess = 2.d0 * sigma0_in
  alpha_guess  = -2.d0!0.5d0 * alpha_in
  fknee_guess  = 2.d0 * fknee_in

  call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)
  seed = int(rand_uni(rng_handle)*100000.d0)
  call initialize_glitchy_noise_mod(seed_in=seed)

  allocate(tod(numsamp), mask(numsamp), ft(nfreq), powspec(nfreq), templates(numsamp))
  mask = .true.
!  mask(50000+1:50000+glength) = .false.

!  do i = 1, numglitch
!     tall = i*1000
!     mask(tall+1:tall+100) = .false.
!  end do

!  mask(1001:1100) = .false.
!  mask(2001:2100) = .false.
!  mask(3001:3100) = .false.
!  mask(4001:4100) = .false.
!  mask(5001:5100) = .false.
!  mask(6001:6100) = .false.
!  mask(7001:7100) = .false.
!  mask(8001:8100) = .false.
!  mask(9001:9100) = .false.
!  mask(9801:9900) = .false.

!  mask(8001:9000) = .false.

!  mask(4001:9000) = .false.
!  mask(34001:39000) = .false.
!  mask(18001:23000) = .false.

!  mask(30000:40000) = .false.
!  mask(50000:60000) = .false.
!  mask(70000:75000) = .false.
!  mask(85000:95000) = .false.
!  mask(numsamp/2:numsamp*2/3) = .false.
  allocate(params(11,numsim))
  params = 0.d0
  call getarg(1, outprefix)
  call int2string(glength, txtnum)
  !  outprefix = trim(outprefix)//'_samp'//txtnum
  !  outprefix = trim(outprefix)//'_nglitch'//txtnum
  !  outprefix = trim(outprefix)//'_glength'//txtnum
  outfile = trim(outprefix)//'_tasks.txt'
  call init_task_list(tasks, outfile, numsim, MPI_COMM_WORLD)
  do while (get_next_task(tasks, s))
     write(*,*) 'bim',int(s,i2b), int(myid, i2b),'= myid    '
!     write(*,*) myid, 'sim', s

     ! Create simulated TOD
     do i = 1, numsamp
        tod(i) = rand_gauss(rng_handle)
     end do
     call fft(tod, ft, 1)
     ft(1) = 0.d0
     do i = 2, nfreq
        nu = ind2freq(i, samprate, nfreq)
        ft(i) = ft(i) * sqrt((sigma0_in**2 * (1.d0 + (nu/fknee_in)**alpha_in)))
     end do
     call fft(tod, ft, -1)
     if(iargc() > 1) call add_stuff(tod, mask, templates)
    
     ! Estimate parameters
     sigma0 = sigma0_guess
     alpha  = alpha_guess
     fknee  = fknee_guess
     sigma0 = get_approximate_sigma0(tod)
     call wall_time(t1)
     call fit_1overf_profile_with_glitches(samprate, tod, mask, sigma0, alpha, &
          & fknee, numit, myid, s, dosample, var, skew, templates)
     call wall_time(t2)
     write(*,*) 'sim',int(s,i2b), int(myid, i2b),'= myid    ', real(sigma0,sp), real(alpha,sp), real(fknee,sp)

!     write(*,fmt='(a,e10.3,a,e10.3)') 'Estimated sigma0 = ', sigma0, ' -- true sigma0 = ', sigma0_in
!     write(*,fmt='(a,e10.3,a,e10.3)') 'Estimated alpha  = ', alpha,  ' -- true alpha  = ', alpha_in
!     write(*,fmt='(a,e10.3,a,e10.3)') 'Estimated fknee  = ', fknee,  ' -- true fknee  = ', fknee_in

     params(1,s) = sigma0
     params(2,s) = alpha
     params(3,s) = fknee
     params(4,s) = numit
     params(5,s) = var(1)
     params(6,s) = var(2)
     params(7,s) = var(3)
     params(8,s) = skew(1)
     params(9,s) = skew(2)
     params(10,s) = skew(3)
     params(11,s) = t2-t1
 end do
  allocate(samleparams(11,numsim))
  samleparams=0.d0
  call mpi_reduce(params, samleparams, size(params), mpi_double_precision, mpi_sum, 0, MPI_COMM_WORLD, ierr)
  if (myid==0) call output_noise(outprefix, samleparams)
  deallocate(tod, mask, ft, powspec, params, samleparams, templates)
  call mpi_finalize(ierr)

contains

  subroutine read_mystery_data(fname, mask, data, samprate)
    implicit none
    character(len=*)          :: fname
    type(hdf_file)            :: hfile
    real(dp)                  :: samprate
    integer(i4b)              :: nsamp, nsim, ext(7)
    integer(i4b), allocatable :: mask_in(:)
    logical(lgt), allocatable :: mask(:)
    real(sp),     allocatable :: data(:,:)
    call open_hdf_file(fname, hfile, "r")
    call get_size_hdf(hfile, "data", ext)
    nsamp = ext(1)
    nsim  = ext(2)
    allocate(mask(nsamp), mask_in(nsamp), data(nsamp,nsim))
    call read_hdf(hfile, "mask", mask_in)
    call read_hdf(hfile, "data", data)
    call read_hdf(hfile, "samprate", samprate)
    call close_hdf_file(hfile)
    mask = mask_in /= 0
    deallocate(mask_in)
  end subroutine

  subroutine add_stuff(tod, mask, tplate)
    implicit none
    real(sp),               intent(inout) :: tod(:)
    real(sp),     optional, intent(inout) :: tplate(:)
    logical(lgt), optional, intent(inout) :: mask(:)
    character(len=512)        :: filename
    type(hdf_file)            :: hfile
    integer(i4b)              :: numsamp
    real(sp),     allocatable :: templates(:,:)
    integer(i4b), allocatable :: mask_in(:)
    integer(i4b) :: i
    numsamp = size(tod)
    call getarg(2, filename)    
    call open_hdf_file(filename, hfile, "r")
    allocate(templates(numsamp,3))
    call read_hdf(hfile, "templates", templates)
    if (present(mask)) then
       allocate(mask_in(numsamp))
       call read_hdf(hfile, "mask", mask_in)
       mask = mask_in /= 0
       deallocate(mask_in)
    end if
    call close_hdf_file(hfile)
!    tod = tod + templates(:,2)
    if (present(tplate)) tplate = templates(:,2)
    deallocate(templates)
  end subroutine add_stuff

  subroutine get_numsamp(numsamp, samprate)
    implicit none
    integer(i4b), intent(out) :: numsamp
    real(dp),     intent(out) :: samprate
    character(len=512)        :: filename
    type(hdf_file)            :: hfile
    integer(i4b)              :: ext(7)
    call getarg(2, filename)    
    call open_hdf_file(filename, hfile, "r")
    call get_size_hdf(hfile, "data", ext)
    numsamp = ext(1)
    call read_hdf(hfile, "samprate", samprate)
    call close_hdf_file(hfile)
  end subroutine get_numsamp

  function get_approximate_sigma0(data) result(sigma0)
    implicit none
    real(sp)                  :: data(:), sigma0
    real(sp),     allocatable :: ps(:)
    complex(spc), allocatable :: ft(:)
    allocate(ps(size(data)/2+1),ft(size(data)/2+1))
    call fft(data,ft,1)
    call extract_powspec(ft,ps)
    sigma0 = mean(dble(ps(size(ps)/2:)))**0.5
   deallocate(ps,ft)
  end function

  subroutine output_noise(outprefix, samleparams)
    implicit none
    character(len=*), intent(in) :: outprefix
    real(dp)                     :: samleparams(:,:)
    character(len=512)           :: outfile
    integer(i4b)                 :: numsim, tumsim, np
    real(dp)                     :: av(4), std(4), tav, tstd
    numsim = size(samleparams(1,:))
    np=0
    if (dosample) np=8 !obs
    tumsim = numsim-np
    tav = 0.d0
    do s = np+1, numsim
       tav  = tav + samleparams(11,s)
    end do
    tav = tav/tumsim
    tstd = 0.d0
    do s = np+1, numsim
       tstd = tstd + (samleparams(11,s)-tav)**2
    end do
    tstd = tstd/tumsim
    tstd = sqrt(tstd)
    write(*,*) ''
    write(*,*) 'time', tav, tstd
    if (.not. dosample) then
       tav = median(samleparams(11,:))
       write(*,*) 'median time', tav
       outfile = trim(outprefix) // '_medtime.dat'
    else
       outfile = trim(outprefix) // '_avtime.dat'
    end if
    open(22, file = outfile)
    write(22,*) glength, tav 
    close(22)    
    av = 0.d0
    do s = 1, numsim
       av  = av + samleparams(1:4,s)
    end do
    av  = av/numsim
    std = 0.d0
    do s = 1, numsim
       std  = std  + (samleparams(1:4,s)-av)**2
    end do
    std = std/numsim
    std = sqrt(std)
    write(*,fmt='(a,4e16.6)') 'av ', av
    write(*,fmt='(a,4e16.6)') 'std', std
    write(*,fmt='(a,3e16.6)') 'sig', (av(1)-sigma0_in)/std(1), (av(2)-alpha_in)/std(2), (av(3)-fknee_in)/std(3)
    outfile = trim(outprefix) // '_log.txt' 
    open(23, file=outfile)
    write(23,*) 'time', tav, tstd
    write(23,*) 'numsim  =', numsim
    write(23,*) 'numsamp =', numsamp
    write(23,*) 'seed    =', seed
    write(23,fmt='(a,3e16.6)') ' guess', sigma0_guess, alpha_guess, fknee_guess
    write(23,*) '         SIGMA_0         ALPHA           F_KNEE          ITERATIONS'
    write(23,fmt='(a,3e16.6)') ' input', sigma0_in, alpha_in, fknee_in
    write(23,fmt='(a,4e16.6)') ' av   ', av
    write(23,fmt='(a,4e16.6)') ' std  ', std
    write(23,fmt='(a,3e16.6)') ' sigma', (av(1)-sigma0_in)/std(1), (av(2)-alpha_in)/std(2), (av(3)-fknee_in)/std(3)   
    open(22, file = trim(outprefix) // '_timing.dat')
    do s = 1, numsim
       write(22,*) samleparams(11,s)
    end do
    close(22)
    if (.not. dosample) then
       outfile = trim(outprefix) // '_sims.dat'
    else
       outfile = trim(outprefix) // '_sampsims.dat'
       samleparams(1,:) = (samleparams(1,:)-sigma0_in)/samleparams(5,:)
       samleparams(2,:) = (samleparams(2,:)-alpha_in)/samleparams(6,:)
       samleparams(3,:) = (samleparams(3,:)-fknee_in)/samleparams(7,:)
    end if
    open(22, file=outfile)!, recl=1024)
    do s = 1, numsim
       write(22,*) samleparams(1,s), samleparams(2,s), samleparams(3,s)!, samleparams(4,s)
    end do
    close(22)
    if (dosample) then
       av = 0.d0
       do s = 1, numsim
          av(1:3) = av(1:3) + samleparams(1:3,s)
       end do
       av(1:3) = av(1:3)/numsim
       std = 0.d0
       do s = 1, numsim
          std(1:3) = std(1:3) + (samleparams(1:3,s)-av(1:3))**2
       end do
       std = std/numsim
       std = sqrt(std)
       write(23,*) 'VAR      SIGMA_0         ALPHA           F_KNEE'
       write(23,fmt='(a,4e16.6)') ' mean ', av(1:3)
       write(23,fmt='(a,4e16.6)') ' var  ', std(1:3)
       write(*,*) 'VAR      SIGMA_0         ALPHA           F_KNEE'
       write(*,fmt='(a,4e16.6)') ' mean ', av(1:3)
       write(*,fmt='(a,4e16.6)') ' var  ', std(1:3)
       outfile = trim(outprefix) // '_var.dat'
       open(22, file=outfile)
       do s = 1, numsim
          write(22,*) samleparams(5,s), samleparams(6,s), samleparams(7,s)
       end do
       close(22)
       
       av = 0.d0
       do s = 1, numsim
          av(1:3) = av(1:3) + samleparams(8:10,s)
       end do
       av(1:3) = av(1:3)/numsim
       std = 0.d0
       do s = 1, numsim
          std(1:3) = std(1:3) + (samleparams(8:10,s)-av(1:3))**2
       end do
       std = std/numsim
       std = sqrt(std)
       write(23,*) 'SKEW      SIGMA_0         ALPHA           F_KNEE'
       write(23,fmt='(a,4e16.6)') ' mean ', av(1:3)
       write(23,fmt='(a,4e16.6)') ' std  ', std(1:3)
       write(*,*) 'SKEW      SIGMA_0         ALPHA           F_KNEE'
       write(*,fmt='(a,4e16.6)') ' mean ', av(1:3)
       write(*,fmt='(a,4e16.6)') ' std  ', std(1:3)
       outfile = trim(outprefix) // '_skew.dat'
       open(22, file=outfile)
       do s = 1, numsim
          write(22,*) samleparams(8,s), samleparams(9,s), samleparams(10,s)
       end do
       close(22)
    end if
    close(23)
    write(*,*)'Written to file = ',trim(outfile)
  end subroutine output_noise

end program noise_tester
