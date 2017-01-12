! noise_stat produces a fits file with a primary array of
! type,module,diode,freq_bin,ces_ind, where type is 0 (mean) or
! 1 (corrected median). The power spectrum is binned linearly
! into nfreq bins, and the mean/med statistic is calculated for
! each ces-diode-bin.
program noise_stat
  use quiet_utils
  use quiet_fileutils
  use quiet_system_mod
  use l2_fileutils
  use quiet_fft_mod
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_mpi_mod
  use quiet_task_mod
  implicit none

  character(len=512) :: lockfile, outfile, odir, parfile
  integer(i4b)       :: i, j, k, l, m, n, ierr, id, nproc, unit, cid
  integer(i4b)       :: nfreq, nces, nmod, ndi
  integer(i4b), parameter :: tmean = 1, tmed = 2, ntype = 2
  real(dp), dimension(:,:,:,:,:), allocatable :: stats, ostats
  real(dp), parameter :: medcorr = 0.6931547180d0
  type(task_list)    :: tasks

  call getarg(1, parfile)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

  call get_parameter(0, parfile, 'OUTPUT_DIR',      par_string=odir)
  call get_parameter(0, parfile, 'NUM_FREQ',        par_int=nfreq)

  lockfile = trim(odir) // "/lock.dat"
  outfile  = trim(odir) // "/corrs.fits"
  call mkdirs(trim(lockfile), .true.)
  call mkdirs(trim(outfile),  .true.)

  call initialize_ces_mod(parfile)
  call initialize_module_mod(parfile)
  call init_task_list(tasks, lockfile, get_num_ces(), MPI_COMM_WORLD)

  nces = get_num_ces()
  nmod = get_num_modules()
  ndi  = get_num_diodes()

  if(id == 0) then
     unit = getlun()
     ierr = 0
     call rm(trim(outfile), noerr=.true.)
     call ftinit(unit, trim(outfile), 1, ierr)
     call ftiimg(unit, -64, 5, (/ nces, nfreq, ndi, nmod, ntype /), ierr)
     call assert(ierr == 0, "Error opening '" // trim(outfile) // "' for writing!")
  end if

  allocate(stats(nces,nfreq,ndi,nmod,ntype))
  stats = 0

  do while(get_next_task(tasks, cid))
     call process_ces(stats, cid)
  end do

  allocate(ostats(nces,nfreq,ndi,nmod,ntype))
  call mpi_reduce(stats, ostats, size(stats), mpi_double_precision, &
   & mpi_sum, 0, mpi_comm_world, ierr)
  deallocate(stats)

  ! Output in fits file
  if(id == 0) then
     call ftpprd(unit, 1, 1, size(ostats), ostats, ierr)
     call ftclos(unit, ierr)
  end if

  deallocate(ostats)
  call free_task_list(tasks)
  call mpi_finalize(ierr)
contains

  ! Read the ces, get the power spectrum, bin it and compute the stats
  ! for each diode
  subroutine process_ces(stats, cid)
    implicit none
    real(dp)        :: stats(:,:,:,:,:)
    integer(i4b)    :: cid, mod, di, bin, i, j, k, n, from, to
    type(quiet_ces) :: ces
    type(module_struct), dimension(:), allocatable :: l2data
    complex(dp), dimension(:), allocatable :: ft
    real(dp),    dimension(:), allocatable :: ps
    ces = get_ces_info(cid)
    write(*,fmt="(i3,a,i5,a,i5,a,i5)") id, " processing ces ", cid, &
     & " run ", ces%run, " seg ", ces%seg
    call l2_read(ces%filename, data=l2data)
    n = size(l2data(1)%time)
    m = n/2+1
    allocate(ft(m), ps(m))
    do i = 1, size(l2data)
       mod = l2data(i)%module_number+1
       do di = 1, ndi
          call fft(l2data(i)%tod(di-1,:), ft, 1)
          call extract_powspec(ft, ps)
          ps = ps / (median(ps(size(ps)/2:))/medcorr)
          do bin = 1, nfreq
             from = (bin-1)*m/nfreq+1
             to   = bin*m/nfreq
             call calc_stats(stats(cid,bin,di,mod,:), ps(from:to))
          end do
       end do
    end do
    deallocate(ft, ps)
    call deallocate_module_struct(l2data)
  end subroutine

  subroutine calc_stats(res, data)
    implicit none
    real(dp)     :: res(:), data(:)
    res(tmean) = mean(data)
    res(tmed)  = median(data)/medcorr
  end subroutine
end program
