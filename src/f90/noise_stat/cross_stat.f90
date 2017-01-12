! This is almost the same as noise_stat, but computes the correlation per
! bin instead of the power spectrum, and also does not provide samples,
! just a total average due to memory considerations.
program cross_stat
  use quiet_utils
  use quiet_fileutils
  use quiet_system_mod
  use l2_fileutils
  use quiet_fft_mod
  use quiet_ces_mod
  use quiet_module_mod
  use quiet_mpi_mod
  use quiet_task_mod
  use quiet_acceptlist_mod
  implicit none

  character(len=512) :: lockfile, outfile, odir, parfile
  integer(i4b)       :: i, j, k, l, m, n, ierr, id, nproc, unit, cid
  integer(i4b)       :: nfreq, nces, nmod, ndi
  integer(i4b), parameter :: tmean = 1, tmed = 2, ntype = 2
  real(dp),     dimension(:,:,:,:,:), allocatable :: stats, ostats
  integer(i4b), dimension(:,:,:,:),   allocatable :: counts, ocounts
  real(dp), parameter :: medcorr = 0.6931547180d0
  type(task_list)    :: tasks

  call getarg(1, parfile)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

  call get_parameter(0, parfile, 'OUTPUT_DIR',      par_string=odir)
  call get_parameter(0, parfile, 'NUM_FREQ',        par_int=nfreq)

  lockfile = trim(odir) // "/lock.dat"
  outfile  = trim(odir) // "/corrstat.fits"
  call mkdirs(trim(lockfile), .true.)
  call mkdirs(trim(outfile),  .true.)

  call initialize_ces_mod(parfile)
  call initialize_module_mod(parfile)
  call initialize_acceptlist_mod(parfile)
  call init_task_list(tasks, lockfile, get_num_ces(), MPI_COMM_WORLD)

  nces = get_num_ces()
  nmod = get_num_modules()
  ndi  = get_num_diodes()

  if(id == 0) then
     unit = getlun()
     ierr = 0
     call rm(trim(outfile), noerr=.true.)
     call ftinit(unit, trim(outfile), 1, ierr)
     call ftiimg(unit, -64, 5, (/ nfreq, ndi, nmod, ndi, nmod /), ierr)
     call assert(ierr == 0, "Error opening '" // trim(outfile) // "' for writing!")
  end if

  allocate(stats(nfreq,ndi,nmod,ndi,nmod), counts(ndi,nmod,ndi,nmod))
  stats  = 0
  counts = 0

  do while(get_next_task(tasks, cid))
     call process_ces(stats, counts, cid)
  end do


  allocate(ostats(nfreq,ndi,nmod,ndi,nmod), ocounts(ndi,nmod,ndi,nmod))
  call mpi_reduce(stats, ostats, size(stats), mpi_double_precision, &
   & mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(counts, ocounts, size(counts), mpi_integer, &
   & mpi_sum, 0, mpi_comm_world, ierr)
  deallocate(stats, counts)

  ! Output in fits file
  if(id == 0) then
     do i = 1, size(ocounts,4)
        do j = 1, size(ocounts,3)
           do l = 1, size(ocounts,2)
              do m = 1, size(ocounts,1)
                 ostats(:,m,l,j,i) = ostats(:,m,l,j,i) / ocounts(m,l,j,i)
              end do
           end do
        end do
     end do
     do i = 1, size(ostats,1); ostats(i,1,1,1,1) = i; end do
     call ftpprd(unit, 1, 1, size(ostats), ostats, ierr)
     call ftclos(unit, ierr)
  end if

  deallocate(ostats, ocounts)
  call free_task_list(tasks)
  call mpi_finalize(ierr)
contains

  ! Read the ces, get the power spectrum, bin it and compute the stats
  ! for each diode
  subroutine process_ces(stats, counts, cid)
    implicit none
    real(dp)        :: stats(:,:,:,:,:)
    integer(i4b)    :: counts(:,:,:,:)
    integer(i4b)    :: cid, mod1, di,di1, mod2, di2, bin, i, i1, i2, j, k, n, nm, from, to
    integer(i4b)    :: mod
    real(dp)        :: corr
    type(quiet_ces) :: ces
    type(module_struct), dimension(:),     allocatable :: l2data
    complex(dp),         dimension(:,:,:), allocatable :: ft
    ces = get_ces_info(cid)
    if(.not. is_accepted(ces%run, ces%seg)) return
    write(*,fmt="(i3,a,i5,a,i5,a,i5)") id, " processing ces ", cid, &
     & " run ", ces%run, " seg ", ces%seg
    call l2_read(ces%filename, data=l2data)
    nm = size(l2data)
    n  = size(l2data(1)%time)
    m  = n/2+1
    ! First get all the fts, to avoid calculating them multiple times
    allocate(ft(m,ndi,nm))
    do i = 1, nm
       do di = 1, ndi
          mod = l2data(i)%module_number+1
          if(.not. is_accepted(ces%run, ces%seg, mod-1, di-1)) cycle
          call fft(l2data(i)%tod(di-1,:), ft(:,di,i), 1)
       end do
    end do
    ! Then calculate all the correlations per bin
    do i1 = 1, nm
       mod1 = l2data(i1)%module_number+1
       do di1 = 1, ndi
          if(.not. is_accepted(ces%run, ces%seg, mod1-1, di1-1)) cycle
          do i2 = 1, nm
             mod2 = l2data(i2)%module_number+1
             do di2 = 1, ndi
                if(.not. is_accepted(ces%run, ces%seg, mod2-1, di2-1)) cycle
                do bin = 1, nfreq
                   from = int(int(bin-1,i8b)*m/nfreq+1,i4b)
                   to   = int(int(bin,i8b)*m/nfreq,i4b)
                   call calc_corr_mean(corr, ft(from:to,di1,i1), ft(from:to,di2,i2))
                   stats(bin,di2,mod2,di1,mod1) = stats(bin,di2,mod2,di1,mod1) + corr
                end do
                counts(di2,mod2,di1,mod1) = counts(di2,mod2,di1,mod1) + 1
             end do
          end do
       end do
    end do
    deallocate(ft)
    call deallocate_module_struct(l2data)
  end subroutine

  subroutine calc_corr_mean(res, ft1, ft2)
    implicit none
    real(dp)     :: res, v1, v2
    complex(dp)  :: ft1(:), ft2(:)
    res = mean(dble(ft1*conjg(ft2)))/sqrt(mean(dble(ft1*conjg(ft1)))*mean(dble(ft2*conjg(ft2))))
  end subroutine

  subroutine calc_corr_med(res, ft1, ft2)
    implicit none
    real(dp)     :: res, v1, v2
    complex(dp)  :: ft1(:), ft2(:)
    v1  = median(dble(ft1*conjg(ft1)))/medcorr
    v2  = median(dble(ft2*conjg(ft2)))/medcorr
    res = corrmed2mean(median(dble(ft1*conjg(ft2)))/sqrt(v1*v2))
  end subroutine


end program
