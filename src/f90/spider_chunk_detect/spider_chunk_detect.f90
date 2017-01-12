program spider_chunk_detect
  use healpix_types
  use quiet_fileutils
  use quiet_pointing_mod
  use quiet_module_mod
  use quiet_shared_output_mod
  use quiet_patch_detect_mod
  use quiet_hdf_mod
  implicit none

  type l1struct
     integer(i4b) :: nsamp, ndi
     real(sp),     allocatable, dimension(:,:) :: point, tod
     integer(i4b), allocatable, dimension(:,:) :: flag
     real(dp),     allocatable, dimension(:)   :: time
  end type l1struct

  real(dp)             :: max_az_dev, max_el_dev, max_dk_dev
  real(dp)             :: rst_el_max_dev
  real(dp)             :: az_stall_lim, el_stall_lim
  real(dp)             :: ces_delay, ces_mindur, cas_mindur, rst_mindur
  real(dp)             :: rst_az_tol, rst_az_stall, rst_az_amp_lim, samprate, t_tot
  integer(i4b)         :: target_dt, chunk_size, glitchtol, i, j, k, samples, n, isun, unit
  integer(i4b)         :: subchunk_size, az_stall_timeout, el_stall_timeout, ierr, nproc, myid
  integer(i4b)         :: ind1, ind2

  integer(i4b)         :: foo, numfile, ext(7), nsamp, ndi, ngood
  logical(lgt)         :: object_detect, ok
  character(len=512)   :: listfile, parfile, l1prefix, output_dir, range_ofilename, filename
  character(len=512)   :: pswitch_file
  real(dp),           dimension(:), allocatable :: pswitch
  character(len=512), dimension(:), allocatable :: filelist
  type(l1struct)       :: data
  type(hdf_file)       :: file

  unit = getlun()
  call getarg(1, parfile)

  ! These are implementation details that probably don't need to be
  ! configurable.
  chunk_size  = 10000 ! 10000 frames = 100 s
  target_dt   = 10    ! 10 ms
  glitchtol   = 300   ! ignore mode /= 3 glitches shorter than this

  call get_parameter(0, parfile, "LEVEL1_FILELIST", par_string=listfile, &
       & desc="List of level1-files to consider, relative to LEVEL1_DIR. Fills the same role" // &
       & " as L1_DATABASE, but does not rely on run/seg/obj classification from Chicago.")
  call get_parameter(0, parfile, "MIN_DUR", par_dp=ces_mindur, &
       & desc="Minimum duration of a ces, in seconds. 100 is sensible.")
  call get_parameter(0, parfile, "OUTPUT_DIR",    par_string=output_dir, &
       & desc="Directory where the output files will be placed.")
  call get_parameter(0, parfile, "CES_DELAY", par_dp=ces_delay, &
       & desc="Number of seconds to cut out at the start of the CES.")

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid,  ierr)

  call initialize_quiet_hdf_mod

  open(unit,file=trim(listfile))
  read(unit,*) numfile
  do i = 1, numfile
     read(unit,*) filename
     
     call open_hdf_file(filename, file, "r")
     call get_size_hdf(file, "tod", ext)
     data%nsamp = ext(1); ndi = ext(2)
     allocate(data%time(data%nsamp), data%tod(data%nsamp,ndi), data%flag(data%nsamp,ndi), &
          & data%point(data%nsamp,3))
     call read_hdf(file, "samprate",   samprate)
     call read_hdf(file, "time",       data%time)
     call read_hdf(file, "orig_point", data%point)
     call read_hdf(file, "tod",        data%tod)
     call read_hdf(file, "flag",       data%flag)     

     ! Eliminate flag = 2 => cosmic ray
     where (data%flag == 2)
        data%flag = 0
     end where

!!$              open(58,file='tod2.dat')
!!$              do k = 1, 5000000
!!$                 write(58,*) k, data%tod(k,1), data%flag(k,1)
!!$              end do
!!$              close(58)
!!$              call mpi_finalize(ierr)
!!$              stop

     ngood = 0
     t_tot = 0.d0
     do j = 1, ndi

        ! Search for good chunks
        ind1 = 1
        open(58,file='tod.dat')
        do while (ind1 < data%nsamp)
           
           ! Skip first CES_DELAY seconds
           if (all(data%flag(ind1:min(ind1+int(ces_delay*samprate),data%nsamp),j) == 0)) then
              ind1 = ind1 + int(ces_delay*samprate)
           else
              ! Drop segment if it's not clean
              ind1 = ind1 + int(ces_delay*samprate)
              cycle
           end if

           ! Find non-flagged region
           do while (data%flag(ind1,j) > 0 .and. ind1 < data%nsamp)
              ind1 = ind1+1
           end do
           ind2 = ind1+1
           do while (data%flag(ind2,j) == 0 .and. ind2 < data%nsamp)
              ind2 = ind2+1
           end do

           ! Drop segment if it's too short
           if (ind2-ind1 < ces_mindur*samprate) then
              ind1 = ind2+1
              cycle
           end if

           ! Output good segment
           ngood = ngood+1
           t_tot = t_tot + int((ind2-ind1)/samprate)/3600.d0
           write(*,*) 'ngood = ', ngood, ind1, ind2, int((ind2-ind1)/samprate), real(t_tot,sp)

           if (ngood == 805) then
              do k = ind1, ind2
                 write(58,*) k, data%tod(k,j)
              end do
           end if
           write(58,*)
!!$              call mpi_finalize(ierr)
!!$              stop
!           end if

           ! Prepare for next chunk
           ind1 = ind2+1

        end do
        write(*,fmt='(a,i5,a,i5,a,f6.2)') 'det = ', j, ' -- ngood = ', ngood, ', good frac = ', &
             & t_tot / (data%nsamp/samprate/3600.d0)
        close(58)

     end do




!     do j = 1+myid, ndi, numprocs

!     end do
     call close_hdf_file(file)

     deallocate(data%time, data%tod, data%flag, data%point)
  end do
  close(unit)

  call cleanup_quiet_hdf_mod

  call mpi_finalize(ierr)


contains

end program spider_chunk_detect
