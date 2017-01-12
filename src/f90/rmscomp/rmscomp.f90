program rmscomp
  use quiet_fileutils
  use quiet_mpi_mod
  use alm_tools
  implicit none

  character(len=512)   :: infile, outfile, param, beamfile
  integer(i4b)         :: i, j, k, l, m
  integer(i4b)         :: nreal, lmax, nside_in, nside_out, seed, myid, nprocs, ierr, root, ordering, npix_out, npix_in, nmaps
  real(dp)             :: fwhm
  real(dp),     dimension(:,:),   allocatable :: rms, outmap, map_in, map_out
  real(dp),     dimension(:,:),   allocatable :: beam_in, beam_out
  real(dp),     dimension(:,:),   pointer     :: weights, pixwin_in, pixwin_out
  complex(dpc), dimension(:,:,:), allocatable :: alms
  type(planck_rng)     :: rng_handle

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  root = 0

  if (iargc() /= 8) then
     if (myid == root) then
        write(*,*) 'Usage: rmscomp [rmsfile] [input beamfile] [output fwhm] [output lmax]'
        write(*,*) '                   [output nside] [seed] [number of realizations] [outfile]'
     end if
     call mpi_finalize(ierr)
     stop
  end if

  ! Read parameters
  call getarg(1,infile)
  call getarg(2,beamfile)
  call getarg(3,param)
  read(param,*) fwhm
  call getarg(4,param)
  read(param,*) lmax
  call getarg(5,param)
  read(param,*) nside_out
  call getarg(6,param)
  read(param,*) seed
  call getarg(7,param)
  read(param,*) nreal
  call getarg(8,outfile)
  npix_out = 12*nside_out**2

  call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)

  ! Read RMS file
  npix_in = getsize_fits(infile, nside=nside_in, nmaps=nmaps)
  call read_map(rms, ordering, infile)
  do i = 1, nmaps
     if (ordering == 2) call convert_nest2ring(nside_in, rms(:,i))
  end do

  ! Set up beams
  allocate(beam_in(0:lmax,nmaps), beam_out(0:lmax,nmaps))
  call read_beam(beamfile, beam_in)
  call gaussbeam(fwhm, lmax, beam_out)
  !do l = 0, lmax
  !   beam_out(l,1) = exp(-0.5d0*l*(l+1.d0) * (fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0)))**2)
  !end do
  call read_pixwin(nside_in,  nmaps, pixwin_in)
  call read_pixwin(nside_out, nmaps, pixwin_out)
  beam_in  = beam_in  * pixwin_in
  beam_out = beam_out * pixwin_out

  allocate(outmap(0:npix_out-1,nmaps), map_in(0:npix_in-1,nmaps), map_out(0:npix_out-1,nmaps), alms(nmaps,0:lmax,0:lmax))
  call read_ringweights(nside_in, nmaps==3, weights)
  outmap = 0.d0
  do i = 1+myid, nreal, nprocs
     write(*,*) 'Iteration = ', i, ' of ', nreal
     do k = 1, nmaps
        do j = 0, npix_in-1
           map_in(j,k) = rms(j,k) * rand_gauss(rng_handle)
        end do
     end do
     if (nmaps == 1) then
        call map2alm(nside_in, lmax, lmax, map_in(:,1), alms, [0.d0,0.d0], weights)
     else
        call map2alm(nside_in, lmax, lmax, map_in, alms, [0.d0,0.d0], weights)
     end if
     do k = 1, nmaps
        do l = 0, lmax
           if (beam_in(l,k) > 1.d-12) then 
              alms(k,l,0:l) = alms(k,l,0:l) / beam_in(l,k) * beam_out(l,k)
           end if
        end do
     end do
     if (nmaps == 1) then
        call alm2map(nside_out, lmax, lmax, alms, map_out(:,1))
     else
        call alm2map(nside_out, lmax, lmax, alms, map_out)
     end if
     outmap = outmap + map_out**2
  end do
  map_out = outmap
  call mpi_reduce(map_out, outmap, size(outmap), MPI_DOUBLE_PRECISION, MPI_SUM, &
        & root, MPI_COMM_WORLD, ierr)

  if (myid == root) then
     outmap = sqrt(outmap / nreal)
     do k = 1, nmaps
        if (ordering == 2) call convert_ring2nest(nside_out, outmap(:,k))
     end do
     call write_map(outmap, ordering, outfile)
  end if

  deallocate(beam_in, beam_out, outmap, map_out, map_in, alms, weights, rms)
  call mpi_finalize(ierr)

end program rmscomp
