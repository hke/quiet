 program modulation
  use modulation_evidence
  use modulation_mcmc
  use modulation_utils
  use modulation_chisq
  implicit none

  include 'mpif.h'

  real(dp)           :: evidence, sigma_n, f_rms, prior_f, q_rms, prior_q, n_rms, prior_n
  real(dp)           :: theta_rms, nullval
  real(dp)           :: theta_s, phi_s, A_s, q_s, n_s
  integer(i4b)       :: lmin, lmax, lfix, nummaps, ordering, unit, nmaps
  integer(i4b)       :: myid, ierr, numprocs, root, i, j, k, l, n, iargc, nside_map, seed
  integer(i4b)       :: temp_i, numpoint, numbin, lmax_f, numsamp, chain, l_pivot, l_mod
  logical(lgt)       :: precomp_pls, anynull
  character(len=1)   :: map_text
  character(len=3)   :: proc_id
  character(len=128) :: paramfile, resultfile, operation, covering, filename, prefix
  character(len=128) :: model, cmbfile, rmsfile, beamfile, paramtext, binfile, method
  character(len=128) :: clfile
  type(planck_rng)   :: rng_handle

  integer(i4b),       allocatable, dimension(:,:)   :: bins
  logical(lgt),       allocatable, dimension(:)     :: sample_bin
  real(dp),           allocatable, dimension(:,:)   :: cmbmap, temp_map, beam
  real(dp),           allocatable, dimension(:,:)   :: fid_spectrum
  real(dp),           allocatable, dimension(:,:)   :: prior_cl
  real(dp),           allocatable, dimension(:)     :: def_value
  real(dp),           allocatable, dimension(:)     :: spectrum, propdens
  real(dp),           pointer,     dimension(:,:)   :: pixwin, pls
  integer(i4b),      dimension(MPI_STATUS_SIZE)     :: status
  character(len=80), dimension(1:180)               :: header


  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

  root = 0
  unit = 25+myid
  chain = myid+1

  if (iargc() /= 1) then
     write(*,*) ''
     write(*,*) 'Usage: modulation [parameter file]'
     write(*,*) ''

     call mpi_finalize(ierr)
     stop
  end if 


  call getarg(1,paramfile)
  call get_parameter(unit, paramfile, 'METHOD',          par_string=method)
  call get_parameter(unit, paramfile, 'LMAX',            par_int=lmax)
  call get_parameter(unit, paramfile, 'LMAX_F',          par_int=lmax_f)
  call get_parameter(unit, paramfile, 'NUM_SAMPLES',     par_int=numsamp)
  call get_parameter(unit, paramfile, 'NMAPS',           par_int=nmaps)
  call get_parameter(unit, paramfile, 'SEED',            par_int=seed)
  call get_parameter(unit, paramfile, 'F_RMS',           par_dp=f_rms)
  call get_parameter(unit, paramfile, 'F_PRIOR',         par_dp=prior_f)
  call get_parameter(unit, paramfile, 'Q_RMS',           par_dp=q_rms)
  call get_parameter(unit, paramfile, 'Q_PRIOR',         par_dp=prior_q)
  call get_parameter(unit, paramfile, 'N_RMS',           par_dp=n_rms)
  call get_parameter(unit, paramfile, 'N_PRIOR',         par_dp=prior_n)
  call get_parameter(unit, paramfile, 'THETA_RMS_DEG',   par_dp=theta_rms)
  call get_parameter(unit, paramfile, 'L_PIVOT',         par_int=l_pivot)
  call get_parameter(unit, paramfile, 'L_MOD',           par_int=l_mod)
  call get_parameter(unit, paramfile, 'SIGMA_N',         par_dp=sigma_n)
  call get_parameter(unit, paramfile, 'MAPFILE',         par_string=cmbfile)
  call get_parameter(unit, paramfile, 'BEAMFILE',        par_string=beamfile)
  call get_parameter(unit, paramfile, 'FIDUCIAL_CL_FILE',par_string=clfile)
  call get_parameter(unit, paramfile, 'PRECOMPUTE_PLS',  par_lgt=precomp_pls)
  call get_parameter(unit, paramfile, 'PREFIX',          par_string=prefix)
  call get_parameter(unit, paramfile, 'SINGLE_POINT_LATITUDE', par_dp=theta_s)
  call get_parameter(unit, paramfile, 'SINGLE_POINT_LONGITUDE', par_dp=phi_s)
  call get_parameter(unit, paramfile, 'SINGLE_POINT_AMPLITUDE', par_dp=A_s)
  call get_parameter(unit, paramfile, 'SINGLE_POINT_Q',         par_dp=Q_s)
  call get_parameter(unit, paramfile, 'SINGLE_POINT_N',         par_dp=n_s)
  theta_rms    = theta_rms * pi / real(180.d0)
  theta_s      = (90.d0-theta_s) * pi/180.d0
  phi_s        = phi_s * pi/180.d0

  call initialize_util_mod(paramfile)

  allocate(temp_map(nmaps,0:npix-1))
  allocate(beam(0:lmax, 1))
  allocate(fid_spectrum(0:lmax, 1))
  allocate(cmbmap(0:npix-1, nmaps))
        
  ! Read the data
  temp_i = getsize_fits(cmbfile, nside=nside_map, ordering=ordering)
  call read_bintab(trim(cmbfile), cmbmap, npix, nmaps, nullval, anynull)
  do i = 1, nmaps
     if (ordering == 2) call convert_nest2ring(nside,cmbmap(:,i))
  end do

  call read_beam(beamfile, lmax, beam(:,1:1))
  call read_beam(clfile, lmax, fid_spectrum(:,1:1))

  call read_pixwin(nside, 1, pixwin)
  do l = 0, lmax
     beam(l,:) = beam(l,:) * pixwin(l,1)
  end do
  deallocate(pixwin)
  

  ! Initialize random number generator
  if (myid == root) then
     
     call rand_init(rng_handle, seed)
     
     do i = 1, numprocs-1
        seed = nint(rand_uni(rng_handle)*1000000.)
        call mpi_send(seed, 1, MPI_INTEGER, i, 98, MPI_COMM_WORLD, ierr)
     end do
     
  else 
     
     call mpi_recv(seed, 1, MPI_INTEGER, root, 98, MPI_COMM_WORLD, status, ierr)
     call rand_init(rng_handle, seed)
     
  end if


  ! Precompute Legendre polynomials
  if (precomp_pls) then
     call precompute_pls(lmax, pls)
  else
     allocate(pls(1,1))
  end if


  ! Do the computations
  if (trim(method) == 'mcmc') then

     if (myid == root) write(*,*) '     Doing the MCMC analysis'

     call generate_chain(prefix, chain, unit, rng_handle, cmbmap, real(beam(:,1),dp), &
          & sigma_n, numsamp, real(fid_spectrum(:,1),dp), l_pivot, l_mod, prior_f, f_rms, &
          & prior_q, q_rms, prior_n, n_rms, theta_rms, pls)

  else if (trim(method) == 'evidence') then

     if (myid == root) write(*,*) '     Computing Bayesian evidence'

     call compute_evidence(chain, unit, rng_handle, cmbmap, &
          & real(beam(:,1),dp), sigma_n, lmax_f, numsamp, &
          & real(fid_spectrum(:,1),dp), l_pivot, l_mod, prior_f, prior_q, &
          & prior_n, pls, evidence)

  else if (trim(method) == 'evaluate_single_point') then

     if (myid == root) write(*,*) '     Evaluate single likelihood point'

     if (myid == root) then
        call evaluate_single_point(cmbmap, real(beam(:,1),dp), &
             & sigma_n, real(fid_spectrum(:,1),dp), l_pivot, l_mod, theta_s, phi_s, A_s, q_s, n_s, pls)
     end if


  else if (trim(method) == 'output_chisq_map') then

     if (myid == root) then
        
        write(*,*) '     Computing chi-squared map'

!        call compute_chisq_map(unit, cmbmap, real(beam(:,1),dp), sigma_n, &
!             & real(fid_spectrum(:,1),dp), l_pivot, pls)

     end if

  else

     write(*,*) 'Unknown method = ', trim(method)

  end if


  ! Clean up
  call mpi_finalize(ierr)

  if (allocated(beam))        deallocate(beam)
  if (allocated(cmbmap))      deallocate(cmbmap)
  if (allocated(sample_bin))  deallocate(sample_bin)
  deallocate(pls)



end program modulation




