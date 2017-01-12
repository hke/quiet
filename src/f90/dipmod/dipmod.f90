 program dipmod
  use dipmod_evidence
  use dipmod_mcmc
  use dipmod_utils
  use dipmod_chisq
  implicit none

  include 'mpif.h'

  real(dp)           :: evidence, sigma_n, f_rms, prior_f, q_rms, prior_q, n_rms, prior_n
  real(dp)           :: theta_rms
  integer(i4b)       :: nside, lmin, lmax, lfix, nummaps, npix, ordering, unit, nmaps
  integer(i4b)       :: myid, ierr, numprocs, root, i, j, k, l, n, iargc, nside_map, seed
  integer(i4b)       :: temp_i, numpoint, numbin, numsamp, chain, numtemp, l_pivot
  integer(i4b)       :: numsim
  logical(lgt)       :: precomp_pls
  character(len=1)   :: map_text
  character(len=3)   :: proc_id
  character(len=128) :: paramfile, maskfile, resultfile, operation, covering, filename
  character(len=128) :: model, cmbfile, rmsfile, beamfile, paramtext, binfile, method
  character(len=128) :: clfile
  type(planck_rng)   :: rng_handle

  integer(i4b),       allocatable, dimension(:,:)   :: bins
  logical(lgt),       allocatable, dimension(:)     :: mask, sample_bin
  real(sp),           allocatable, dimension(:,:)   :: cmbmap, temp_map, beam, templates
  real(sp),           allocatable, dimension(:,:)   :: fid_spectrum
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
  unit = getlun()
  chain = myid+1

  if (iargc() /= 1) then
     write(*,*) ''
     write(*,*) 'Usage: dipmod [parameter file]'
     write(*,*) ''

     call mpi_finalize(ierr)
     stop
  end if 


  call getarg(1,paramfile)
  call get_parameter(unit, paramfile, 'METHOD',          par_string=method)
  call get_parameter(unit, paramfile, 'LMAX',            par_int=lmax)
  call get_parameter(unit, paramfile, 'NUM_SAMPLES',     par_int=numsamp)
  call get_parameter(unit, paramfile, 'NSIDE',           par_int=nside)
  call get_parameter(unit, paramfile, 'NMAPS',           par_int=nmaps)
  call get_parameter(unit, paramfile, 'SEED',            par_int=seed)
  call get_parameter(unit, paramfile, 'NUMSIM',          par_int=numsim)
  call get_parameter(unit, paramfile, 'F_RMS',           par_dp=f_rms)
  call get_parameter(unit, paramfile, 'F_PRIOR',         par_dp=prior_f)
  call get_parameter(unit, paramfile, 'Q_RMS',           par_dp=q_rms)
  call get_parameter(unit, paramfile, 'Q_PRIOR',         par_dp=prior_q)
  call get_parameter(unit, paramfile, 'N_RMS',           par_dp=n_rms)
  call get_parameter(unit, paramfile, 'N_PRIOR',         par_dp=prior_n)
  call get_parameter(unit, paramfile, 'THETA_RMS_DEG',   par_dp=theta_rms)
  call get_parameter(unit, paramfile, 'L_PIVOT',         par_int=l_pivot)
  call get_parameter(unit, paramfile, 'SIGMA_N',         par_dp=sigma_n)
  call get_parameter(unit, paramfile, 'MASKFILE',        par_string=maskfile)
  call get_parameter(unit, paramfile, 'MAPFILE',         par_string=cmbfile)
  call get_parameter(unit, paramfile, 'BEAMFILE',        par_string=beamfile)
  call get_parameter(unit, paramfile, 'FIDUCIAL_CL_FILE',par_string=clfile)
  call get_parameter(unit, paramfile, 'NUMTEMPLATES',    par_int=numtemp)
  call get_parameter(unit, paramfile, 'PRECOMPUTE_PLS',  par_lgt=precomp_pls)               
  npix         = 12*nside**2
  theta_rms    = theta_rms * pi / real(180.d0)
  

  ! Read mask file and convert to ring if necessary
  allocate(mask(0:npix-1))
!  call get_maskparam(unit, maskfile, nside, ordering)
!  call read_mask(unit, maskfile, mask)
  
!  if (ordering == 2) then
!     call convert_nest2ring_mask(mask)
!  end if
        
  allocate(temp_map(nmaps,0:npix-1))
  allocate(beam(0:lmax, 1))
  allocate(fid_spectrum(0:lmax, 1))
  allocate(cmbmap(0:npix-1, nmaps))
  allocate(templates(0:npix-1, numtemp))
        
  ! Read the data
  call read_map(cmbmap, ordering, cmbfile)
  do i = 1, nmaps
     if (ordering == 2) call convert_nest2ring(nside,cmbmap(:,i))
  end do
  
  do i = 1, numtemp
     call int2string(i, map_text)
     paramtext = 'TEMPLATE' // map_text
     call get_parameter(unit, paramfile, trim(paramtext), par_string=filename)
     ! HKE: COMMENTED OUT TO MAKE THINGS COMPILE RIGHT NOW. MUST BE FIXED!!!
     !call read_map(templates(:,i:i), ordering, filename)
     if (ordering == 2) call convert_nest2ring(nside, templates(:,i))
  end do
           
  ! HKE: COMMENTED OUT TO MAKE THINGS COMPILE RIGHT NOW. MUST BE FIXED!!!
  !call read_beam(beamfile, lmax, beam(:,1:1))
  !call read_beam(clfile, lmax, fid_spectrum(:,1:1))

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
     call precompute_pls(mask, lmax, pls)
  else
     allocate(pls(1,1))
  end if


  if (myid == root) then

     ! Do the computations
     if (trim(method) == 'mcmc') then

        if (myid == root) write(*,*) '     Doing the MCMC analysis'
        
        call generate_chain(0, unit, rng_handle, cmbmap, templates, real(beam(:,1),dp), &
             & mask, sigma_n, numsamp, real(fid_spectrum(:,1),dp), l_pivot, prior_f, f_rms, &
             & prior_q, q_rms, prior_n, n_rms, theta_rms, pls)
        
     else if (trim(method) == 'evidence') then
        
        if (myid == root) write(*,*) '     Computing Bayesian evidence'
        
        call compute_evidence(0, unit, rng_handle, cmbmap, templates, &
             & real(beam(:,1),dp), mask, sigma_n, numsamp, &
             & real(fid_spectrum(:,1),dp), l_pivot, prior_f, prior_q, &
             & prior_n, pls, evidence)
        
     else if (trim(method) == 'output_chisq_map') then
        
        if (myid == root) then
           
           write(*,*) '     Computing chi-squared map'
           
           call compute_chisq_map(unit, cmbmap, templates, real(beam(:,1),dp), mask, sigma_n, &
                & real(fid_spectrum(:,1),dp), l_pivot, pls)
           
        end if
        
     else
        
        write(*,*) 'Unknown method = ', trim(method)
        
     end if
     
  end if

  do i = myid+1, numsim, numprocs

     call generate_simulation(rng_handle, beam(:,1), sigma_n, fid_spectrum(:,1), cmbmap(:,1))

     ! Do the computations
     if (trim(method) == 'mcmc') then

        if (myid == root) write(*,*) '     Doing the MCMC analysis'
        
        call generate_chain(i, unit, rng_handle, cmbmap, templates, real(beam(:,1),dp), &
             & mask, sigma_n, numsamp, real(fid_spectrum(:,1),dp), l_pivot, prior_f, f_rms, &
             & prior_q, q_rms, prior_n, n_rms, theta_rms, pls)
        
     else if (trim(method) == 'evidence') then
        
        if (myid == root) write(*,*) '     Computing Bayesian evidence'
        
        call compute_evidence(i, unit, rng_handle, cmbmap, templates, &
             & real(beam(:,1),dp), mask, sigma_n, numsamp, &
             & real(fid_spectrum(:,1),dp), l_pivot, prior_f, prior_q, &
             & prior_n, pls, evidence)
        
     else if (trim(method) == 'output_chisq_map') then
        
        if (myid == root) then
           
           write(*,*) '     Computing chi-squared map'
           
           call compute_chisq_map(unit, cmbmap, templates, real(beam(:,1),dp), mask, sigma_n, &
                & real(fid_spectrum(:,1),dp), l_pivot, pls)
           
        end if
        
     else
        
        write(*,*) 'Unknown method = ', trim(method)
        
     end if


  end do




  ! Clean up
  call mpi_finalize(ierr)

  if (allocated(beam))        deallocate(beam)
  if (allocated(cmbmap))      deallocate(cmbmap)
  if (allocated(sample_bin))  deallocate(sample_bin)
  deallocate(pls)



end program dipmod




