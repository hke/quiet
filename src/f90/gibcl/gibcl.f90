!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  gibcl, a Gibbs sampler for CMB analysis.             !
!  June 2011 Unni Fuskeland                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program gibcl

  use healpix_types
  use pix_tools
  use fitstools
  use quiet_fileutils
  use math_tools
  use alm_tools
  use quiet_postutils
  use quiet_utils
  use rngmod
  use quiet_healpix_mod
  use gibcl_utils
  use gibcl_noise_mod


  implicit none

  include "mpif.h"

!! declare variables
  integer(i4b)       :: iargc
  integer(i4b)       :: myid, numprocs, ierr, root, seed
  character (len=128) :: seed_chr, paramfile

  !external random
  type(planck_rng)     :: rng_handle

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0

  ! Initalise random generator
  call getarg(2,seed_chr) ; read (unit=seed_chr,fmt=*) seed
  call rand_init(rng_handle,seed+myid)

  ! Initialize modules
  call getarg(1, paramfile)
  call initialize_noise_mod(paramfile)

  if (iargc() /= 2) then
     write(*,*) ''
     write(*,*) 'Usage: mpirun -n [number] gibcl [parameter file] [seed]'
     write(*,*) ''
  else if (myid == root) then
     write(*,*) '-------------------- Cl gibbs sampler starting  --------------------'
     write(*,*)
  end if
  call main

  ! And exit MPI
  call mpi_finalize(ierr)
  if (myid == root) then 
     write(*,*)
     write(*,*) '-------------------- Cl gibbs sampler finished  ----------'
  end if

contains

  subroutine main
    implicit none

!! declare variables
    integer(i4b)                                :: unit,n_gibbs, n, i,l, m, npix, nside, ordering, nmaps, lmax, indexmax
    character(len=256)                          :: paramfile, clfile, beamfile, wmapfile, teoriclfile
    character(len=256)                          :: simmapoutfile, smapoutfile, preClfile, filename,pre
    character(len=5)                            :: filn
    character(len=3)                            :: chr_myid
    complex(dpc), allocatable, dimension(:,:,:) ::s_Tlm,d_Tlm
    real(dp), allocatable, dimension(:)         :: s
    real(dp), allocatable, dimension(:,:)       :: d_map, s_map
    real(dp), pointer    , dimension(:,:)       :: wmap_map
    real(dp), allocatable, dimension(:,:)       :: cls_teori, bp, Cl
    real(dp)                                    :: beam_fwhm, sigma_noise
    real(dp), pointer,     dimension(:,:)       :: pixwin
    logical(lgt)                                :: simulated_data, chatty

    chatty = .false.

!! Read parameters from the parameter file

    unit = getlun()
    call getarg(1, paramfile)
    if (myid == root) write(*,*) "Reading from ", trim(paramfile)
    
    call get_parameter(unit, paramfile,  'SIMULATED_DATA', par_lgt=simulated_data)
    call get_parameter(unit, paramfile,  'CLFILE',      par_string=clfile)
    call get_parameter(unit, paramfile,  'BEAMFILE',    par_string=beamfile)
    call get_parameter(unit, paramfile,  'WMAPFILE',    par_string=wmapfile)
    call get_parameter(unit, paramfile,  'BEAM_FWHM',   par_dp=beam_fwhm)
    call get_parameter(unit, paramfile,  'SIMMAPOUTFILE',   par_string=simmapoutfile)
    call get_parameter(unit, paramfile,  'SMAPOUTFILE',   par_string=smapoutfile)
    call get_parameter(unit, paramfile,  'TEORICLFILE',   par_string=teoriclfile)
    call get_parameter(unit, paramfile,  'PRECLFILE',   par_string=preClfile)
    call get_parameter(unit, paramfile,  'PRE',         par_string=pre)
    call get_parameter(unit, paramfile,  'LMAX',        par_int=lmax)
    call get_parameter(unit, paramfile,  'NSIDE',       par_int=nside)
    call get_parameter(unit, paramfile,  'SIGMA_NOISE', par_dp=sigma_noise)
    call get_parameter(unit, paramfile,  'N_GIBBS',     par_int=n_gibbs)
    indexmax = lmax*(lmax+2)

    !get the file structure (pre+myid)
    call int2string(myid,chr_myid)
    pre=trim(pre)//trim(chr_myid)//'/'
    if (myid == root) write(*,*) "Output directory: ",trim(pre)

!! Read powspec, maps and organize data

    !Read in powerspectrum
    if (myid == root) write(*,*) "Loading in file ", trim(clfile)
    allocate(cls_teori(0:lmax,1))   
    call read_powspec(clfile, cls_teori)
   
    !Get beam
    allocate(bp(0:lmax,1))
    if(simulated_data) then
       !if simulated data, generate beam from beam_fwhm
       if (myid == root) write(*,*) "Generating beam from FWHM= ", beam_fwhm
       call lag_beam(beam_fwhm,bp(:,1))
    else
       !if wmap data, read in beam from file
       if (myid == root)write(*,*) "Loading in beam from file ", trim(beamfile)
       call read_beam(beamfile, bp)
    end if

    !Get pixwin and multiply with beam
    call read_pixwin(nside, 1, pixwin)
    bp = bp * pixwin
    deallocate(pixwin)

    !Either make a simulation or read in data
    npix = 12*nside**2
    if (myid == root) write(*,*) "nside = ",nside," lmax = ",lmax,"npix= ",npix
    allocate(d_map(0:npix-1,1))
    if(simulated_data) then
       if (myid == root) write(*,*) "Simulating data"
       call simulate_data(simulated_data,nside,lmax,cls_teori(0:lmax,1),sigma_noise,d_map(0:npix-1,1), bp(0:lmax,1),pre)
    else if(.not.simulated_data .and. nside.lt.512) then
       if (myid == root) write(*,*) "Reading wmap data map from file ",trim(wmapfile)
       allocate(wmap_map(0:npix-1,1:3))
       call read_map(wmap_map,ordering,wmapfile,nside=nside)
       if (ordering == 2) call convert_nest2ring(nside, wmap_map(:,1))
       d_map(:,1)=wmap_map(:,1)
       deallocate(wmap_map)
       call simulate_data(simulated_data,nside,lmax,cls_teori(0:lmax,1),sigma_noise,d_map(0:npix-1,1), bp(0:lmax,1),pre)
    end if


!!do the gibbs sampler

    if (myid == root) write(*,*) "doing the gibbs sampler with ",n_gibbs," steps"

    !initialise Cl
    allocate(Cl(0:lmax,1))
    Cl = 0.d0
    do l=2,lmax
       Cl(l,1) = 2000.*2*pi/ ( l*(l+1) )
    end do
    filename=trim(preClfile)//'00000.txt'
    if (myid == root) write (*,*) "writing ",n_gibbs," Cl_out  to ",trim(teoriclfile)
    call plot_cls(Cl,trim(pre)//filename)

    allocate(s(0:indexmax))

    do n = 1,n_gibbs
       !sample from the Gaussian distribution
       if (chatty.and.myid == root) write(*,*) "cl2s"
       call cg_cl2s_sample(sqrt(Cl(:,1)),bp(:,1),d_map(:,1),s)

       !sample from the inverse gamma function
       if (chatty.and.myid == root) write(*,*) "s2cl"
       call s2cl_sample(Cl(:,1),s)

       !write out cls for each step
       call int2string(n,filn)
       filename=trim(preClfile)//trim(filn)//'.txt'
       call plot_cls(Cl,trim(pre)//filename)
    end do



!! Write out stuff to file

    !write out cls_teori *l*(l+1)/2pi
    if (myid == root) write (*,*) "writing cls_teori to ",trim(teoriclfile)
    call plot_cls(cls_teori ,trim(pre)//teoriclfile)

    !write out data map (either the simulation or real data)
    if (myid == root) write (*,*) "writing d-map to ",trim(simmapoutfile)
    call write_map(d_map, ordering=ring, trim(pre)//trim(simmapoutfile))

    !write out the last signal cmb map
    if (myid == root) write (*,*) "writing s-map to ",trim(smapoutfile)
    allocate(s_map(0:npix-1,1),s_Tlm(1,0:lmax,0:lmax))
    call real_multiplication(bp(:,1),s)   !beamsmooth map
    call alm_real2complex(s, s_Tlm)
    call alm2map(nside,lmax,lmax,s_Tlm,s_map(:,1))
    call write_map(s_map, ordering=ring, trim(pre)//trim(smapoutfile))


!! Clean up

    deallocate(cls_teori)
    deallocate(bp)
    deallocate(d_map)    
    deallocate(Cl)
    deallocate(s)
    deallocate(s_map)    
    deallocate(s_Tlm)

   end subroutine main



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine cg_cl2s_sample(sqroot_Cl,bp,d_map,s)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! sample s from the distribution P(s | Cl,d)
     implicit none

     real(dp), dimension(0:), intent(in)              :: sqroot_Cl, bp
     real(dp), dimension(0:), intent(out)             :: s
     real(dp), dimension(0:), intent(in)              :: d_map

     integer(i4b)                                   :: nside, npix,i,j,l,m, ii,lmax, indexmax
     real(dp), allocatable, dimension(:)            :: rhs

     lmax  = size(sqroot_Cl)-1
     indexmax=lmax*(lmax+2)
     npix  = size(d_map,1)
     nside = nint(sqrt(real(npix,dp)/12.d0))

     allocate(rhs(0:indexmax))

     ! Set up right-hand side
     call generate_rhs(d_map, rhs, bp, sqroot_Cl)

     !do the CG with preconditioning
     s = 0.d0   !initialise
     call cg_precon(s,rhs,bp,sqroot_Cl)

     !wrap up results

     call real_multiplication(sqroot_Cl,s)

     deallocate(rhs)

   end subroutine cg_cl2s_sample


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine cg_precon(x,rhs,bp,sqroot_Cl)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !do the conjugate gradients search with preconditioner

     implicit none
     real(dp), dimension(0:), intent(in)              :: sqroot_Cl,bp,rhs
     real(dp), dimension(0:), intent(inout)           :: x

     !for cg
     real(dp), allocatable, dimension(:)            :: Avec,Mvec,r,d,s,q
     integer(i4b)                                   :: cg_steps_max,lmax,npix,indexmax,l,i
     real(dp)                                       :: err, alpha, beta, delta_new, delta_old, delta_0
     logical(lgt)  :: chatty

     chatty = .false.

     !do the cg, Mat*Ax=Mat*rhs

     lmax = size(sqroot_Cl)-1
     indexmax=lmax*(lmax+2)

     allocate(r(0:indexmax))
     allocate(d(0:indexmax))
     allocate(q(0:indexmax))
     allocate(s(0:indexmax))
     allocate(Avec(0:indexmax))
     allocate(Mvec(0:indexmax))

     cg_steps_max=4000
     err = 1.e-6
     i = 0

     call compute_Avec(x ,Avec, bp, sqroot_Cl)
     r = rhs - Avec
     call compute_Mvec(r,Mvec, bp, sqroot_Cl)
     d = Mvec
     delta_new = dot_product(r,d)
     delta_0 = delta_new
     do while((i.lt.cg_steps_max) .and. (delta_new.gt.(err**2 *delta_0))  )
        call  compute_Avec(d ,Avec, bp, sqroot_Cl)
        q = Avec
        alpha = delta_new/dot_product(d,q)
        x = x + alpha*d
        if (modulo(i+1,50).eq.0) then
           call compute_Avec(x ,Avec, bp, sqroot_Cl)
           r = rhs - Avec
        else
           r = r - alpha*q
        end if
        call compute_Mvec(r,Mvec, bp, sqroot_Cl)
        s = Mvec
        delta_old = delta_new
        delta_new = dot_product(r,s)
        beta = delta_new/delta_old
        d = s + beta*d
        if(chatty.and.myid==0) write(*,*)i,sqrt(delta_new/delta_0), x(6)
        i=i+1
     end do

     if(myid==0) write(*,*)i ,"gibbs steps. err= ",err 

     deallocate(r)
     deallocate(d)
     deallocate(q)
     deallocate(s)
     deallocate(Avec)
     deallocate(Mvec)

   end subroutine cg_precon


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_Mvec(x, Mx, bp, sqroot_Cl)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !calculates the preconditioner matrix multiplied by a vector

     implicit none
     real(dp), dimension(0:), intent(in)               :: x
     real(dp), dimension(0:), intent(out)              :: Mx
     real(dp), dimension(0:), intent(in)               :: sqroot_Cl,bp

     real(dp)     :: Nl
     integer(i4b) :: i, l, lmax,m
     real(dp), allocatable, dimension(:)               :: lvec

     Nl   = get_N_l()
     lmax = size(sqroot_Cl)-1
     allocate(lvec(0:lmax))

     Mx = x
     lvec = 1.d0/ ( 1.d0 + (bp(l)**2 * sqroot_Cl(l)**2 / Nl) )
     call real_multiplication(lvec,Mx)
     deallocate(lvec)


   end subroutine compute_Mvec

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_Avec(x, Ax, bp, sqroot_Cl)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !the matrix A is too large to keep in memory, so we have to calculate the matrix multiplication A*x 
     !and just store the vector each time we need A.

     implicit none
     real(dp), dimension(0:), intent(in)               :: x
     real(dp), dimension(0:), intent(out)              :: Ax
     real(dp), dimension(0:), intent(in)               :: sqroot_Cl,bp

     integer(i4b)                                      :: l, i, m, lmax
     real(dp), allocatable, dimension(:)               :: lvec

     lmax = size(sqroot_Cl)-1
     allocate(lvec(0:lmax))
     Ax = x

     lvec = bp * sqroot_Cl
     call real_multiplication(lvec , Ax)
     call multiply_by_invN(.false., Ax)
     call real_multiplication(lvec , Ax)

     Ax = Ax + x
     deallocate(lvec)
    

   end subroutine compute_Avec

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine generate_rhs(d_map, rhs, bp, sqroot_Cl)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! set ut the right hand side of the equation

     implicit none
     real(dp), dimension(0:), intent(in)              :: sqroot_Cl, bp
     real(dp), dimension(0:), intent(in)              :: d_map
     real(dp), dimension(0:), intent(out)             :: rhs

     integer(i4b)                                   :: i,j,l,m, ii,lmax
     real(dp), allocatable, dimension(:)            :: alms,bpcl

     lmax  = size(sqroot_Cl)-1
     allocate(alms(0:size(rhs)-1))
     allocate(bpcl(0:lmax))
     bpcl = bp * sqroot_Cl

     ! Initialize rhs
     rhs = 0.d0

     ! Add mean-field term
     call map2alm_real(d_map, alms)
     call multiply_by_invN(.false., alms)
     call real_multiplication(bpcl , alms)
     rhs = rhs + alms

     ! Add first fluctuation term
     do i = 0, size(rhs,1)-1
        rhs(i) = rhs(i) + rand_gauss(rng_handle)
     end do

     ! Add second fluctuation term
     do i = 0, size(alms)-1
        alms(i) = rand_gauss(rng_handle)
     end do
     call multiply_by_invN(.true., alms)
     call real_multiplication(bpcl , alms)
     rhs = rhs + alms

     deallocate(alms)
     deallocate(bpcl)

   end subroutine generate_rhs


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine s2cl_sample(Cl,s)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! sample Cl from the distribution P(Cl | s,d)

     implicit none
     real(dp), dimension(0:), intent(in)         :: s
     real(dp), dimension(0:), intent(out)        :: Cl

     real(dp), allocatable, dimension(:)      :: sigma_l
     real(dp)                                 :: z, sum_alm
     integer(i4b)                             :: l,i,m,lmax

     lmax = size(Cl)-1
     allocate(sigma_l(0:lmax))

     !transform s to sigma_l
     sigma_l = 0.d0
     i = 4
     do l = 2,lmax
        do m = -l, l
           sigma_l(l) = sigma_l(l)+ s(i)**2
           i            = i+1
        end do
        sigma_l(l) = sigma_l(l) / real(2*l+1,dp)
     end do

     ! do the sampling and get Cl
     Cl= 0.d0
     do l = 2,lmax
        z=0
        do i = 0, 2*l-1
           z = z + (rand_gauss(rng_handle)**2)
        end do
        Cl(l) =  real(2*l+1,dp) * sigma_l(l)  / z
     end do

     deallocate(sigma_l)

   end subroutine s2cl_sample



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine simulate_data(simulated_data,nside,lmax,cls,sigma_noise,map,bp,pre)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     implicit none

     logical(lgt),        intent(in)                    :: simulated_data
     integer(i4b),        intent(in)                    :: lmax,nside
     real(dp)   ,         intent(in)                    :: sigma_noise
     real(dp)   ,         intent(in),    dimension(0:)  :: bp, cls
     character(len=256)  ,intent(in)                    :: pre
     real(dp)   ,         intent(inout), dimension(0:)  :: map

     complex(dpc), allocatable, dimension(:,:,:)    :: alm_Tlm
     real(dp),     allocatable, dimension(:,:)      :: sigma_l     
     integer(i4b)                                   :: npix, i, l, m
     character(len=256)                             :: outfile

     npix=size(map)
     allocate(alm_Tlm(1,0:lmax,0:lmax)) 

     !simulate signal
     if(simulated_data) then
        allocate(sigma_l(0:lmax,1))
        alm_Tlm(1,0:1,:) = cmplx(0.d0,0.d0)
        do l = 2,lmax
           alm_Tlm(1,l,0)= sqrt(cls(l))*cmplx(rand_gauss(rng_handle),0.d0)
           do m= 1,l
              alm_Tlm(1,l,m)= (sqrt(cls(l)/2.d0))*cmplx(rand_gauss(rng_handle),rand_gauss(rng_handle))
           end do
           alm_Tlm(1,l,:) = bp(l) * alm_Tlm(1,l,:)
        end do
        
        sigma_l = 0.d0
        do l = 2,lmax
           sigma_l(l,1) = abs(alm_Tlm(1,l,0))**2
           do m=1,l
              sigma_l(l,1) = sigma_l(l,1)+ 2*abs(alm_Tlm(1,l,m))**2
           end do
           sigma_l(l,1) = sigma_l(l,1) / (2*l+1)
        end do
        outfile=trim(pre)//'sigma_l_outfile.txt'
        call plot_cls(sigma_l,outfile)
        deallocate(sigma_l)

     !if wmap data
     else if(.not. simulated_data) then
        !remove monopol and dipole
        call remove_monopol_dipole(map)
        !beamsmooth
        call map2alm(nside,lmax,lmax,map,alm_Tlm,[0.d0,0.d0],get_hpix_ringweights(nside))
        do l = 0,lmax
           alm_Tlm(1,l,:) = bp(l) * alm_Tlm(1,l,:)
        end do
     end if

     call alm2map(nside, lmax, lmax, alm_Tlm, map)
     deallocate(alm_Tlm)

     !simulate noise, add white noise to map to create d = As + n
     do i = 0, npix-1
        map(i) = map(i) + (sigma_noise *rand_gauss(rng_handle))
     end do

   end subroutine simulate_data




end program gibcl

