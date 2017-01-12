program sindex
  use healpix_types
  use quiet_lx_mod
  use pix_tools
  use fitstools
  use quiet_fileutils
  use math_tools
  use alm_tools
  use quiet_postutils
  use quiet_utils
!  use quiet_mapfile_mod
  use quasi_newton_mod
  use sindex_maxlike_mod
  use newton_minimizing_mod
  use powell_mod
  use rngmod


  implicit none

  include "mpif.h"

!  type banddata
!     real(dp),     allocatable, dimension(:,:,:)     :: map
!     real(dp),     allocatable, dimension(:,:,:,:,:) :: cov
!     real(dp),     allocatable, dimension(:)         :: freq
!     integer(i4b), allocatable, dimension(:)         :: pixels, map2mask
!     integer(i4b)    :: n, npix, nside, nmaps, ordering, pol, num_bands
!  end type banddata

  integer(i4b)       :: iargc
  integer(i4b)       :: unit, myid, numprocs, ierr, root, nregions, startreg, ii
  character(len=30)  :: kommando, nregions_char
  type(planck_rng)   :: rng

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0
!  unit        = 42+myid

  call rand_init(rng,1347402)

  if (iargc() == 0) then
     call give_user_info
  else if (myid == root) then
     write(*,*) '-------------------- QUIET spectral index calculator  --------------------'
     write(*,*)
  end if

  ! Get name of main kommando
  call getarg(1,kommando)

  if (kommando == 'index') then
     call getarg(3,nregions_char) ; read (unit=nregions_char,fmt=*) nregions
     if(iargc() /= 4) then 
        startreg=1
     else
        call getarg(4,nregions_char) ; read (unit=nregions_char,fmt=*) startreg
     end if
     if (myid==0) write(*,*) 'Find spectral index from two maps. ',nregions,' regions'

     do ii = startreg+myid, nregions, numprocs
        call spectral_index(ii)
     end do
  else
     call give_user_info
  end if

  ! And exit
  call mpi_finalize(ierr)
  if (myid == root) then 
     write(*,*)
     write(*,*) '-------------------- QUIET spectral index calculator finished  ----------'
  end if
  
contains



  !-----------------------------------------------------------------------------------------------
  ! subroutine spectral_index
  !-----------------------------------------------------------------------------------------------

  subroutine spectral_index(region)
    implicit none
    
    integer(i4b),       intent(in) :: region 
!    type(banddata)                 :: data
    type(banddata), allocatable, dimension(:) :: data
    character(len=256)             :: parfile, directory, subdir, outprefix, outdirprefix, outfile, outpostfilename
    character(len=5)               :: region_char
    integer(i4b)                   :: i, j, jj, k, n, p, comp, maskcomp, lstop, nband=2, subreg, regpix,v
    integer(i4b)                   :: nmaps, nside, npix, ordering, ierr, n_bootstrap
    integer(i4b)                   :: nmaps_mask, nmaps_rms
    integer(i4b)                   :: areanside, areanpix, num, max_nside, ns, ns_max, max_npix, teller
    real(dp)                       :: nullval, fwhm, uncertlimit
    logical(lgt)                   :: dogrid, dogridbm, domaxlike, doscatter, doratio, doalpha
    real(dp)                       :: healnan=-1.6375d30, root=0, m(2), factor(2), invN(2,2)

    real(dp),     allocatable, dimension(:,:)   :: outmap
    real(dp),     allocatable, dimension(:,:)   :: spectral, offset, uncert, fg, chi2
    integer(i4b), allocatable, dimension(:)     :: reg_pixels

    unit=42
    ! Get name of parameter file
    if (iargc() < 3) then
       write(*,*) 'Usage: mpirun -n N sindex index parfile nregions [startreg]'
       call give_user_info
    else 
       call getarg(2, parfile)
    end if

    ! Get parameters
    call get_parameter(unit, parfile,  'DIRECTORY',           par_string=directory)
    call get_parameter(unit, parfile,  'OUTPREFIX',           par_string=outprefix)
    call get_parameter(unit, parfile,  'INCLUDE_MAXLIKE',     par_lgt=domaxlike)
    call get_parameter(unit, parfile,  'INCLUDE_SCATTERPLOT', par_lgt=doscatter)
    call get_parameter(unit, parfile,  'INCLUDE_NAIVE_RATIO', par_lgt=doratio)
    call get_parameter(unit, parfile,  'GRID',                par_lgt=dogrid)
    call get_parameter(unit, parfile,  'GRID_BM',             par_lgt=dogridbm)
    call get_parameter(unit, parfile,  'ALPHA',               par_lgt=doalpha)
    call get_parameter(unit, parfile,  'N_BOOTSTRAP',         par_int=n_bootstrap)

    if (dogrid) write(*,*) "Running with grid method"

    call sindex_int2string(region, region_char)
    directory=trim(directory) //trim(region_char)//'/'
    outdirprefix=trim(directory) // trim(outprefix)
    write(*,*) trim(outdirprefix)

    if (file_exists(trim(directory) //"regions.txt")) then
       ! Read how many dirs in each strip
       unit = getlun()
       open(unit,file=trim(directory)//"regions.txt",action="read",status="old")
       read(unit,*) subreg
       close(unit)

       ! Read input maps, covs and freqs (and convert to nest if necessary)
       allocate(data(1:subreg))    
       do i=1,subreg
          call sindex_int2string(i,region_char)
          subdir = trim(directory)//"dir"//trim(region_char)          
          call read_banddata(subdir, data(i), myid, domaxlike)
       end do
    else
       subreg = 1
       allocate(data(1:subreg))    
       call read_banddata(directory, data(1), myid, domaxlike)
    end if


  
!!$    ! Write input maps
!!$    outfile = trim(outprefix) // '_input_map1.fits'
!!$    call write_map(data%map(:,:,1), data%pixels, data%nside, data%ordering, outfile)
!!$    outfile = trim(outprefix) // '_input_map2.fits'
!!$    call write_map(data%map(:,:,2), data%pixels, data%nside, data%ordering, outfile)


    !find number of pixels in (the largest) region (the strip):
    regpix=0
    do i=1,subreg
       regpix=regpix + data(i)%n
    end do

    !find the pixels-array for (the largest) region (the strip)
    allocate(reg_pixels(regpix))
    jj = 1
    do i = 1, subreg
       do j = 1, data(i)%n
          reg_pixels(jj) = data(i)%pixels(j)
          jj=jj+1
       end do
    end do


    !!!!!!!! MAXLIKE !!!!!!!!!!
    if (domaxlike) then

       nmaps=data(1)%nmaps
       allocate(spectral(regpix,nmaps))
       allocate(uncert(regpix,nmaps))

       comp=nmaps
       if (data(1)%pol==2) comp=nmaps+1

       allocate(offset(regpix,comp))
       allocate(fg(regpix,comp))
       allocate(chi2(regpix,comp))
       chi2=healnan
       call maxlikeplot(data, spectral(1,:), uncert(1,:), offset(:,:), fg(:,:), dogrid, dogridbm, directory, outdirprefix, doalpha)

       do i = 1, regpix
          spectral(i,:) = spectral(1,:)
          uncert(i,:) = uncert(1,:)
       end do

       ! Write maxlike spectral index map and uncertainty and offset to file
       outfile = trim(outdirprefix) // '_spectral_index_maxlike.fits'
       call write_map(spectral, reg_pixels, data(1)%nside, data(1)%ordering, outfile)
       outfile = trim(outdirprefix) // '_uncert_maxlike.fits'
       call write_map(uncert, reg_pixels, data(1)%nside, data(1)%ordering, outfile)

       if (doalpha) then
          outfile = "_maxlike_alpha_spectral_index.txt"
          outpostfilename = "maxlike"
          call calculate_spectral_alpha(outdirprefix,  outfile, outpostfilename, reg_pixels, data(1)%nside, data(1)%ordering)
       end if

       if (data(1)%pol==1) then
          outfile = trim(outdirprefix) // '_offset_maxlike.fits'
          call write_map(offset, reg_pixels, data(1)%nside, data(1)%ordering, outfile)
          outfile = trim(outdirprefix) // '_fg.fits'
          call write_map(fg, reg_pixels, data(1)%nside, data(1)%ordering, outfile)

       else if (data(1)%pol==2) then
          allocate(outmap(regpix,nmaps))
          outmap=healnan

          outfile = trim(outdirprefix) // '_offset_maxlike.fits'
          outmap(:,2:3) =offset(:,1:2)
          call write_map(outmap, reg_pixels, data(1)%nside, data(1)%ordering, outfile)
          outfile = trim(outdirprefix) // '_QUoffset_maxlike.fits'
          outmap(:,2:3) =offset(:,3:4)
          call write_map(outmap, reg_pixels, data(1)%nside, data(1)%ordering, outfile)

          outfile = trim(outdirprefix) // '_fg.fits'
          outmap(:,2:3) =fg(:,1:2)
          call write_map(outmap, reg_pixels, data(1)%nside, data(1)%ordering, outfile)
          outfile = trim(outdirprefix) // '_QUfg.fits'
          outmap(:,2:3) =fg(:,3:4)
          call write_map(outmap, reg_pixels, data(1)%nside, data(1)%ordering, outfile)

          !make a chisquare map
          m(1) =0.d0
          chi2=0.d0
          teller=1
          do i =1,subreg
             n=data(i)%n
             do j=1,n
                do v=1,nfreq
                   m(2) =offset(teller,3)
                   factor(1) = data(i)%map(j,2,v) - fg(teller,3)*ant2thermo(data(i)%freq(v))*(data(i)%freq(v)/data(i)%freq(1))**spectral(teller,1) - m(v)
                   chi2(teller,2) =chi2(teller,2)+ factor(1)**2 /data(i)%cov(j,j,v)

                   m(2) =offset(teller,4)
                   factor(2) = data(i)%map(j,3,v) - fg(teller,4)*ant2thermo(data(i)%freq(v))*(data(i)%freq(v)/data(i)%freq(1))**spectral(teller,1) - m(v)
                   chi2(teller,3) =chi2(teller,3)+ factor(2)**2 / data(i)%cov(j+n,j+n,v)

                   call get_invNelements(i,j,v,invN)
                   chi2(teller,1) =chi2(teller,1) + dot_product(factor, matmul(invN,factor))
!                   call dump_matrix(invN)
!                   write(*,*)data(i)%cov(j,j,v),"=",invN(1,1), data(i)%cov(j+n,j+n,v),"=",invN(2,2)
                end do
                teller=teller+1
             end do
          end do

          chi2(:,1)=chi2(:,2)+chi2(:,3)

          outfile = trim(outdirprefix) // '_QUchisquare.fits'
          outmap(:,1:3) =chi2(:,1:3)
          call write_map(outmap, reg_pixels, data(1)%nside, data(1)%ordering, outfile)

          call terminate_maxlike_mod

       end if
       
       ! Clean up
       if (allocated(spectral)) deallocate(spectral)
       if (allocated(uncert))   deallocate(uncert)
       if (allocated(offset))   deallocate(offset)
       if (allocated(fg))       deallocate(fg)
       if (allocated(chi2))       deallocate(chi2)
       if (allocated(outmap))   deallocate(outmap)

    end if  !end domaxlike


    !!!!!!!!!! SCATTER !!!!!!!!!!!
    if (doscatter) then

       allocate(spectral(regpix,data(1)%nmaps))
       allocate(uncert(regpix,data(1)%nmaps))
       call scatterplot(data, regpix, spectral(1,:), uncert(1,:), n_bootstrap, trim(outdirprefix), doalpha)
       do j = 1, regpix
          spectral(j,:) = spectral(1,:)
          uncert(j,:) = uncert(1,:)
       end do

       ! Write scatter spectral index map and uncertainty and offset to file
       outfile = trim(outdirprefix) // '_spectral_index_scatter.fits'
       call write_map(spectral, reg_pixels, data(1)%nside, data(1)%ordering, outfile)
       outfile = trim(outdirprefix) // '_uncert_scatter.fits'
       call write_map(uncert, reg_pixels, data(1)%nside, data(1)%ordering, outfile)

       if (doalpha) then
          outfile = "_scatter_alpha_spectral_index.txt"
          outpostfilename = "scatter"
          call calculate_spectral_alpha(outdirprefix,  outfile, outpostfilename, reg_pixels, data(1)%nside, data(1)%ordering)
       end if

       ! Clean up
       if (allocated(spectral)) deallocate(spectral)
       if (allocated(uncert)) deallocate(uncert)

    end if  !end scatter


    !!!!!!!!!! NAIVE RATIO !!!!!!!!!!
    if (doratio) then
       do i = 1, subreg
          call sindex_int2string(i,region_char)
          subdir = trim(directory)//"dir"//trim(region_char) //'/'// trim(outprefix)
          call ratioplot(myid, data(i), trim(subdir))
       end do
    end if

    ! Clean up
    do i = 1, subreg
       call free_banddata(data(i))
    end do
    deallocate(data)
    deallocate(reg_pixels)

  end subroutine spectral_index


  subroutine calculate_spectral_alpha(outdirprefix, filename, outpostfilename, pixels, nside, ordering)
    implicit none
    character(len=256) , intent(in)              :: outdirprefix, filename, outpostfilename
    integer(i4b), intent(in)                     :: nside, ordering
    integer(i4b), dimension(:), intent(in)       :: pixels

    character(len=256)                           :: outfile
    integer(i4b)                                 :: unit, i, ai(1:18)
    real(dp)                                     :: betas(1:18), uncerts(1:18)
    real(dp)                                     :: spectral_index, u_syst, u_stat
    real(dp), allocatable, dimension(:,:)        :: map

    unit = getlun()
    open(unit,file=trim(outdirprefix)//trim(filename),action="read",status="old")
    do i=1,18 
       read(unit, *) ai(i), betas(i), uncerts(i)
    end do
    close(unit)

    spectral_index = sum(betas/(uncerts**2))/ sum(1./(uncerts**2)) 
    u_stat = MINVAL(uncerts)

    !write map with spectral indices
    allocate(map(1:size(pixels),1 ))
    map(:,1) = spectral_index 
    outfile = trim(outdirprefix) //'_spectral_index_' // trim(outpostfilename) //'_alphamean.fits'
    call write_map(map, pixels, nside, ordering, outfile)
    deallocate(map)

    !write map with uncertainties
    allocate(map(1:size(pixels),1 ))
    map(:,1) = u_stat
    outfile = trim(outdirprefix) //'_uncert_' // trim(outpostfilename) //'_alphamean.fits'
    call write_map(map, pixels, nside, ordering, outfile)

    !write file with value of spectral index in region
    unit = getlun()
    outfile = trim(outdirprefix) //'_value_beta_' // trim(outpostfilename) //'.txt'
    open(unit,file=trim(outfile),action="write")
    write(unit,*) spectral_index
    close(unit)

    deallocate(map)

  end subroutine calculate_spectral_alpha


  !---------------------------------------------------------------------
  ! Find spectral index from real plot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine maxlikeplot(data, spectral, uncert, offset, fg, dogrid, dogridbm, directory, outdirprefix, doalpha)
    implicit none

    type(banddata), dimension(:),  intent(in)      :: data 
    real(dp),        dimension(:),     intent(out) :: spectral, uncert
    real(dp),        dimension(:,:),   intent(out) :: fg, offset
    logical(lgt)    ,intent(in)                    :: dogrid, dogridbm
    character(len=256) , intent(in)              :: directory, outdirprefix
    logical(lgt)    ,intent(in), optional        :: doalpha

    character(len=256)                           :: outfile
    integer(i4b)                                 :: pol, ai
    integer(i4b)                                 :: i, k, unit
    real(dp)                                     :: healnan=-1.6375d30, a_spectral, a_uncert
    real(dp),        allocatable, dimension(:,:) :: fgmat, offmat

    ! Initialize output
    spectral   = healnan
    uncert     = healnan
    fg         = healnan
    offset     = healnan

    pol = data(1)%pol

    if (pol==1) then
       k=1
       write(*,*) k
       call maxlike(k,data,spectral(1),offset,uncert(1),fg,dogrid, dogridbm, directory)
    else if (pol==2) then

       ! do the pol spectral index for each alpha
       if (doalpha==.true.) then
          outfile = trim(outdirprefix) // '_maxlike_alpha_spectral_index.txt'
          open(43, file=trim(outfile))
          do ai=0,85,5
             write(*,*)"Calling maxlike with a alpha=",ai
             call maxlike(ai,data,a_spectral,offset(:,3:4),a_uncert,fg(:,3:4),dogrid, dogridbm, directory)
             ! Write spectral index to one file
             write(43,*) ai, a_spectral, a_uncert
          end do
          close(43)
          fg = healnan
          offset = healnan
       end if

       do k = 2,3
          write(*,*) k, "Q or U"
          allocate(offmat(size(offset,1),1), fgmat(size(fg,1),1))
          call maxlike(k,data,spectral(k),offmat,uncert(k),fgmat,dogrid, dogridbm, directory)
          offset(:,k-1)=offmat(:,1)
          fg(:,k-1)=fgmat(:,1)
          deallocate(offmat, fgmat)
       end do
       k=4
       write(*,*) k, "Q+U"
       call maxlike(k,data,spectral(1),offset(:,3:4),uncert(1),fg(:,3:4),dogrid, dogridbm, directory)
    end if


  end subroutine maxlikeplot


  !---------------------------------------------------------------------
  ! Find spectral index from maximum likelihood
  !----------------------------------------------------------------------

  subroutine maxlike(k, data, spectral, offset, uncert, fg, dogrid, dogridbm, directory)
    implicit none

    type(banddata), dimension(:),     intent(in)      :: data 
    real(dp),                             intent(out) :: spectral, uncert
    integer(i4b),                         intent(in)  :: k
    real(dp),             dimension(:,:), intent(out) :: offset,fg
    logical(lgt)    ,intent(in)                       :: dogrid, dogridbm
    character(len=256) , intent(in)                   :: directory

    integer(i4b)                                 :: i, j, v, unit, num, int, l, i1,i2, nlik
    real(dp)                                     :: healnan=-1.6375d30
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:)   :: A
    real(dp),        allocatable, dimension(:,:) :: y, cov_params
    logical(lgt)  :: chatty
real(dp) :: foo(100)

    !variables for dfpmin
    INTEGER(I4B)              :: iter
    REAL(dp)                  :: gtol
    REAL(dp)                  :: fret
    REAL(dp), ALLOCATABLE, DIMENSION(:)    :: p, m
    integer(i4b)              :: ierr

    chatty = .false.

    !initialise p and m (spectral index and offset)
    allocate(p(1))
    p(1) = -1.d0
    allocate(m(1:2))       
    gtol = 1d-3
    if (chatty) write(*,*) "Start values ",p

    !important: initialise the maxlike module
    call setup_maxlike_mod(data,k)


    !do the search for spectral index.
    if (dogrid) then
       if (dogridbm) then
          call gridmin_bm(p, uncert,directory)
       else
          call gridmin(p, uncert, directory)
       end if
    else

       !       call dfpmin(p,gtol,iter,fret,func,dfunc,ierr)
       !an optional search method
       write(*,*) "Searching for chi square minimum with powell"
       call powell(p,func,ierr,iter,1.d-6)
       !       call newton_min(p,gtol,1.d0,iter,dfunc,ddfunc,ierr, 18.d0)
       foo(1)=func(p)
       if (k==4 ) then
          open(43, file=trim(directory)//"/chisq_4.txt")
          write(43,*)foo(1)/(2.*sum(data(:)%n) )
          close(43)
       else if (k==2) then
          open(43, file=trim(directory)//"/chisq_2.txt")
          write(43,*)foo(1)/(sum(data(:)%n) )
          close(43)
       else if (k==3) then
          open(43, file=trim(directory)//"/chisq_3.txt")
          write(43,*)foo(1)/(sum(data(:)%n) )
          close(43)
       end if
       if (chatty )write(*,*) "newton min needed ",iter," iterations"

       ! Calculate uncertainties
       allocate(cov_params(1,1))
       call calculate_cov_params2(p,cov_params)    !program 2 is only uncertainty in beta, not beta and m
       write(*,*)"uncert cov2",sqrt(cov_params(1,1))
       uncert = sqrt(cov_params(1,1))

       if (allocated(cov_params)) deallocate(cov_params)
       call calculate_cov_params_full(p,uncert)
       write(*,*)"uncert cov_full",uncert
    end if  !end else(dogrid)


    ! Present values
    spectral = p(1)
    call get_region_A(p,offset,fg)

    write(*,*) 'Spectral index', spectral, "+-",uncert
    if (k>1) write(*,*) 'Region 1 offset value: ', offset(1,1)

    ! Clean up
    if (k /=4) call terminate_maxlike_mod

    deallocate(p,m)


  end subroutine maxlike



  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scatterplot(data, regpix, spectral, uncert, n_bootstrap, outprefix, doalpha)
    implicit none

    type(banddata), dimension(:),  intent(in)    :: data 
    real(dp),          dimension(:), intent(out) :: spectral, uncert
    character(len=*), intent(in)                 :: outprefix
    integer(i4b), intent(in)                     :: regpix, n_bootstrap
    logical(lgt)    ,intent(in), optional        :: doalpha

    integer(i4b)                                 :: pol
    integer(i4b)                                 :: n, nb, i, j, k, jj, ai, v
    real(dp)                                     :: healnan=-1.6375d30
    character(len=5)                             :: filnum, ai_char
    character(len=256)                           :: outfile 
    real(dp),      allocatable, dimension(:,:)   :: polmap, map2
    real(dp),      allocatable, dimension(:,:,:) :: map
    real(dp)                                     :: a_spectral, a_uncert, mean(2), tmpspectral,tmpuncert, notused

!   if (index/=4677 .and. index/=4701 .and. index/=4697 .and. index/=4344) then
!   if (index/=4865) then
!       spectral = healnan
!       return
!    end if

    pol = data(1)%pol

    !make a big map for each region (strip)
    allocate(map(regpix, data(1)%nmaps, data(1)%num_bands))
    jj = 1
    do i = 1, size(data)
       do j = 1, data(i)%n
          map(jj,:,:) = data(i)%map(j,:,:)
          jj=jj+1
       end do
    end do

    ! Find num of existing pixels for given area
    n = 0
    do i = 1, regpix
       if ( healok(map(i,pol,1)) .and. healok(map(i,pol,2)) ) then
          n = n+1
       end if
    end do

    ! Return if too few pixels
    if (n<1) then
       spectral = healnan
       return
    end if

    !use all pixels in bootstrap if n_bootstrap=-1
    if (n_bootstrap .eq. -1) then
       nb=n
    else
       nb=n_bootstrap
    endif

    ! Calculate scale factor for Q and U separately
    if (pol==2) then

       ! do the pol spectral index for each alpha
       if (doalpha==.true.) then
          allocate(map2(size(map,1),data(1)%num_bands))
          outfile = trim(outprefix) // '_scatter_alpha_spectral_index.txt'
          open(43, file=trim(outfile))
          do ai=0,85,5
             do v=1,2
                map2(:,v)=map(:,2,v)*cos(2.d0*ai*pi/180.) + map(:,3,v)*sin(2.d0*ai*pi/180.)
             end do
             call scattercalc(map2(:,1), map2(:,2), n, a_spectral, data(1)%freq(1:2), notused, a_uncert)
             call bootstrap_scatter(map2(:,1), map2(:,2), n, nb, tmpspectral, data(1)%freq(1:2), tmpuncert)
             a_uncert=sqrt(a_uncert**2+tmpuncert**2)

             ! Write spectral index to one file
             write(43,*) ai, a_spectral, a_uncert

             ! Write scatter plot to file
             call int2string(ai,ai_char)
             outfile = trim(outprefix) // '_scatter_a' //trim(ai_char) // '.txt'
             open(42, file=trim(outfile))
             do j=1,2
                mean(j)= sum(map2(:,j))/size(map2(:,j))
             end do
             do i = 1, regpix
                write(42,*) map2(i,1)-mean(1), map2(i,2)-mean(2)
             end do
             close(42)
          end do
          close(43)
       end if


       do k=2,3
          call scattercalc(map(:,k,1), map(:,k,2), n, spectral(k), data(1)%freq(1:2), notused, uncert(k))
          call bootstrap_scatter(map(:,k,1), map(:,k,2), n, nb, tmpspectral, data(1)%freq(1:2), tmpuncert)
          uncert(k)=sqrt(uncert(k)**2+tmpuncert**2)
       end do

       allocate(polmap(pol*regpix, 2))
       polmap = reshape(map(:,2:3,:), [pol*regpix,2])
       ! Calculate scale factor for T or P separately
       n=n*pol
       call scattercalc(polmap(:,1), polmap(:,2), n, spectral(1), data(1)%freq(1:2), notused, uncert(1))
       call bootstrap_scatter(polmap(:,1), polmap(:,2), n, nb, tmpspectral, data(1)%freq(1:2), tmpuncert)
       uncert(1)=sqrt(uncert(1)**2+tmpuncert**2)
       deallocate(polmap)  


       ! Write scatter plot to file
       outfile = trim(outprefix) // '_scatter_2.txt'
       open(42, file=trim(outfile))
       outfile = trim(outprefix) // '_scatter_3.txt'
       open(43, file=trim(outfile))
       outfile = trim(outprefix) // '_scatter_4.txt'
       open(44, file=trim(outfile))
       do i = 1, regpix
          do k = 2,3
             write(44,*) map(i,k,1), map(i,k,2)
             if(k==2) then;     write(42,*) map(i,k,1), map(i,k,2)
             elseif(k==3) then; write(43,*) map(i,k,1), map(i,k,2) ;endif
          end do
       end do
       close(42)
       close(43)
       close(44)
       write(*,*) 'Scatter plot written to file = ', trim(outfile)
       
    else if (pol==1) then
       call scattercalc(map(:,1,1), map(:,1,2), n, spectral(1), data(1)%freq(1:2),notused, uncert(1))
       call bootstrap_scatter(map(:,1,1), map(:,1,2), n, nb, tmpspectral, data(1)%freq(1:2), tmpuncert)
       uncert(1)=sqrt(uncert(1)**2+tmpuncert**2)

       ! Write scatter plot to file
       outfile = trim(outprefix) // '_scatter_1.txt'
       unit = getlun()
       open(unit, file=trim(outfile))
       do i = 1, regpix
          do k = 1,1
             write(unit,*) map(i,k,1), map(i,k,2)
          end do
       end do
       close(unit)
       write(*,*) 'Scatter plot written to file = ', trim(outfile)

    end if


    deallocate(map)

  end subroutine scatterplot


  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------
  subroutine ellipse_scatter(map1, map2, n, spectral, freq, uncert_tot, uncert_stat, uncert_syst)
    ! calculates the spectral index by fitting an ellipse to the scatter-data
    implicit none

    real(dp),          dimension(:), intent(in)  :: map1, map2, freq
    integer(i4b),                    intent(in)  :: n
    real(dp),                        intent(out) :: spectral
    real(dp),          optional,     intent(out) :: uncert_stat, uncert_syst, uncert_tot

    integer(i4b)                                 :: i, j
    real(dp)                                     :: healnan=-1.6375d30 
    real(dp)                                     :: mu(2), a, phi, dphi, solutions, dtheta, sigma_stat
    real(dp), allocatable, dimension(:,:)        :: map, C, Vn, invsq_W, L, u
    real(dp), allocatable, dimension(:)          :: eigen

    allocate(map(2,n))
    allocate(C(2,2))
    !rearrange maps
    map(1,:) = map1
    map(2,:) = map2
    !find mean
    do i=1,2
       mu(i) = sum(map(i,:)) /n
    end do
    !calculate the empirical cov matrix
    do i=1,2
       do j=1,2
          C(i,j) = sum((map(i,:)-mu(i)) * (map(j,:)-mu(j)) ) / n
       end do
    end do

    !eigendecompose and find the slope parameter and spectral index
    allocate(eigen(2), Vn(2,2))
    call get_eigen_decomposition(34, C, eigen, Vn)
!call dump_matrix(C)
!call dump_matrix(Vn)
    a = Vn(2,2)/Vn(2,1)
    spectral = log( a*(ant2thermo(freq(1))/ant2thermo(freq(2))) ) / log(freq(2)/freq(1))


    !diagonalize C
    allocate(invsq_W(2,2), L(2,2), u(2,2))
    invsq_W = 0.d0
    do i=1,2
       invsq_W(i,i) = 1.d0/sqrt(eigen(i))
    end do
!    call cholesky_decompose(C, L)

    !find the uncertainty
!    dtheta = sqrt( (1.d0/(n-2.d0)) * 1.d0/eigen(1)  )   !where the data is highly correlated
    dtheta = sqrt( (1.d0/(n-2.d0)) * (eigen(1)+eigen(2))/((eigen(1)-eigen(2))**2)  ) 
    sigma_stat = (1.d0+a**2) * dtheta
    uncert_stat = sigma_stat / a / log(freq(2)/freq(1)) 
    uncert_tot = uncert_stat

    write(*,*) "ellipse spectral index=", spectral, uncert_tot

!!$    u=matmul(invsq_W,transpose(Vn))
!!$!call dump_matrix(u)
!!$    phi=0.d0
!!$    dphi=1.0d0
!!$    do i=0,90
!!$       phi=i*dphi
!!$       solution=(u(1,1)*cos(phi)+u(2,1)*sin(phi))**2 + (u(1,2)*cos(phi) +u(2,2)*sin(phi))**2
!!$!       write(*,*) solution
!!$       if (solution.gt.1.) write(*,*) "l askfjlaiji"
!!$    end do
!!$!    stop

    deallocate(map, C, eigen, Vn, invsq_W, L)
  end subroutine ellipse_scatter


  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------
  subroutine scattercalc(map1, map2, n, spectral, freq, uncert_tot, uncert_stat, uncert_syst)
    ! least squares method when there is noise in both y and x
    ! called the effective variance method.
    implicit none

    real(dp),          dimension(:), intent(in)  :: map1, map2, freq
    integer(i4b),                    intent(in)  :: n
    real(dp),                        intent(out) :: spectral
    real(dp),          optional,     intent(out) :: uncert_stat, uncert_syst, uncert_tot

    integer(i4b)                                 :: i, j
    real(dp)                                     :: healnan=-1.6375d30 
    real(dp)                                     :: m, dm, Vx, Vy, Cxy, dV, nn, aa,b, dtheta, theta, dp2, sigma_syst, sigma_stat, m_syst(2)
    real(dp)                                     :: mu1, mu2


    !y=map2, x=map1   y = mx + b
    nn = n*1.d0
    Vx  = sum(map1*map1)/nn - sum(map1)/nn*sum(map1)/nn
    Vy  = sum(map2*map2)/nn - sum(map2)/nn*sum(map2)/nn 
    Cxy = sum(map1*map2)/nn - sum(map1)/nn*sum(map2)/nn
    dV  = Vx-Vy
    aa  = (Vy - Vx) / ( 2*Cxy )

!    m = ( (Vy-Vx) + sqrt(4*Cxy*Cxy + dV*dV) ) / (2.d0 * Cxy)
!    m = Cxy / Vx  !the method without noise in both y and x

    m = aa + (Cxy/abs(Cxy)) * sqrt(1+aa*aa)
    b = sum(map2)/nn - m * sum(map1)/nn
    
    spectral = log( m*(ant2thermo(freq(1))/ant2thermo(freq(2))) ) / log(freq(2)/freq(1))
    if(m.lt.0.d0) spectral = - log( -m*(ant2thermo(freq(1))/ant2thermo(freq(2))) ) / log(freq(2)/freq(1))  !if m is negative, log(negative number) is healnan!

!write(*,*) "antthermo,slfjiweoi", ant2thermo(freq(1))/ant2thermo(freq(2)), m,b
    if (present(uncert_stat) .or. present(uncert_tot)) then
       dtheta = sqrt( (1.d0/(nn-2.d0)) * (Vx+Vy)/ (dV*dV + 4.d0*Cxy*Cxy) ) 
       sigma_stat = sqrt((1+m**2)**2 * dtheta**2)
       if (present(uncert_stat)) uncert_stat = sqrt(sigma_stat**2 / m**2 / log(freq(2)/freq(1))**2 )
    end if

    if (present(uncert_syst) .or. present(uncert_tot)) then
       ! Divide points into two components, above and below the best-fit line
    
       do j = 1, 2
          nn = 0; Vx = 0; Vy = 0; Cxy = 0; mu1 = 0; mu2 = 0
          do i = 1, n
             if ((j == 1 .and. map2(i) > b + m*map1(i)) .or. (j == 2 .and. map2(i) < b + m*map1(i))) then
                nn  = nn  + 1
                mu1 = mu1 + map1(i)
                mu2 = mu2 + map2(i)
                Vx  = Vx  + map1(i)*map1(i)
                Vy  = Vy  + map2(i)*map2(i)
                Cxy = Cxy + map1(i)*map2(i)
             end if
          end do
          mu1 = mu1 / nn
          mu2 = mu2 / nn
          Vx  = Vx / nn - mu1**2
          Vy  = Vy / nn - mu2**2
          Cxy = Cxy / nn - mu1*mu2
          dV  = Vx-Vy
          aa  = (Vy - Vx) / ( 2*Cxy )
          m_syst(j) = aa + (Cxy/abs(Cxy)) * sqrt(1+aa*aa)
       end do
       sigma_syst = 0.5 * abs(m_syst(1)-m_syst(2))
       if (present(uncert_syst)) uncert_syst = sqrt(sigma_syst**2 / m**2 / log(freq(2)/freq(1))**2 )
    end if

    if (present(uncert_tot)) uncert_tot = sqrt((sigma_stat**2 + sigma_syst**2) / m**2 / log(freq(2)/freq(1))**2 )

  end subroutine scattercalc

  subroutine bootstrap_scatter(map1, map2, n, bn, spectral, freq, uncert_tot, uncert_stat, uncert_syst)
    implicit none

    real(dp),          dimension(:), intent(in)  :: map1, map2, freq
    integer(i4b),                    intent(in)  :: n, bn
    real(dp),                        intent(out) :: spectral
    real(dp),          optional,     intent(out) :: uncert_stat, uncert_syst, uncert_tot

    integer(i4b)                                 :: unit, i, j, p, loops
    real(dp)                                     :: healnan=-1.6375d30 
    real(dp), dimension(:), allocatable          :: bmap1, bmap2, beta
    real(dp)                                     :: r, betamean, betastd
!    type(planck_rng)                             :: rng

    ! Initialize random generator
!    call rand_init(rng,234702)

    loops = 10000
    allocate(beta(loops))
    allocate(bmap1(bn), bmap2(bn))
    unit = getlun()
    open(unit, file="bootstrap.txt")  
    i=1
    do while (i.le.loops)
       do j=1,bn
          r = rand_uni(rng)
          p= (r * n) +1
          bmap1(j)=map1(p)
          bmap2(j)=map2(p)
       end do
       call scattercalc(bmap1, bmap2, bn, beta(i), freq)!, uncert_tot)
       if(beta(i).eq.healnan) then
          write(*,*) "beta is healnan"
!          cycle
       end if
       write(unit,*) beta(i), i!, uncert_tot
       i=i+1
    end do
    close(unit)
    betamean=sum(beta)/loops
    betastd =sqrt(sum((beta-betamean)**2) / loops )
    write(*,*)"bootstrap_scatter, beta", betamean, " +- ", betastd, "n_bootstrap= ",bn
    uncert_tot = betastd

    deallocate(bmap1,bmap2,beta)

  end subroutine bootstrap_scatter

  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scattercalc2(map1, map2, n, spectral, freq, uncert)
    !ordinary least squares method (from squires), do not use!
    implicit none

    real(dp),          dimension(:), intent(in)  :: map1, map2, freq
    integer(i4b),                    intent(in)  :: n
    real(dp),                        intent(out) :: spectral
    real(dp),          optional,     intent(out) :: uncert

    real(dp)                                     :: healnan=-1.6375d30, E, D, m, F, dm

    !y=map2, x=map1
    D = sum(map1*map1) - (sum(map1))**2/n*1.d0
    E = sum(map1*map2) - (sum(map1)*sum(map2))/n*1.d0
    F = sum(map2*map2) - (sum(map2)**2)/n*1.d0

    m = E / D

    spectral = log( m*(ant2thermo(freq(1))/ant2thermo(freq(2))) ) / log(freq(2)/freq(1))

    dm = sqrt( (D*F-E*E)/(D*D)/(n*1.d0-2.d0) )
    uncert = dm / m


    write(*,*) "spectral**",m, spectral, uncert


  end subroutine scattercalc2

  !-----------------------------------------------------------------------------------------------
  ! subroutine ratioplot
  !-----------------------------------------------------------------------------------------------

  subroutine ratioplot(myid, data, outprefix)
    implicit none

    type(banddata)                  :: data
    integer(i4b),     intent(in)    :: myid
    character(len=*), intent(in)    :: outprefix

    real(dp),          allocatable, dimension(:)   :: ratio, kart1, kart2
    real(dp),          allocatable, dimension(:,:) :: outmap
    character(len=256)                             :: outfile
    integer(i4b)                                   :: i, n
    real(dp)                                       :: freq1, freq2
    real(dp)                                       :: healnan=-1.6375d30, root=0

    n=data%n
    freq1 = data%freq(1)
    freq2 = data%freq(2)

    ! Extract maps
    allocate(kart1(data%n))
    allocate(kart2(data%n))
    if (data%pol==1) then 
       kart1 = data%map(:, 1, 1)
       kart2 = data%map(:, 1, 2)
    else if (data%pol==2) then
       kart1 = sqrt(data%map(:, 2, 1)**2 + data%map(:, 3, 1)**2)
       kart2 = sqrt(data%map(:, 2, 2)**2 + data%map(:, 3, 2)**2)
    end if

    
    ! Calculate naive ratio
!    if (myid==0) then 
       allocate(ratio(data%n))
       ratio = kart1/kart2
       write(*,*) log(sum(ratio)/n)/log(freq1/freq2), '= spectral index from average ratio'
       ratio = log(ratio)/log(freq1/freq2)
       write(*,*) sum(ratio)/n, '= average spectral index'
       ! Writing spectral index map to file
       allocate(outmap(0:data%npix-1,data%nmaps))
       outmap = healnan
       outmap(data%pixels,1) = ratio
       outfile = trim(outprefix) // '_spectral_index_from_ratio.fits'
       call write_map(outmap, data%ordering, outfile)
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Writing ratio-based mask to file
!       outmap = 0.d0
!       do i = 1, n
!          if (ratio(i)>-3.7d0 .and. ratio(i)<-2.3d0) then
!             outmap(data%pixels(i),:) = 1.d0
!          end if
!       end do
!       outfile = trim(outprefix) // '_mask_from_ratio.fits' 
!       call write_map(outmap, data%ordering, outfile)
!       write(*,*) '* Mask from ratio written to file = ', trim(outfile)
       deallocate(ratio)    
       deallocate(kart1)    
       deallocate(kart2)    
       deallocate(outmap)
!    end if

  end subroutine ratioplot

  !-----------------------------------------------------------------------------------------------
  ! subroutine read_banddata
  !-----------------------------------------------------------------------------------------------

  subroutine read_banddata(dir, data, myid, domaxlike)
    implicit none

    character(len=*), intent(in)    :: dir
    type(banddata),   intent(out)   :: data
    integer(i4b),     intent(in)    :: myid
    character(len=512)                                  :: infofile, mapname, covname, icovname, notinuse2
    character(len=512), dimension(:),       allocatable :: mapnames, covnames
    character(len=5)                                    :: i_char
    real(dp),           dimension(:,:),     allocatable :: kart, mat, Vn
    real(dp),           dimension(:),       allocatable :: eigen
    integer(i4b),       dimension(:),       allocatable :: pixels
    integer(i4b)                                        :: unit, i, nside, ordering, n, notinuse, ncomp
    logical(lgt)                                        :: inv, domaxlike
real(dp)      :: theta, phi

    ! Read infofile
    unit = getlun()
    infofile = trim(dir) // '/info.txt' 
    open(unit, file=trim(infofile))
    read(unit, *) data%num_bands, notinuse, notinuse2, data%pol
    allocate(data%freq(data%num_bands))
    allocate(mapnames(data%num_bands))
    allocate(covnames(data%num_bands))
    do i = 1, data%num_bands
       read(unit, *) data%freq(i), mapname, covname, icovname
       if (myid==0) write(*,*) 'Freqency',i,'is',data%freq(i)
       mapnames(i) = trim(dir) // "/" // trim(mapname)
       covnames(i) = trim(dir) // "/" // trim(covname)
    end do
    close(unit)

    ! Checking polariation
    if (data%pol==1) then
       data%nmaps = 1
       if (myid==0) write(*,*) 'Running at temperature data' 
    else if (data%pol==2) then
       data%nmaps = 3
       if (myid==0) write(*,*) 'Running at polarisation data' 
    else
       write(*,*) data%pol, '= pol. Unknown number. Quiting'
       stop
    end if

    ! Read the first map to get the pixels etc.
    if (myid==0) write(*,*) trim(mapnames(1))
    call read_map(kart, data%pixels, data%nside, data%ordering, trim(mapnames(1)))
    data%npix  = 12*data%nside**2
    data%n     = size(kart, 1)
    allocate(data%map(data%n, data%nmaps, data%num_bands))
    data%map(:,:,1) = kart(:,1:data%nmaps)
    deallocate(kart)

    ! Then all the rest
    do i = 2, data%num_bands
       if (myid==0) write(*,*) trim(mapnames(i))
       call read_map(kart, pixels, nside, ordering, trim(mapnames(i)))
       call assert(size(kart,1) == data%n  .and. nside == data%nside &
            & .and. ordering == data%ordering .and. all(pixels == data%pixels), &
            & "Error in read bandata: Map " // trim(itoa(i)) // " is inconsistent!")
       data%map(:,:,i) = kart(:,1:data%nmaps)
       deallocate(kart, pixels)
    end do
    deallocate(mapnames)


    ! If input maps are ringed, then convert to nest
    if (data%ordering ==1) then 
       if (myid==0) write(*,*) 'Changing ordering from ring to nest'
       call convert_ring2nest_sparse(data%nside, data%pixels)
       data%ordering = 2
    end if

!!$!!find pixels less than -60muK
!!$    do i=1,data%n
!!$       do notinuse=2,3
!!$          if (data%map(i,notinuse,1) > 70. .or. data%map(i,notinuse,2) > 70.) then
!!$             call pix2ang_nest(nside, data%pixels(i), theta, phi)
!!$             theta=theta*RAD2DEG; phi = phi * RAD2DEG
!!$             write(*,*)"theta, phi",i,theta,phi
!!$          end if
!!$       end do
!!$    end do
!!$!!

!    !change data with gaussian distributed data
!    do i=1,data%n
!       do notinuse=2,3
!          data%map(i,notinuse,1) = rand_gauss(rng)*40.d0
!          data%map(i,notinuse,2) = rand_gauss(rng)*40.d0
!       end do
!    end do

    ! And bulid map2mask
    allocate(data%map2mask(0:data%npix-1))
    data%map2mask = -1
    do i = 1, data%n
       data%map2mask(data%pixels(i)) = i
    end do

    ! Write some more info to screen
    if (myid==0) then 
       write(*,*) data%n,'= n'
       write(*,*) data%ordering,'= ordering'
       write(*,*) data%nside,'= nside'
    end if

    ! Read covmatrices only if maxlike method
    if(domaxlike) then
       !    allocate(data%cov(data%n,data%pol,data%n,data%pol, data%num_bands))
       allocate(data%cov(data%n*data%pol,data%n*data%pol,data%num_bands))
       do i = 1, data%num_bands
          if (myid==0) write(*,*) trim(covnames(i))
          call read_covmat(mat, pixels, nside, ordering, ncomp, covnames(i))
          call assert(size(mat,1) == data%n*data%pol .and. size(mat,2) == data%n*data%pol & 
               & .and. nside == data%nside .and. ordering == data%ordering &
               & .and. all(pixels == data%pixels), &
               & "Error in read bandata: Cov " // trim(itoa(i)) // " is inconsistent!")
          data%cov(:,:,i) = mat
          
          !write eigenvalues to file
          allocate(eigen(data%n*data%pol), Vn(data%n*data%pol,data%n*data%pol))
          call get_eigen_decomposition(34, mat, eigen, Vn)
          call sindex_int2string(i,i_char)
          call dump_matrix(eigen,trim(dir) // '/eigen_band'// trim(i_char) //'.txt')
          deallocate(eigen, Vn)

          deallocate(mat, pixels)

       end do
       deallocate(covnames)
    end if


  end subroutine read_banddata

!!$  subroutine free_banddata(data)
!!$    implicit none
!!$    type(banddata) :: data
!!$    if(allocated(data%map))   deallocate(data%map)
!!$    if(allocated(data%cov))   deallocate(data%cov)
!!$    if(allocated(data%pixels))   deallocate(data%pixels)
!!$    if(allocated(data%map2mask)) deallocate(data%map2mask)
!!$    if(allocated(data%freq))     deallocate(data%freq)
!!$
!!$  end subroutine

  !-----------------------------------------------------------------------------------------------
  ! subroutine give_user_info
  !-----------------------------------------------------------------------------------------------

  subroutine give_user_info
    implicit none

    if (myid == root) then
       write(*,*) 'Usage: mpirun -n N sindex index parfile nregions'
    end if
    call mpi_finalize(ierr)
    stop

  end subroutine give_user_info

  !---------------------------------------------------------------------
  ! Find spectral index from maximum likelihood, grid method
  !----------------------------------------------------------------------

  subroutine maxlike_grid(map1, map2, icov, freq, spectral, offset, uncert, n, fg, polbeta, poloff, hpindex)
    implicit none

    real(dp),             dimension(:),   intent(in)  :: map1, map2, freq
    real(dp),             dimension(:,:,:), intent(in)  :: icov
    real(dp),                             intent(out) :: spectral, offset, uncert
    integer(i4b),                         intent(in)  :: n
    integer(i4b), optional,               intent(in)  :: hpindex
    real(dp),     optional, dimension(:), intent(out) :: fg
    real(dp),     optional,               intent(in)  :: polbeta, poloff

    integer(i4b)                                 :: i, j, p, v, unit, num, int, k, l, i1,i2
    integer(i4b)                                 :: nfreq, ngrid   
    real(dp)                                     :: healnan=-1.6375d30, tall, beta, finebeta, fineoffset
    real(dp)                                     :: betamin, betamax, offmin, offmax, freq0
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    integer(i4b),    allocatable, dimension(:)   :: in2red
    real(dp),        allocatable, dimension(:)   :: pos, m, A, Ateller, a2t, prob, chi_vector
    real(dp),        allocatable, dimension(:,:) :: grid, y, nevn
    logical(lgt)  :: chatty

    chatty = .false.

    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if
    num = size(map1)
    if (chatty) write(*,*) num,'=num', n ! size(fg),'= size(fg)', n

    write(*,*) "Calculating the chi square with grid"

    ! Put data into place
    nfreq=size(freq)
    freq0 = freq(1)
    allocate(a2t(nfreq))
    do v = 1, nfreq
       a2t(v)=ant2thermo(freq(v))
    end do
    allocate(m(nfreq))
    m(1)=0.d0 
    allocate(y(n,nfreq))
    allocate(in2red(n))
    j=0
    do i = 1, num
       if ( healok(map1(i)) .and. healok(map2(i))) then
          j = j+1
          y(j,1) = map1(i)
          y(j,2) = map2(i)
          in2red(j) = i
       end if
    end do

    ! Calculate course grid = chi_square
    ngrid   = 100
    betamin = -10.d0
    betamax = 10.d0
    offmin  = -100.d0
    offmax  = 100.d0
    allocate(grid(ngrid, ngrid))
    allocate(A(n))
    allocate(Ateller(n))
    allocate(nevn(n,n))
    allocate(chi_vector(n))
    allocate(pos(nfreq))

    do k = 1, 2

       ! If you want an additional loop, loop from 1 to 3, uncomment this and set the one below for k=3.
!!$       if (k == 2) then ! Calculate fine grid = chi_square
!!$          finebeta = 1.d0
!!$          fineoffset = 10.d0
!!$          beta = (betamax-betamin)/(ngrid-1)*(pos(2)-1) + betamin
!!$          m(2) = (offmax-offmin)/(ngrid-1)*(pos(1)-1) + offmin
!!$          betamin = beta -finebeta
!!$          betamax = beta +finebeta
!!$          offmin  = m(2) -fineoffset
!!$          offmax  = m(2) +fineoffset
!!$          if (chatty) write(*,*) betamin, 'betamin, betamax', betamax
!!$          if (chatty) write(*,*) offmin, 'offmin, offmax', offmax
!!$       end if

       if (k == 2) then ! Calculate grid 3 based on confidence levels on chi_square
            ! For a Gaussian distribution, the 1, 2 and 3 sigma (68%, 95% and 99.7%) confidence regions 
            ! correspond to where chi_square = -2 lnL increases by 2.3, 6.17 and 11.8 from its minimum value.
          do j = pos(2),ngrid      ! loop over beta to find where grid = 3sigma
             if ( grid(pos(1),j) >= 11.8 + minval(grid) ) exit
          end do
          do i = pos(1),ngrid      ! loop over m to find where grid = 3 sigma 
             if ( grid(i,pos(2)) >= 11.8 + minval(grid) ) exit
          end do
          finebeta = (betamax-betamin)/(ngrid-1)*(j - pos(2))
          fineoffset = (offmax-offmin)/(ngrid-1)*(i - pos(1))
          beta = (betamax-betamin)/(ngrid-1)*(pos(2)-1) + betamin
          m(2) = (offmax-offmin)/(ngrid-1)*(pos(1)-1) + offmin
          betamin = beta -finebeta
          betamax = beta +finebeta
          offmin  = m(2) -fineoffset
          offmax  = m(2) +fineoffset
          if (chatty) write(*,*) grid(pos(1),j), 'grid(pos(1),j), grid(i,pos(2))', grid(i,pos(2))
          if (chatty) write(*,*) j- pos(2), 'number grid beta, number grid m', i - pos(1)
          if (chatty) write(*,*) betamin, 'betamin, betamax', betamax
          if (chatty) write(*,*) offmin, 'offmin, offmax', offmax       
       end if

       grid= 0.d0
       do j = 1, ngrid
          beta = (betamax-betamin)/(ngrid-1)*(j-1) + betamin
          do i = 1, ngrid
             m(2) = (offmax-offmin)/(ngrid-1)*(i-1) + offmin
             Ateller = 0.d0
             nevn = 0.d0
             do v = 1, nfreq
                Ateller(:) = Ateller(:) + matmul( (y(:,v)-m(v))*a2t(v)*(freq(v)/freq0)**beta , icov(:,:,v) )
                nevn(:,:) = nevn(:,:) + ( ((a2t(v)**2)*(freq(v)/freq0)**(2.d0*beta)) * icov(:,:,v) )
             end do

             ! A = Ateller/nevn -> nevn * A = Ateller 
             call solve_system_real(nevn,A,Ateller)
             do v = 1, nfreq
                chi_vector = y(:,v) - A(:)*a2t(v)*(freq(v)/freq0)**beta - m(v)
                grid(i,j) = grid(i,j) + dot_product( chi_vector, matmul( icov(:,:,v) , chi_vector ) )
             end do
          end do
       end do

       !write out a 2 dim contour plot of chisquare values
       if(chatty) then
          open(55,file="chisquare.txt")
          do j = 1, ngrid
             do i = 1, ngrid
                beta = (betamax-betamin)/(ngrid-1)*(j-1) + betamin
                m(2) = (offmax-offmin)/(ngrid-1)*(i-1) + offmin
                write(55,fmt="(3e15.7)") beta,m(2),grid(i,j)
             end do
             write(55,*)
          end do
          close(55)
       end if

       pos=minloc(grid)
       if (chatty) write(*,*) 'positionon grid ',k,': ', pos(1), pos(2)
       if (pos(1)==1 .or. pos(1)==ngrid .or. pos(2)==1 .or. pos(2)==ngrid) then
          spectral = healnan
          offset   = healnan
          uncert   = healnan
          exit
       end if
    end do


    ! Check that we didn't not hit the grid
    if (pos(1)==1 .or. pos(1)==ngrid .or. pos(2)==1 .or. pos(2)==ngrid) then
       uncert   = healnan
       spectral = healnan
       offset   = healnan
       if (present(fg)) fg = healnan
    else
       ! Find spectral index and offset
       beta = (betamax-betamin)/(ngrid-1)*(pos(2)-1) + betamin
       spectral = beta
       m(2) = (offmax-offmin)/(ngrid-1)*(pos(1)-1) + offmin
       offset = m(2)
       if (chatty) write(*,*) 'Spectral index', spectral
       if (chatty) write(*,*) 'Offset', offset
       ! Calculate uncertainty
       allocate(prob(ngrid))
       tall = minval(grid)
       grid = exp(-(grid-tall)/2.d0)
       prob=sum(grid,1)
       prob=prob/sum(prob)
       prob = -2.d0*log(prob)
       if (chatty) write(*,*) 'prob', prob(pos(2)), minval(prob), pos(2)
       do i = pos(2),ngrid
          if (prob(i) >= minval(prob) + 1.d0) exit
       end do
       int=i
       do i = pos(2),1, -1
          if (prob(i) >= minval(prob) + 1.d0) exit
       end do
       if (int > ngrid .or. i < 1) then
          uncert = healnan
       else
          uncert = (betamax-betamin)/(ngrid-1)*(int-i)/2.d0
       end if
       if (uncert==0.d0) uncert=healnan
       if (chatty) write(*,*) 'delta beta', uncert, int,i
       if (chatty) write(*,*)
       deallocate(prob)
       ! Calculate foreground
       if (present(fg)) then
          if (present(polbeta) .and. present(poloff)) then
             beta = polbeta    ! obs changing to pol(=Q+U) values
             m(2) = poloff     ! obs changing to pol(=Q+U) values
          end if
          Ateller = 0.d0
          nevn = 0.d0
          do v = 1, nfreq
             Ateller(:) = Ateller(:) + matmul( (y(:,v)-m(v))*a2t(v)*(freq(v)/freq0)**beta , icov(:,:,v) )
             nevn(:,:) = nevn(:,:) + ( ((a2t(v)**2)*(freq(v)/freq0)**(2.d0*beta)) * icov(:,:,v) )
          end do
          ! A = Ateller/nevn -> nevn*A=Ateller
          call solve_system_real(nevn,A,Ateller)
          A = A*ant2thermo(freq0)
          fg(in2red)=A(1:size(fg))       !!this is somewhat fishy, polbeta and poloff is not set?
          if (chatty) write(*,*) 'fg', fg(in2red(1)), fg(1)
       end if
    end if

    if (present(hpindex)) then
       ! Write contour  plot to file
       call int2string(hpindex, filnum)
       filename='contour_'//filnum//'.txt'
       open(42, file=trim(filename))
       write(42,*) ngrid, ngrid
       do j = 1, ngrid
          beta = (betamax-betamin)/(ngrid-1)*(j-1) + betamin
          do i = 1, ngrid
             m(2) = (offmax-offmin)/(ngrid-1)*(i-1) + offmin
             write(42,*) m(2), beta, grid(i,j)
          end do
       end do
       close(42)
       write(*,*) 'written to file = ', trim(filename)
    end if

    
    ! Clean up
    deallocate(A)
    deallocate(Ateller)
    deallocate(nevn)
    deallocate(chi_vector)
    deallocate(pos)
    deallocate(grid)
    deallocate(y)
    deallocate(m)
    deallocate(a2t)
 
  end subroutine maxlike_grid

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sindex_int2string(region, region_char)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer(i4b),       intent(in)    :: region 
    character(len=5),   intent(out)   :: region_char

    character(len=1)    :: char_1
    character(len=2)    :: char_2
    character(len=3)    :: char_3
    character(len=4)    :: char_4
    character(len=5)    :: char_5

    if (region.LT.10) then
       call int2string(region,char_1)
       region_char = char_1
    else if (region.GT.9.AND.region.LT.100) then
       call int2string(region,char_2)
       region_char = char_2
    else if (region.GT.99.AND.region.LT.1000) then
       call int2string(region,char_3)
       region_char = char_3
    else if (region.GT.999.AND.region.LT.10000) then
       call int2string(region,char_4)
       region_char = char_4
    else if (region.GT.9999.AND.region.LT.100000) then
       call int2string(region,char_5)
       region_char = char_5
    end if

  end subroutine sindex_int2string



end program sindex
