program mdfit
  use quiet_utils
  use fitstools
  use powell_mod
  use rngmod
  use quiet_mapfile_mod
  use pix_tools
  !use quiet_mpi_mod
  implicit none
  
!  include 'mpif.h'
  
  integer(i4b)       :: myid, ierr, numprocs, root, i, j, k, l, n, iargc
  integer(i4b)       :: nside, nband, nscatter, npix, ordering, unit, niter, seed, nstep, ntemp, numval, iter
  integer(i4b)       :: nside_map, npix_map, p, q, myloc(1), npos, nsim
  logical(lgt)       :: anynull, enforce_positivity, exist, ok, debug
  real(dp)           :: vec(3), nullval, d, b_lim, mcmc_rms, s, prior(2), theta, phi, mu_res, sigma_res, subsum, sigma_n
  character(len=256) :: paramfile, filename, outfile, maskfile
  character(len=8000) :: md_input
  character(len=2)   :: itext, jtext
  real(dp),     allocatable, dimension(:)     :: md, md_prop, md_fix, md_init, map_pos, nu
  real(dp),     allocatable, dimension(:,:)   :: T, T_map, offset, md_tot, maps, rms, mask, buffer, outmap, betamap, maps_full
  real(dp),     allocatable, dimension(:,:,:) :: spar
  integer(i4b), allocatable, dimension(:,:)   :: sind
  integer(i4b), allocatable, dimension(:)     :: pix, pixlist, i2x
  type(planck_rng) :: rng_handle
  real(dp), allocatable, dimension(:,:) :: A, A_red, M, M_red, invN, A_pos
  real(dp), allocatable, dimension(:)   :: b, x, x_red, res, b_pos

!  call mpi_init(ierr)
!  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
!  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)

!  root = 0
  unit = getlun()

  if (iargc() /= 1) then
     write(*,*) ''
     write(*,*) 'Usage: mdfit [parameter file]'
     write(*,*) ''
     
!     call mpi_finalize(ierr)
     stop
  end if
 
  call getarg(1,paramfile)
  call get_parameter(unit, paramfile, 'DEBUG',                par_lgt=debug)
  call get_parameter(unit, paramfile, 'NSIDE',                par_int=nside)
  call get_parameter(unit, paramfile, 'NSIDE_MAP',            par_int=nside_map)
  call get_parameter(unit, paramfile, 'NBAND',                par_int=nband)
  call get_parameter(unit, paramfile, 'NSCATTER',             par_int=nscatter)
  call get_parameter(unit, paramfile, 'OUTFILE',              par_string=outfile)
  call get_parameter(unit, paramfile, 'ENFORCE_POSITIVITY',   par_lgt=enforce_positivity)
  call get_parameter(unit, paramfile, 'MIN_GAL_LATITUDE',     par_dp=b_lim)
  call get_parameter(unit, paramfile, 'NUM_BOOTSTRAP',        par_int=nsim)
  call get_parameter(unit, paramfile, 'SEED',                 par_int=seed)
  call get_parameter(unit, paramfile, 'MCMC_RMS',             par_dp=mcmc_rms)
  call get_parameter(unit, paramfile, 'BETA_PRIOR_LOW',       par_dp=prior(1))
  call get_parameter(unit, paramfile, 'BETA_PRIOR_HIGH',      par_dp=prior(2))
  call get_parameter(unit, paramfile, 'MASKFILE',             par_string=maskfile)
  call get_parameter(unit, paramfile, 'SCALE',                par_dp=s)
  call get_parameter(unit, paramfile, 'NOISE_CUTOFF',         par_dp=sigma_n)
  b_lim = cos((90.d0-b_lim)*pi/180.d0); if (b_lim < 1e-8) b_lim = 0.d0
  npix         = 12*nside**2
  npix_map     = 12*nside_map**2

  call rand_init(rng_handle, seed)
  !call initialize_random_seeds(MPI_COMM_WORLD, seed, rng_handle)

  ntemp = 4
  allocate(mask(0:npix_map-1,1))
  if (trim(maskfile) /= 'none') then
     j = getsize_fits(maskfile, ordering=ordering)
     call read_bintab(maskfile, mask(:,1:1), npix_map, 1, nullval, anynull)
     if (ordering == 1) call convert_ring2nest(nside_map,mask(:,1))
  else
     mask = 1.d0
  end if
  
  allocate(spar(0:npix-1,2,nscatter), sind(nscatter,2), md_init(nband*ntemp), md(nband*ntemp), md_prop(nband*ntemp), T(0:npix-1,ntemp), md_tot(nband*ntemp,nsim), md_fix(nband*ntemp), A_pos(10000, nband*ntemp), b_pos(10000))
  allocate(maps(0:npix_map-1,nband), rms(0:npix_map-1,nband), T_map(0:npix_map-1,ntemp), nu(nband))
  allocate(maps_full(0:npix_map-1,nband))
  allocate(buffer(0:npix_map-1,1), map_pos(0:npix_map-1))
  allocate(pixlist(0:npix_map-1))
  md_init = 0.d0
  A_pos   = 0.d0

  ! Precompute templates
  T(:,1) = 1.d0
  do i = 0, npix-1
     call pix2vec_nest(nside, i, vec)
     T(i,2:4) = vec
  end do

  if (.false.) then
     q = npix_map / npix
     do k = 1, 4
        do i = 0, npix-1
           T_map(i*q:(i+1)*q-1,k) = T(i,k)
        end do
     end do
  else
     T_map(:,1) = 1.d0
     do i = 0, npix_map-1
        call pix2vec_nest(nside_map, i, vec)
        T_map(i,2:4) = vec
     end do
  end if

  ! Read band data
  npos = 0
  do i = 1, nband
     call int2string(i,itext)
     call get_parameter(unit, paramfile, 'MAP' // itext, par_string=filename)
     j = getsize_fits(filename, nside=nside_map, ordering=ordering)
     call read_bintab(filename, maps(:,i:i), npix_map, 1, nullval, anynull)
     if (ordering == 1) call convert_ring2nest(nside_map,maps(:,i))
     maps_full(:,i) = maps(:,i)

     call get_parameter(unit, paramfile, 'FREQ' // itext, par_dp=nu(i))
     
     if (allocated(mask)) then
        where (mask(:,1) < 0.5d0) 
           maps(:,i) = -1.6375d30
        end where
     end if

!!$     if (i == 2) then
!!$        write(*,*) 'Warning: subtracting 100 uK'
!!$        where (maps(:,i) /= -1.6375d30) 
!!$           maps(:,i) = maps(:,i) - 100.d0 * T_map(:,1) + 30.d0*T_map(:,2) - 60.d0*T_map(:,3) + 90.d0 * T_map(:,4)
!!$           !maps(:,i) = maps(:,i) - 30.d0*T_map(:,2)
!!$        end where
!!$     else
!!$        where (maps(:,i) /= -1.6375d30) 
!!$           maps(:,i) = maps(:,i) + 50.d0 * T_map(:,1) + 10.d0*T_map(:,2) + 20.d0*T_map(:,3) - 15.d0 * T_map(:,4)
!!$           !maps(:,i) = maps(:,i) - 30.d0*T_map(:,2)
!!$        end where
!!$     end if

     call get_parameter(unit, paramfile, 'RMS' // itext, par_string=filename)
     j = getsize_fits(filename, nside=nside_map, ordering=ordering)
     call read_bintab(filename, rms(:,i:i), npix_map, 1, nullval, anynull)
     if (ordering == 1) call convert_ring2nest(nside_map,rms(:,i))
     
     call get_parameter(unit, paramfile, 'MDFIX' // itext, par_string=md_input)
     read(md_input,*) md_fix((i-1)*ntemp+1:i*ntemp)

     ! Locate the coldest pixels, separated by at least 10 degrees. Allow for 
     ! some negative values created by noise, thresholded by user
     where (maps(:,i) /= -1.6375d30) 
        map_pos = maps(:,i) - sigma_n * rms(:,i)
     elsewhere
        map_pos = -1.6375d30
     end where
     do while (any(map_pos /= -1.6375d30))
        myloc = minloc(map_pos, map_pos /= -1.6375d30)
        k     = myloc(1)-1
        npos  = npos+1
        b_pos(npos) = map_pos(k)
        A_pos(npos,(i-1)*ntemp+1) = 1.d0
        call pix2vec_nest(nside_map, k, vec)
        A_pos(npos,(i-1)*ntemp+2:(i-1)*ntemp+4) = vec
        call query_disc(nside_map, vec, 10.d0*pi/180.d0, pixlist, k, nest=1)
        map_pos(pixlist(0:k-1)) = -1.6375d30
     end do
  end do
  prior = (nu(2)/nu(1))**prior
  prior(1) = -1.d30
  prior(2) =  1.d30


  k = count(md_fix == 1.d30)
  allocate(i2x(k))
  k = 1
  do i = 1, nband
     do j = 1, ntemp
        if (md_fix((i-1)*ntemp+j) == 1.d30) then
           i2x(k) = (i-1)*ntemp+j
           k      = k+1
        end if
     end do
  end do

  md_tot = 0.d0
  md     = 0.d0
  iter = 0
  !do iter = 1, 5
  ok = .false.
  do while (.not. ok .and. iter < 10)
     iter = iter+1

     ! Read the scatter plot data
     write(*,*) 'Evaluating scatter plots'
     do i = 1, nscatter
        call int2string(i,itext)
        call get_parameter(unit, paramfile, 'XBAND' // itext, par_int=sind(i,1))
        call get_parameter(unit, paramfile, 'YBAND' // itext, par_int=sind(i,2))
        
        q = npix_map / npix
        spar = -1.6375d30
        do j = 0, npix-1
           if (abs(T(j,4)) < b_lim) cycle
           if (mod(j,10) == 0 .and. debug) write(*,*) 'Computing scatter plot ', j, ' of ', npix
           call scattercalc(maps(j*q:(j+1)*q-1,sind(i,2)), maps(j*q:(j+1)*q-1,sind(i,1)), &
                & mask(j*q:(j+1)*q-1,1), q, spar(j,1,i), spar(j,2,i))

!           if (abs(phi) < 0.1d0 .and. abs(theta-0.5d0*pi) < 0.1d0) then
           if (.false.) then
              call pix2ang_nest(nside, j, theta, phi)
              write(*,*) 'l = ', RAD2DEG * (phi)
              write(*,*) 'b = ', RAD2DEG * (0.5d0*pi-theta)
              open(58,file='scatter.dat')
              write(*,*) spar(j,:,i), log(spar(j,1,i))/log(nu(2)/nu(1))
              do k = 0, q-1
                 if (maps(j*q+k,sind(i,2)) /= -1.6375d30) then
                    write(58,*) maps(j*q+k,sind(i,1)), maps(j*q+k,sind(i,2))
                 end if
              end do
              close(58)
              stop
           end if
           if (j == 2) stop
        end do
     end do

     open(58,file='slope.dat')
     do i = 0, npix-1
        if (spar(i,1,1) /= -1.6375d30) write(58,*) i, spar(i,:,1)
     end do
     close(58)
     stop

     ! Find number of valid pixels
     if (iter == 1) then
        numval = 0
        do i = 0, npix-1 
           if (any(spar(i,1,:) == -1.6375d30) .or. &
                & spar(i,1,1) < prior(1) .or. spar(i,1,1) > prior(2)) then
              spar(i,1,:) = -1.6375d30
           else
              numval = numval+1
           end if
        end do
        allocate(pix(numval))
        j = 1
        do i = 0, npix-1 
           if (any(spar(i,1,:) == -1.6375d30)) cycle
           pix(j) = i
           j      = j+1
        end do
     end if

     allocate(betamap(0:npix-1,1))
     do i = 0, npix-1
        if (spar(i,1,1) /= -1.6375d30) then
           betamap(i,1) = spar(i,1,1) !log(spar(i,1,1)) / log(nu(2)/nu(1))
        else
           betamap(i,1) = -1.6375d30
        end if
     end do
     if (debug) write(*,*) 'beta = ', mean(betamap(pix,1)), sqrt(variance(betamap(pix,1)))
     call int2string(iter, itext)
     call write_map(betamap, 2, 'beta'//itext//'.fits')
     deallocate(betamap)
     
     ! Set up basic equation set
     if (iter > 1) deallocate(A, x, b, M, res)
     allocate(A(numval*nscatter,nband*ntemp), x(nband*ntemp), b(numval*nscatter))
     allocate(M(nband*ntemp, nband*ntemp), res(numval*nscatter))
     A = 0.d0; b = 0.d0
     do i = 1, nscatter
        do j = 1, numval
           do k = 1, ntemp
              A((i-1)*numval+j,(sind(i,1)-1)*ntemp+k) = -spar(pix(j),1,i) * T(pix(j),k)
              A((i-1)*numval+j,(sind(i,2)-1)*ntemp+k) =                     T(pix(j),k)
           end do
           b((i-1)*numval+j) = spar(pix(j),2,i)
        end do
     end do
     M   = matmul(transpose(A), A)
     call invert_matrix(M)
     x   = matmul(M, matmul(transpose(A), b))
     res = matmul(A,x)-b

     do j = 1, nscatter
        mu_res    = mean(res((j-1)*numval+1:j*numval))
        sigma_res = sqrt(variance(res((j-1)*numval+1:j*numval)))
        do i = 1, numval
           ! Cut at 4 sigma
           if (abs(res((j-1)*numval+i)-mu_res)/sigma_res > 4.d0) then
              spar(pix(i),:,j) = -1.6375d30
              A((j-1)*numval+i,:) = 0.d0
              b((j-1)*numval+i)   = 0.d0
           end if
        end do
     end do

     allocate(betamap(0:npix-1,1))
     do j = 1, nscatter
        call int2string(j, jtext)
        betamap = -1.6375d30
        betamap(pix,1) = res((j-1)*numval+1:j*numval)
        call write_map(betamap, 2, 'res_scatter'//jtext//'_'//itext//'.fits')
        open(58,file='res_hist'//jtext//'_'//itext//'.dat')
        do k = 1, numval
           write(58,*) k, res((j-1)*numval+k)
        end do
        close(58)
     end do
     deallocate(betamap)

     call int2string(iter, itext)
     do j = 1, nscatter
        call int2string(j, jtext)
        open(unit,file='lines'//itext//'_'//jtext//'.dat')
        do i = 0, npix-1
           if (spar(i,1,j) /= -1.6375d30) then
              do k = -100, 100
                 write(unit,*) s*real(k,dp), spar(i,1,1)*s*k+spar(i,2,1)
              end do
              write(unit,*)
           end if
        end do
        close(unit)
        open(unit,file='res_sq'//itext//'_'//jtext//'.dat')
        do k = -100, 100
           subsum = 0.d0
           do l = 0, npix-1
              if (spar(l,1,j) /= -1.6375d30) subsum = subsum + spar(l,1,j)*s*k+spar(l,2,j)
           end do
           mu_res = subsum/numval
           subsum = 0.d0
           do l = 0, npix-1
              if (spar(l,1,j) /= -1.6375d30) subsum = subsum + (spar(l,1,j)*s*k+spar(l,2,j)-mu_res)**2
           end do
           write(unit,*) s*real(k,dp), sqrt(subsum/numval)
        end do
        close(unit)
     end do
     
     close(unit)

     if (enforce_positivity) then
        write(*,*) 'Optimizing with positivity constraint'
        md = 0.d0
        call mcmc_search(rng_handle, md)
        if (debug) write(*,*) 'md mcmc = ', real(md,sp)
        if (debug) write(*,*) 'dist mcmc = ', dist(md)
        !call powell(x_red, dist, ierr)  
     else
        md = x
     end if

     md_tot(:,1) = md_tot(:,1) + md
     do i = 1, nband
        do j = 1, ntemp
           where (maps(:,i) /= -1.6375d30)
              maps(:,i) = maps(:,i) - md((i-1)*ntemp+j) * T_map(:,j)
           end where
           maps_full(:,i) = maps_full(:,i) - md((i-1)*ntemp+j) * T_map(:,j)
        end do
     end do

     do i = 1, nband
        write(*,fmt='(i6,a,4f16.3)') i, ' -- ', md_tot((i-1)*ntemp+1:i*ntemp,1)
        write(unit,fmt='(4f16.8)') md_tot((i-1)*ntemp+1:i*ntemp,1)
     end do

     ok = all(abs(md / md_tot(:,1)) < 1.d-2) .or. all(2.d0*md == md_tot(:,1))

     do i = 1, nband
        call int2string(i,itext)
        filename = 'mdcorrmap'//itext//'.fits'
        call write_map(maps(:,i:i), 2, filename)
        filename = 'mdcorrmap_full'//itext//'.fits'
        call write_map(maps_full(:,i:i), 2, filename)
     end do

  end do

  !deallocate(A, x, b, M)
  !deallocate(md, T, spar, sind)
  
contains

  subroutine grid_search(md)
    implicit none

    real(dp), dimension(1:), intent(inout) :: md
    
    integer(i4b) :: i, j, numbin
    real(dp)     :: d0, x0(2)
    real(dp), allocatable, dimension(:)   :: p
    real(dp), allocatable, dimension(:,:) :: d, limits, x
    
    numbin = 1000

    allocate(limits(2,2), d(numbin,numbin), x(numbin,2), p(size(md)))
    limits(1,1) = -0
    limits(2,1) =  2d7
    limits(1,2) =  4d6
    limits(2,2) =  6d6
    do i = 1, numbin
       do j = 1, 2
          x(i,j) = limits(1,j) + (limits(2,j)-limits(1,j)) * (i-1) / real(numbin-1,dp)
       end do
    end do
    
    p = md
    x0 = 0
    d0 = 1.d30
!    open(68,file='grid.dat')
    do i = 1, numbin
       !if (mod(i,100) == 0) write(*,*) i
       p(1) = x(i,1)
       do j = 1, numbin
          p(5) = x(j,2)
          d(i,j) = dist(p)
          if (d(i,j)<d0) then
             d0    = d(i,j)
             x0(1) = p(1)
             x0(2) = p(5)
          end if
!          write(68,*) real(p(1),sp), real(p(5),sp), real(d(i,j),sp)
       end do
    end do
!    close(68)
!    stop

    md(1) = x0(1)
    md(5) = x0(2)
  end subroutine grid_search

  subroutine mcmc_search(handle, md)
    implicit none

    type(planck_rng), intent(inout) :: handle
    real(dp), dimension(1:), intent(inout) :: md

    integer(i4b) :: i, j, k, n, m
    real(dp) :: d, d0
    real(dp), allocatable, dimension(:) :: p, p0, accept, reject
    real(dp), allocatable, dimension(:), save :: dp

    n = size(md)
    m = 1000000
    allocate(p(n), p0(n), accept(n), reject(n))
    if (.not. allocated(dp)) then
       allocate(dp(n))
       if (iter == 1) then
          dp = mcmc_rms
       else
          dp = 1.d0 !mcmc_rms
       end if
    end if
    
!!$    p  = 0.d0
!!$    open(58,file='d.dat')
!!$    do i = -200000, -150000
!!$       p(n) = i
!!$       write(58,*) i, dist(p)
!!$    end do
!!$    close(58)
!!$    stop

    p0 = md
    d0 = dist(p0)
    do while (d0 == 1.d30)
       do j = 1, n, 4
          if (md_fix(j) == 1.d30) then
             p0(j) = p0(j) - s
             d0    = dist(p0)
          end if
       end do
    end do

    p  = p0
    accept = 0.d0
    reject = 0.d0
    do i = 1, m
       do j = 1, n
          if (md_fix(j) == 1.d30) then
             p(j) = p0(j) + dp(j) * rand_gauss(rng_handle)
          else
             cycle
          end if
          d = dist(p)
       
          if (d < d0)  then
             p0 = p
             d0 = d
             accept(j) = accept(j) + 1.d0
          else
             reject(j) = reject(j) + 1.d0
          end if
       end do
       if (debug .and. mod(i,10000)==0) then
          write(*,*) 'p0 = ', real(p0,sp)
          write(*,*) 'd0 = ', real(d0,sp)
       end if
       if (mod(i,1000) == 0) then
          if (all(accept == 0)) exit
          accept(j) = 0; reject(j) = 0
       end if
    end do
    md = p0
    !md(2:4) = 0.d0
    !md(6:8) = 0.d0

!    write(*,*) 'final dp = ', real(dp,sp)

    deallocate(p, p0)

  end subroutine mcmc_search

  function dist(md)
    use healpix_types
    implicit none
    real(dp), dimension(:), intent(in), optional :: md
    real(dp)                                     :: dist

    integer(i4b) :: i

    if (enforce_positivity) then
       do i = 1, npos
          if (b_pos(i)-sum(A_pos(i,:)*(md_tot(:,1)+md)) < 0.d0) then
!!$             write(*,*) 'A   = ', real(A_pos(i,:),sp)
!!$             write(*,*) 'b   = ', real(b(i),sp)
!!$             write(*,*) 'md  = ', real(md_tot(:,1)+md,sp)
!!$             write(*,*) 'tot = ', real(b_pos(i)-sum(A_pos(i,:)*(md_tot(:,1)+md)),sp)
!!$             write(*,*)
             dist = 1.d30
             return
          end if
       end do
    end if

    dist = sum((matmul(A,md)-b)**2)

  end function dist

  subroutine scattercalc(map1, map2, mask, n, slope, offset)
    implicit none

    real(dp),          dimension(:), intent(in)  :: map1, map2
    real(dp),          dimension(:), intent(in)  :: mask
    integer(i4b),                    intent(in)  :: n
    real(dp),                        intent(out) :: slope, offset

    integer(i4b)                                 :: i, j, p, k, numpairs
    real(dp)                                     :: healnan=-1.6375d30, tall, abugmed, bbugmed
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:)   :: y, beta, vector, ymap, xmap, a, b, abug
    real(dp),        allocatable, dimension(:,:) :: x, matrix2, h
    real(dp)                                     :: vec(2), amed, bmed, xmin, xmax, xlimit
   
    allocate(ymap(n))
    allocate(xmap(n))
    xmap=0.d0
    ymap=0.d0
    p=0
    do i = 1, size(map1)
       if ( healok(map1(i)) .and. healok(map2(i)) .and. mask(i)==1.d0) then
          p = p+1
          ymap(p) = map1(i)
          xmap(p) = map2(i)
       end if
    end do
    if (p < n/4) then
       deallocate(ymap, xmap)
       slope  = -1.6375d30
       offset = -1.6375d30
       return
    end if

 
!!$    allocate(y(p))
!!$    allocate(beta(2))
!!$    allocate(x(p,2))
!!$    allocate(matrix2(2,2))
!!$    matrix2(1,1) = p*1.d0
!!$    matrix2(1,2) = -sum(xmap)
!!$    matrix2(2,1) = matrix2(1,2)
!!$    matrix2(2,2) = sum(xmap*xmap)
!!$    matrix2 = matrix2/(p*sum(xmap*xmap)-(sum(xmap))**2)
!!$    vec(1) = sum(xmap*ymap)
!!$    vec(2) = sum(ymap)
!!$    beta = matmul(matrix2, vec) 

!!$    numpairs = p*(p-1)/2
    allocate(a(numpairs))
!!$    !allocate(abug(numpairs))
    allocate(b(p))
!!$    k = 0
!!$    do i = 1, p
!!$       do j = i+1, p	
!!$ 	     k = k+1
!!$	     a(k) = (ymap(i)-ymap(j))/(xmap(i)-xmap(j))
!!$!	  abug(k) = abs(ymap(i)-ymap(j))/abs(xmap(i)-xmap(j))
!!$!	  if (xmap(i)>xmap(j)) then
!!$!	     a(k) = (ymap(i)-ymap(j))/(xmap(i)-xmap(j))
!!$!	  else
!!$!	     a(k) = (ymap(j)-ymap(i))/(xmap(j)-xmap(i))
!!$!	  end if   
!!$       end do
!!$    end do   
!!$    if (k /= numpairs) then
!!$       write(*,*) 'wrong again!'
!!$       stop
!!$    end if   
!!$    amed = median(a)
!!$!    abugmed = median(abug)
!!$    b = ymap - amed*xmap
!!$    bmed = median(b)
!!$    slope  = amed
!!$    offset = bmed
!!$!    b = ymap - abugmed*xmap
!!$!    bbugmed = median(b)

    slope = sum((xmap(1:p)-mean(xmap(1:p)))*(ymap(1:p)-mean(ymap(1:p)))**2) / sum((xmap(1:p)-mean(xmap(1:p)))**2 * (ymap(1:p)-mean(ymap(1:p))))
    b = ymap - slope*xmap
    offset = median(b)

    open(58,file='scatter.dat')
    do i = 1, p
       write(58,*) xmap(i), ymap(i)
    end do
    close(58)
    write(*,*) slope, offset

    !stop

   ! Clean up
    if (allocated(ymap)) deallocate(ymap)
    if (allocated(xmap)) deallocate(xmap)
    if (allocated(a))    deallocate(a)
    if (allocated(b))    deallocate(b)
!    deallocate(y)
!    deallocate(x)
!    deallocate(beta)
!    deallocate(matrix2)

  end subroutine scattercalc

  
end program mdfit




