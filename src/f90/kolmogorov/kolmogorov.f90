program kolmogorov
   use quiet_utils
   use ziggurat
   use rngmod
   use math_tools
   use sort_utils
   use quiet_fileutils
   use quiet_healpix_mod
   use spline_1d_mod
   use quiet_mpi_mod
   implicit none

   !call multikolmo
   !call gurzadyan
   call arnold

contains

  subroutine multikolmo
    implicit none
    character(len=512)                 :: arg, mapfile, psfile, beamfile, noise_mapfile
    character(len=512)                 :: mapdir, maskfile
    integer(i4b)                       :: i, j, k, l, nsamp, nneigh, order, nside, npix
    integer(i4b)                       :: onpix, nmap, lmax
    integer(i4b)                       :: order2, nside2, nstep, ierr, myid, nproc
    integer(i4b)                       :: avg_npix, nhist, h, unit, err, ncorrect
    type(zig_rng)                      :: rng
    type(spline_type)                  :: corr, kspline
    logical(lgt)                       :: anynull, domask
    real(dp)                           :: dev, m, v, rad, pte, gurz, theta, phi
    real(dp)                           :: nullval
    real(dp),     allocatable          :: maps(:,:), noise(:), ps(:)
    real(dp),     allocatable          :: nmap_tmp(:,:), ps_tmp(:,:), beam(:,:)
    real(dp),     allocatable          :: probs(:,:), mask_tmp(:,:)
    real(dp),     pointer              :: pixwin(:,:)
    integer(i4b), allocatable          :: hists(:,:,:), thists(:,:,:)
    integer(i4b), allocatable, target  :: neighs(:)
    logical(lgt), allocatable          :: mask(:)
    integer(i4b),              pointer :: pix(:)
    domask = .false.
    call getarg(1, mapdir)
    call getarg(2, psfile)
    call getarg(3, beamfile)
    call getarg(4, noise_mapfile)
    call getarg(5, arg); read(arg,*) rad
    call getarg(6, arg); read(arg,*) nsamp
    call getarg(7, arg); read(arg,*) nmap
    call getarg(8, arg); read(arg,*) nhist
    if(iargc() >= 9) then
       call getarg(9, maskfile)
       domask = .true.
    end if

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    call dset(level=1, id=myid)
    call zig_init(rng, 18698265+myid)

    ! Read maps
    i    = getsize_fits(trim(mapdir) // "/sim000.fits", nside=nside, ordering=order)
    npix = 12*nside**2
    allocate(maps(0:npix-1, nmap),noise(0:npix-1), mask(0:npix-1))
    do i = 1, nmap
       call read_bintab(trim(mapdir) // "/sim" // trim(itoa(i-1,3)) // ".fits", &
        & maps(:,i:i), npix, 1, nullval, anynull)
       call dmem("Read map " // trim(itoa(i-1)))
    end do
    call read_sparse_as_full(nmap_tmp, nside2, order2, noise_mapfile)
    call assert(nside==nside2 .and. order == order2, "Inconsistent nside or ordering")
    noise = nmap_tmp(:,1)
    deallocate(nmap_tmp)
    if(domask) then
       call read_sparse_as_full(mask_tmp, nside2, order2, maskfile)
       call assert(nside==nside2 .and. order == order2, "Inconsistent nside or ordering")
       mask = mask_tmp(:,1) /= 0
       deallocate(mask_tmp)
    end if

    ! Set up power spectrum
    call read_powspec_ascii(psfile, ps_tmp)
    lmax = ubound(ps_tmp,1)
    allocate(beam(0:lmax,1),ps(0:lmax))
    beam = 0
    call read_beam(beamfile, beam)
    pixwin => get_hpix_pixwin(nside)
    do i = 0, min(size(beam,1),size(pixwin,1))-1
       beam(i,1) = beam(i,1) * pixwin(i+1,1)
    end do
    ps = ps_tmp(0:lmax,1) * beam(0:lmax,1)**2
    deallocate(beam, ps_tmp)

    ! Build a correlation function from the power spectrum, and spline it
    ! 50k is overkill for now, but if the beam changes, this could change.
    call build_corrspline(ps, 50000, corr)
    if(myid == 0) then
       open(40,file="corrfun.txt")
       do i = 1, 50000
          m = (i-1)*pi/(50000-1)
          write(40,*) m*RAD2DEG, splint(corr, m)
       end do
    end if
    call dmem("Precomputed correlations")

    ! Set up kolmogorov correction, to account for our finite number
    ! of samples
    rad = rad*DEG2RAD
    avg_npix = nint(npix*(pi*rad**2)/(4*pi))
    ncorrect = 5000000
    call build_ks_correction(rng, avg_npix, ncorrect, kspline)
    if(myid == 0) then
       open(40,file="kscorr.txt")
       do i = 1, 10000
          m = (i-1)*1.0/(10000-1)
          write(40,*) m, splint(kspline, m)
       end do
    end if
    call dmem("Precomputed ks correction")

    ! Sample random pixels
    allocate(neighs(npix), probs(nmap,2),hists(nhist,nmap,2))
    hists = 0
    i     = 0
    do while(i < nsamp)
       j = modulo(zig_int(rng),npix)
       call pix2ang(nside, order, j, theta, phi)
       if(domask) then
          ! Make sure no pixel falls withing the mask
          call query_disc(nside, pix2vec(nside, order, j), rad, neighs, nneigh, nest=1)
          if(any(.not. mask(neighs(:nneigh)))) cycle
       else
          if(abs(pi/2-theta) < 30*DEG2RAD) cycle
       end if
       i = i+1
       if(modulo(i-1,nproc) /= myid) cycle
       call query_disc(nside, pix2vec(nside, order, j), rad, neighs, nneigh, nest=1)
       pix => neighs(:nneigh)
       ! Now go through all the simulations, getting a probability for each
!call plot_circle(maps(pix,1:1), corr, noise(pix), pix, nside, order, "disk.fits")
!stop
       call calc_ks_prob(maps(pix,:), corr, noise(pix),pix, nside, order, probs(:,1), rng)
       call calc_ks_prob_gurz(maps(pix,:), probs(:,2))
       ! Correct ks statistic
       do k = 1, nmap
          probs(k,1) = splint(kspline, probs(k,1))
          probs(k,2) = splint(kspline, probs(k,2))
       end do
       ! Build histogram
       do k = 1, nmap
          do l = 1, 2
             h = min(floor(probs(k,l)*nhist)+1,nhist)
             hists(h,k,l) = hists(h,k,l) + 1
          end do
       end do
       call dmem("Step " // trim(itoa(i)))
    end do
    deallocate(neighs, maps, noise, probs)
    call dmem("Analysed data")

    ! Collect all results
    allocate(thists(nhist,nmap,2))
    thists = hists
    call mpi_reduce(thists, hists, size(hists), mpi_integer, mpi_sum, 0, mpi_comm_world, err)
    deallocate(thists)
    call dmem("Reduce")

    if(myid == 0) then
       ! And output
       unit = getlun()
       open(unit,file="hist.txt")
       do i = 1, nhist
          write(unit,'(9e15.7)') (i-0.5)/nhist, real(hists(i,1,1),dp), &
           & mean(real(hists(i,2:,1),dp)), sqrt(variance(real(hists(i,2:,1),dp))), &
           & quantile(real(hists(i,2:,1),dp),[0.025d0,0.16d0,0.5d0,0.84d0,0.975d0])
       end do
       close(unit)
       open(unit,file="gist.txt")
       do i = 1, nhist
          write(unit,'(9e15.7)') (i-0.5)/nhist, real(hists(i,1,2),dp), &
           & mean(real(hists(i,2:,2),dp)), sqrt(variance(real(hists(i,2:,2),dp))), &
           & quantile(real(hists(i,2:,2),dp),[0.025d0,0.16d0,0.5d0,0.84d0,0.975d0])
       end do
       close(unit)

       ! Also output the full histograms
       open(unit,file="hist_full.txt")
       do i = 1, nhist
          write(unit,'(e15.7)',advance="no") (i-0.5)/nhist
          do j = 1, nmap
             write(unit,'(i8)',advance="no") hists(i,j,1)
          end do
          write(unit,*)
       end do
       close(unit)

       open(unit,file="gist_full.txt")
       do i = 1, nhist
          write(unit,'(e15.7)',advance="no") (i-0.5)/nhist
          do j = 1, nmap
             write(unit,'(i8)',advance="no") hists(i,j,2)
          end do
          write(unit,*)
       end do
       close(unit)
    end if

    write(stderr,*) myid, "done"
    deallocate(hists)
    call mpi_finalize(ierr)

  end subroutine

  subroutine arnold
    implicit none
    character(len=10000)           :: line
    character(len=64), allocatable :: str(:)
    integer(i4b)                   :: i, n
    real(dp)                       :: mval
    real(dp),          allocatable :: arr(:)
    call getarg(1, line); read(line,*) mval
    do
       read(*,'(a)') line
       n = num_tokens(line, ", ")
       allocate(str(n), arr(n))
       call get_tokens(line, ", ", str)
       do i = 1, n; read(str(i),*) arr(i); end do
       call ksp_uni(arr, mval)
       deallocate(str, arr)
    end do
  end subroutine

  ! For arnold
  subroutine ksp_uni(samples, mval)
    implicit none
    real(dp)     :: samples(:), res, dev, x(size(samples)), mval
    integer(i4b) :: j
    x    = samples
    call quicksort_real(x)
    dev = 0
write(*,'(4e15.7)') 0d0, 0d0, 0d0, 0d0
    do j = 1, size(x)
       dev = max(dev, abs(min(1d0,x(j)/mval)-real(j-1,dp)/size(x)))
       dev = max(dev, abs(min(1d0,x(j)/mval)-real(j,dp)/size(x)))
write(*,'(4e15.7)') x(j), real(j-1,dp)/size(x), min(1d0,x(j)/mval), dev
write(*,'(4e15.7)') x(j), real(j,dp)/size(x), min(1d0,x(j)/mval), dev
    end do
write(*,'(4e15.7)') 1d10, 1d0, 1d0, dev
    res = size(samples)**0.5*dev
write(*,*) res
write(*,*) ks_p(res)
  end subroutine

  subroutine calc_ks_prob(vals, corr, noise, pix, nside, order, prob, rng)
    implicit none
    real(dp)          :: vals(:,:), noise(:), prob(:)
    type(spline_type) :: corr
    integer(i4b)      :: pix(:), nside, order, n, k, i, j, err, nmap
    integer(i4b), allocatable :: piv(:)
    real(dp),     allocatable :: mat(:,:), pos(:,:), x(:,:), chol(:,:)
type(zig_rng) :: rng
    n    = size(pix)
    nmap = size(vals,2)
    allocate(mat(n,n),pos(3,n),x(n,nmap),chol(n,n), piv(n))
    ! Build equation system
    do k = 1, n
       pos(:,k) = pix2vec(nside, order, pix(k))
    end do
    call make_cmbmat(corr, pos, mat)
    do j = 1, size(mat,1)
       mat(j,j) = mat(j,j) + noise(j)**2 ! Add the white noise
    end do
    call cholesky_decompose(mat, chol)

    ! And solve
    err = 0
    x   = vals
    call dgesv(n, nmap, chol, n, piv, x, n, err)
    call assert(err == 0, "Error in solving equation system in calc_ks_prob")

    ! And extract the probability of each
    do k = 1, nmap
       prob(k) = ks_gauss(x(:,k), n)
    end do
    deallocate(mat,pos,x,chol,piv)
  end subroutine

  subroutine calc_ks_prob_gurz(vals, prob)
    implicit none
    real(dp)          :: vals(:,:), prob(:)
    integer(i4b)      :: i
    do i = 1, size(vals,2)
       prob(i) = ks_gauss((vals(:,i)-mean(vals(:,i)))/sqrt(variance(vals(:,i))), size(vals,1))
    end do
  end subroutine

  function calc_ks_gurz(vals) result(pte)
    implicit none
    real(dp) :: vals(:), pte
    pte = ks_gauss((vals-mean(vals))/sqrt(variance(vals)), size(vals))
  end function

  subroutine make_cmbmat(corr, pos, mat)
    implicit none
    type(spline_type) :: corr
    real(dp)     :: mat(:,:), pos(:,:)
    integer(i4b) :: i, j, n, nside, order
    real(dp), allocatable :: x(:)
    n = size(pos,2)
    allocate(x(size(mat)))
    do i = 1, n
       x((i-1)*n+i) = 0
       do j = i+1, n
          call angdist(pos(:,i), pos(:,j), x((i-1)*n+j))
          x((j-1)*n+i) = x((i-1)*n+j)
       end do
    end do
    do i = 1, size(x)
       x(i) = splint(corr, x(i))
    end do
    mat = reshape(x,[size(mat,1),size(mat,2)])
    deallocate(x)
  end subroutine

  subroutine build_corrspline(ps, n, corr)
    implicit none
    real(dp)          :: ps(:)
    integer(i4b)      :: n, i
    type(spline_type) :: corr
    real(dp), allocatable :: x(:), y(:)
    allocate(x(n),y(n))
    x = (irange(n)-1)*pi/(n-1)
    ! Beware: If we have some extra white noise, we will have
    ! a sharp spike at entry 0, which will be hard to spline.
    call cl2corr(ps, y, x)
    call spline(corr, x, y, boundary=[0d0,0d0], regular=.true.)
    deallocate(x,y)
  end subroutine


  ! Make a cumulative distribution function of the results of
  ! the KS test on data that actually follows the distribution.
  ! Then spline this.
  subroutine build_ks_correction(rng, nsamp, nreal, kspline)
    implicit none
    type(zig_rng)     :: rng
    integer(i4b)      :: nsamp, nreal, i, j, myid, nproc, ierr
    type(spline_type) :: kspline
    real(dp), allocatable :: x(:), y(:), samps(:), ps(:), ps_tmp(:)
    allocate(samps(nsamp),ps(nreal),ps_tmp(nreal))
    call mpi_comm_rank(MPI_COMM_WORLD, myid,  ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    ps_tmp = 0
    do i = 1+myid, nreal, nproc
       do j = 1, nsamp
          samps(j) = zig_gauss(rng)
       end do
       ps_tmp(i) = ks_gauss(samps, nsamp)
    end do
    call mpi_allreduce(ps_tmp, ps, size(ps), mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    call quicksort_real(ps)
    allocate(x(nreal+2),y(nreal+2))
    x(1)         = 0;  y(1)       = 0
    x(nreal+2)   = 1;  y(nreal+2) = 1
    x(2:nreal+1) = ps
    y(2:nreal+1) = real(irange(nreal),dp)/(nreal+1)
    call spline(kspline, x, y, linear=.true.)
    deallocate(x,y,samps,ps,ps_tmp)
  end subroutine

  function ks_gauss(samples, neff) result(res)
    implicit none
    real(dp)     :: samples(:), res, dev, x(size(samples))
    integer(i4b) :: j, neff
    x = samples
    call quicksort_real(x)
    dev = 0
    do j = 1, size(x)
       dev = max(dev, abs(pnorm(x(j),0d0,1d0)-real(j-1,dp)/size(x)))
       dev = max(dev, abs(pnorm(x(j),0d0,1d0)-real(j,dp)/size(x)))
    end do
    res = ks_prob(dev, neff)
  end function

  function pnorm(x,m,v)
    implicit none
    real(dp) :: x, pnorm, m, v
    pnorm = 0.5d0*(1+erf((x-m)/sqrt(2d0*v)))
  end function

  function ks_prob(maxdev, nobs) result(pte)
    implicit none
    real(dp)     :: maxdev, pte, nroot
    integer(i4b) :: nobs
    nroot = nobs**0.5d0
    pte = ks_p(nroot*maxdev)
    !pte = ks_p((nroot + 0.12 + 0.11/nroot)*maxdev)
  end function

  function ks_p(x) result(p)
    implicit none
    real(dp) :: x, p
    real(dp), parameter :: A = pi**2/8, B = sqrt(2*pi)
    if(x < 0.07302) then
       p = 0
    else if(x < 1.081) then
       p = B/(x*exp(A/x**2))
    else
       p = 1-2*exp(-2*x**2)
    end if
  end function

  function ks_q(x) result(q)
    implicit none
    real(dp) :: x, q
    real(dp), parameter :: A = pi**2/8, B = sqrt(2*pi)
    if(x < 1.081) then
       q = 1 - ks_p(x)
    else
       q = 2*exp(-2*x**2)
    end if
  end function

  function ks_p2(x) result(p)
    implicit none
    real(dp) :: x, p
    real(dp), parameter :: A = pi**2/8, B = sqrt(2*pi)
    integer(i4b) :: i
    if(x < 0.07302) then
       p = 0
    else if(x < 1.081) then
       p = 0
       do i = 1, 10
          p = p + exp(-(2*i-1)**2*A/x**2)
       end do
       p = p*B/x
    else
       p = 0
       do i = 1, 10
          p = p + (-1)**(i-1)*exp(-2*i**2*x**2)
       end do
       p = 1-2*p
    end if
  end function

  function ks_q2(x) result(q)
    implicit none
    real(dp) :: x, q
    real(dp), parameter :: A = pi**2/8, B = sqrt(2*pi)
    integer(i4b) :: i
    if(x < 1.081) then
       q = 1 - ks_p2(x)
    else
       q = 0
       do i = 1, 10
          q = q + (-1)**(i-1)*exp(-2*i**2*x**2)
       end do
       q = 2*q
    end if
  end function

  subroutine kstest
    implicit none
    character(len=512) :: arg
    integer(i4b)       :: i, j, n, ns
    type(zig_rng)      :: rng
    type(spline_type)  :: kspline
    real(dp)           :: y
    real(dp), allocatable :: x(:)
    call getarg(1,arg); read(arg,*) n
    call getarg(2,arg); read(arg,*) ns
    call zig_init(rng, 11)
    call build_ks_correction(rng, n, ns, kspline)
    allocate(x(n))
    do j = 1, ns
       do i = 1, n
          x(i) = zig_gauss(rng)
       end do
       y = splint(kspline, ks_gauss(x,n))
       write(*,*) y
    end do
  end subroutine





  ! Try to repeat gurzadyan's analysis
  subroutine gurzadyan
    implicit none
    character(len=512)                 :: arg, mapfile, psfile, beamfile, noise_mapfile
    integer(i4b)                       :: i, j, k, nsamp, nneigh, order, nside, npix, onpix
    integer(i4b)                       :: order2, nside2, nstep, ierr, myid, nproc
    integer(i4b)                       :: avg_npix, donoise
    type(zig_rng)                      :: rng
    type(spline_type)                  :: corr, kspline
    real(dp)                           :: dev, m, v, rad, noise, pte, gurz, theta, phi
    real(dp),     allocatable          :: map(:,:), nmap(:,:), ps(:,:), beam(:,:)
    real(dp),     allocatable          :: ptemap(:,:), gurzmap(:,:), tmpmap(:,:)
    real(dp),     pointer              :: pixwin(:,:)
    integer(i4b), allocatable, target  :: neighs(:)
    integer(i4b),              pointer :: pix(:)
    call getarg(1, mapfile)
    call getarg(2, psfile)
    call getarg(3, beamfile)
    call getarg(4, noise_mapfile)
    call getarg(5, arg); read(arg,*) rad
    call getarg(6, arg); read(arg,*) nsamp
    call getarg(7, arg); read(arg,*) donoise
    call zig_init(rng, 186982652)
    call read_sparse_as_full(map,  nside,  order,  mapfile)
    call read_sparse_as_full(nmap, nside2, order2, noise_mapfile)
    call assert(nside==nside2 .and. order == order2, "Inconsistent nside or ordering")
    call read_powspec_ascii(psfile, ps)
    allocate(beam(0:ubound(ps,1),1))
    beam = 0
    call read_beam(beamfile, beam)
    ! Apply pixel window. Should make sure simulation also has pixel window.
    pixwin => get_hpix_pixwin(nside)
    do i = 0, min(size(beam,1),size(pixwin,1))-1
       beam(i,1) = beam(i,1) * pixwin(i+1,1)
    end do

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
!call zig_init(rng,1+myid)

    npix = 12*nside**2

    ! Add some noise to map (skip this if we have a real map)
    if(donoise > 0) then
       do i = 0, size(map,1)-1
          map(i,1) = map(i,1) + zig_gauss(rng)*nmap(i,1)
       end do
    end if
    !nmap = nmap*1.07
    ! Apply beam to power spectrum. This should match the beam of the map
    ps(:,1) = ps(:,1) * beam(:,1)**2
    ! Build a correlation function from the power spectrum, and spline it
    ! 50k is overkill for now, but if the beam changes, this could change.
    call build_corrspline(ps(:,1), 50000, corr)

    !! Sample random pixels
    !rad = rad*DEG2RAD
    !allocate(neighs(npix))
    !! Set up kolmogorov correction, to account for our finite number
    !! of samples
    !avg_npix = nint(npix*(pi*rad**2)/(4*pi))
    !call build_ks_correction(rng, avg_npix, 10*nsamp, kspline)
    !i = 0
    !do while(i < nsamp)
    !   j = modulo(zig_int(rng),npix)
    !   call pix2ang(nside, order, j, theta, phi)
    !   if(abs(pi/2-theta) < 30*DEG2RAD) cycle
    !   i = i+1
    !   if(modulo(i-1,nproc) /= myid) cycle
    !   call query_disc(nside, pix2vec(nside, order, j), rad, neighs, nneigh, nest=order-1)
    !   pix => neighs(:nneigh)
    !   pte  = splint(kspline, calc_ks_pte (map(pix,1), corr, nmap(pix,1), pix, nside, order))
    !   gurz = splint(kspline, calc_ks_gurz(map(pix,1)))
    !   write(*,'(i8,2f10.7)') i, pte, gurz
    !end do
    !deallocate(neighs)

    ! Make pte-maps. Not sure why this should be useful. They will look
    ! like noise if everything is OK. Better with something with more contrast
    nstep = nint(rad)
    onpix = npix/4**nstep
    allocate(ptemap(0:onpix-1,1), gurzmap(0:onpix-1,1), pix(4**nstep))
    ! must be in nest
    call set_ordering(nest, order, map)
    call set_ordering(nest, order, nmap)
    ptemap  = 0
    gurzmap = 0
    ! Set up kolmogorov correction, to account for our finite number
    ! of samples
    call build_ks_correction(rng, 4**nstep, 10*onpix, kspline)
    do i = myid, onpix-1, nproc
       !call pix2ang(nside/2**nstep, nest, i, theta, phi)
       !if(abs(pi/2-theta) < pi/6) cycle
       pix  = irange(i*4**nstep,(i+1)*4**nstep-1)
       pte  = splint(kspline, calc_ks_pte(map(pix,1),corr,nmap(pix,1),pix,nside,nest))
       gurz = splint(kspline, calc_ks_gurz(map(pix,1)))
       ptemap (i,1) = pte  !qnorm(pte)
       gurzmap(i,1) = gurz !qnorm(gurz)
       write(*,'(i8,2f16.12,2e15.7)') i, pte, gurz, ptemap(i,1), gurzmap(i,1)
    end do
    allocate(tmpmap(0:onpix-1,1))
    call mpi_reduce(ptemap,  tmpmap, size(ptemap), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    ptemap = tmpmap
    call mpi_reduce(gurzmap, tmpmap, size(ptemap), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
    gurzmap = tmpmap
    deallocate(tmpmap)
    if(myid == 0) then
       call write_map(ptemap,  nest, "ptemap.fits")
       call write_map(gurzmap, nest, "gurzmap.fits")
    end if
    deallocate(ptemap, gurzmap, pix)

    write(stderr,*) myid, "done"
    call mpi_finalize(ierr)
  end subroutine

  function calc_ks_pte(vals, corr, noise, pix, nside, order) result(pte)
    implicit none
    real(dp)          :: vals(:), noise(:), pte
    type(spline_type) :: corr
    integer(i4b)      :: pix(:), nside, order, n, k, j, neff
    real(dp), allocatable :: mat(:,:), pos(:,:), x(:), ichol(:,:)
    n = size(pix)
    allocate(mat(n,n),pos(3,n),x(n),ichol(n,n))
    do k = 1, n
       pos(:,k) = pix2vec(nside, order, pix(k))
    end do
    call make_cmbmat(corr, pos, mat)
    do j = 1, size(mat,1)
       mat(j,j) = mat(j,j) + noise(j)**2 ! Add the white noise
    end do
    call cholesky_decompose(mat, ichol)
    call solve_linear_system(ichol, x, vals)
    neff=n
    !call eigen_pow(mat, -0.5d0, ichol, neff)
    !x = matmul(ichol, vals)
!call eigen_pow(mat,0.5d0, ichol)
!vals = matmul(ichol, x)
    pte = ks_gauss(x, neff)
    deallocate(mat,pos,x,ichol)
  end function

  subroutine plot_circle(vals, corr, noise, pix, nside, order, oname)
    implicit none
    character(len=*)  :: oname
    real(dp)          :: vals(:,:), noise(:)
    type(spline_type) :: corr
    integer(i4b)      :: pix(:), nside, order, n, k, i, j, err, nmap
    integer(i4b), allocatable :: piv(:)
    real(dp),     allocatable :: mat(:,:), pos(:,:), x(:,:), chol(:,:), omap(:,:)
    n    = size(pix)
    nmap = size(vals,2)
    allocate(mat(n,n),pos(3,n),x(n,nmap),chol(n,n), piv(n))
    ! Build equation system
    do k = 1, n
       pos(:,k) = pix2vec(nside, order, pix(k))
    end do
    call make_cmbmat(corr, pos, mat)
    do j = 1, size(mat,1)
       mat(j,j) = mat(j,j) + noise(j)**2 ! Add the white noise
    end do
    call cholesky_decompose(mat, chol)

    ! And solve
    err = 0
    x   = vals
    call dgesv(n, nmap, chol, n, piv, x, n, err)
    call assert(err == 0, "Error in solving equation system in calc_ks_prob")

    ! Plot the circle before and after
    allocate(omap(n,2))
    omap(:,1) = vals(:,1)
    omap(:,2) = x(:,1)
    call write_sparse_map(omap, pix, nside, order, oname)

    deallocate(mat,pos,x,chol,piv,omap)
  end subroutine

end program
