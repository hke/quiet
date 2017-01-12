!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a module to sindex                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module sindex_maxlike_mod
  use quiet_utils
  use quiet_fileutils
  use rngmod
  implicit none

!  public :: setup2_maxlike_mod, setup_maxlike_mod, terminate_maxlike_mod, calculate_A, calculate_cov_params, func, dfunc, ddfunc

  type banddata
     real(dp),     allocatable, dimension(:,:,:)     :: map, Pmap
     real(dp),     allocatable, dimension(:,:,:)     :: cov, Pcov
     real(dp),     allocatable, dimension(:,:)       :: P
     real(dp),     allocatable, dimension(:)         :: freq
     integer(i4b), allocatable, dimension(:)         :: pixels, map2mask
     integer(i4b)    :: n, npix, nside, nmaps, ordering, pol, num_bands
  end type banddata

!  private
  type(banddata), allocatable, dimension(:)  :: d
  real(dp), allocatable, dimension(:)        :: a2t 
  real(dp), allocatable, dimension(:)        :: freq
  integer(i4b)                               :: nfreq, k, subreg, pol, mnum
  real(dp)                                   :: healnan=-1.6375d30
  logical(lgt)                               :: chatty=.false.
contains



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setup_maxlike_mod(in_data, in_k)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type(banddata), dimension(1:), intent(in) :: in_data
    integer(i4b), intent(in)                 :: in_k
    integer(i4b)                             :: v, i, npix, j,q, ntot, u,w,n2, i1,i2, n_new, teller
    !
    real(dp), allocatable,dimension(:)       :: eigen
    real(dp), allocatable,dimension(:,:)     :: Vn
    type(planck_rng)     :: rng_handle



    subreg = size(in_data)
    nfreq=size(in_data(1)%freq)
    allocate(freq(nfreq))
    freq = in_data(1)%freq
    pol=in_data(1)%pol
    k = in_k
    mnum=1
    if (k==4) mnum=2

!
!! the component (k) covariance matrix is read in and inverted. 
!  NB in this routine it is only the inverted covariance matrix which is stored.
!
    allocate(d(1:subreg))
!    data = in_data
    do i=1,subreg
       ntot = in_data(i)%n
       if (k==4) ntot = in_data(i)%n*pol

       allocate(d(i)%cov(ntot,ntot,nfreq))
       allocate(d(i)%map(ntot,1,nfreq))
!       allocate(eigen(ntot), Vn(ntot,ntot))
       d(i)%n = in_data(i)%n

       do v=1,nfreq
          if (k==1) then
             d(i)%cov(:,:,v) = in_data(i)%cov(:,:,v)
             d(i)%map(:,1,v) = in_data(i)%map(:,k,v)
          else if (k==4) then
!!$             n2 = d(i)%n
!!$             d(i)%cov(:,:,v) =0.d0
!!$             d(i)%cov(1:n2,1:n2,v) = in_data(i)%cov(1:n2,1:n2,v)        !QQ
!!$             d(i)%cov(n2+1:ntot,n2+1:ntot,v) = in_data(i)%cov(n2+1:ntot,n2+1:ntot,v)   !UU
!!$             do j=1, n2
!!$                d(i)%cov(n2+j,j,v) = 0.0d0   !QU
!!$                d(i)%cov(j,n2+j,v) = 0.0d0
!!$             end do
!             do i1=1,ntot
!                do i2=1,ntot
!                   d(i)%cov(i1,i2,v) = in_data(i)%cov(i1,i2,v)
!                end do
!             end do
             d(i)%cov(:,:,v) = in_data(i)%cov(:,:,v)                  !get the whole cov matrix
             d(i)%map(1:in_data(i)%n,1,v) =in_data(i)%map(:,2,v)      !stack the maps
             d(i)%map(in_data(i)%n+1:ntot,1,v) =in_data(i)%map(:,3,v)
          else if(k==2) then
             d(i)%cov(:,:,v) = in_data(i)%cov(1:ntot,1:ntot,v)        !only the QQ part of matrix
             d(i)%map(:,1,v) = in_data(i)%map(:,k,v)
          else if(k==3) then
             d(i)%cov(:,:,v) = in_data(i)%cov(ntot+1:2*ntot,ntot+1:2*ntot,v)   ! only the UU part
             d(i)%map(:,1,v) = in_data(i)%map(:,k,v)
          else  ! this is the alpha-case, alpha=0..90
             d(i)%cov(:,:,v) = in_data(i)%cov(1:ntot,1:ntot,v)*(cos(2.d0*k*pi/180.))**2 + &
                             & in_data(i)%cov(ntot+1:2*ntot,ntot+1:2*ntot,v)*(sin(2.d0*k*pi/180.))**2 +&
                             & 2* in_data(i)%cov(1:ntot,ntot+1:2*ntot,v)*sin(2.d0*k*pi/180.)*cos(2.d0*k*pi/180.)
!!!             if (v==1) then
!!!                d(i)%map(:,1,v) = in_data(i)%map(:,2,v)*cos(2.d0*(k+3.d0)*pi/180.) + in_data(i)%map(:,3,v)*sin(2.d0*(k+3.d0)*pi/180.)
!!!             else
                d(i)%map(:,1,v) = in_data(i)%map(:,2,v)*cos(2.d0*k*pi/180.) + in_data(i)%map(:,3,v)*sin(2.d0*k*pi/180.)
!!!             end if
          end if

!!$          !This is with change of basisGet the eigenmodes wanted
!!$          if (v==1) then
!!$             call get_eigen_decomposition(34, d(i)%cov(:,:,1), eigen, Vn)
!!$             teller=0
!!$             do j=1,ntot
!!$!!!                if (eigen(j)< 0.d0) Vn(:,j)=0.d0
!!$                if (eigen(j)< 1.d-4) then
!!$                   teller=teller+1
!!$                else
!!$                   exit
!!$                end if
!!$             end do
!!$             write(*,*) "antall moder under terskel = ",teller
!!$             n_new=ntot-teller
!!$             allocate(d(i)%P(ntot,n_new))
!!$             d(i)%P = Vn(:,teller+1:ntot)
!!$             allocate(d(i)%Pcov(n_new, n_new,nfreq), d(i)%Pmap(n_new,1,nfreq))
!!$          end if
!!$          d(i)%Pcov(:,:,v) = matmul(transpose(d(i)%P), matmul(d(i)%cov(:,:,v),d(i)%P) )
!!$          d(i)%Pmap(:,1,v) = matmul(transpose(d(i)%P),d(i)%map(:,1,v) )
!!$!!!          d(i)%cov(:,:,v) = matmul(transpose(Vn), matmul(d(i)%cov(:,:,v),Vn) )
!!$!!!          d(i)%map(:,1,v) = matmul(transpose(Vn),d(i)%map(:,1,v) )
!!$!!!          call invert_singular_matrix(d(i)%cov(:,:,v),1.d-12)
!!$          call invert_singular_matrix(d(i)%Pcov(:,:,v),1.d-12)
!!$       end do !end v
!!$       d(i)%npix = d(i)%n
!!$       d(i)%n = n_new
!!$       if (k==4) d(i)%n = n_new / 2  !should check that this is ok!
!!$       deallocate(eigen, Vn)
          !!! end this is with change of basis

          !!! this is without the change of basis
          d(i)%npix = d(i)%n
          if (v==1) allocate(d(i)%P(ntot,ntot),d(i)%Pcov(ntot, ntot,nfreq), d(i)%Pmap(ntot,1,nfreq))
          d(i)%Pcov(:,:,v) = d(i)%cov(:,:,v)
          d(i)%Pmap(:,:,v) = d(i)%map(:,:,v)
          d(i)%P=0.d0
          do j=1,ntot
             d(i)%P(j,j)=1.d0
          end do
          call invert_singular_matrix(d(i)%Pcov(:,:,v),1.d-20)  !-12 or -20
       end do !end v
       !!! end this is without the change of basis

       deallocate(d(i)%map, d(i)%cov)
    end do !end i subreg

    if(.not.allocated(a2t)) then
       allocate(a2t(nfreq))
       do v = 1, nfreq
          a2t(v)=ant2thermo(freq(v))
       end do
    end if    

    if (chatty) then
       i=1
       do v=1,2
          write(*,*)"covdiag",v, d(i)%cov(1,1,v), d(i)%cov(2,2,v), d(i)%cov(3,3,v), d(i)%cov(4,4,v),d(i)%cov(5,5,v), d(i)%cov(6,6,v), d(i)%cov(7,7,v)
          write(*,*)"map",v, d(i)%map(1,1,v), d(i)%map(2,1,v), d(i)%map(3,1,v), d(i)%map(4,1,v),d(i)%map(5,1,v), d(i)%map(6,1,v), d(i)%map(7,1,v)
       end do
    end if

  end subroutine setup_maxlike_mod


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine terminate_maxlike_mod
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer(i4b)               :: i

    do i = 1, size(d)
       call free_banddata(d(i))
    end do
    deallocate(a2t,freq,d)
  end subroutine terminate_maxlike_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine free_banddata(data2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    type(banddata) :: data2
    if(allocated(data2%map))   deallocate(data2%map)
    if(allocated(data2%Pmap))   deallocate(data2%Pmap)
    if(allocated(data2%cov))   deallocate(data2%cov)
    if(allocated(data2%Pcov))   deallocate(data2%Pcov)
    if(allocated(data2%P))   deallocate(data2%P)
    if(allocated(data2%pixels))   deallocate(data2%pixels)
    if(allocated(data2%map2mask)) deallocate(data2%map2mask)
    if(allocated(data2%freq))     deallocate(data2%freq)

  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_Am(p,A,m, reg)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    REAL(dp),        dimension(1:), INTENT(IN)    :: p
    real(dp),        dimension(1:), intent(out)   :: A
    real(dp),        dimension(1:), intent(out)   :: m
    integer(i4b), intent(in)                     :: reg
    real(dp),        allocatable, dimension(:)   :: W1, W2, W, x
    real(dp),        allocatable, dimension(:,:) :: Penvector,envector, S1, Q2, S2, Q1, Q
    integer(i4b)                                 :: v, num,n, i
    real(dp)                                     :: beta

    if (chatty) write(*,*) "calculating Am"
    
    m  = healnan
    beta = p(1)
    A = 0.d0

    i=reg
    n = size(d(i)%Pmap,1)

    allocate(Q1(1:n,1:n))
    allocate(W1(1:n))
    Q1 = 0.d0
    W1 = 0.d0
    do v = 1, nfreq
       Q1 = Q1 + ( (a2t(v)**2)*((freq(v)/freq(1))**(2.d0*beta)) * d(i)%Pcov(:,:,v) )
       W1 = W1 + ( (a2t(v)*(freq(v)/freq(1))**beta) * matmul( d(i)%Pcov(:,:,v), d(i)%Pmap(:,1,v) ) )
    end do
    v=2  
    allocate(envector(1:d(i)%npix*mnum,mnum))
    allocate(Penvector(1:n,mnum))
    allocate(Q2(1:n,mnum))
    allocate(S1(1:n,mnum))
    allocate(S2(mnum,mnum))
    allocate(W2(mnum))

    if (k==4) then
       envector = 0.d0
       envector(1:n/2,1) =1.d0
       envector(n/2+1:n,2) =1.d0
    else
       envector(:,1) = 1.d0
    end if
    Penvector = matmul(transpose(d(i)%P),envector)   !!!P!!!
    deallocate(envector)

    Q2  = matmul(a2t(v) * ((freq(v)/freq(1))**beta)* d(i)%Pcov(:,:,v) , Penvector)
    S2  = matmul(transpose(Penvector), matmul( d(i)%Pcov(:,:,v), Penvector) )
    W2  = matmul(transpose(Penvector), matmul( d(i)%Pcov(:,:,v), d(i)%Pmap(:,1,v) ) )
    S1 = Q2    

    if (chatty) write(*,*) "n=", n, "mnum=",mnum
    !put the matrix elements into the big matrix
    allocate(Q(1:n+mnum,1:n+mnum))
    allocate(W(1:n+mnum))
    Q(1  :n ,1  :n ) = Q1(1  :n ,1  :n )
    Q(1  :n ,n+1:n+mnum )   = S1
    Q(n+1:n+mnum , 1:n ) = transpose(Q2)
    Q(n+1:n+mnum ,n+1:n+mnum ) = S2
    W( 1 : n ) = W1
    W(n+1:n+mnum)  = W2

!!$    !write out stuff
!!$    write(*,*) "Q1"
!!$    call dump_matrix(Q1,fmt="(e16.8)")
!!$    write(*,*) "S1"
!!$    call dump_matrix(S1,fmt="(e16.8)")
!!$    write(*,*) "W1"
!!$    call dump_matrix(W1,fmt="(e16.8)")
!!$    write(*,*) "Q2"
!!$    call dump_matrix(Q2,fmt="(e16.8)")
!!$    write(*,*) "S2"
!!$    call dump_matrix(S2,fmt="(e16.8)")
!!$    write(*,*) "W2"
!!$    call dump_matrix(W2,fmt="(e16.8)")
!!$    write(*,*) "Q"
!!$    call dump_matrix(Q,fmt="(e16.8)")
!!$    write(*,*) "W"
!!$    call dump_matrix(W,fmt="(e16.8)")
    !end write out

    deallocate(Q1,Q2)
    deallocate(S1,S2)
    deallocate(W1,W2)

    !solve the system of linear equations
    allocate(x(1:n+mnum))
    call solve_linear_system(Q,x,W)

    A(1:n) = x(1:n)
    m(1) = x(n+1)
    if (k==4) m(2) = x(n+mnum)

!!$    write(*,*) "x"
!!$    call dump_matrix(x,fmt="(e16.8)")

    deallocate(x,Q,W)
    deallocate(Penvector)


  end subroutine calculate_Am

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_Am0(p,A,m, reg)  !m=0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !routine where m=0, it just calculates the amplitude, A
    implicit none
    REAL(dp),        dimension(1:), INTENT(IN)    :: p
    real(dp),        dimension(1:), intent(out)   :: A
    real(dp),        dimension(1:), intent(out)   :: m
    integer(i4b), intent(in)                     :: reg
    real(dp),        allocatable, dimension(:)   :: W1, W2, W, x
    real(dp),        allocatable, dimension(:,:) :: S1, Q2, S2, Q1, Q
    integer(i4b)                                 :: v, num,n, i
    real(dp)                                     :: beta

    if (chatty) write(*,*) "calculating Am"
    
    m  = 0.d0
    beta = p(1)
    A = 0.d0

    i=reg
    n = size(d(i)%Pmap,1)

    allocate(Q1(1:n,1:n))
    allocate(W1(1:n))
    Q1 = 0.d0
    W1 = 0.d0
    do v = 1, nfreq
       Q1 = Q1 + ( (a2t(v)**2)*((freq(v)/freq(1))**(2.d0*beta)) * d(i)%Pcov(:,:,v) )
       W1 = W1 + ( (a2t(v)*(freq(v)/freq(1))**beta) * matmul( d(i)%Pcov(:,:,v), d(i)%Pmap(:,1,v) ) )
    end do

    if (chatty) write(*,*) "n=", n, "mnum=",mnum
    !put the matrix elements into the big matrix

    !solve the system of linear equations
    allocate(x(1:n))
    call solve_linear_system(Q1,x,W1)
    A(1:n) = x(1:n)

    deallocate(x,Q1,W1)


  end subroutine calculate_Am0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_Am3(p,A,m, reg)  !m=0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !routine where m is a fixed value (set equal the invalue of m), it just calculates the amplitude, A
    ! NB only run for k=2 or 3  (not 4!!)
    implicit none
    REAL(dp),        dimension(1:), INTENT(IN)    :: p
    real(dp),        dimension(1:), intent(out)   :: A
    real(dp),        dimension(1:,1:), intent(inout) :: m
    integer(i4b), intent(in)                     :: reg
    real(dp),        allocatable, dimension(:)   :: W1, W2, W, x
    real(dp),        allocatable, dimension(:,:) :: S1, Q2, S2, Q1, Q
    integer(i4b)                                 :: v, num,n, i
    real(dp)                                     :: beta

    if (chatty) write(*,*) "calculating Am"
    
    beta = p(1)
    A = 0.d0

    i=reg
    n = size(d(i)%Pmap,1)

    allocate(Q1(1:n,1:n))
    allocate(W1(1:n))
    Q1 = 0.d0
    W1 = 0.d0
    do v = 1, nfreq
       Q1 = Q1 + ( (a2t(v)**2)*((freq(v)/freq(1))**(2.d0*beta)) * d(i)%Pcov(:,:,v) )
       W1 = W1 + ( (a2t(v)*(freq(v)/freq(1))**beta) * matmul( d(i)%Pcov(:,:,v), ( d(i)%Pmap(:,1,v) -m(1,v) ) )  )
    end do

    if (chatty) write(*,*) "n=", n, "mnum=",mnum
    !put the matrix elements into the big matrix

    !solve the system of linear equations
    allocate(x(1:n))
    call solve_linear_system(Q1,x,W1)
    A(1:n) = x(1:n)

    deallocate(x,Q1,W1)


  end subroutine calculate_Am3



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gridmin(p, uncert, directory)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    real(dp), dimension(1:), intent(out)   :: p
    real(dp), intent(out)                  :: uncert
    character(len=256) , intent(in)        :: directory
    character(len=256)                     :: outfile
    character(len=1)                       :: k_char
    real(dp), allocatable, dimension(:,:)  :: chisq, beta, grid, grid2
    real(dp), allocatable, dimension(:)    :: dbeta, temp
    integer(i4b)                           :: ngrid, i, pos(1), l, j, loops, np
    real(dp)                               :: betamin, betamax, finebeta, finebeta2, mval
    real(dp)                               :: slett, norm_factor

    ngrid = 30
    loops = 2
    allocate(chisq(ngrid,loops))
    allocate(beta(ngrid,loops))
    betamin = -5.d0! -3.7d0
    betamax = -1.0d0! -1.5d0


    do l = 1, loops

       if (l==1) then
       else if (l==3) then
          mval = minval(chisq(:,l-1))
          do j = pos(1),ngrid      ! loop over beta to find where grid = 3sigma=11.8  2sigma=6.17 1sigma=2.3
             if ( chisq(j,l-1) >= 3. + mval ) exit
          end do
          finebeta = (betamax-betamin)/(ngrid-1.d0)*(j - pos(1))
          betamax = beta(pos(1),l-1) + finebeta
          betamin = beta(pos(1),l-1) - finebeta
       else
          finebeta =(betamax-betamin)/5.
          betamax = beta(pos(1),l-1) + finebeta
          betamin = beta(pos(1),l-1) - finebeta
       end if

       ! write out file with gnuplot-format
       if(l==loops) then
          call int2string(k,k_char)
          outfile = trim(directory) // 'liks_' // k_char // '.txt'
          open(42, file=trim(outfile))
!          write(42,*) ngrid
       end if
       do i = 1,ngrid
          beta(i,l) = (betamax-betamin)/(ngrid-1.d0)*(i-1) + betamin
          chisq(i,l) = func( [beta(i,l)] )
          if(l==loops )write(42,*) beta(i,l),chisq(i,l)
       end do
       if(l==loops) close(42)
       pos=minloc(chisq(:,l))
       if (pos(1)==1 .or. pos(1)==ngrid) then
          p = healnan
          write(*,*) "NB: gridmin hit the border of the grid"
          cycle  !exit loop
       end if
    end do

    !reorganize matrix
write(*,*)"the beta at chisq max:",beta(pos(1),l-1)
    np = ngrid*loops
    allocate(dbeta(np))
    allocate(grid(np,2), grid2(np,2))
    grid(:,1)=reshape( beta (:,:),[np])
    grid(:,2)=reshape( chisq(:,:),[np])

    !sort the matrix according to value of beta (grid(:,1))
    do i=1,np
       pos = minloc(grid(:,1))
       grid2(i,:)=grid(pos(1),:)
       grid(pos(1),1)=99.
    end do
    grid=grid2
    deallocate(grid2)

    !find the step size, dbeta
    dbeta(1)= (grid(2,1)-grid(1,1) )/2.d0
    do i=2,np-1
       dbeta(i)=( grid(i+1,1)-grid(i-1,1) )/2.d0
    end do
    dbeta(np)=(grid(np,1)-grid(np-1,1) )/2.d0

    !evaluate likelihood from chisquare and normalize
    grid(:,2)=grid(:,2)-minval(grid(:,2))
    grid(:,2)=exp( - grid(:,2)/2.d0 )
    norm_factor=sum(grid(:,2)*dbeta)
    grid(:,2)=grid(:,2)/norm_factor

    !calculate spectral index and its standard deviation
    p(1) = sum(grid(:,1)*grid(:,2)* dbeta  )
    uncert = sum( ((grid(:,1)-p(1))**2) *grid(:,2) *dbeta  )
    uncert = sqrt(uncert)

!write(*,*)"AA"
!call dump_matrix(grid)
    write(*,*) "spectral index from grid =", p(1), "+-",uncert

    deallocate(chisq, beta, dbeta, grid)

  end subroutine gridmin

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gridmin_bm(p, uncert, directory)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    real(dp), dimension(1:), intent(out)   :: p
    real(dp), intent(out)                  :: uncert
    character(len=256) , intent(in)        :: directory
    character(len=256)                     :: outfile
    character(len=1)                       :: k_char
    real(dp), allocatable, dimension(:,:,:):: chisq
    real(dp), allocatable, dimension(:,:)  :: beta, off, grid, grid2, m, lik
    real(dp), allocatable, dimension(:)    :: dbeta, doff, temp, vec, Pb, Pm
    integer(i4b)                           :: ngrid, i, pos(2), l, j, loops
    real(dp)                               :: betamin, betamax, finebeta, finebeta2, mval
    real(dp)                               :: slett, norm_factor, mmin, mmax, uncert_b, uncert_m


    ngrid = 30
    loops = 1
    allocate(chisq(ngrid,ngrid,loops))
    allocate(beta(ngrid,loops))
    allocate(off(ngrid,loops))
    allocate(m(1,2))
    betamin = -6.0d0
    betamax = 0.0d0
    mmin = -10.
    mmax = 10.

    do l = 1, loops

       if (l==1) then
       else
          finebeta =(betamax-betamin)/5.
          betamax = beta(pos(1),l-1) + finebeta
          betamin = beta(pos(1),l-1) - finebeta
       end if

       ! write out file with gnuplot-format
       if(l==1) then
          call int2string(k,k_char)
          outfile = trim(directory) // 'liks_' // k_char // '.txt'
          open(42, file=trim(outfile))
          write(42,*) ngrid, ngrid
       end if
       do i = 1,ngrid
          beta(i,l) = (betamax-betamin)/(ngrid-1.d0)*(i-1) + betamin
          m=0.d0
          do j=1,ngrid
             off(j,l) = (mmax-mmin)/(ngrid-1.d0)*(j-1) +mmin
             m(1,2) = off(j,l)
             chisq(i,j,l) =func_bm( [beta(i,l)] , m )
             if(l==1 )write(42,*) beta(i,l), m(1,2), chisq(i,j,l)
          end do
       end do
       if(l==1) close(42)
       pos=minloc(chisq(:,:,l))
       if (pos(1)==1 .or. pos(1)==ngrid .or.pos(2)==1 .or. pos(2)==ngrid) then
          p = healnan
          write(*,*) "NB: gridmin hit the border of the grid"
          cycle  !exit loop
       end if
    end do
    l=loops

    !reorganize matrix
write(*,*)"the beta at chisq max:",beta(pos(1),l)
    allocate(dbeta(ngrid))
    allocate(doff(ngrid))
    allocate(lik(ngrid,ngrid))
    allocate(vec(ngrid))
    allocate(Pb(ngrid),Pm(ngrid))

    !find the step size, dbeta
    dbeta(1)= (beta(2,l)-beta(1,l) )/2.d0
    doff(1) = (off(2,l) -off(1,l)  )/2.d0
    do i=2,ngrid-1
       dbeta(i)=( beta(i+1,l)-beta(i-1,l) )/2.d0
       doff(i) =( off(i+1,l) -off(i-1,l)  )/2.d0
    end do
    dbeta(ngrid)=(beta(ngrid,l)-beta(ngrid-1,l) )/2.d0
    doff(ngrid)=(off(ngrid,l)-off(ngrid-1,l) )/2.d0

    !evaluate likelihood from chisquare and normalize
    chisq(:,:,l)=chisq(:,:,l)-minval(chisq(:,:,l))
    lik(:,:)=exp( - chisq(:,:,l)/2.d0 )
    do i =1,ngrid
       vec(i)=sum(lik(i,:)*doff  )
    end do
    norm_factor=sum(vec*dbeta)
    lik(:,:)=lik(:,:)/norm_factor
    write(*,*)"norm_factor",norm_factor

    !marginalize
    do i=1,ngrid
       Pb(i)=sum(lik(i,:)*doff  )
       Pm(i)=sum(lik(:,i)*dbeta )
    end do

    !calculate spectral index and its standard deviation
    p(1) = sum(Pb* beta(:,l)*dbeta )
    m(1,2)=sum(Pm* off(:,l)* doff )
write(*,*)"b,m",p(1),m(1,2)

    uncert_b = sum( ((beta(:,l)-p(1))**2) *Pb *dbeta  )
    uncert_b = sqrt(uncert_b)
    uncert_m = sum( ((off(:,l)-m(1,2))**2)*Pm *dbeta  )
    uncert_m = sqrt(uncert_m)
write(*,*)"uncert",uncert_b,uncert_m    
    uncert=uncert_b

!write(*,*)"AA"
!call dump_matrix(grid)
    write(*,*) "spectral index from grid =", p(1), "+-",uncert

    deallocate(chisq, beta, off,dbeta,m, doff, lik ,vec, Pb, Pm)

  end subroutine gridmin_bm


  !!!!!!!!!!!!!!!!!!!!
  FUNCTION func_bm(p,m)
  !!!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(1:), INTENT(IN),optional :: p
    real(dp), dimension(:,:)          :: m
    REAL(dp)                           :: func_bm
    integer(i4b)                       :: v, n, i, ntot
    real(dp),        allocatable, dimension(:)   :: A, chi_vector, chi
    real(dp),        allocatable, dimension(:,:) :: Pm_vector, m_vector
    real(dp)                                     :: beta
integer, save :: it = 0

    allocate(chi(1:subreg))
chi=0.d0
    do i = 1, subreg
       n = d(i)%n
       ntot = d(i)%n
       if (k==4) ntot = d(i)%n*pol

       allocate(chi_vector(ntot))
       allocate(A(ntot))
       allocate(Pm_vector(ntot,nfreq))
       allocate(m_vector(d(i)%npix*mnum,nfreq))

       beta = p(1)
       m_vector = 0.d0 
       m_vector(1:n,2)= m(1,2)
       if (k==4)  m_vector(n+1:ntot,2) = m(2,2)
       Pm_vector = matmul(transpose(d(i)%P),m_vector)
       deallocate(m_vector)

       call calculate_Am3(p,A,m(:,:),i)

       chi(i)=0.d0
       do v = 1, nfreq
          chi_vector = d(i)%Pmap(:,1,v) - A(:)*a2t(v)*(freq(v)/freq(1))**beta - Pm_vector(:,v)
          chi(i) = chi(i) + dot_product( chi_vector, matmul( d(i)%Pcov(:,:,v) , chi_vector ) )
       end do

!!$       !write out stuff
!!$       write(*,*) "chi_vector"
!!$       call dump_matrix(chi_vector,fmt="(e16.8)")  
!!$       write(*,*) "Pm_vector"
!!$       call dump_matrix(Pm_vector,fmt="(e16.8)")  
!!$       write(*,*) "m"
!!$       call dump_matrix(m,fmt="(e16.8)")  
!!$       write(*,*) "A"
!!$       call dump_matrix(A,fmt="(e16.8)")  
!!$       write(*,*) "chi"
!!$       call dump_matrix(chi,fmt="(e16.8)")  
!!$       !end write out stuff

!write(*,*) m
       deallocate(chi_vector, A, Pm_vector)
    end do
!write(*,*)"chi",k,chi
    func_bm = SUM(chi)
!it=it+1
!write(*,'(a6,i6,'//itoa(2)//'e15.7)') 'func_bm', it, p, func_bm

    deallocate(chi)

  END FUNCTION func_bm


  !!!!!!!!!!!!!!!!!!!!
  FUNCTION func(p)
  !!!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(1:), INTENT(IN),optional :: p
    REAL(dp)                           :: func
    integer(i4b)                       :: v, n, i, ntot, pix
    real(dp),        allocatable, dimension(:)   :: A, chi_vector, chi
    real(dp),        allocatable, dimension(:,:) :: m, Pm_vector,m_vector
    real(dp)                                     :: beta
integer, save :: it = 0

    allocate(chi(1:subreg))
chi=0.d0
pix=0
    do i = 1, subreg
       n = d(i)%n
       ntot = d(i)%n
       if (k==4) ntot = d(i)%n*pol

       allocate(chi_vector(ntot))
       allocate(m(mnum,nfreq))
       allocate(A(ntot))
       allocate(m_vector(d(i)%npix*mnum,nfreq))
       allocate(Pm_vector(ntot,nfreq))

       m  = 0.d0
       beta = p(1)

       call calculate_Am(p,A,m(:,2),i)
       m_vector = 0.d0 
       m_vector(1:n,2)= m(1,2)
       if (k==4)  m_vector(n+1:ntot,2) = m(2,2)
       Pm_vector = matmul(transpose(d(i)%P),m_vector)
       deallocate(m_vector)

       chi(i)=0.d0
       do v = 1, nfreq
          chi_vector = d(i)%Pmap(:,1,v) - A(:)*a2t(v)*(freq(v)/freq(1))**beta - Pm_vector(:,v)
          chi(i) = chi(i) + dot_product( chi_vector, matmul( d(i)%Pcov(:,:,v) , chi_vector ) )
       end do

!!$       !write out stuff
!!$       write(*,*) "chi_vector"
!!$       call dump_matrix(chi_vector,fmt="(e16.8)")  
!!$       write(*,*) "Pm_vector"
!!$       call dump_matrix(Pm_vector,fmt="(e16.8)")  
!!$       write(*,*) "m"
!!$       call dump_matrix(m,fmt="(e16.8)")  
!!$       write(*,*) "A"
!!$       call dump_matrix(A,fmt="(e16.8)")  
!!$       write(*,*) "chi"
!!$       call dump_matrix(chi,fmt="(e16.8)")  
!!$       !end write out stuff

!write(*,*) m
       deallocate(chi_vector, m, A, Pm_vector)
       pix=pix+d(i)%n
    end do
!write(*,*)"chi",k,chi
    func = SUM(chi)
it=it+1
write(*,'(a6,i6,'//itoa(3)//'e15.7)') 'func', it, p, func, func/real(mnum*pix,dp)

    deallocate(chi)

  END FUNCTION func

  !!!!!!!!!!!!!!!!!!!
  FUNCTION dfunc(p)
  !!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(1:), INTENT(IN) :: p
    REAL(dp), DIMENSION(size(p))       :: dfunc
    integer(i4b)                       :: v, n, i, ntot
    real(dp),        allocatable, dimension(:)   :: A, chi_vector, vec1, minusto, dbeta, dm, dm2
    real(dp),        allocatable, dimension(:,:) :: m, Pm_vector, m_vector
    real(dp)                                     :: beta
integer, save :: it


!    allocate(dm(1:subreg))
!    allocate(dm2(1:subreg))    
    allocate(dbeta(1:subreg))
!    dm = 0.d0
!    dm2 = 0.d0
    dbeta=0.d0
    do i = 1, subreg
       n = d(i)%n
       ntot = d(i)%n
       if (k==4) ntot = d(i)%n*pol

       allocate(vec1(ntot))
!       allocate(minusto(ntot))
       allocate(chi_vector(ntot))
       allocate(m(mnum,nfreq))
       allocate(A(ntot))
       allocate(m_vector(d(i)%npix*mnum,nfreq))
       allocate(Pm_vector(ntot,nfreq))

       m  = 0.d0
       beta = p(1)
!       minusto(:) = -2.d0
       call calculate_Am(p,A,m(:,2),i)
       m_vector = 0.d0 
       m_vector(1:n,2)= m(1,2)
       if (k==4)  m_vector(n+1:ntot,2) = m(2,2)
       Pm_vector = matmul(transpose(d(i)%P),m_vector)
       deallocate(m_vector)

       do v = 1, nfreq
          chi_vector = d(i)%Pmap(:,1,v) - A(:)*a2t(v)*((freq(v)/freq(1))**beta) - Pm_vector(:,v)
          vec1 = -2.d0 * A(:) *a2t(v)*((freq(v)/freq(1))**beta) * log(freq(v)/freq(1))
          dbeta(i) = dbeta(i) + dot_product( vec1, matmul( d(i)%Pcov(:,:,v) , chi_vector ) )
       end do

       v = nfreq
!       dm(i) = dot_product( minusto, matmul( d(i)%Pcov(:,:,v) , chi_vector ) )  !only the v=2 comp for dm!
       deallocate(m, vec1, chi_vector, A, Pm_vector)
!       deallocate(minusto)
    end do

    dfunc = SUM(dbeta)
    deallocate(dbeta)
!    deallocate(dm,dm2)

!it=it+1
!write(*,'(a6,i6,'//itoa(2)//'e15.7)') 'dfunc', it, p, dfunc

  END FUNCTION dfunc



  !!!!!!!!!!!!!!!!!!
  function ddfunc(p)
  !!!!!!!!!!!!!!!!!!
    !NB! this is acutally the inverse of d2chisq!    
    implicit none
    real(dp), dimension(1:), intent(in) :: p
    real(dp),allocatable, dimension(:,:) :: mat
    real(dp), dimension(size(p),size(p)) :: ddfunc

!!$    if(k==4) then
!!$       allocate(mat(3,3))
!!$    else
!!$       allocate(mat(2,2))
!!$    end if
    allocate(mat(1,1))
    call calculate_cov_params2(p,mat)
    ddfunc(1,1) = mat(1,1)
    deallocate(mat)

  end function ddfunc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_cov_params2(p,cov_params)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !routine that calculates a matrix with the second derivatives of the chisquare
    implicit none
    REAL(dp),        DIMENSION(1:), INTENT(IN)    :: p
    real(dp),        dimension(1:,1:), intent(out) :: cov_params

    real(dp),        allocatable, dimension(:)   :: A, chi2_vector, vec1, envector
    real(dp),        allocatable, dimension(:)   :: d2b2, d2m2, d2bm, d2bw, d2w2, d2mw
    real(dp)                                     :: beta
    integer(i4b)                                 :: v, ntot, n, i
    real(dp),        allocatable, dimension(:,:) :: m, Pm_vector, m_vector

    allocate(d2b2(1:subreg))
    d2b2 = 0.d0

    do i = 1, subreg
       n = d(i)%n
       ntot = d(i)%n
       if (k==4) ntot = d(i)%n*pol

       allocate(chi2_vector(ntot))
       allocate(vec1(ntot))
       allocate(m(mnum,nfreq))
       allocate(A(ntot))
       allocate(m_vector(d(i)%npix*mnum, nfreq))
       allocate(Pm_vector(ntot, nfreq))
       m  = 0.d0
       beta = p(1)
       call calculate_Am(p,A,m(:,2),i)
       m_vector = 0.d0 
       m_vector(1:n,2)= m(1,2)
       if (k==4)  m_vector(n+1:ntot,2) = m(2,2)
       Pm_vector = matmul(transpose(d(i)%P),m_vector)

       do v = 1, nfreq
          chi2_vector = d(i)%Pmap(:,1,v) - 2*A(:)*a2t(v)*((freq(v)/freq(1))**beta) - Pm_vector(:,v)
          vec1 = -2.d0 * A(:) *a2t(v)*((freq(v)/freq(1))**beta) * ( (log(freq(v)/freq(1)))**2 )
          d2b2(i) = d2b2(i) + dot_product( vec1 , matmul( d(i)%Pcov(:,:,v) , chi2_vector ) )
       end do
       deallocate(A,m,chi2_vector,vec1,Pm_vector,m_vector)
    end do  !end loop subreg

    cov_params(1,1) = SUM(d2b2)
    !since the matrix should contain -lnL instead of chisq, we must divide by two!
    cov_params = cov_params/2.d0

    cov_params(1,1) = 1.d0 /  cov_params(1,1)   !invert it

    deallocate(d2b2)

!write(*,*) "B"
!call dump_matrix(cov_params,fmt='(e15.7)')


  end subroutine calculate_cov_params2


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calculate_cov_params_full(p,uncert)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !routine that calculates a matrix with all the second derivatives of the chisquare
    implicit none
    REAL(dp),        DIMENSION(:), INTENT(IN)    :: p
    REAL(dp),        INTENT(out)                 :: uncert

    real(dp),        allocatable, dimension(:,:) :: cov_params
    real(dp),        allocatable, dimension(:)   :: A, chi2_vector, vec1, Penvector,envector, cov_vector
    real(dp),        allocatable, dimension(:)   :: d2b2, d2bm, d2bw, d2Ab
    real(dp),        allocatable, dimension(:,:) :: d2A2, d2m2,d2w2,d2Am, d2Aw
    real(dp)                                     :: beta, factor,factor2
    integer(i4b)                                 :: v, n, i, ntot, totpix, t, regpix
    real(dp),        allocatable, dimension(:,:) :: m, Pm_vector, m_vector

    regpix = sum(d(:)%n)
    totpix = regpix*mnum
    write(*,*)"Number of pixels in region",regpix, totpix

    ! we call m(Q)= m and m(U)= w 
    allocate(d2bm(subreg))
    allocate(d2bw(subreg))

    allocate(d2w2(subreg,subreg))
    allocate(d2m2(subreg,subreg))

    allocate(d2b2(subreg))


    allocate(d2Am(subreg,totpix))
    allocate(d2Aw(subreg,totpix))

    allocate(d2A2(totpix,totpix))
    allocate(d2Ab(totpix))

    d2bw = 0.d0
    d2b2 = 0.d0
    d2m2 = 0.d0
    d2w2 = 0.d0
    d2Am = 0.d0
    d2Aw = 0.d0
    d2bm = 0.d0
    d2A2 = 0.d0
    d2Ab = 0.d0

    t=0
    do i = 1, subreg
       n = d(i)%n
       ntot = d(i)%n
       if (k==4) ntot = d(i)%n*pol

       allocate(chi2_vector(ntot))
       allocate(vec1(ntot))
       allocate(envector(d(i)%npix*mnum))
       allocate(Penvector(ntot))
       allocate(m(mnum,nfreq))
       allocate(A(ntot))
       allocate(m_vector(d(i)%npix*mnum, nfreq))
       allocate(Pm_vector(ntot, nfreq))
       m  = 0.d0
       beta = p(1)
       call calculate_Am(p,A,m(:,2),i)
       envector(:) = 1.d0
       m_vector = 0.d0 
       m_vector(1:n,2)= m(1,2)
       if (k==4)  m_vector(n+1:ntot,2) = m(2,2)
       Penvector = matmul(transpose(d(i)%P),envector)
       Pm_vector = matmul(transpose(d(i)%P),m_vector)
       deallocate(envector,m_vector)

       do v = 1, nfreq
          factor = a2t(v) *(freq(v)/freq(1))**beta
          chi2_vector = d(i)%Pmap(:,1,v) - 2*A(:)*factor - Pm_vector(:,v)
          vec1 = -2.d0 * A(:) *factor * log(freq(v)/freq(1))
          d2b2(i) = d2b2(i) + dot_product( log(freq(v)/freq(1)) *vec1 , matmul( d(i)%Pcov(:,:,v) , chi2_vector ) )
          d2A2(t+1:t+ntot,t+1:t+ntot) = d2A2(t+1:t+ntot,t+1:t+ntot) + 2.d0*factor**2 * d(i)%Pcov(:,:,v)
          factor2= -2.d0 *factor * log(freq(v)/freq(1))
          d2Ab(t+1:t+ntot) = d2Ab(t+1:t+ntot) + factor2* matmul( d(i)%Pcov(:,:,v) , chi2_vector )
       end do

       v = nfreq      !only the v=2 comp for derivation with respect to m!
       allocate(cov_vector(ntot))
       cov_vector = matmul( d(i)%Pcov(:,:,v) , Penvector )

       d2bm(i) =       dot_product( -vec1(1:n) ,   cov_vector(1:n) )
       d2m2(i,i) = 2.d0* dot_product( Penvector(1:n), cov_vector(1:n) ) 
       d2Am(i,t+1:t+n) = cov_vector(1:n) *a2t(v)*(freq(v)/freq(1))**beta 
       if (k==4) then
          d2bw(i) =       dot_product( -vec1(n+1:ntot) ,  cov_vector(n+1:ntot) )
          d2w2(i,i) = 2.d0* dot_product( Penvector(n+1:ntot),cov_vector(n+1:ntot) )
          d2Aw(i,t+1:t+n) = cov_vector(n+1:ntot) *a2t(v)*(freq(v)/freq(1))**beta 
       end if

       deallocate(A,m,Pm_vector,chi2_vector,vec1,Penvector, cov_vector)
       t=t+ntot
    end do !end subreg


    allocate(cov_params(totpix+subreg*mnum+1,totpix+subreg*mnum+1))

    cov_params(1:totpix,1:totpix) = d2A2
    cov_params(totpix+1:totpix+subreg,1:regpix) = d2Am
    cov_params(1:regpix,totpix+1:totpix+subreg) = transpose(d2Am)

    cov_params(totpix+1:totpix+subreg,totpix+1:totpix+subreg) = d2m2

    cov_params(1:totpix,totpix+1+mnum*subreg) = d2Ab
    cov_params(totpix+1+mnum*subreg,1:totpix) = d2Ab
    cov_params(totpix+1+mnum*subreg,totpix+1:totpix+subreg) = d2bm
    cov_params(totpix+1:totpix+subreg,totpix+1+mnum*subreg) = d2bm

    cov_params(totpix+1+mnum*subreg,totpix+1+mnum*subreg) = sum(d2b2)

    if (k==4) then
       cov_params(totpix+subreg+1:totpix+2*subreg,regpix+1:totpix) = d2Aw
       cov_params(regpix+1:totpix,totpix+subreg+1:totpix+2*subreg) = transpose(d2Aw)

       cov_params(totpix+1+subreg:totpix+subreg*2,totpix+subreg+1:totpix+subreg*2) = d2w2
       cov_params(totpix+1+mnum*subreg,totpix+subreg+1:totpix+mnum*subreg) = d2bw
       cov_params(totpix+subreg+1:totpix+mnum*subreg,totpix+1+mnum*subreg) = d2bw
    end if


    !since the matrix should contain -lnL instead of chisq, we must divide by two!
    cov_params = cov_params/2.d0

!write(*,*) "B"
!call dump_matrix(cov_params,fmt='(e15.7)')

    deallocate(d2bm,d2bw,d2m2,d2w2,d2b2,d2Am,d2Aw,d2A2,d2Ab)
    write(*,*)"uncert before invert",cov_params(totpix+1+mnum*subreg,totpix+1+mnum*subreg)

    call invert_singular_matrix(cov_params,1.d-20) !when using invert_matrix() I get some Nan´s, so must use singular
!call dump_matrix(cov_params,fmt='(e15.7)')
    write(*,*)"uncert after invert",cov_params(totpix+1+mnum*subreg,totpix+1+mnum*subreg)

    uncert = sqrt(cov_params(totpix+1+mnum*subreg,totpix+1+mnum*subreg))

    deallocate(cov_params)

  end subroutine calculate_cov_params_full


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_region_A(p,offset,fg)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    REAL(dp),     dimension(1:),    INTENT(IN)    :: p
    real(dp),     dimension(1:,1:), intent(out)   :: offset
    real(dp),     dimension(1:,1:), intent(out)   :: fg
    integer(i4b)                                  :: i, lower, upper, n!, teller, v, j
    real(dp),     allocatable, dimension(:)       :: A, m
!    real(dp)                                      :: factor

    allocate(m(2))
    lower = 0
    upper = 0

    if (k==4) then
       do i = 1, subreg
          n = d(i)%n
          lower = upper +1
          upper = upper +n
          allocate(A(1:2*n))
          call calculate_Am(p,A,m,i)
          fg(lower:upper,1) = A(1:n)  
          fg(lower:upper,2) = A(n+1:2*n)  
          offset(lower:upper,1) = m(1)
          offset(lower:upper,2) = m(2)
          write(*,*) "offset for Q+U, subreg=",i, m(1), m(2) 
          deallocate(A)
       end do
    else
       do i = 1, subreg
          lower = upper +1
          upper = upper+d(i)%n
          call calculate_Am(p,fg(lower:upper,1),m,i)
          offset(lower:upper,1) = m(1)
       end do
    end if
 
!!$    !make a chisquare map
!!$    m(1) =0.d0
!!$    chi2=0.d0
!!$    teller=1
!!$    do i =1,subreg
!!$       n=d(i)%n
!!$       do j=1,n
!!$          m(2) =offset(teller,1)
!!$          do v=1,nfreq
!!$             factor = d(i)%Pmap(j,1,v) - fg(teller,1)*a2t(v)*(freq(v)/freq(1))**p(1) - m(v)
!!$             chi2(teller,1) =chi2(teller,1)+ factor**2 !* d(i)%Pcov(j,j,v)
!!$             if(k==4) then
!!$                m(2) =offset(teller,2)
!!$                factor = d(i)%Pmap(j+n,1,v) - fg(teller,2)*a2t(v)*(freq(v)/freq(1))**p(1) - m(v)
!!$                chi2(teller,2) =chi2(teller,2)+ factor**2 !* d(i)%Pcov(j+n,j+n,v)
!!$             end if
!!$          end do
!!$          teller=teller+1
!!$       end do
!!$    end do

    !must convert to thermal temperature
    fg = fg*a2t(1)



    deallocate(m)

  end subroutine get_region_A


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_invNelements(i, pixe,freqv,invN)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer(i4b),                    INTENT(IN)   :: i, pixe, freqv
    real(dp),     dimension(1:,1:), intent(out)   :: invN

    invN(1,1) = d(i)%Pcov(pixe,pixe,freqv)
    invN(2,2) = d(i)%Pcov(pixe+d(i)%n,pixe+d(i)%n,freqv)
    invN(1,2) = d(i)%Pcov(pixe,pixe+d(i)%n,freqv)
    invN(2,1) = invN(1,2)

  end subroutine get_invNelements



  !!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_chisquare2(p)
  !!!!!!!!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(1:), INTENT(IN):: p
    REAL(dp)                           :: get_chisquare2
    integer(i4b)                       :: v, n, i, ntot
    real(dp)                           :: A(4), d_vec(4,2), Te(4,2), chisquares(4)
    real(dp)                           :: m(2)
    real(dp)                           :: beta, alpha

    beta = p(1)
    alpha=a2t(2)*(freq(2)/freq(1))**beta
    m = [1,1]
    m = (1.d0 - alpha/2.d0) *m
    A(1) = - alpha * (1.d0-alpha/2.d0)/(1+alpha*alpha)
    A(2) = (1.d0 + alpha + alpha*alpha/2.d0) / (1.d0+alpha*alpha)
    A(3) = A(2)
    A(4) = A(1)

    d_vec(:,1)=[0,1,1,0]
    d_vec(:,2)=[0,2,2,0]
    Te(:,1) = [1,1,0,0]
    Te(:,2) = [0,0,1,1]    
    do i=1,4
       v=1
       chisquares(i) = (d_vec(i,v)-A(i) )  **2 
       v=2
       chisquares(i) =chisquares(i)+ (d_vec(i,v)-A(i)*alpha -  dot_product(m,Te(i,:) ) ) **2 
    end do
    get_chisquare2 = sum(chisquares)
write(*,*)"chisq2=",get_chisquare2  

  END FUNCTION get_chisquare2

  !!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION get_chisquare3(p)
  !!!!!!!!!!!!!!!!!!!!!!!!!
    USE healpix_types
    IMPLICIT NONE
    REAL(dp), DIMENSION(1:), INTENT(IN):: p
    REAL(dp)                           :: get_chisquare3
    integer(i4b)                       :: v, n, i, ntot
    real(dp)                           :: A(4), d_vec(4,2), Te(4,2), chisquares(4)
    real(dp)                           :: m(2)
    real(dp)                           :: beta, alpha

    beta = p(1)
    alpha=a2t(2)*(freq(2)/freq(1))**beta
    m = [1,1]
    m = (1.d0 - alpha/2.d0) *m *2.d0
    A(1:4) = 1.d0

    d_vec(:,1)=[1,1,1,1]
    d_vec(:,2)=[2,2,2,2]
    Te(:,1) = [1,1,0,0]
    Te(:,2) = [0,0,1,1]    
    do i=1,4
       v=1
       chisquares(i) = (d_vec(i,v)-A(i) )  **2 
       v=2
       chisquares(i) =chisquares(i)+ (d_vec(i,v)-A(i)*alpha -  dot_product(m,Te(i,:) ) ) **2 
    end do
    get_chisquare3 = sum(chisquares)
write(*,*)"chisq2=",get_chisquare3

  END FUNCTION get_chisquare3


end module sindex_maxlike_mod
