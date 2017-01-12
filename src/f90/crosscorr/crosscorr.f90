program crosscorr
  use quiet_utils
  use quiet_fileutils
  use alm_tools
  implicit none

    integer(i4b) :: i, j, k, l, m, q, nside, npix, nmaps, lmax, num_fg_comp, ordering, nside_corr
    real(dp)     :: zbounds(2), n, corr, mu1, mu2, sigma1, sigma2
    complex(dpc) :: corr_cmplx, auto1, auto2
    real(dp), allocatable, dimension(:) :: mu, sigma
    real(dp), allocatable, dimension(:,:) :: weights, mask, corr_map, inmap
    real(dp), allocatable, dimension(:,:,:) :: maps
    complex(dpc), allocatable, dimension(:,:,:,:) :: alms
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmb
    character(len=512) :: filename, prefix
    character(len=2) :: chain_text
    character(len=5) :: iter_text
    character(len=2) :: c1_text, c2_text

    if (iargc() < 7) then
       write(*,*) 'Usage: crosscorr [lmax] [nside] [nside_corr] [maskfile] [prefix] '
       write(*,*) '                      [map1] [map2] [map3] ...'
       stop
    end if

    call getarg(1,prefix)
    read(prefix,*) lmax

    call getarg(2,prefix)
    read(prefix,*) nside
    npix = 12*nside**2
    nmaps       = 1
    num_fg_comp = 2

    call getarg(3,prefix)
    read(prefix,*) nside_corr

    call getarg(4, filename)
!    allocate(mask(0:npix-1,1), inmap(0:npix-1,1))
    write(*,*) 'Reading mask = ', trim(filename)
    call read_map(mask, ordering, filename)
    if (ordering == 2) call convert_nest2ring(nside, mask(:,1))

    call getarg(5, prefix)

    num_fg_comp = iargc()-5
    allocate(maps(0:npix-1,1,num_fg_comp))
    do i = 1, num_fg_comp
       call getarg(5+i, filename)
       write(*,*) 'Reading map = ', trim(filename)
       call read_map(inmap, ordering, filename)
       maps(:,:,i) = inmap
       if (ordering == 2) call convert_nest2ring(nside, maps(:,1,i))
       deallocate(inmap)
    end do

    allocate(alms(nmaps,0:lmax,0:lmax, num_fg_comp), weights(1:2*nside, nmaps))
    weights = 1.d0
    zbounds = 0.d0

    ! Compute mean and standard deviation
    allocate(mu(num_fg_comp), sigma(num_fg_comp))
    mu = 0.d0
    n  = sum(mask(:,1))
    do j = 1, num_fg_comp
       mu(j) = sum(mask(:,1) * maps(:,1,j)) / n
    end do

    sigma = 0.d0
    do j = 1, num_fg_comp
       sigma(j) = sqrt(sum(mask(:,1)*(maps(:,1,j)-mu(j))**2) / (n-1))
    end do
    
    ! Compute alms
    alms = 0.d0
    do i = 1, num_fg_comp
       write(*,*) 'Computing alms for map ', i
       if (nmaps == 1) then
          call map2alm(nside, lmax, lmax, (maps(:,1,i)-mu(i))*mask(:,1), alms(:,:,:,i), zbounds, weights)
       else
          call map2alm(nside, lmax, lmax, (maps(:,:,i)-mu(i))*mask,      alms(:,:,:,i), zbounds, weights)
       end if
    end do
    
    call convert_ring2nest(nside, mask(:,1))
    do i = 1, num_fg_comp
       call convert_ring2nest(nside, maps(:,1,i))
    end do
    
    ! Output real-space correlation coefficients
    filename = trim(prefix) // '_fg_pix_corr.dat'
    open(58,file=trim(filename), recl=1024)
    do i = 1, num_fg_comp
       do j = 1, i
          corr = sum(mask(:,1) * (maps(:,1,i)-mu(i)) * (maps(:,1,j)-mu(j))) / (n-1)
          corr = corr / (sigma(i)*sigma(j))
          write(58,fmt='(f6.3,a)',advance="no") corr, ' '
       end do
       write(58,*) 
    end do

    ! Output cross spectra
    do i = 1, num_fg_comp
       call int2string(i,c1_text)
       do j = i+1, num_fg_comp
          call int2string(j,c2_text)
          filename = trim(prefix) // '_fg_cross_corr_' // c1_text // '_' // c2_text // '.dat'
          open(58,file=trim(filename))
          do l = 2, lmax
             corr_cmplx = alms(1,l,0,i)*alms(1,l,0,j)
             auto1      = alms(1,l,0,i)*alms(1,l,0,i)
             auto2      = alms(1,l,0,j)*alms(1,l,0,j)
             do m = 1, l
                corr_cmplx = corr_cmplx + alms(1,l,m,i) * conjg(alms(1,l,m,j)) + alms(1,l,m,j) * conjg(alms(1,l,m,i))
                auto1      = auto1      + 2.d0 * alms(1,l,m,i) * conjg(alms(1,l,m,i))
                auto2      = auto2      + 2.d0 * alms(1,l,m,j) * conjg(alms(1,l,m,j))
             end do
             write(58,*) l, real(corr_cmplx,sp) / sqrt(real(auto1,dp)*real(auto2,dp))
          end do
          close(58)
       end do
    end do

    ! Output cross-maps
    allocate(corr_map(0:12*nside_corr**2-1,1))
    q = (nside / nside_corr)**2
    do i = 1, num_fg_comp
       call int2string(i,c1_text)
       do j = i+1, num_fg_comp
          call int2string(j,c2_text)
          filename = trim(prefix) // '_fg_cross_map_' // c1_text // '_' // c2_text // '.fits'

          corr_map = -1.6375d30
          do k = 0, 12*nside_corr**2-1
             n = sum(mask(k*q:(k+1)*q-1,1))
             if (n < 3) cycle
             mu1    = sum(mask(k*q:(k+1)*q-1,1) * maps(k*q:(k+1)*q-1,1,i)) / n
             mu2    = sum(mask(k*q:(k+1)*q-1,1) * maps(k*q:(k+1)*q-1,1,j)) / n
             sigma1 = sqrt(sum(mask(k*q:(k+1)*q-1,1) * (maps(k*q:(k+1)*q-1,1,i)-mu1)**2) / (n-1))
             sigma2 = sqrt(sum(mask(k*q:(k+1)*q-1,1) * (maps(k*q:(k+1)*q-1,1,j)-mu2)**2) / (n-1))
             corr_map(k,1) = sum(mask(k*q:(k+1)*q-1,1) * (maps(k*q:(k+1)*q-1,1,i)-mu1) * (maps(k*q:(k+1)*q-1,1,j)-mu2)) / (n-1) / (sigma1*sigma2)
          end do
          call convert_nest2ring(nside_corr, corr_map(:,1))
          call write_map(corr_map, 1, filename)             
       end do
    end do

    deallocate(alms, mu, sigma, weights, mask, corr_map)

  end program crosscorr
