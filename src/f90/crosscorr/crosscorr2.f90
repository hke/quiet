program crosscorr
  use quiet_utils
  use quiet_fileutils
  use alm_tools
  implicit none

    integer(i4b) :: i, j, k, l, m, q, nside, npix, nmaps, lmax, num_fg_comp, ordering, nside_corr
    real(dp)     :: zbounds(2), n, corr, mu1, mu2, sigma1, sigma2, sigma0
    complex(dpc) :: corr_cmplx, auto1, auto2
    real(dp), allocatable, dimension(:) :: mu, sigma, threshold
    real(dp), allocatable, dimension(:,:) :: weights, mask, corr_map, inmap, mask2
    real(dp), allocatable, dimension(:,:,:) :: maps
    complex(dpc), allocatable, dimension(:,:,:,:) :: alms
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmb
    character(len=512) :: filename, prefix, tmp
    character(len=6), allocatable, dimension(:) :: labels
    character(len=2) :: chain_text, itext, jtext
    character(len=5) :: iter_text
    character(len=2) :: c1_text, c2_text

    if (iargc() < 7) then
       write(*,*) 'Usage: crosscorr [lmax] [nside] [nside_corr] [maskfile] [prefix] '
       write(*,*) '                      [map1] [label1] [threshold1] [map2] [label2] [threshold2] ...'
       stop
    end if

    call getarg(1,prefix)
    read(prefix,*) lmax

    call getarg(2,prefix)
    read(prefix,*) nside
    npix = 12*nside**2
    nmaps       = 1

    call getarg(3,prefix)
    read(prefix,*) nside_corr

    call getarg(4, filename)
!    allocate(mask(0:npix-1,1), inmap(0:npix-1,1))
    write(*,*) 'Reading mask = ', trim(filename)
    call read_map(mask, ordering, filename)
    if (ordering == 2) call convert_nest2ring(nside, mask(:,1))
    allocate(mask2(0:npix-1,1))

    call getarg(5, prefix)

    sigma0      = 70
    num_fg_comp = (iargc()-5)/3
    allocate(maps(0:npix-1,1,num_fg_comp), threshold(num_fg_comp), labels(num_fg_comp))
    do i = 1, num_fg_comp
       call getarg(5+(3*i)-2, filename)
       call getarg(5+(3*i)-1, labels(i))
       call getarg(5+(3*i),   tmp)
       read(tmp,*) threshold(i)
       write(*,*) 'Reading map = ', trim(filename)
       call read_map(inmap, ordering, filename)
       maps(:,:,i) = inmap
       if (ordering == 2) call convert_nest2ring(nside, maps(:,1,i))
       deallocate(inmap)
       if (threshold(i) /= 1.d30) then
          where (maps(:,1,i) < threshold(i))
             maps(:,1,i) = 0.d0
          end where
       end if
    end do

    ! Output real-space correlation coefficients
    filename = trim(prefix) // '_fg_pix_corr.dat'
    open(58,file=trim(filename), recl=1024) 
    filename = trim(prefix) // '_temp_amp.dat'
    open(60,file=trim(filename), recl=1024) 
    do i = 1, num_fg_comp
       call int2string(i, itext)
       do j = 1, i-1
          call int2string(j, jtext)
          mask2 = mask
          where (maps(:,1,i) == 0.d0 .or. maps(:,1,j) == 0.d0)
             mask2(:,1) = 0.d0
          end where

          n      = sum(mask2(:,1))
          mu1    = sum(mask2(:,1) * maps(:,1,i)) / n
          mu2    = sum(mask2(:,1) * maps(:,1,j)) / n
          sigma1 = sqrt(sum(mask2(:,1)*(maps(:,1,i)-mu1)**2) / (n-1))
          sigma2 = sqrt(sum(mask2(:,1)*(maps(:,1,j)-mu2)**2) / (n-1))

          corr = sum(mask2(:,1) * (maps(:,1,i)-mu1) * (maps(:,1,j)-mu2)) / (n-1)
          corr = corr / (sigma1*sigma2)
          write(58,fmt='(f6.3,a)',advance="no") corr, ' '

!!$          open(59,file=trim(prefix)//'_scatter_'//trim(labels(i))//'_'//trim(labels(j))//'.dat')
!!$          do k = 0, npix-1
!!$             if (mask2(k,1) > 0.5d0) then
!!$                write(59,*) maps(k,1,i), maps(k,1,j)
!!$             end if
!!$          end do
!!$          close(59)

          corr = sum((maps(:,1,i)*maps(:,1,j)*mask2(:,1)/sigma0**2 / sigma1 / sigma2)) / &
               & sum((maps(:,1,i)*maps(:,1,i)*mask2(:,1)/sigma0**2 / sigma1 / sigma1))

          write(60,fmt='(f6.3,a)',advance="no") corr, ' '

       end do
       write(58,*) labels(i)
       write(60,*) labels(i)
    end do
    close(58)
    close(60)


  end program crosscorr
