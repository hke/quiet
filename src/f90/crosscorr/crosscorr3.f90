program crosscorr
  use quiet_utils
  use quiet_fileutils
  use alm_tools
  implicit none

    integer(i4b) :: i, j, k, l, m, q, nside, npix, nmaps, lmax, num_fg_comp, ordering, nside_corr
    real(dp)     :: zbounds(2), n, corr, mu1, mu2, sigma1, sigma2
    complex(dpc) :: corr_cmplx, auto1, auto2
    real(dp), allocatable, dimension(:) :: mu, sigma, threshold
    real(dp), allocatable, dimension(:,:) :: weights, mask, corr_map, inmap, mask2
    real(dp), allocatable, dimension(:,:,:) :: maps
    complex(dpc), allocatable, dimension(:,:,:,:) :: alms
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmb
    character(len=512) :: filename, prefix, tmp, outfile
    character(len=2) :: chain_text
    character(len=5) :: iter_text
    character(len=2) :: c1_text, c2_text

    if (iargc() < 6) then
       write(*,*) 'Usage: crosscorr [maskfile] [outfile] [map1] [threshold1] [map2] [threshold2] ...'
       stop
    end if

    call getarg(1, filename)
    write(*,*) 'Reading mask = ', trim(filename)
    call read_map(mask, ordering, filename)
    npix  = size(mask,1)
    nmaps = size(mask,2)
    nside = npix2nside(npix)
    if (ordering == 2) then
       do i = 1, nmaps
          call convert_nest2ring(nside, mask(:,i))
       end do
    end if
    allocate(mask2(0:npix-1,nmaps))

    call getarg(2, outfile)

    num_fg_comp = (iargc()-2)/2
    allocate(maps(0:npix-1,nmaps,num_fg_comp), threshold(num_fg_comp))
    do i = 1, num_fg_comp
       call getarg(2+(2*i)-1, filename)
       call getarg(2+(2*i),   tmp)
       read(tmp,*) threshold(i)
       write(*,*) 'Reading map = ', trim(filename)
       call read_map(inmap, ordering, filename)
       maps(:,:,i) = inmap
       if (ordering == 2) then
          do j = 1, nmaps
             call convert_nest2ring(nside, maps(:,j,i))
          end do
       end if
       deallocate(inmap)
       if (threshold(i) /= 1.d30) then
          where (maps(:,:,i) < threshold(i))
             maps(:,:,i) = 0.d0
          end where
       end if
    end do

    ! Output template coefficients
    open(60,file=trim(outfile), recl=1024) 
    do i = 1, num_fg_comp
       do j = 1, i-1
          mask2 = mask
          where (maps(:,:,i) == 0.d0 .or. maps(:,:,j) == 0.d0)
             mask2(:,:) = 0.d0
          end where

          corr = sum(maps(:,:,i)*maps(:,:,j)*mask2(:,:)) / sum(maps(:,:,i)*maps(:,:,i)*mask2(:,:))
          write(*,*) sum(maps(:,:,i)), sum(maps(:,:,j)), sum(mask2)

          write(60,fmt='(f10.6,a)',advance="no") corr, ' '

       end do
       write(60,*) 
    end do
    close(60)


  end program crosscorr
