program crosscorr_pix
  use quiet_utils
  use quiet_fileutils
  use alm_tools
  implicit none

    integer(i4b) :: i, j, k, l, m, q, nside, npix, nmaps, lmax, num_fg_comp, ordering, nside_corr, n, bf(2), nx, ny
    real(dp)     :: zbounds(2), mu1, mu2, sigma1, sigma2, unit, x0, y0, dx, dy
    complex(dpc) :: corr_cmplx, auto1, auto2
    real(dp), allocatable, dimension(:) :: mu, sigma, x,y
    real(dp), allocatable, dimension(:,:) :: weights, mask, map1, map2, corr, chisq, chisqmap
    complex(dpc), allocatable, dimension(:,:,:,:) :: alms
    complex(dpc), allocatable, dimension(:,:,:)   :: alms_cmb
    logical(lgt) :: exist
    character(len=512) :: filename, prefix, temp, file1, file2, filelist, maskfile, chisqfile
    character(len=2) :: chain_text
    character(len=5) :: iter_text
    character(len=2) :: c1_text, c2_text

    if (iargc() < 7) then
       write(*,*) 'Usage: crosscorr_pix [nside] [ngrid] [X gridsize] [Y gridsize] [filelist] [mask] [prefix]'
       stop
    end if
    unit = 58

    call getarg(1,temp)
    read(temp,*) nside
    npix = 12*nside**2

    call getarg(2,temp)
    read(temp,*) n

    call getarg(3,temp)
    read(temp,*) dx

    call getarg(4,temp)
    read(temp,*) dy

    call getarg(5,filelist)
    call getarg(6,maskfile)
    call getarg(7,prefix)

    write(*,*) 'Reading mask = ', trim(maskfile)
    call read_map(mask, ordering, maskfile)
    if (ordering == 2) call convert_nest2ring(nside, mask(:,1))

    allocate(corr(-n:n,-n:n), chisq(-n:n,-n:n), x(-n:n), y(-n:n))
    corr = 2
    chisq = 1.d30
    do i = -n, n
       x(i) = i*dx
       y(i) = i*dy
    end do

    open(unit,file=trim(filelist))
    do while (.true.) 
       read(unit,*,end=10) x0, y0, file1, file2, chisqfile
       inquire(file=trim(file1), exist=exist)
       if (.not. exist) cycle
       inquire(file=trim(file2), exist=exist)
       if (.not. exist) cycle
       inquire(file=trim(chisqfile), exist=exist)
       if (.not. exist) cycle
       write(*,*) 'Processing ', x0, y0
       call read_map(map1, ordering, file1)
       if (ordering == 2) call convert_nest2ring(nside, map1(:,1))
       call read_map(map2, ordering, file2)
       if (ordering == 2) call convert_nest2ring(nside, map2(:,1))
       call read_map(chisqmap, ordering, chisqfile)
       if (ordering == 2) call convert_nest2ring(nside, chisqmap(:,1))
       
       i = nint(x0/dx)
       j = nint(y0/dy)

       ! Cross-correlation
       m          = sum(mask(:,1))
       mu1        = sum(map1(:,1)*mask(:,1)) / m
       sigma1     = sqrt(sum((map1(:,1)-mu1)**2*mask(:,1)) / (m-1))
       mu2        = sum(map2(:,1)*mask(:,1)) / m
       sigma2     = sqrt(sum((map2(:,1)-mu2)**2*mask(:,1)) / (m-1))
       corr(i,j) = sum((map1(:,1)-mu1)*(map2(:,1)-mu2)*mask(:,1))/(sigma1*sigma2)/(m-1)

       ! Chisq
       chisq(i,j) = sum(chisqmap(:,1)*mask(:,1)) / m
       write(*,*) corr(i,j), chisq(i,j)

       deallocate(map1,map2,chisqmap)
    end do
10  close(unit)

    
    nx = 0
    do i = -n, n
       if (any(corr(i,:) /= 2.d0)) nx = nx+1
    end do

    ny = 0
    do i = -n, n
       if (any(corr(:,i) /= 2.d0)) ny = ny+1
    end do

    write(*,*) nx, ny

    ! Output correlation coefficients
    filename = trim(prefix) // '_corr.dat'
    open(unit,file=trim(filename), recl=1024)
    write(unit,*) nx, ny
    do i = -n, n
       if (any(corr(i,:) /= 2.d0)) then
          do j = -n, n
             if (any(corr(:,j) /= 2.d0)) then
                write(unit,*) x(i), y(j), corr(i,j), chisq(i,j)
             end if
          end do
       end if
    end do
    close(unit)

    ! Find best-fit points
    corr = abs(corr)
    bf    = minloc(corr)-n-1
    write(*,fmt='(a,2f8.3)') 'Minimum correlation     = ', x(bf(1)), y(bf(2))
    write(*,fmt='(a,2f8.3)') 'Lowest correlation      = ', corr(bf(1), bf(2))

    bf    = minloc(chisq)-n-1
    write(*,fmt='(a,2f8.3)') 'Minimum chisq           = ', x(bf(1)), y(bf(2))
    write(*,fmt='(a,2f8.3)') 'Lowest chisq            = ', chisq(bf(1), bf(2))

    
  end program crosscorr_pix
