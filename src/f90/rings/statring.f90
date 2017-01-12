program statring
  use healpix_types
  use pix_tools
  use fitstools
  use quiet_fileutils
  use math_tools
  use udgrade_nr
  use quiet_utils
  implicit none

  include "mpif.h"

  integer(i4b)       :: i, j, p, bin, numbins, numsims, maxpix
  integer(i4b)       :: nside, npix, ordering, nmaps, unit, npix_in, ordering_in, nmaps_in
  character(len=256) :: parfile, filename, filelist, outdir, prefix
  real(dp)           :: minverdi, maxverdi, gap, nullval, pixval
  logical(lgt)       :: anynull
  real(dp), allocatable, dimension(:)     :: av, sigma
  real(dp), allocatable, dimension(:,:)   :: hist, maxvals
  real(dp), allocatable, dimension(:,:,:) :: maps

  if (iargc() == 0) then
!     call print_usage
  else 
     write(*,*) '-------------------- Concentric rings statistics  --------------------'
     write(*,*)
  end if

  ! Read parameters
  call getarg(1,parfile)
  unit=13
  call get_parameter(unit, parfile, 'FILELIST',           par_string=filelist)
  call get_parameter(unit, parfile, 'OUTDIR',             par_string=outdir)
  call get_parameter(unit, parfile, 'PREFIX',             par_string=prefix)
  call get_parameter(unit, parfile, 'NUMSIMS',            par_int   =numsims)
  call get_parameter(unit, parfile, 'NUMBINS',            par_int   =numbins)
  write(*,*) 'filelist = ', trim(filelist)
  write(*,*) 'outdir = ', trim(outdir)
  write(*,*) 'prefix = ', trim(prefix)
  write(*,*) 'numsims = ', numsims
  write(*,*) 'numbins = ', numbins


  ! Read wmap
  open(42,file=trim(filelist))
  read(42,'(a)') filename
  npix = getsize_fits(trim(filename), nmaps=nmaps, ordering=ordering, nside=nside)
  write(*,*) nmaps, ordering, nside, npix
  allocate(maps(0:npix-1,1,0:numsims))
  call read_bintab(filename, maps(:,:,0),  npix, nmaps, nullval, anynull)
 write(*,*) nmaps, ordering, nside, npix
  allocate(av(1:numbins))
  allocate(sigma(1:numbins))
  allocate(hist(0:numsims, 1:numbins))
  allocate(maxvals(0:numsims,2))
  maxvals = 0.d0
  hist = 0.d0
  av = 0.d0
  sigma = 0.d0

  ! Read sims
  do i = 1, numsims  
     read(42,'(a)') filename
     call read_bintab(trim(filename), maps(:,:,i),  npix, nmaps, nullval, anynull)
 !    if (npix_in /= npix .or. nmaps_in /=1) then
 !       write(*,*)npix_in, '=npix,',nmaps_in,'=nmaps. Quiting'
 !       stop
 !    end if
  end do
  close(42)

  maxverdi = 0.d0
  minverdi = 1d30
  do i = 0, numsims
     do p = 0, npix-1
        pixval = maps(p,1,i)
        if (pixval /=-1.6375d30) then
           if (pixval > maxverdi) then
              maxverdi = pixval
              maxpix = p
           end if
           if (pixval < minverdi) minverdi = pixval
        end if
     end do
     maxvals(i,1) = maxverdi
     maxvals(i,2) = maxpix
     maxverdi = 0.d0
  end do
  maxverdi = maxval(maxvals(:,1))
!  minverdi=minval(maps(:,:,1:numsims))
!  maxverdi=maxval(maps(:,:,1:numsims))
  gap=maxverdi-minverdi
  write(*,*) 'min max gap', minverdi, maxverdi, gap
  do i = 0, numsims
     do p = 0, npix-1
        pixval = maps(p,1,i)
        if (pixval /=-1.6375d30) then
           bin = min(int((maps(p,1,i)-minverdi)/gap*numbins)+1,numbins)
           hist(i,bin) = hist(i,bin) + 1.d0
        end if
     end do
  end do

  do bin = 1, numbins
     av(bin) = sum(hist(1:numsims,bin))/real(numsims,dp)
     hist(1:numsims,bin) = (hist(1:numsims,bin)-av(bin))**2
     sigma(bin) =  sqrt(sum(hist(1:numsims,bin))/real(numsims-1,dp))
  end do

!  do bin = 1, numbins
!     write(*,*) bin, (bin-0.5)/30 *gap+minverdi
!  end do


  filename = trim(outdir) // '/' // trim(prefix) // '_stat.txt'
  open(43, file=trim(filename))
  do bin = 1, numbins
     write(43,*) (bin-0.5)/30 *gap+minverdi, av(bin)
  end do
  write(43,*)
  do bin = 1, numbins
     write(43,*) (bin-0.5)/30 *gap+minverdi, av(bin)-sigma(bin)
  end do
  do bin = numbins,1,-1
     write(43,*) (bin-0.5)/30 *gap+minverdi, av(bin)+sigma(bin)
  end do
  write(43,*)
  do bin = 1, numbins
     write(43,*) (bin-0.5)/30 *gap+minverdi, av(bin)-2.d0*sigma(bin)
  end do
  do bin = numbins, 1, -1
     write(43,*) (bin-0.5)/30 *gap+minverdi, av(bin)+2.d0*sigma(bin)
  end do
  write(43,*)
  do bin = 1, numbins
     write(43,*) (bin-0.5)/30 *gap+minverdi, hist(0,bin)
  end do
  close(43)

  filename = trim(outdir) // '/' // trim(prefix) // '_hist.txt'
  open(43, file=trim(filename))
  do i = 0, numsims
     write(43,*) i, maxvals(i,1), maxvals(i,2)
  end do
  close(43)


end program statring
