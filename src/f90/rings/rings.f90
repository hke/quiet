program rings
  use healpix_types
  use pix_tools
  use fitstools
  use quiet_fileutils
  use quiet_mapfile_mod
  use math_tools
  use udgrade_nr
  use quiet_utils
  implicit none

  include "mpif.h"

  integer(i4b)       :: iargc, i, j, k, q, n
  integer(i4b)       :: unit, myid, numprocs, ierr, root, numbin, nside_center, npix_center, minbin
  integer(i4b)       :: nside, npix, ordering, nmaps, nummaps, numrings, count2, mycount2
  character(len=30)  :: kommando
  character(len=4)   :: real_text
  character(len=8)   :: i_text
  character(len=256) :: mapname, maskname, parfile, filename, filelist, outdir
  logical(lgt)       :: output_histograms
  real(dp)           :: threshold, binsize, mu, l, b, scale
  real(dp),              dimension(3)     :: vec
  real(dp), allocatable, dimension(:)     :: rms, avg, chisq
  !real(dp), allocatable, dimension(:,:)   :: mask, map
  real(dp), pointer,     dimension(:,:)   :: mask, map, pos, mypos, pos_out, myrms, cov, pte, chisq_out, rms_buff
  real(dp), pointer,     dimension(:,:,:) :: totrms

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0
  unit        = 42+myid

  if (iargc() == 0 .and. myid == root) then
     call print_usage
  else if (myid == root) then
     write(*,*) '-------------------- Concentric rings searcher  --------------------'
     write(*,*)
  end if

  ! Read parameters
  call getarg(1,parfile)

  call get_parameter(unit, parfile, 'CMB_SCALE_FACTOR',   par_dp=scale)
  call get_parameter(unit, parfile, 'MASKFILE',           par_string=maskname)
  call get_parameter(unit, parfile, 'RMS_THRESHOLD',      par_dp=threshold)
  call get_parameter(unit, parfile, 'BINSIZE',            par_dp=binsize)
  call get_parameter(unit, parfile, 'NUMBIN',             par_int=numbin)
  call get_parameter(unit, parfile, 'MINBIN',             par_int=minbin)
  call get_parameter(unit, parfile, 'NSIDE_CENTER',       par_int=nside_center)
  call get_parameter(unit, parfile, 'NUM_REALIZATIONS',   par_int=nummaps)
  call get_parameter(unit, parfile, 'FILELIST',           par_string=filelist)
  call get_parameter(unit, parfile, 'OUTDIR',             par_string=outdir)
  call get_parameter(unit, parfile, 'OUTPUT_HISTOGRAMS',  par_lgt=output_histograms)
  npix_center = 12*nside_center**2

  call read_map_dp(maskname, nside, ordering, nmaps, mask) 
  npix = 12*nside**2
  if (ordering == 2) call convert_nest2ring(nside, mask(:,1))

  ! Allocate work arrays
  allocate(rms(numbin))
  allocate(avg(numbin))
  allocate(mypos(0:npix_center-1,1))
  allocate(pos(0:npix_center-1,1))
  allocate(pos_out(0:npix-1,1))
  allocate(myrms(numbin,0:npix_center-1))
  allocate(rms_buff(numbin,0:npix_center-1))
  if (myid == root) allocate(totrms(numbin,0:npix_center-1,nummaps))

  ! Loop over pixel positions
  if (myid == root) open(28,file='ring_count.dat')
  if (myid == root) open(29,file='ring_count_pluss2.dat')
  open(85,file=trim(filelist))
  do k = 1, nummaps

     if (myid == root) write(*,*) 'Map = ', k

     ! Read in data
     call int2string(k, real_text)
     read(85,'(a)') mapname
     call read_map_dp(mapname, nside, ordering, nmaps, map)
     if (ordering == 2) call convert_nest2ring(nside, map(:,1))
     map = map * scale

     mycount2    = 0
     mypos       = 0
     myrms       = 0.d0
     do i = 0+myid, npix_center-1, numprocs
        if (mod(i,100) == 0) write(*,*) i, npix_center
        call pix2vec_ring(nside_center, i, vec)
!        call ang2vec(53.d0*pi/180.d0, 105.04d0*pi/180.d0, vec)
!        call ang2pix_ring(nside_center, 53.d0*pi/180.d0, 105.04d0*pi/180.d0, j)
!        call pix2vec_ring(nside_center, j, vec)
        ! Penrose January 2011
!        call ang2vec(122.25*pi/180.d0, 271.96d0*pi/180.d0, vec)
        call compute_rms(nside, map(:,1), mask(:,1), vec, rms, binsize, numbin)
        myrms(:,i) = rms
        mu = sum(rms(minbin:numbin)) / real(numbin-minbin+1,dp)
        if (count(mu-rms(minbin:numbin) > threshold) >= 2) mycount2 = mycount2 + 1
        mypos(i,1) = count(mu-rms(minbin:numbin) > threshold)
        if (sum(rms) == 0) mypos(i,1) = -1.6375d30
        if (output_histograms) then
           call int2string(i,i_text)
           filename = trim(outdir) // '/rms_real' // real_text // '_pix' // i_text // '.dat'
           open(58,file=filename)
           do j = 1, numbin
              write(58,*) (j-1)*binsize, rms(j)
              write(58,*)     j*binsize, rms(j)
           end do
           close(58)
        end if
!        stop
     end do

     ! Collect results from all processors
     call mpi_reduce(mycount2, count2, 1, mpi_integer, mpi_sum, root, MPI_COMM_WORLD, ierr)
     if (myid == root) write(29,*) k, count2

     call mpi_reduce(mypos, pos, size(pos), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
     if (myid == root) then
        filename = trim(outdir) // '/pos_real' // real_text // '.fits'
        call udgrade_ring(pos, nside_center, pos_out, nside)
        where (mask(:,1:1) == 0.d0) pos_out = -1.6375d30
        call write_map(pos_out, 1, filename)
     end if

     call mpi_reduce(myrms, rms_buff, size(myrms), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)

     ! Output mean and rms per pixel
     if (myid == root) then
        totrms(:,:,k) = rms_buff
        call output_mean_and_rms(k, totrms(:,:,k))
     end if

     deallocate(map)

  end do
  close(28)
  deallocate(pos_out)

  if (myid == root) then

     call search_multiple_rings_gp(totrms, 12.d0, 3)
!     call search_multiple_rings(totrms, 12.d0, 3)
!     call compute_wavelet_stat('MH_w1', [0.5d0, -1.d0, 0.5d0], totrms, minbin, nside_center) ! MH filter of width 1
!     call compute_wavelet_stat('MH_w2', [0.25d0, 0.25d0, -0.5d0, -0.5d0, 0.25d0, 0.25d0], totrms, minbin, nside_center) ! MH filter of width 2
!     call compute_wavelet_stat('wedge', [0.5d0, -0.125d0, -0.75d0, -0.125d0, 0.5d0], totrms, minbin, nside_center) ! Wedge filter
!     call compute_chisq_stat(totrms, nside_center)
  end if

  ! Clean up
  deallocate(rms, mask)

  call mpi_finalize(ierr)

contains

  subroutine compute_rms(nside, map, mask, vec, rms, binsize, numbins)
    implicit none

    real(dp),     dimension(0:), intent(in)  :: map
    real(dp),     dimension(0:), intent(in)  :: mask
    real(dp),     dimension(3),  intent(in)  :: vec
    real(dp),     dimension(1:), intent(out) :: rms
    real(dp),                    intent(in)  :: binsize
    integer(i4b),                intent(in)  :: nside, numbins

    integer(i4b)                                    :: i, j, pix, bin, numpix
    real(dp)                                        :: radius, av, bs
    real(dp),     allocatable, dimension(:,:)       :: mask_out
    integer(i4b), allocatable, dimension(:)         :: pixlist, numinring
    integer(i4b), allocatable, dimension(:,:)       :: pixrings
    real(dp),     allocatable, dimension(:,:), save :: pixcenters

    if (.not. allocated(pixcenters)) then
       allocate(pixcenters(0:12*nside**2-1, 3))
       do i = 0, 12*nside**2-1
          call pix2vec_ring(nside, i, pixcenters(i,:)) 
       end do
    end if

    bs     = binsize*pi/180.d0 ! in radians
    radius = bs*numbins
    allocate(pixlist(0:nside**2-1))
    call query_disc(nside, vec, radius, pixlist, numpix)

!    write(*,*) vec
!    if (myid == 0 .and. vec(3) < 0.1d0) then
!       write(*,*) vec
!       allocate(mask_out(0:12*nside**2-1,1))
!       mask_out = -1.6375d30
!       mask_out(pixlist(0:numpix-1),1) = mask(pixlist(0:numpix-1))
!       call write_map('mask.fits', mask_out, 1)
!       call mpi_finalize(ierr)
!       stop
!       write(*,*) sum(mask(pixlist(0:numpix-1))), numpix
!    end if
    if (sum(mask(pixlist(0:numpix-1))) < 0.6 * numpix) then
       rms = 0.d0
       return
    end if

    allocate(pixrings(numbins,numpix))
    allocate(numinring(numbins))
    numinring = 0
    do i = 0, numpix-1
       if (mask(pixlist(i)) > 0) then
          bin = min(floor(acos(sum(vec*pixcenters(pixlist(i),:)))/bs)+1,numbins)
          numinring(bin) = numinring(bin) + 1
          pixrings(bin, numinring(bin)) = pixlist(i)
       end if
    end do

    rms = 0.d0
    do bin = 1, numbins
       if (numinring(bin) < 2) then
          rms(bin) = 0.d0
          cycle
       end if
       av = 0.d0
       do i = 1, numinring(bin)
          pix = pixrings(bin, i)
          av = av + map(pix)*mask(pix)
       end do
       av = av/numinring(bin)
       do i = 1, numinring(bin)
          pix = pixrings(bin, i)
          rms(bin) = rms(bin) + (map(pix)*mask(pix)-av)**2
       end do
       rms(bin) = rms(bin)/(numinring(bin)-1)
       rms(bin) = sqrt(rms(bin))
    end do
  end subroutine compute_rms


  subroutine output_mean_and_rms(k, totrms)
    implicit none

    integer(i4b),                   intent(in) :: k
    real(dp),     dimension(1:,0:), intent(in) :: totrms

    integer(i4b) :: j, n, q, numbin, npix
    real(dp), allocatable, dimension(:) :: avg, rms

    numbin = size(totrms(:,0))
    npix   = size(totrms(1,:))

    allocate(avg(numbin), rms(numbin))

    n = count(totrms(numbin,:) > 0)
    do j = 1, numbin
       n      = 0
       avg(j) = 0.d0
       do q = 0, npix-1
          if (totrms(j,q) > 0.d0) then
             avg(j) = avg(j) + totrms(j,q)
             n      = n + 1
          end if
       end do
       avg(j) = avg(j) / real(n,dp)
           
       rms(j) = 0.d0
       do q = 0, npix-1
          if (totrms(j,q) > 0.d0) then
             rms(j) = rms(j) + (totrms(j,q)-avg(j))**2
          end if
       end do
       rms(j) = sqrt(rms(j) / (n-1))
    end do
    filename = trim(outdir) // '/mean_rms_real' // real_text // '.dat'
    open(unit,file=trim(filename))
    do j = 1, numbin
       write(unit,*) real((j-0.5)*binsize,sp), real(avg(j),sp)
    end do
    write(unit,*) 
    do j = 1, numbin
       write(unit,*) real((j-0.5)*binsize,sp), real(avg(j)-rms(j),sp)
    end do
    do j = numbin, 1, -1
       write(unit,*) real((j-0.5)*binsize,sp), real(avg(j)+rms(j),sp)
    end do
    write(unit,*) 
    do j = 1, numbin
       write(unit,*) real((j-0.5)*binsize,sp), real(avg(j)-2*rms(j),sp)
    end do
    do j = numbin, 1, -1
       write(unit,*) real((j-0.5)*binsize,sp), real(avg(j)+2*rms(j),sp)
    end do
    close(unit)

    deallocate(avg, rms)

  end subroutine output_mean_and_rms

  subroutine compute_chisq_stat(totrms,nside)
    implicit none

    real(dp), dimension(1:,0:,1:), intent(inout) :: totrms
    integer(i4b),                  intent(in)    :: nside

    integer(i4b) :: i, j, k, l, q, p, nabolist(1:8), numnabo
    real(dp),     allocatable, dimension(:,:)       ::  res2

    allocate(pos_out(0:npix_center-1,1))
    allocate(pte(0:npix_center-1,1))
    allocate(cov(numbin,numbin))
    allocate(chisq(nummaps))
    allocate(chisq_out(0:npix_center-1,nummaps))
    allocate(res2(0:npix_center-1,1))
    
    do i = 0, npix_center-1
       do j = 1, nummaps
          totrms(:,i,j) = totrms(:,i,j) - sum(totrms(:,i,j)) / real(numbin,dp)
       end do
    end do
    
    do i = 0, npix_center-1
       if (any(totrms(:,i,:) == 0.d0)) then
          pos_out(i,1) = -1.6375d30
          pte(i,1)     = -1.6375d30
       else
          avg = 0.d0
          do q = 1, numbin
             avg(q) = sum(totrms(q,i,:)) / real(nummaps,dp)
          end do
          do j = 1, numbin
             do k = 1, numbin
                cov(j,k) = sum((totrms(j,i,:)-avg(j))*(totrms(k,i,:)-avg(k)))/real(nummaps-1,dp)
             end do
          end do
          call invert_singular_matrix(cov, 1.d-12)
          do j = 1, nummaps
             chisq(j) = sum((totrms(:,i,j)-avg) * matmul(cov, (totrms(:,i,j)-avg)))
          end do
          chisq_out(i,:) = chisq(:)
          pte(i,1)     = count(chisq(2:nummaps) > chisq(1)) / real(nummaps-1,dp)
       end if
    end do
    
    open(95,file='chisq.dat')
    do j = 1, nummaps
       do k = 0, npix_center-1
          if (chisq_out(k,j) /= 0.d0) then
             write(95,*) real(chisq_out(k,j),sp)
          end if
       end do
    end do
    close(95)
   
 
    call convert_ring2nest(nside, chisq_out) 
    do i = 1, nummaps
       res2 = 0.d0
       do p = 0, npix_center-1
          if (chisq_out(p,i)== 0.d0) chisq_out(p,i) = -1.6375d30
          res2(p,1) = chisq_out(p,i)
       end do
       do p = 0, npix_center-1
          call neighbours_nest(nside, p, nabolist, numnabo)
          do j=1,numnabo
             if (chisq_out(nabolist(j),i)==-1.6375d30) res2(p,1) = -1.6375d30
          end do
       end do
       call convert_nest2ring(nside, res2)

       call int2string(i, real_text)
       filename = trim(outdir) // '/chisq_real' // real_text // '.fits'
       call write_map(res2, 1, filename)
       write(*,*) 'Real = ', i, ', max chi pos = ', maxloc(chisq_out(:,i))-1, real(chisq_out(maxloc(chisq_out(:,i))-1,i),sp), real(maxval(chisq_out(:,i)),sp)
    end do
    call write_map(pte, 1, 'pte_real1.fits')
    deallocate(pos_out)
    
  end subroutine compute_chisq_stat


  subroutine compute_wavelet_stat(prefix, filter, rms, minbin, nside)
    implicit none

    character(len=*),              intent(in)    :: prefix
    real(dp), dimension(1:),       intent(in)    :: filter
    real(dp), dimension(1:,0:,1:), intent(in)    :: rms
    integer(i4b),                  intent(in)    :: minbin, nside

    integer(i4b)          :: numbin, npix, nummaps, n, i, j, p, bin, nabolist(1:8), numnabo
    real(dp)              :: svar
    character(len=4)      :: mapnr
    character(len=512)    :: filename
    real(dp), dimension(2):: pos
    real(dp),     allocatable, dimension(:,:)       :: res, res2
    real(dp),     allocatable, dimension(:)         :: filbins

    numbin  = size(rms(:,0,1))
    npix    = size(rms(1,:,1))
    nummaps = size(rms(1,0,:))
    n       = size(filter)
    
    allocate(res(0:npix-1, 1))
    allocate(res2(0:npix-1, 1))
    allocate(filbins(1:numbin))

    do i = 1, nummaps
       res = 0.d0
       do p =0, npix-1
          filbins = 0.d0
          do bin = minbin, numbin -n +1
             if (rms(bin,p,i)==0.d0) then
                res(p,1)=-1.6375d30
                exit
             else
                do j=1,n
                   filbins(bin)= filbins(bin) + filter(j)*rms(bin+j-1,p,i)
                end do
              end if
           end do
! Enten leite etter maxval
!           res(p,1) = maxval(filbins)
! Eller sum
           where(filbins < 12.d0) filbins = 0.d0
           res(p,1) = sum(filbins)
           if (res(p,1)==0.d0) res(p,1)=-1.6375d30
        end do

          open(39,file='filtered.dat')
          do bin = 1, numbin
             write(39,*) rms(bin,p,i), filbins(bin)
          end do
          close(39)
!          if (maxval(filbins) > 20d0 .and. trim(prefix) == 'wedge') then
!             write(*,*) trim(prefix)
!             stop
!          end if

       end do
       
        call convert_ring2nest(nside, res) 
        do p = 0, npix-1
           res2(p,1) = res(p,1)
           call neighbours_nest(nside, p, nabolist, numnabo)
           do j=1,numnabo
              if (res(nabolist(j),1)==-1.6375d30) res2(p,1) = -1.6375d30
           end do
        end do
!        call convert_nest2ring(nside, res)
        call convert_nest2ring(nside, res2)

       call int2string(i, mapnr)
       filename = trim(outdir) // '/' // trim(prefix) // '_map' // mapnr // '.fits'
       write(*,*) 'nan = ', count(res == -1.6375d30)
       write(*,*) 'nan2 = ', count(res2 == -1.6375d30)
 !      call write_map(filename, res, 1)
 !      filename = trim(outdir) // '/' // trim(prefix) // '2_map' // mapnr // '.fits'
       call write_map(res2, 1, filename)
       write(*,*) 'Writing ', trim(filename)
       pos=maxloc(res)
       write(*,*) 'Maximum in pix no:', pos(1)-1
       write(*,*) ''
    end do

    deallocate(res)

  end subroutine compute_wavelet_stat


  subroutine search_multiple_rings_gp(rms, threshold, maxwidth)
    implicit none

    real(dp),     dimension(1:,0:,1:), intent(in) :: rms
    real(dp),                          intent(in) :: threshold
    integer(i4b),                      intent(in) :: maxwidth

    integer(i4b) :: numbin, npix, nummaps, n
    integer(i4b) :: numring, width
    integer(i4b) :: i, j, p
    logical(lgt) :: ring
    real(dp)     :: mu
    character(len=4)   :: mapnr
    character(len=128) :: filename
    integer(i4b), dimension(2) :: pos
    real(dp), allocatable, dimension(:,:) :: res

    numbin  = size(rms(:,0,1))
    npix    = size(rms(1,:,1))
    nummaps = size(rms(1,0,:))

    allocate(res(0:npix-1,1))

    do i = 1, nummaps
       do p = 0, npix-1
          
          if (any(rms(:,p,i) == 0.d0)) then
             res(p,1) = -1.6375d30
          else
             
             numrings = 0
             ring     = .false.
             mu       = sum(rms(:,p,i)) / numbin
             do j = 5, 32
                if (.not. ring .and. mu-rms(j,p,i) > threshold) then
                   ! New ring begins
                   ring  = .true.
                   width = 1
                end if
                if (ring) then
                   if (mu-rms(j,p,i) < threshold) then
                      ! Existing ring closes
                      ring     = .false.
                      numrings = numrings + 1
                   else
                      width = width + 1
                   end if
                end if
             end do
             res(p,1) = numrings
             
          end if
          
       end do
       
       call int2string(i, mapnr)
       filename = trim(outdir) // '/ringcount_map' // mapnr // '.fits'
       call write_map(res, 1, filename)
       
       pos=maxloc(res)
       write(*,*) 'Maximum number of rings in pix no:', pos(1)-1, res(pos(1)-1,1)
       
    end do
    
    deallocate(res)
    
  end subroutine search_multiple_rings_gp


  subroutine search_multiple_rings(rms, threshold, maxwidth)
    implicit none

    real(dp),     dimension(1:,0:,1:), intent(in) :: rms
    real(dp),                          intent(in) :: threshold
    integer(i4b),                      intent(in) :: maxwidth

    integer(i4b) :: numbin, npix, nummaps, n
    integer(i4b) :: numring, width
    integer(i4b) :: i, j, p
    logical(lgt) :: ring
    character(len=4)   :: mapnr
    character(len=128) :: filename
    integer(i4b), dimension(2) :: pos
    real(dp), allocatable, dimension(:,:) :: res

    numbin  = size(rms(:,0,1))
    npix    = size(rms(1,:,1))
    nummaps = size(rms(1,0,:))

    allocate(res(0:npix-1,1))

    do i = 1, nummaps
       do p = 0, npix-1

          if (any(rms(:,p,i) == 0.d0)) then
             res(p,1) = -1.6375d30
          else

             numrings = 0
             ring     = .false.
             do j = 2, numbin
                if (.not. ring .and. rms(j-1,p,i)-rms(j,p,i) > threshold) then
                   ! New ring begins
                   ring  = .true.
                   width = 1
                end if
                if (ring) then
                   if (width <= maxwidth .and. (rms(j,p,i)-rms(j-1,p,i) > threshold)) then
                      ! Existing ring closes
                      ring     = .false.
                      numrings = numrings + 1
                   else if (width > maxwidth) then
                      ring = .false.
                   else
                      width = width + 1
                   end if
                end if
             end do
             res(p,1) = numrings
             
          end if

       end do

       call int2string(i, mapnr)
       filename = trim(outdir) // '/ringcount_map' // mapnr // '.fits'
       call write_map(res, 1, filename)

       pos=maxloc(res)
       write(*,*) 'Maximum number of rings in pix no:', pos(1)-1, res(pos(1)-1,1)

    end do

    deallocate(res)

  end subroutine search_multiple_rings


  subroutine print_usage

    write(*,*) 'Usage: mpirun -n X rings param.txt'

  end subroutine print_usage

  subroutine read_map_dp(filename, nside, ordering, nmaps, map)
    implicit none

    character(len=*),                           intent(in)   :: filename
    integer(i4b),                               intent(out)  :: nside, ordering, nmaps
    real(dp),           pointer, dimension(:,:)              :: map

    integer(i4b)    :: i, temp
    real(dp)        :: nullval
    logical(lgt)    :: anynull

    temp = getsize_fits(filename, nside=nside, ordering=ordering, nmaps=nmaps)
    allocate(map(0:12*nside**2-1,nmaps))
    call read_bintab(trim(filename), map, 12*nside**2, nmaps, nullval, anynull)
  end subroutine read_map_dp

end program rings
