program sindex
  use healpix_types
  use l2_fileutils
  use pix_tools
  use fitstools
  use quiet_fileutils
  use math_tools
  use alm_tools
  use quiet_postutils
  use quiet_utils

  implicit none

  include "mpif.h"

  integer(i4b)       :: iargc
  integer(i4b)       :: unit, myid, numprocs, ierr, root
  character(len=30)  :: kommando

  ! Initialize MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0
  unit        = 42+myid

  if (iargc() == 0) then
     call give_user_info
  else if (myid == root) then
     write(*,*) '-------------------- QUIET spectral index calculator  --------------------'
     write(*,*)
  end if

  ! Get name of main kommando
  call getarg(1,kommando)

  if (kommando == 'index') then
     write(*,*) 'Find spectral index from two maps'
     call spectral_index(unit)
  else   if (kommando == 'leakage') then
     write(*,*) 'Estimate I->QU leakage parameters'
     call leakage(unit)
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
  ! subroutine leakage
  !-----------------------------------------------------------------------------------------------

  subroutine leakage(unit)
    implicit none
    
    integer(i4b),       intent(in) :: unit
    
    character(len=256)             :: tempfile, listfile, infreq1, infreq2, outprefix, outfile, rmsfile1, rmsfile2
    character(len=256)             :: beamfile1, beamfile2, maskfile, infwhm, innside, inradius, incutradius, insindex
    integer(i4b)                   :: i, j, k, n, comp, pol, nlist, p
    integer(i4b)                   :: nmaps, nside, npix, ordering, ierr
    integer(i4b)                   :: numdiodes
    integer(i4b)                   :: areanside, areanpix, num
    real(dp)                       :: nullval, freq1, freq2, fwhm, a2t, radius, cutradius, freq_corr, specindex
    logical(lgt)                   :: anynull, mismatch, inv
    real(dp)                       :: healnan=-1.6375d30, root=0
  
    real(dp),       pointer,     dimension(:,:)   :: tmap, map2, mask, rms1, rms2!, beam1, beam2
    real(dp),       allocatable, dimension(:)     :: leak, uncert
!    integer(i4b),   allocatable, dimension(:)     :: map2heal, antimap2heal
    character(len=256), pointer, dimension(:)     :: maplist
    character(len=3)                              :: filnum
    real(dp),                    dimension(3)     :: gc_vec
    integer(i4b),   allocatable, dimension(: )    :: listpix

    ! Get parameters
    if (iargc() /= 12) then
       write(*,*) 'leakage takes 11 parameters'
       call give_user_info
    else 
       call getarg(2, tempfile)
       call getarg(3, listfile)
       call getarg(4, beamfile1)
       call getarg(5, beamfile2)
       call getarg(6, infreq1)
       call getarg(7, infreq2)
       call getarg(8, insindex)
       call getarg(9, infwhm)
       call getarg(10, inradius)
       call getarg(11, incutradius)
!       call getarg(10, rmsfile1)
!       call getarg(11, rmsfile2)
       call getarg(12, outprefix)
    end if
    
    read(infreq1,*) freq1
    read(infreq2,*) freq2
    read(infwhm,*)  fwhm
    read(inradius,*)  radius
    read(incutradius,*)  cutradius
    read(insindex,*)  specindex
    if (myid==root) write(*,*) 'Estimating leakage from ', trim(infreq1), ' GHz map and ', trim(infreq2), ' GHz map'
    if (myid==root) write(*,*) 'using area with radius ', trim(inradius), ' degrees.'
    write(*,*) 'freq1', freq1
    write(*,*) 'freq2', freq2
    write(*,*) 'fwhm ', fwhm
    write(*,*) 'radius ', radius
    write(*,*) 'cutradius ', cutradius
    write(*,*) 'spectral index ', specindex
    cutradius = (180-cutradius)*pi/180.d0   ! radius in radians
    gc_vec = [-1,0,0]                       ! because we want to remove what's OUTSIDE the disk
    freq_corr= (freq1/freq2)**specindex        
    write(*,*) 'freq correction ', freq_corr

    a2t = ant2thermo(freq2)
   
    ! Read maps and mask and convert to nest if necassary
    call read_mapfile(tempfile,  npix, nside, nmaps, ordering, tmap,  nest=.true.)
!    call read_mapfile(maskfile,  npix, nside, nmaps, ordering, mask,  nest=.true., check_nside=nside)
!    call read_mapfile(rmsfile1,  npix, nside, nmaps, ordering, rms1,  nest=.true., check_nside=nside)
!    call read_mapfile(rmsfile2,  npix, nside, nmaps, ordering, rms2,  nest=.true., check_nside=nside)

    ! Always use temperature data
    pol = 1
    nmaps = 1

    ! Write maps
    outfile = trim(outprefix) // '_input_tmap.fits'
    call write_map(trim(outfile), tmap, ordering)
    write(*,*) '* Map written to file = ', trim(outfile)

    ! Put missing pixels to zero
    where (abs(-1.6375d30 - tmap) < 1d25) tmap=0.d0
    ! Convolve temperature map with common gaussian beam
    if (myid==root) write(*,fmt='(a,f5.1,a)') 'Smoothing to common gaussian beam size of ',fwhm,' arc minutes.'
    call beam_convolve(tmap, ordering, nside, nmaps, beamfile1, fwhm)
    ! Restore unhit pixels
    where (tmap == 0.d0) tmap = healnan

    ! Write maps
    outfile = trim(outprefix) // '_beam_convolved_tmap.fits'
    call write_map(trim(outfile), tmap, ordering)
    write(*,*) '* Map written to file = ', trim(outfile)

    ! Read filelist
    call read_filelist(unit, listfile, numdiodes, maplist)
    allocate(leak(numdiodes))
    allocate(uncert(numdiodes))
    allocate(listpix(0:npix-1))
    ! Loop over files
    do i = 1+myid, numdiodes, numprocs
       write(*,*) 'Processing file no.', i, 'out of',numdiodes
       call read_mapfile(trim(maplist(i)), npix, nside, nmaps, ordering, map2, nest=.true., check_nside=nside)
       where (abs(healnan - map2) > 1d25) map2=map2*a2t*1000.d0
       ! Write map
       call int2string(i, filnum)
       outfile = trim(outprefix) // '_input_map' // filnum // '.fits'
       call write_map(trim(outfile), map2, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Find pixels outside a cutradius deg radius around gc center
       call query_disc(nside, gc_vec, cutradius, listpix, nlist, nest=1)
       do p=0,nlist-1
          map2(listpix(p), :)=0.d0
       end do
        ! Write map
       call int2string(i, filnum)
       outfile = trim(outprefix) // '_input_map' // filnum // '_cut.fits'
       call write_map(trim(outfile), map2, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Beam convolve
       where (abs(healnan - map2) < 1d25) map2=0.d0
       call beam_convolve(map2, ordering, nside, 1, beamfile2, fwhm)
       where (map2 == 0.d0) map2 = healnan
       ! Write map
       outfile = trim(outprefix) // '_beam_convolved_map' // filnum // '.fits'
       call write_map(trim(outfile), map2, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
       call estimate_leakage(tmap(:,1), map2(:,1), nside, radius, leak(i), uncert(i))
       write(*,*) i, 'leakage=', leak(i)*100.d0
       write(*,*)
    end do

    ! Write leakage coeffisients to file
    outfile = trim(outprefix) // '_leakage.txt'
    open(unit, file=trim(outfile))
    do i = 1, numdiodes
       write(unit,*) 'i2qu_leak  0  0  -1  -1 ', (i-1)/4, modulo(i-1,4), leak(i)*100.d0*freq_corr
       write(*,*) 'i2qu_leak  0  0  -1  -1 ', (i-1)/4, modulo(i-1,4), leak(i)*100.d0*freq_corr
    end do    
    close(unit)
    write(*,*) '* Leakage coeffisients written to file = ', trim(outfile)

    ! Clean up
    deallocate(leak)
    deallocate(uncert)
    deallocate(listpix)

  end subroutine leakage

  !-----------------------------------------------------------------------------------------------
  ! subroutine spectral_index
  !-----------------------------------------------------------------------------------------------

  subroutine spectral_index(unit)
    implicit none
    
    integer(i4b),       intent(in) :: unit
    
    character(len=256)             :: infile1, infile2, infreq1, infreq2, outprefix, outfile, rmsfile1, rmsfile2
    character(len=256)             :: beamfile1, beamfile2, maskfile, infwhm, innside, inuncert
    integer(i4b)                   :: i, j, k, n, p, comp, pol, maskcomp, lstop, nband=2
    integer(i4b)                   :: nmaps, nside, npix, ordering, ierr
    integer(i4b)                   :: nmaps_mask, nmaps_rms
    integer(i4b)                   :: areanside, areanpix, num, max_nside, ns, ns_max, max_npix
    real(dp)                       :: nullval, freq1, freq2, fwhm, uncertlimit
    logical(lgt)                   :: anynull, mismatch, inv, multi
    real(dp)                       :: healnan=-1.6375d30, root=0
  
    real(dp),     pointer,     dimension(:,:)   :: map1, map2, rms1, rms2
    real(dp),     allocatable, dimension(:,:)   :: outmap, areamap, uncert, multimap, multiuncert
    real(dp),     allocatable, dimension(:,:)   :: samlekart, samlekart2, samleuncert, samleuncert2 
    real(dp),     allocatable, dimension(:,:)   :: offset, fg, samleoffset, samlefg
    real(dp),     allocatable, dimension(:)     :: ratio, kart1, kart2, var1, var2
    real(dp),     allocatable, dimension(:,:,:) :: kart
    integer(i4b), allocatable, dimension(:)     :: map2heal, antimap2heal, locmask
    integer(i4b), pointer,     dimension(:)     :: mask


    ! Get parameters
    if (iargc() /= 13 .and. iargc() /= 14) then
       write(*,*) 'index takes 12(13) parameters'
       call give_user_info
    else 
       call getarg(2, infile1)
       call getarg(3, infile2)
       call getarg(4, beamfile1)
       call getarg(5, beamfile2)
       call getarg(6, infreq1)
       call getarg(7, infreq2)
       call getarg(8, infwhm)
       call getarg(9, maskfile)
       call getarg(10, rmsfile1)
       call getarg(11, rmsfile2)
       call getarg(12, innside)
       call getarg(13, outprefix)
       if (iargc()==14) then
          multi = .true.
          call getarg(14, inuncert)
          read(inuncert,*) uncertlimit
          write(*,*) 'Multiresolution: uncertlimit =', uncertlimit
       else
          multi = .false.
       end if
    end if

    read(infreq1,*) freq1
    read(infreq2,*) freq2
    read(infwhm,*)  fwhm
    read(innside,*) max_nside
    ns_max = nint(log(real(max_nside,dp))/log(2.d0))+1
    if (max_nside==0) ns_max=0
    if (2**(ns_max-1) /= max_nside) then
       write(*,*) 'nside =',max_nside,' not a factor of 2. Quitng'
       stop
    end if
    if (myid==root) write(*,*) 'Finding spectral index from ', trim(infreq1), ' GHz map and ', trim(infreq2), ' GHz map'
    if (myid==0) then 
       write(*,*) 'freq1', freq1
       write(*,*) 'freq2', freq2
       write(*,*) 'fwhm', fwhm
       write(*,*) 'nside large pixels', max_nside
    end if
  
    ! Read maps and convert to nest if necassary
    call read_mapfile(infile1,  npix, nside, nmaps, ordering, map1, nest=.true.)
    call read_mapfile(infile2,  npix, nside, nmaps, ordering, map2, nest=.true., check_nside=nside, check_nmaps=nmaps)
    ! Checking polariation
    if (nmaps==1) then
       pol = 1
       if (myid==0) write(*,*) 'Running at temperature data' 
    else if (nmaps==3) then
       pol = 2
       if (myid==0) write(*,*) 'Running at polarisation data' 
    else
       write(*,*) nmaps, '= nmaps. Unknown number. Quiting'
    end if
    ! Read mask and convert to nest if necassary
    call read_and_check_mask(myid, maskfile, mask, nmaps, ordering, npix, pol)
    ! Read rms maps and convert to nest if necassary
    call read_mapfile(rmsfile1, npix, nside, nmaps_rms,  ordering, rms1, nest=.true., check_nside=nside)
    if (nmaps_rms<nmaps) then
       if (myid==root) write(*,*) 'We dont have smoothed rms for right component. Quiting'
       call mpi_finalize(ierr)
       stop
    end if    
    call read_mapfile(rmsfile2, npix, nside, nmaps_rms,  ordering, rms2, nest=.true., check_nside=nside)
    if (nmaps_rms<nmaps) then
       if (myid==root) write(*,*) 'We dont have smoothed rms for right component. Quiting'
       call mpi_finalize(ierr)
       stop
    end if    

    ! Write input maps
    if (myid==0) then 
       outfile = trim(outprefix) // '_input_map1.fits'
       call write_map(trim(outfile), map1, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_input_map2.fits'
       call write_map(trim(outfile), map2, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! Check for common pixels
    n = count(healok(map1(:, nmaps)) .and. healok(map2(:, nmaps)) .and. mask == 1)
    if (myid==root) write(*, *) 'The two input maps have ', n, 'common pixels'
    if (n==0) return

    ! Save map2heal
    j=0
    k=0
    allocate(map2heal(n))
    allocate(antimap2heal(npix-n))
    if (pol == 1) then
       do i = 0, npix-1
          if (healok(map1(i, nmaps)) .and. healok(map2(i, nmaps)) .and. mask(i)==1) then
             j = j+1
             map2heal(j) = i
          else
             k = k+1
             antimap2heal(k) = i
             mask(i) = 0
          end if
       end do
    else if (pol==2) then
       do i = 0, npix-1
          if (healok(map1(i, 2)) .and. healok(map2(i, 2)) .and. healok(map1(i, 3)) .and. healok(map2(i, 3)) .and. mask(i)==1) then     !!!!OBS only for pol
             j = j+1
             map2heal(j) = i
          else
             k = k+1
             antimap2heal(k) = i
             mask(i) = 0
          end if
       end do
    end if
    ! Put missing pixels to zero
    do j = 1, nmaps
       do i = 0, npix-1
          if (healnot(map1(i, j))) map1(i, j) = 0.d0
          if (healnot(map2(i, j))) map2(i, j) = 0.d0
       end do
    end do

!!$    ! Write maps
!!$    outfile = trim(outprefix) // '_input_map1_after.fits'
!!$    call write_map(trim(outfile), map1, ordering)
!!$    write(*,*) '* Map written to file = ', trim(outfile)
!!$    outfile = trim(outprefix) // '_input_map2_after.fits'
!!$    call write_map(trim(outfile), map2, ordering)
!!$    write(*,*) '* Map written to file = ', trim(outfile)

    ! Convolve both maps with common gaussian beam (skipping if filename is 'none')
    if (beamfile1 =='none') then
       if (myid==root) write(*,fmt='(a)') 'NB: No smoothing of beams since beam filename1 is none.'
    else
       if (myid==root) write(*,fmt='(a,f5.1,a)') 'Smoothing to common gaussian beam size of ',fwhm,' arc minutes.'
       call beam_convolve(map1, ordering, nside, nmaps, beamfile1, fwhm, onlypol=.true.)
       call beam_convolve(map2, ordering, nside, nmaps, beamfile2, fwhm, onlypol=.true.)
    end if

    ! Restoring missing pixels
    map1(antimap2heal,:) = healnan
    map2(antimap2heal,:) = healnan
    deallocate(antimap2heal)

    ! Write beam_convolved maps
    if (myid==0) then 
       outfile = trim(outprefix) // '_beam_conv_map1.fits'
       call write_map(trim(outfile), map1, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_beam_conv_map2.fits'
       call write_map(trim(outfile), map2, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
    end if

    ! loop over nsides of large pixels
    do ns = 0, 0!ns_max
       if (all(mask==0)) exit
       if (allocated(kart1)) deallocate(kart1)
       if (allocated(kart2)) deallocate(kart2)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(areamap)) deallocate(areamap)
       if (allocated(uncert)) deallocate(uncert)
       ! Divide into healpix areas
       areanside=max_nside/(2**ns)
       if (myid==root) write(*,fmt='(a,i3,a)') 'Scatterplot: Choosing regions to be nside =',areanside,' healpix pixels.'
       areanpix = 12*areanside**2
if (max_nside==0) areanpix=1
       num = (npix/areanpix)   ! pixels per large healpix
!       num = (nside/areanside)**2   ! pixels per large healpix
       allocate(areamap(0:areanpix-1, nmaps))
       allocate(uncert(0:areanpix-1, nmaps))
       areamap = healnan
       uncert  = healnan
       allocate(kart1(pol*num))
       allocate(kart2(pol*num))
       allocate(locmask(num))
       do i = 0+myid, areanpix-1, numprocs
          locmask  = mask(i*num:(i+1)*num-1)
          if (pol==1) then 
             kart1 = map1(i*num:(i+1)*num-1, 1)
             kart2 = map2(i*num:(i+1)*num-1, 1)
          else if (pol==2) then
             kart1(1:num)       = map1(i*num:(i+1)*num-1, 2)
             kart1(num+1:2*num) = map1(i*num:(i+1)*num-1, 3)
             kart2(1:num)       = map2(i*num:(i+1)*num-1, 2)
             kart2(num+1:2*num) = map2(i*num:(i+1)*num-1, 3)
          end if
!          call scatterplot(kart1, kart2, locmask, pol, i, areamap(i, :), freq1, freq2, uncert(i, :) )
          mask(i*num:(i+1)*num-1) = locmask
       end do
    end do
    ! Collect to one root map
    allocate(samlekart(0:areanpix-1, nmaps))
    allocate(samleuncert(0:areanpix-1, nmaps))
    where(areamap == healnan) areamap = 0.d0
    where(uncert == healnan)  uncert  = 0.d0
    call mpi_reduce(areamap, samlekart,   size(areamap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(uncert,  samleuncert, size(areamap), mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
write(*,*) 'A3 scatter areamap samlekart', size(areamap), size(samlekart)

    ! Write spectral index map and uncertainty to file
    if (myid==0) then 
       where(samlekart == 0.d0)   samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       outfile = trim(outprefix) // '_spectral_index_from_scatter.fits'
       call write_map(trim(outfile), samlekart, ordering)
       write(*,*) '* Spectral index map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_spectral_index_uncertainty_from_scatter.fits'
       call write_map(trim(outfile), samleuncert, ordering)
       write(*,*) '* Uncertainty in spectral index written to file = ', trim(outfile)
    end if
    
    ! Clean up
    if (allocated(kart1)) deallocate(kart1)
    if (allocated(kart2)) deallocate(kart2)
    if (allocated(areamap)) deallocate(areamap)
    if (allocated(uncert)) deallocate(uncert)
    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samleuncert)) deallocate(samleuncert)

    ! Do the same for maxlike
    ! loop over nsides of large pixels
    if (max_nside==0) then
       max_npix=1
    else
       max_npix=12*max_nside**2
    end if
    allocate(multimap(0:max_npix-1, nmaps))
    allocate(samlekart(0:max_npix-1, nmaps))
    allocate(samlekart2(0:max_npix-1, nmaps))
    allocate(multiuncert(0:max_npix-1, nmaps))
    allocate(samleuncert(0:max_npix-1, nmaps))
    allocate(samleuncert2(0:max_npix-1, nmaps))
    multimap=healnan
    samlekart2=healnan
    multiuncert=healnan
    samleuncert2=healnan
    ! loop over nsides of large pixels
    if (multi) then
       lstop = ns_max
    else
       lstop = 0
    end if
    do ns = 0, lstop
       if (all(mask==0)) exit
       if (allocated(kart)) deallocate(kart)
       if (allocated(var1)) deallocate(var1)
       if (allocated(var2)) deallocate(var2)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(areamap)) deallocate(areamap)
       if (allocated(uncert)) deallocate(uncert)
       if (allocated(offset)) deallocate(offset)
       if (allocated(fg)) deallocate(fg)
       ! Divide into healpix areas
       areanside=max_nside/(2**ns)
       if (myid==root) write(*,fmt='(a,i3,a)') 'Maxlike: Choosing regions to be nside =',areanside,' healpix pixels.'
       areanpix = 12*areanside**2
       if (areanside==0) areanpix=1
       num = (npix/areanpix)   ! pixels per large healpix
       allocate(areamap(0:areanpix-1, nmaps))
       allocate(uncert(0:areanpix-1, nmaps))
       allocate(offset(0:areanpix-1, nmaps))
       allocate(fg(0:npix-1, nmaps))
       allocate(kart(num, nmaps, nband))
       allocate(var1(pol*num))
       allocate(var2(pol*num))
       allocate(locmask(num))
       areamap = healnan
       uncert  = healnan
       offset  = healnan
       uncert  = healnan
       fg      = healnan
       multimap = healnan
       samlekart = healnan
       do i = 0+myid, areanpix-1, numprocs
          locmask  = mask(i*num:(i+1)*num-1)
          kart(:,:,1) = map1(i*num:(i+1)*num-1,:)
          kart(:,:,2) = map2(i*num:(i+1)*num-1,:)
          if (pol==1) then 
             var1  = rms1(i*num:(i+1)*num-1, 1)
             var2  = rms2(i*num:(i+1)*num-1, 1)
          else if (pol==2) then
             var1(1:num)        = rms1(i*num:(i+1)*num-1, 2)
             var1(num+1:2*num)  = rms1(i*num:(i+1)*num-1, 3)
             var2(1:num)        = rms2(i*num:(i+1)*num-1, 2)
             var2(num+1:2*num)  = rms2(i*num:(i+1)*num-1, 3)
          end if
          call maxlikeplot(kart, var1, var2, locmask, pol, i, freq1, freq2, areamap(i,:), offset(i,:), uncert(i,:), fg(i*num:(i+1)*num-1, :))
          if (myid==root .and. max_nside==0) write(*,*) 'Maxlike: spectral index =',areamap(i,1),' Uncert =',uncert(i,1)
          if (multi) then
             do j = i*(4**ns), i*(4**ns) + 4**ns -1
                if (healnot(multimap(j,1)) .and. abs(uncert(i,1))<uncertlimit) then 
                   multimap(j,:) = areamap(i,:)
                   multiuncert(j,:) = uncert(i,:)
                end if
             end do
             mask(i*num:(i+1)*num-1) = locmask 
          end if
       end do
       where(multimap == healnan) multimap = 0.d0
       where(multiuncert == healnan) multiuncert = 0.d0
       call mpi_allreduce(multimap,samlekart, size(multimap),mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
       call mpi_allreduce(multiuncert,samleuncert, size(multiuncert),mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
       where(samlekart   == 0.d0) samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       where(samlekart2 == healnan) samlekart2 = samlekart         
       where(samleuncert2 == healnan) samleuncert2 = samleuncert
    end do

    ! Collect to one root map
    allocate(samleoffset(0:areanpix-1, nmaps))
    allocate(samlefg(0:npix-1, nmaps))
    where(areamap == healnan) areamap = 0.d0
    where(uncert == healnan)  uncert  = 0.d0
    where(offset == healnan)  offset  = 0.d0
    where(fg == healnan)      fg  = 0.d0
    if (.not. multi) then
       call mpi_reduce(areamap,samlekart2,   size(areamap),mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
       call mpi_reduce(uncert,  samleuncert2, size(uncert),mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    end if

    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samleuncert)) deallocate(samleuncert)
    allocate(samlekart(0:npix-1, nmaps))
    allocate(samleuncert(0:npix-1, nmaps))
    samlekart = healnan
    samleuncert=healnan
    if (ns_max == 0) then
       ns = nint(log(real(nside,dp))/log(2.d0)) +3
    else
       ns = nint(log(real(nside,dp))/log(2.d0)) - nint(log(real(max_nside,dp))/log(2.d0))  
    end if
    ns = 4**ns
    do k = 1, nmaps
       do i = 0, npix-1
          if (mask(i)==1) then
             samlekart(i,k) = samlekart2(i/ns, k)
             samleuncert(i,k) = samleuncert2(i/ns, k)
          end if
       end do
    end do

    call mpi_reduce(offset,  samleoffset, size(offset),  mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(fg,      samlefg,     size(fg),      mpi_double_precision, mpi_sum, root, MPI_COMM_WORLD, ierr)


    ! Write maxlike spectral index map and uncertainty and offset to file
    if (myid==0) then 
       where(samlekart   == 0.d0) samlekart   = healnan
       where(samleuncert == 0.d0) samleuncert = healnan
       where(samleoffset == 0.d0) samleoffset = healnan
       where(samlefg     == 0.d0) samlefg = healnan
       outfile = trim(outprefix) // '_spectral_index_from_maxlike.fits'
       call write_map(trim(outfile), samlekart, ordering)
       write(*,*) '* Spectral index map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_offset.fits'
       call write_map(trim(outfile), samleoffset, ordering)
       write(*,*) '* Offset map written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_uncert.fits'
       call write_map(trim(outfile), samleuncert, ordering)
       write(*,*) '* Spectral index uncert written to file = ', trim(outfile)
       outfile = trim(outprefix) // '_fg.fits'
       call write_map(trim(outfile), samlefg, ordering)
       write(*,*) '* foreground written to file = ', trim(outfile)
       ! Write histogram to file
       outfile = trim(outprefix) // '_histogram.txt'
       open(unit, file=trim(outfile))
       do i = 0, areanpix-1
          if (healok(samleuncert(i,1))) write(unit, *) (areamap(i,1)-(-2.7))/samleuncert(i,1)
       end do
       do j = 1, nmaps
          do i = 0, areanpix-1
             if (healok(samleuncert(i,j))) then 
                offset(i,j) = (samlekart(i,j)-(-2.7))/samleuncert(i,j)
             else
                offset(i,j) = healnan
             end if
          end do
       end do
       outfile = trim(outprefix) // '_delta.fits'
       call write_map(trim(outfile), offset, ordering)
       write(*,*) '* Delta spectral index written to file = ', trim(outfile)
    end if

    ! Clean up
    if (allocated(kart1)) deallocate(kart1)
    if (allocated(kart2)) deallocate(kart2)
    if (allocated(var1)) deallocate(var1)
    if (allocated(var2)) deallocate(var2)
    if (allocated(locmask)) deallocate(locmask)
    if (allocated(areamap)) deallocate(areamap)
    if (allocated(uncert)) deallocate(uncert)
    if (allocated(offset)) deallocate(offset)
    if (allocated(fg)) deallocate(fg)
    if (allocated(samlekart)) deallocate(samlekart)
    if (allocated(samlekart2)) deallocate(samlekart2)
    if (allocated(samleuncert)) deallocate(samleuncert)
    if (allocated(samleuncert2)) deallocate(samleuncert2)
    if (allocated(samleoffset)) deallocate(samleoffset)
    if (allocated(samlefg)) deallocate(samlefg)
    if (allocated(multimap)) deallocate(multimap)
    if (allocated(multiuncert)) deallocate(multiuncert)


    ! Extract maps
    allocate(kart1(n))
    allocate(kart2(n))
    if (pol==1) then 
       kart1 = map1(map2heal, 1)
       kart2 = map2(map2heal, 1)
    else if (pol==2) then
       kart1 = sqrt(map1(map2heal, 2)**2 + map1(map2heal,3)**2)
       kart2 = sqrt(map2(map2heal, 2)**2 + map2(map2heal,3)**2)
    end if
 
    if (myid==0) then 
       ! Write scatter plot to file
       outfile = trim(outprefix) // '_scatter_plot.txt'
       open(unit, file=trim(outfile))
       if (pol==1) then 
          do i = 1, n
             write(unit,*) kart2(i), kart1(i)
          end do
       else if (pol==2) then
          do i = 1, n
             write(unit,*) map2(map2heal(i), 2),  map1(map2heal(i), 2)
          end do
          write(unit,*)
          do i = 1, n
             write(unit,*) map2(map2heal(i), 3),  map1(map2heal(i), 3)
          end do
       end if
       close(unit)
       write(*,*) 'Scatter plot written to file = ', trim(outfile)
    end if
    ! Starting to clean up
    deallocate(map1)
    deallocate(map2)
    
    ! Calculate naive ratio
    if (myid==0) then 
       allocate(ratio(n))
       ratio = kart1/kart2
       write(*,*) log(sum(ratio)/n)/log(freq1/freq2), '= spectral index from average ratio'
       ratio = log(ratio)/log(freq1/freq2)
       write(*,*) sum(ratio)/n, '= average spectral index'
       ! Writing spectral index map to file
       allocate(outmap(0:npix-1,nmaps))
       outmap = healnan
       outmap(map2heal,1) = ratio
       outfile = trim(outprefix) // '_spectral_index_from_ratio.fits'
       call write_map(trim(outfile), outmap, ordering)
       write(*,*) '* Map written to file = ', trim(outfile)
       ! Writing ratio-based mask to file
       outmap = 0.d0
       do i = 1, n
          if (ratio(i)>-3.7d0 .and. ratio(i)<-2.3d0) then
             outmap(map2heal(i),:) = 1.d0
          end if
       end do
       outfile = trim(outprefix) // '_mask_from_ratio.fits'
       call write_map(trim(outfile), outmap, ordering)
       write(*,*) '* Mask from ratio written to file = ', trim(outfile)
       deallocate(ratio)    
       deallocate(outmap)
    end if

    ! Clean up
    deallocate(map2heal)
    deallocate(kart1)
    deallocate(kart2)
    deallocate(mask)
    deallocate(rms1)
    deallocate(rms2)

  end subroutine spectral_index

  !---------------------------------------------------------------------
  ! Find spectral index from real plot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine maxlikeplot(map, sig1, sig2, mask, pol, index, freq1, freq2, spectral, offset, uncert, fg)

    real(dp),        dimension(:,:,:), intent(in)  :: map
    real(dp),        dimension(:),     intent(in)  :: sig1, sig2
    integer(i4b),    dimension(:),     intent(in)  :: mask
    integer(i4b),                      intent(in)  :: pol, index
    real(dp),        dimension(:),     intent(out) :: spectral, offset, uncert
    real(dp),        dimension(:,:),   intent(out) :: fg
    real(dp),                          intent(in)  :: freq1, freq2

    integer(i4b)                                 :: i, j, k, n, p, v, unit, num   
    integer(i4b)                                 :: nfreq, ngrid   
    real(dp)                                     :: healnan=-1.6375d30, tall, nevn, beta
    real(dp)                                     :: betamin, betamax, offmin, offmax, freq0
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:,:) :: polmap

!    if (index/=4677 .and. index/=4701 .and. index/=4697 .and. index/=4344 .and. index/=4333) then
!    if (index/=4701 .and. index/=4524) then
!    if (index/=4865) then
!       spectral = healnan
!       return
!    end if

    num = size(map(:,1,1))
    ! Initialize output
    spectral   = healnan
    offset     = healnan
    uncert     = healnan

    ! Find num of existing pixels for given area
    n = 0
    do i = 1, num
       if ( healok(map(i,pol,1)) .and. healok(map(i,pol,2)) ) then
          n = n+1
       end if
    end do
    ! Return if too few pixels
!    if (n<num/2) return
    if (n<1) return
    write(*,*) n,'=n', num,'=num' , index,'=hpindex'

    ! else proceed
    if (pol==1) then
       call maxlike(map(:,1,1), map(:,1,2), sig1(:), sig2(:), freq1, freq2, spectral(1), offset(1), uncert(1), n, fg=fg(:,1))!,hpindex=index)
    else if (pol==2) then
       do k=2,3
          write(*,*) k
          call maxlike(map(:,k,1), map(:,k,2), sig1(1+num*(k-2):num*(k-1)), sig2(1+num*(k-2):num*(k-1)), freq1, freq2, spectral(k), offset(k), uncert(k), n, fg=fg(:,k), polbeta=spectral(1), poloff=offset(1))
       end do
       allocate(polmap(pol*num,2))
       polmap(1:num,:)       = map(:,2,:)
       polmap(num+1:2*num,:) = map(:,3,:)
       call maxlike(polmap(:,1), polmap(:,2), sig1(:), sig2(:), freq1, freq2, spectral(1), offset(1), uncert(1), n*pol)!, hpindex=index)
       deallocate(polmap)
    end if

  end subroutine maxlikeplot

  !---------------------------------------------------------------------
  ! Find spectral index from maximum likelihood
  !----------------------------------------------------------------------

  subroutine maxlike(map1, map2, sig1, sig2, freq1, freq2, spectral, offset, uncert, n, fg, polbeta, poloff, hpindex)

    real(dp),               dimension(:), intent(in)  :: map1, map2, sig1, sig2
    real(dp),                             intent(out) :: spectral, offset, uncert
    real(dp),                             intent(in)  :: freq1, freq2
    integer(i4b),                         intent(in)  :: n
    integer(i4b), optional,               intent(in)  :: hpindex
    real(dp),     optional, dimension(:), intent(out) :: fg
    real(dp),     optional,               intent(in)  :: polbeta, poloff

    integer(i4b)                                 :: i, j, p, v, unit, num, int, k
    integer(i4b)                                 :: nfreq, ngrid   
    real(dp)                                     :: healnan=-1.6375d30, tall, beta, finebeta, fineoffset
    real(dp)                                     :: betamin, betamax, offmin, offmax, freq0
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    integer(i4b),    allocatable, dimension(:)   :: in2red
    real(dp),        allocatable, dimension(:)   :: pos, m, freq, A, a2t, prob, nevn
    real(dp),        allocatable, dimension(:,:) :: grid, y, sigsq
    logical(lgt)  :: chatty

    chatty = .false.

    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if
    num = size(map1)
    if (chatty) write(*,*) num,'=num', n ! size(fg),'= size(fg)', n

    ! Put data into place
    nfreq=2
    allocate(freq(nfreq))
    freq(1) = freq1
    freq(2) = freq2
    freq0 = 23.d0
    allocate(a2t(nfreq))
    do v = 1, nfreq
       a2t(v)=ant2thermo(freq(v))
    end do
    allocate(m(nfreq))
    m(1)=0.d0 
    allocate(y(n,nfreq))
    allocate(sigsq(n,nfreq))
    allocate(in2red(n))
    j=0
    do i = 1, num
       if ( healok(map1(i)) .and. healok(map2(i))) then
          j = j+1
          y(j,1) = map1(i)
          y(j,2) = map2(i)
          sigsq(j,1) = sig1(i)**2
          sigsq(j,2) = sig2(i)**2
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
    allocate(nevn(n))
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
             A = 0.d0
             nevn = 0.d0
             do v = 1, nfreq
                A(:) = A(:) + (y(:,v)-m(v))*a2t(v)/sigsq(:,v)*(freq(v)/freq0)**beta
                nevn(:) = nevn(:) + ((a2t(v)**2)*(freq(v)/freq0)**(2.d0*beta))/sigsq(:,v)
             end do
             A = A/nevn           
             do v = 1, nfreq
                do p = 1, n ,16
                   grid(i,j) = grid(i,j) +((y(p,v)-A(p)*a2t(v)*(freq(v)/freq0)**beta-m(v))**2)/sigsq(p,v)
                end do
             end do
          end do
       end do
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
          A = 0.d0
          nevn = 0.d0
          do v = 1, nfreq
             A(:) = A(:) + (y(:,v)-m(v))*a2t(v)/sigsq(:,v)*(freq(v)/freq0)**beta
             nevn(:) = nevn(:) + ((a2t(v)**2)*(freq(v)/freq0)**(2.d0*beta))/sigsq(:,v)
          end do
          A = A/nevn           
          A = A*ant2thermo(freq0)
          fg(in2red)=A(1:size(fg))
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
    deallocate(nevn)
    deallocate(pos)
    deallocate(grid)
    deallocate(y)
    deallocate(sigsq)
    deallocate(m)
    deallocate(freq)
    deallocate(a2t)
 
  end subroutine maxlike

  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scatterplot(map1, map2, mask, pol, index, spectral, freq1, freq2, uncert)

    real(dp),          dimension(:), intent(in)  :: map1, map2
    integer(i4b),      dimension(:), intent(in)  :: mask
    integer(i4b),                    intent(in)  :: pol, index
    real(dp),          dimension(:), intent(out) :: spectral, uncert
    real(dp),                        intent(in)  :: freq1, freq2

    integer(i4b)                                 :: lmax, n, l, i, j, k, unit, num   
    real(dp)                                     :: healnan=-1.6375d30, tall
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:)   :: y, beta, vector
    real(dp),        allocatable, dimension(:,:) :: x, matrix2, h
    real(dp),                     dimension(3)   :: slope
    integer(i4b),    allocatable, dimension(:)   :: polmask

!   if (index/=4677 .and. index/=4701 .and. index/=4697 .and. index/=4344) then
!   if (index/=4865) then
!       spectral = healnan
!       return
!    end if

    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if

    num = size(map1)/pol

    ! Find num of existing pixels for given area
    n = 0
    do i = 1, num
!       if ( map1(i)/=healnot .and. map2(i)/=healnot .and. map1(i)>200.d0) then
       if ( healok(map1(i)) .and. healok(map2(i))) then
          n = n+1
       end if
    end do

    ! Return if too few pixels
    if (n<num/2) then
       spectral = healnan
       return
    end if

    ! Calculate scale factor for Q and U separately
    if (pol==2) then
       call scattercalc(map1(1:num), map2(1:num), mask, n, index, spectral(2), freq1, freq2)
       call scattercalc(map1(num+1:2*num), map2(num+1:2*num), mask, n, index, spectral(3), freq1, freq2)
       allocate(polmask(num*pol))
       polmask(1:num)       = mask
       polmask(num+1:2*num) = mask
       ! Calculate scale factor for T or P separately
       n=n*pol
       call scattercalc(map1, map2, polmask, n, index, spectral(1), freq1, freq2)
       deallocate(polmask)    
    else if (pol==1) then
       call scattercalc(map1, map2, mask, n, index, spectral(1), freq1, freq2)
    end if
      
        return


   ! Write scatter plot to file
    call int2string(index, filnum)
    filename='scatter_'//filnum//'.txt'
    open(42, file=trim(filename))
    do i = 1, size(map1)/2
       if ( healok(map1(i)) .and. healok(map2(i))) write(42,*) map2(i), map1(i)
    end do
    write(42,*)
    do i = size(map1)/2 +1, size(map1)
       if ( healok(map1(i)) .and. healok(map2(i))) write(42,*) map2(i), map1(i)
    end do
    write(42,*)
    write(42,*) 0.d0, beta(2)
    write(42,*) -100.d0, -100*beta(1) +beta(2)
    write(42,*) 150.d0, 150*beta(1) +beta(2)
    close(42)
    write(*,*) 'written to file = ', trim(filename)


     ! Write scatter plot to file
    call int2string(index, filnum)
    filename='scatter_'//filnum//'.txt'
    open(42, file=trim(filename))
    do i = 1, n
       write(42,*) y(i), x(i,1)
    end do
    close(42)
    write(*,*) 'written to file = ', trim(filename)

  end subroutine scatterplot

  !---------------------------------------------------------------------
  ! Find spectral index from scatterplot of map1 vs map2
  !----------------------------------------------------------------------

  subroutine scattercalc(map1, map2, mask, n, index, spectral, freq1, freq2, uncert)

    real(dp),          dimension(:), intent(in)  :: map1, map2
    integer(i4b),      dimension(:), intent(in)  :: mask
    integer(i4b),                    intent(in)  :: n, index
    real(dp),                        intent(out) :: spectral
    real(dp),                        intent(in)  :: freq1, freq2
    real(dp),          optional,     intent(out) :: uncert

    integer(i4b)                                 :: i, j, p
    real(dp)                                     :: healnan=-1.6375d30, tall
    character(len=5)                             :: filnum
    character(len=256)                           :: filename 
    real(dp),        allocatable, dimension(:)   :: y, beta, vector, kmap, qmap
    real(dp),        allocatable, dimension(:,:) :: x, matrix2, h
    real(dp)                                     :: slope(3), vec(2)
   
    allocate(y(n))
    allocate(kmap(n))
    allocate(qmap(n))
    allocate(beta(2))
    allocate(x(n,2))
    allocate(matrix2(2,2))
    p=0
    do i = 1, size(map1)
       !       if ( map1(i)/=healnot .and. map2(i)/=healnot .and. map1(i)>200.d0) then
       if ( healok(map1(i)) .and. healok(map2(i)) .and. mask(i)==1.d0) then
          p = p+1
          kmap(p) = map1(i)
          qmap(p) = map2(i)
       end if
    end do

    matrix2(1,1) = p*1.d0
    matrix2(1,2) = -sum(qmap)
    matrix2(2,1) = matrix2(1,2)
    matrix2(2,2) = sum(qmap*qmap)
    matrix2 = matrix2/(p*sum(qmap*qmap)-(sum(qmap))**2)

    vec(1) = sum(qmap*kmap)
    vec(2) = sum(kmap)

    beta = matmul(matrix2, vec) 
    spectral = log(beta(1))/log(freq1/freq2)
    slope(1)    = beta(1)
 !   write(*,*) beta(1), spectral, index, n
 !   write(*,*)
    if (present(uncert)) then
       ! Uncertainty
       allocate(h(n,n))
       allocate(vector(n))
       h = matmul(x, matmul(matrix2, transpose(x)))     !!!!!!!!!!!!!!!!
       h = get_identity(n) - h
       vector = matmul(h,y)                              !!!!!!!!!!!!!!!!!!1
       tall =sum( y*vector)
       matrix2 = matrix2*tall/(n-2) !sjekk n-2
       uncert = matrix2(1,1)
       deallocate(h)
       deallocate(vector)
    end if

   ! Clean up
    deallocate(y)
    deallocate(x)
    deallocate(beta)
    deallocate(matrix2)


  end subroutine scattercalc

  !---------------------------------------------------------------------
  ! Estimate leakage map1 vs map2
  !----------------------------------------------------------------------

  subroutine estimate_leakage(map1, map2, nside, radius, leakage, uncert)

    real(dp),          dimension(:), intent(in)  :: map1, map2
    integer(i4b),                    intent(in)  :: nside
    real(dp),                        intent(in)  :: radius
    real(dp),                        intent(out) :: leakage, uncert

    integer(i4b)                                 :: i, j, p, n, num, nlist, pix   
    real(dp)                                     :: healnan=-1.6375d30, teller, nevner, radradius

    real(dp),        allocatable, dimension(:)   :: quiet, wmap
    real(dp),                     dimension(3)   :: gc_vec
    integer(i4b),    allocatable, dimension(:)   :: listpix

    ! Check that input maps are of same size
    if (size(map1) /= size(map2)) then
       write(*,*) size(map1), size(map2), 'Map1 and map2 not of same size. Quiting'
       stop    
    end if

    ! Find pixels in a 1deg radius around gc center
    gc_vec = [1,0,0]
    radradius = radius*pi/180.d0   ! radius in radians
    allocate(listpix(0:nside**2))
    call query_disc(nside, gc_vec, radradius, listpix, nlist, nest=1) 
    
    ! Find num of existing pixels for given area
    n = 0
    do p=0,nlist-1
       pix=listpix(p)
       if ( healok(map1(pix)) .and. healok(map2(pix))) then
          n = n+1
       end if
    end do
    ! Return if too few pixels
    if (n==0) then
       leakage = healnan
       return
    end if
    write(*,*) n, '= number of pixels'

    ! Calculate leakage parameter
    allocate(quiet(n))
    allocate(wmap(n))
    j=0
    do p=0,nlist-1
       pix=listpix(p)
       if ( healok(map1(pix)) .and. healok(map2(pix))) then
          j = j+1
          wmap(j)  = map1(pix)
          quiet(j) = map2(pix)
       end if
    end do

    teller = 0.d0
    nevner = 0.d0
    do i = 1,n
       teller = teller + quiet(i)*wmap(i)
       nevner = nevner + wmap(i)*wmap(i)
    end do
    leakage = teller/nevner

    uncert = healnan

   ! Clean up
    deallocate(quiet)
    deallocate(wmap)
    deallocate(listpix)

  end subroutine 

  !---------------------------------------------------------------------
  ! Unconvolve map by beamfile, and then convolve with gaussian
  !----------------------------------------------------------------------

  subroutine beam_convolve(map, ordering, nside, nmaps, beamfile, fwhm, onlypol)

    real(dp),          dimension(0:,:), intent(inout) :: map
    character(len=256),                 intent(in)    :: beamfile
    real(dp),                           intent(in)    :: fwhm
    integer(i4b),                       intent(in)    :: nside, nmaps, ordering
    logical(lgt),      optional,        intent(in)    :: onlypol

    complex(dpc),     allocatable, dimension(:,:,:)  :: alm
    real(dp),         allocatable, dimension(:,:)    :: dw8
    real(dp),                      dimension(2)      :: z
    real(dp),         allocatable, dimension(:, :) :: beam
    integer(i4b)      :: lmax, n, l, i, j, m, unit   
    real(dp)          :: tall, sigma
    real(dp),     allocatable, dimension(:)     :: power
    logical(lgt) :: useonlypol

    useonlypol = .false.
    if (present(onlypol)) then
       if (onlypol) useonlypol = .true.
    end if

    ! Read beam file
    if (myid==0) write(*,*) 'Reading from ', trim(beamfile)
    lmax = min(3*nside, 649)
    allocate(beam(0:lmax, 1:nmaps))
    call read_beam(beamfile, beam)
    
    ! Convert map from nest to ringed
    if (ordering ==2) then                       ! map2alm requieres ringed
       call convert_nest2ring(nside ,map)
    end if

!!$    ! output beams
!!$    sigma=fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0))
!!$    open(13, file='totbeam_quiet.txt')
!!$    do l = 0, lmax
!!$       write(13,*) l, 1.d0/beam(l, 2)* exp(-0.5d0*l*(l+1)*sigma**2)
!!$    end do
!!$    close(13)
!!$
!!$    open(13, file='beam_quiet_Q.txt')
!!$    do l = 0, lmax
!!$       write(13,*) l, beam(l, 2)
!!$    end do
!!$    close(13)

    ! Finding alm's for given skymap
    allocate(dw8(1:2*nside, 1:nmaps))
    allocate(alm(1:nmaps, 0:lmax, 0:lmax))
    dw8 = 1.d0
    z = 0.d0
    if (nmaps==1) then
       call map2alm(nside, lmax, lmax, map(:,1), alm, z, dw8)
    else
       call map2alm(nside, lmax, lmax, map, alm, z, dw8)
    end if
    deallocate(dw8)

!    ! and power spectrum
!    allocate(power(0:lmax))
!    power = 0.d0
!    do i = 0, lmax
!       do m = 0, i
!          if (m == 0) then
!             power(i) = power(i) + abs(alm(2, i, m))**2             
!          else
!             power(i) = power(i) + 2*abs(alm(2, i, m))**2
!          end if
!          !write(*,*) i,m,alm(npol,i,m)       
!       end do
!       power(i) = power(i)/(2.d0*i+1.d0)*i*(i+1.d0)/2.d0*pi
!    end do
!    open(13, file='power_before_quiet.txt')
!    do l = 0, lmax
!       write(13,*) l, power(l)
!    end do
!    close(13)

    ! Modifying beam
    sigma=fwhm*pi/180.d0/60.d0/sqrt(8.d0*log(2.d0))
    do i = 1, nmaps
       do l = 0, lmax
          alm(i, l, :) = alm(i, l, :)/beam(l, i) * exp(-0.5d0*l*(l+1.d0)*sigma**2)
       end do
    end do
    deallocate(beam)

!    ! and power spectrum
!    power = 0.d0
!    do i = 0, lmax
!       do m = 0, i
!          if (m == 0) then
!             power(i) = power(i) + abs(alm(2, i, m))**2             
!          else
!             power(i) = power(i) + 2*abs(alm(2, i, m))**2
!          end if
!    !write(*,*) i,m,alm(npol,i,m)       
!       end do
!       power(i) = power(i)/(2.d0*i+1.d0)*i*(i+1.d0)/2.d0*pi
!    end do
!    open(13, file='power_after_quiet.txt')
!    do l = 0, lmax
!       write(13,*) l, power(l)
!    end do
!    close(13)

    ! Back to map again
    if (nmaps==1) then
       call alm2map(nside, lmax, lmax, alm, map(:,1))
    else
       call alm2map(nside, lmax, lmax, alm, map)
    end if

    deallocate(alm)

    ! Convert map from ringed to nest
    call convert_ring2nest(nside ,map)

    ! put pol=sqrt(Q^2+U^2) in map(:,1)
    if (useonlypol) map(:,1) = sqrt(map(:,2)**2 + map(:,3)**2)
    
  end subroutine beam_convolve

  !-----------------------------------------------------------------------------------------------
  ! subroutine read_mapfile
  !-----------------------------------------------------------------------------------------------

  subroutine read_mapfile(infile, npix, nside, nmaps, ordering, map, nest, check_nside, check_nmaps, check_ordering)
    implicit none

    character(len=*),           intent(in)    :: infile
    real(dp),   dimension(:,:), pointer       :: map
    integer(i4b),               intent(out)   :: npix, nmaps, nside, ordering
    logical(lgt),     optional, intent(in)    :: nest
    integer(i4b),     optional, intent(in)    :: check_nmaps, check_nside, check_ordering

    integer(i4b)                   :: root=0
    real(dp)                       :: nullval
    logical(lgt)                   :: anynull, mismatch
   
    ! Read size of map
    if (myid==root) write(*,*) 'Reading from file ', trim(infile)
    npix = getsize_fits(infile, nmaps=nmaps, ordering=ordering, nside=nside)
    if (myid==root) write(*,*) nside, '= nside,', nmaps, '= nmaps,' ,ordering,'= ordering' 
    
    ! Checking that reference map and this map is of same size
    mismatch = .false.
    if (present(check_nside)) then    
       if (nside /= check_nside) mismatch = .true.
    end if
    if (present(check_nmaps)) then    
       if (nmaps /= check_nmaps) mismatch = .true.
    end if
    if (present(check_ordering)) then    
       if (ordering /= check_ordering) mismatch = .true.
    end if
    if (mismatch) then
       if (myid==root) write(*,*) trim(infile),'is not of same size as last map/reference map. Quiting'
       call mpi_finalize(ierr)
       stop
    end if
    
    ! Read map  and convert to nest if necassary
    allocate(map(0:npix-1,nmaps))
    call read_bintab(infile, map, npix, nmaps, nullval, anynull)
    if (present(nest) .and. nest) then
       if (ordering == 1) then                       ! we want nested
          write(*,*) 'Converting from ring to nest'
          call convert_ring2nest(nside ,map)
          ordering = 2
       end if
    end if

  end subroutine read_mapfile

  !-----------------------------------------------------------------------------------------------
  ! subroutine give_user_info
  !-----------------------------------------------------------------------------------------------

  subroutine give_user_info
    implicit none

    if (myid == root) then
       write(*,*) "Usage: 'sindex index parameters', where parameters are:" 
       write(*,*) '(maps, beams amd mask are fits-files. freqs, fwhm and nside are numbers)' 
       write(*,*) 'map1 map2 beam1 beam2 freq1 freq2 fwhm mask rms1 rms2 nside outprefix (uncertlimit)' 
    end if
    call mpi_finalize(ierr)
    stop

  end subroutine give_user_info



end program sindex
