! small programs that are useful

program sindex_utils
  use quiet_pixspace_mod
  use quiet_fileutils
  use rngmod
  use quiet_utils
  use math_tools
   use quiet_ces_mod
   use quiet_utils
   use quiet_hdf_mod
   use quiet_fft_mod
   use ziggurat
   use quiet_fft_mod
   use quiet_lx_mod
   use quiet_gain_mod
   use quiet_pointing_mod
   use quiet_target_mod
   use powell_mod
   use quiet_gain_mod
   use quiet_pointing_mod2
   use quiet_covfile_mod
   use quiet_mapfile_mod

  use sindex_maxlike_mod

  implicit none
  character(len=512) :: command

  call getarg(1, command)
  select case(command)
     case("makescatterfile"); call command_makescatterfile
     case("block2cov"); call command_block2cov
     case("planck2cov"); call command_planck2cov
     case("block2rms"); call command_block2rms
     case("makeWrms"); call command_makeWrms
     case("displaycov"); call command_displaycov
     case("adddiagcov"); call command_adddiagcov
     case("chisquare"); call command_chisquare
     case("a2t"); call command_a2t
     case("noisemask"); call command_noisemask
     case("gauss"); call command_gauss
     case("meanvalue"); call command_meanvalue
     case("regiongraph"); call command_regiongraph
     case("finalregion"); call command_finalregion
     case("normscatter"); call command_normscatter
     case("minmaxscatter"); call command_minmaxscatter
     case("makemask"); call command_makemask
     case("medianfilter"); call command_medianfilter
     case("makeQUmap"); call command_makeQUmap
     case("makebetamap"); call command_makebetamap
     case("removeregions"); call command_removeregions
     case("removepixels"); call command_removepixels
     case("removepts"); call command_removepts
     case("check_m"); call command_check_m
!     case("ringregions"); call command_ringregions
     case("nsidemap"); call command_nsidemap
     case("fitsmake1map"); call command_fitsmake1map
     case("readwritemap"); call command_readwritemap
     case("covmat"); call command_covmat
     case("check"); call command_check2
     case("addcovnoise"); call command_addcovnoise
     case("simulate_rms"); call command_simulate_rms
     case("simulate_map"); call command_simulate_map
     case("simulate"); call command_simulate
     case("histogram"); call command_histogram
     case("plothistogram"); call command_plothistogram
     case("makedatfile"); call command_makedatfile
     case("findlatitude"); call command_findlatitude
     case("calculate_rms_value"); call command_calculate_rms_value
     case("calculate_sensitivity"); call command_calculate_sensitivity
     case default
        write(*,*) "Unrecognized command : '" // trim(command) // "'"
!        call command_help
  end select

contains


!  subroutine command_covmat2covmatrix
!    implicit none
!
!    integer(i4b),                       intent(in) :: unit, ordering, polarization
!    character(len=*),                   intent(in) :: infile, outfile
!    logical(lgt),                       intent(in) :: inv
!    real(dp),           dimension(:,:), intent(in) :: matrix
!    character(len=3)                               :: option
!
!
!    write(*,*) "this converts a matrix from format: n / matrix"
!    write(*,*) "to: n / ordering / pol / matrix / inv (scalapost format)"
!    
!    if(iargc() /= 4) then
!       write(*,*) "Wrong number of arguments!"
!!       call command_help                                                                     !  
!       return
!    end if
!    call getarg(2, infile)
!    call getarg(3, outfile)
!    call getarg(4, option)
!
!    if (option=="inv") then
!       inv=.true.
!    else
!       inv=.false.   
!    end if
!
!    ! polarization = 1  => I only 
!    ! polarization = 2  => QU only   
!    ! polarization = 3  => IQU 
!    polarization=2
!    
!    call read_covmat(matrix, pixels, nside, order, ifname)
!    unit=42
!    call write_covmatrix(unit, outfile, ordering, polarization, matrix, inv)
!
!  end subroutine command_covmat2covmatrix

  subroutine command_block2cov
    implicit none
    integer(i4b)                                   :: unit, nside, ordering, polarization, n, i
    character(len=512)                             :: infile, outfile          
    logical(lgt)                                   :: inv
    integer(i4b), dimension(:), allocatable        :: pixels
    real(dp),     dimension(:,:), allocatable      :: matrix
    real(dp),     dimension(:,:,:), allocatable    :: inmap

    character(len=3)                               :: option

    write(*,*) "this converts a block QU noise map to a covariance matrix"

    if(iargc() /= 4) then          
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return          
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, option)
                                                                               
    if (option=="inv") then
       inv=.true.
    else
       inv=.false. 
    end if

    ! polarization = 1  => I only          
    ! polarization = 2  => QU only
    ! polarization = 3  => IQU
    polarization=2                                   
    write(*,*)polarization, trim(infile)
    call read_map(inmap, pixels, nside, ordering, trim(infile))

    write(*,*)"A"
!    call read_map(inmap, ordering, trim(infile))
    write(*,*)"B"
    n=size(inmap,1)
    write(*,*)n, size(pixels),inmap(1,1,1)
    allocate(matrix(n*2,n*2))
    write(*,*)n
    matrix=0d0
    do i=1,n
       matrix(i,i)    = inmap(i,1,1)
       matrix(n+i,n+i)= inmap(i,2,2)
       matrix(i,n+i)  = inmap(i,1,2)
       matrix(n+i,i)  = inmap(i,2,1)
    end do
    write(*,*)matrix(1,1), matrix(1000,1000)
!    call read_covmat(matrix, pixels, nside, order, ifname)
    unit=42
    write(*,*)matrix(1,n+1), matrix(n+1,1), unit
    call write_covmatrix(unit, trim(outfile), ordering, polarization, matrix, inv)            
  end subroutine command_block2cov

  subroutine command_planck2cov
    implicit none
    integer(i4b)                                   :: unit, nside, ordering, polarization, n, i
    character(len=512)                             :: infile, outfile          
    logical(lgt)                                   :: inv
    integer(i4b), dimension(:), allocatable        :: pixels
    real(dp),     dimension(:,:), allocatable      :: matrix
    real(dp),     dimension(:,:), allocatable      :: inmap, rmsmap

    character(len=3)                               :: option

    write(*,*) "this converts a block QU noise map to a covariance matrix"

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help                                                                                             
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, option)

    if (option=="inv") then
       inv=.true.
    else
       inv=.false.
    end if

    ! polarization = 1  => I only     
    ! polarization = 2  => QU only      
    ! polarization = 3  => IQU
    polarization=2
    write(*,*)polarization, trim(infile)
    call read_map(inmap, pixels, nside, ordering, trim(infile))
    !planck map is in Kelvin, need muK
    inmap=inmap*1e12

    write(*,*)"A"
    n=size(inmap,1)
    write(*,*)n, size(pixels)
    allocate(matrix(n*2,n*2))
    write(*,*)"B",n
    matrix=0d0
    do i=1,n
       matrix(i,i)    = 1d0/inmap(i,8)  !QQ
       matrix(n+i,n+i)= 1d0/inmap(i,10) !UU
       matrix(i,n+i)  = 1d0/inmap(i,9)  !QU
       matrix(n+i,i)  = 1d0/inmap(i,9)  !QU
    end do
    inv=.true.
    write(*,*)matrix(1,1), matrix(1000,1000)
!    call read_covmat(matrix, pixels, nside, order, ifname)                                                           
    unit=42
    write(*,*)matrix(1,n+1), matrix(n+1,1), unit
    call write_covmatrix(unit, trim(outfile), ordering, polarization, matrix, inv)

!optional, write a rms map too
    allocate(rmsmap(n,1:3))
    rmsmap(:,1) = sqrt(inmap(:,5)) *1e6
    rmsmap(:,2) = sqrt(inmap(:,8)) *1e6
    rmsmap(:,3) = sqrt(inmap(:,10))*1e6
    call write_map(rmsmap,pixels, nside, ordering, "outrms.fits")

  end subroutine command_planck2cov


  subroutine command_block2rms
    implicit none
    integer(i4b)                                   :: unit, nside, ordering
    character(len=512)                             :: infile, outfile          
    integer(i4b), dimension(:), allocatable        :: pixels
    real(dp),     dimension(:,:), allocatable      :: rmsmap
    real(dp),     dimension(:,:,:), allocatable    :: inmap
    real(dp)           :: healnan=-1.6375d30

    write(*,*) "this converts a block QU noise map to a noise rms map (only selects out two components, Q and U and take the square root of the elements. The rms (out) map is standard deviations while the block (covariance matrix) is variances."

    if(iargc() /= 3) then          
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return          
    end if
    call getarg(2, infile)
    call getarg(3, outfile)

    call read_map(inmap, pixels, nside, ordering, trim(infile))
    allocate(rmsmap(size(inmap,1),1:3))
    rmsmap(:,1)= healnan
    rmsmap(:,2)= sqrt(inmap(:,1,1))
    rmsmap(:,3)= sqrt(inmap(:,2,2))

    call write_map(rmsmap,pixels, nside,ordering,trim(outfile))

  end subroutine command_block2rms

  subroutine command_makeWrms
    implicit none
    integer(i4b)                                   :: unit, nside, ordering
    character(len=512)                             :: outfile          
    integer(i4b), dimension(:), allocatable        :: pixels1, pixels2, pixels3
    real(dp),     dimension(:,:), allocatable      :: map1,map2,map3,rmsmap
    real(dp)           :: healnan=-1.6375d30

    write(*,*) "this converts a block QU noise map to a noise rms map (only selects out two components, Q and U and take the square root of the elements. The rms (out) map is standard deviations while the block (covariance matrix) is variances."

    if(iargc() /= 2) then          
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return          
    end if
    call getarg(2, outfile)

    call read_map(map1, pixels1, nside, ordering, "wmap9_rms_n512_W1.fits")
    call read_map(map2, pixels2, nside, ordering, "wmap9_rms_n512_W2.fits")
    call read_map(map3, pixels3, nside, ordering, "wmap9_rms_n512_W3.fits")
    allocate(rmsmap(size(map1,1),1:3))
    rmsmap(:,1)= healnan
    rmsmap(:,2)= sqrt( (map1(:,2)+map2(:,2)+map3(:,2))/9.d0 )
    rmsmap(:,3)= sqrt( (map1(:,3)+map2(:,3)+map3(:,3))/9.d0 )

    call write_map(rmsmap,pixels1, nside,ordering,trim(outfile))

  end subroutine command_makeWrms


  subroutine command_adddiagcov
    implicit none
    integer(i4b)                                 :: npix,i,nside,ordering,polarization, min,max
    integer(i8b)                                 :: n
    character(len=512)                           :: covfile1,covfile2,covfile3,outcovfile
    logical(lgt)                                 :: inv
    real(dp),     dimension(:,:), allocatable    :: cov1,cov2,cov3
    integer(i4b), dimension(:), allocatable      :: pixels
    character(len=4)                             :: intext

    write(*,*) "this adds 3 diagonal cov matrices and scales them"

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
!       call command_help                                                          
       return
    end if
    call getarg(2, covfile1)
    call getarg(3, covfile2)
    call getarg(4, covfile3)
    call getarg(5, outcovfile)

    call read_covmatrix(42, trim(covfile1), ordering, polarization, cov1, inv, n)
    call read_covmatrix(43, trim(covfile2), ordering, polarization, cov2, inv, n)
    call read_covmatrix(44, trim(covfile3), ordering, polarization, cov3, inv, n)
    npix=n/2
    do i=1,npix
       cov1(i,i)=cov1(i,i)+cov2(i,i)+cov3(i,i)
       cov1(i,i)=cov1(i,i)*0.3333333333333333333d0*0.3333333333333333333d0
       cov1(i,npix+i)=cov1(i,npix+i)+cov2(i,npix+i)+cov3(i,npix+i)
       cov1(i,npix+i)=cov1(i,npix+i)*0.3333333333333333333d0*0.3333333333333333333d0
       cov1(npix+i,i)=cov1(npix+i,i)+cov2(npix+i,i)+cov3(npix+i,i)
       cov1(npix+i,i)=cov1(npix+i,i)*0.3333333333333333333d0*0.3333333333333333333d0
       cov1(npix+i,npix+i)=cov1(npix+i,npix+i)+cov2(npix+i,npix+i)+cov3(npix+i,npix+i)
       cov1(npix+i,npix+i)=cov1(npix+i,npix+i)*0.3333333333333333333d0*0.3333333333333333333d0
    end do

    call write_covmatrix(45, trim(outcovfile), ordering, polarization, cov1, inv)
    deallocate(cov1,cov2,cov3)
    
  end subroutine command_adddiagcov


  subroutine command_displaycov
    implicit none
    integer(i4b)                                 :: unit, nside, ordering, polarization, min,max
    integer(i8b)                                 :: n
    character(len=512)                           :: covfile
    logical(lgt)                                 :: inv
    real(dp),     dimension(:,:), allocatable    :: matrix
    integer(i4b), dimension(:), allocatable      :: pixels
    character(len=4)                             :: intext

    write(*,*) "this displays (parts of) covariance matrix"

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, covfile)
    call getarg(3, intext)
    read(intext,*) min
    call getarg(4, intext)
    read(intext,*) max

    unit=42
    call read_covmatrix(unit, trim(covfile), ordering, polarization, matrix, inv, n)
    call dump_matrix(matrix(min:max, min:max))
    write(*,*)"inverse?:", inv
    write(*,*)"QU component, n/2+i"
    call dump_matrix(matrix(n/2+min:n/2+max, min:max))


  end subroutine command_displaycov

  subroutine command_chisquare
    implicit none
    integer(i4b)                                 :: unit, nside, ordering, polarization, npix, i, index
    integer(i8b)                                 :: n
    character(len=512)                           :: covfile, wmapfile, quietfile, outfile, maskfile
    logical(lgt)                                 :: inv
    real(dp),     dimension(:,:), allocatable    :: matrix, wmapmap, quietmap, mask
    integer(i4b), dimension(:), allocatable      :: pixels, pixels2, pixels3
    character(len=4)                             :: intext
    real(dp)                                     :: chi2, alpha
    real(dp),     dimension(:), allocatable      :: map

    write(*,*) "this computes a chi_square and outputs to file"

    if(iargc() /= 5 .and. iargc() /= 6) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, covfile)
    call getarg(3, wmapfile)
    call getarg(4, quietfile)
    call getarg(5, outfile)

    unit=42
    call read_covmatrix(unit, trim(covfile), ordering, polarization, matrix, inv, n)
    call read_map(wmapmap, pixels, nside, ordering, trim(wmapfile))
    call read_map(quietmap, pixels2, nside, ordering, trim(quietfile))

    write(*,*)"inverse?:", inv
    write(*,*)"sizes",n,size(matrix,1), size(wmapmap,1), size(quietmap,1)

    !mask out pixels in this mask-file
    if (iargc() == 6) then
       call getarg(6, maskfile)
       call read_map(mask, pixels3, nside, ordering, trim(maskfile))    
       write(*,*)"will apply mask file to data",size(pixels3),"pixels masked out"
       do i=1,size(pixels3)
          index=pixels3(i)
!          write(*,*)i,index
          where(pixels == index)
             wmapmap(:,2)=-1.637500000000E-30
             wmapmap(:,3)=-1.637500000000E-30
             quietmap(:,2)=1.637500000000E-30
             quietmap(:,3)=1.637500000000E-30
          end where
       end do
       deallocate(mask)
    endif

    !calculate the chisquare
    npix=size(wmapmap,1)
    open(unit,file=trim(outfile),action="write")

    alpha=1.0d0
    do i=1,6000
       alpha=alpha+0.0001d0
       allocate(map(1:2*npix))
       map(1:npix)=wmapmap(:,2)-alpha*quietmap(:,2)
       map(npix+1:2*npix)=wmapmap(:,3)-alpha*quietmap(:,3)
!       if(i==1)write(*,*)map
       chi2= dot_product(map,matmul(matrix,map))
!       chi2= dot_product(map(1:npix),matmul(matrix(1:npix,1:npix),map(1:npix))) only Q
       deallocate(map)
       write(unit,*)alpha,chi2
    end do
    close(unit)
    deallocate(wmapmap,quietmap,matrix)

  end subroutine command_chisquare

  subroutine command_a2t
    implicit none
    integer(i4b)                                 :: unit
    real(dp)                                     :: freq1,freq2, beta, factor

    write(*,*) "program a2t"

    if(iargc() /= 1) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if

    !for Q-band
!    freq1=40.51d0
!    freq2=43.1d0
!    beta=-3d0

    !for W-band
    freq1=94.2d0
    freq2=94.5d0
    beta=1.59d0

    !for Q-band gb
!    freq1=40.51d0
!    freq2=43.1d0
!    beta=-2.93d0
   
    factor=ant2thermo(freq1)/ant2thermo(freq2) * ((freq1/freq2)**beta)
    write(*,*)"the factor is:",factor
    write(*,*)"ant2thermo 1:",ant2thermo(freq1)
    write(*,*)"ant2thermo 2:",ant2thermo(freq2)
    write(*,*)"freq 1:",freq1
    write(*,*)"freq 2:",freq2
    write(*,*)"beta:",beta

  end subroutine command_a2t


  subroutine command_noisemask
    implicit none
    integer(i4b)                                 :: unit, nside, ordering, polarization, npix, i, index
    integer(i8b)                                 :: n
    character(len=512)                           :: inmaskfile, wmapcovfile, quietcovfile, outmaskfile
    logical(lgt)                                 :: inv
    real(dp),     dimension(:,:), allocatable    :: wmapcov, quietcov, mask
    integer(i4b), dimension(:), allocatable      :: pixels
    character(len=4)                             :: intext
    real(dp)                                     :: chi2, alpha

    write(*,*) "program noisemask"

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, inmaskfile)
    call getarg(3, wmapcovfile)
    call getarg(4, quietcovfile)
    call getarg(5, outmaskfile)

    unit=42
    call read_covmatrix(unit, trim(wmapcovfile), ordering, polarization, wmapcov,inv, n)
    call read_covmatrix(unit, trim(quietcovfile), ordering, polarization, quietcov,inv, n)
    call read_map(mask, pixels, nside, ordering, trim(inmaskfile))
    npix=size(mask,1)

    !make a mask that contains the pixels we want to remove
    do i=1,npix
       if (quietcov(i,i) > 1.d0* wmapcov(i,i).OR.(quietcov(npix+i,npix+i) > 1.d0* wmapcov(npix+i,npix+i)))then
          !pixel is bad, so we keep it.
       else
          !pixel is ok, so we remove it from the mask
          mask(i,1)=-1.6375d30
       end if
    end do
    mask(:,2)=mask(:,1)
    mask(:,3)=mask(:,1)

    call write_map(mask,pixels, nside,ordering,trim(outmaskfile))
    deallocate(wmapcov,quietcov,mask)
    
  end subroutine command_noisemask

  subroutine command_meanvalue
    implicit none

    character(len=512) :: infile
    character(len=10)  :: option, intext
    integer(i4b)       :: min, max, unit, i, j, ri, regions
    real(dp)           :: b, u, spectral_index, u_stat, u_inv
    real(dp),          dimension(:), allocatable :: betas, uncerts
    real(dp)           :: healnan=-1.6375d30
!    logical (lgt)      :: 

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, intext)
    read(intext,*) min
    call getarg(4, intext)
    read(intext,*) max

    regions=max!80!24 !32

    allocate(betas(regions), uncerts(regions))
    betas=0.d0
    uncerts=0.d0
    open(42,file=trim(infile),action="read",status="old")
    do i=1,regions
       read(42, *) ri, b, u
!       write(*,*) ri, b, u
       betas(ri) = b
       uncerts(ri) = u
!       if (ri == regions) exit
    end do
    close(42)

!    spectral_index = sum(betas)/ real(size(betas),dp)

    spectral_index = (sum(betas(min:max)/(uncerts(min:max)**2))) / (sum(1./(uncerts(min:max)**2))) !/ real(size(betas),dp)
    u_stat = (sum(uncerts(min:max)))/ real((max-min+1),dp)
    u_inv  = sum(1.d0/uncerts(min:max)**2)
    u_inv  = sqrt(1.d0/u_inv)

    write(*,*) "Mean spectral index (invers varians vekting) from regions", min, " to", max
    write(*,*) spectral_index, "+-", u_stat, "mean uncertainty"
    write(*,*) spectral_index, "+-", u_inv, "inv weigh uncertainty"

  end subroutine command_meanvalue

  subroutine command_regiongraph
    implicit none

    character(len=512) :: prefix, mapfile, uncertfile, outfile, prefilename, filename
    character(len=10)  :: option, intext
    character(len=1)   :: i_char_1
    character(len=2)   :: i_char_2
    integer(i4b)       :: unit, nside, order, npix, i, j, imax, imin, comp
    real(dp),     dimension(:,:), allocatable :: map, uncert
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp)           :: healnan=-1.6375d30
    logical (lgt)      :: nops, gal

    if(iargc() /= 6) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, option)
    call getarg(3, intext)
    read(intext,*) comp
    call getarg(4, mapfile)
    call getarg(5, uncertfile)
    call getarg(6, outfile)

    nops=.false.
    gal =.false.
    imin=1
    imax=32
    if(option=="nops") then
       nops=.true.
    else if(option=="gal") then
       gal=.true.
    else if(option=="sim") then
       imax=80
    else if(option=="t") then
       imax=24
    end if

    if (comp==1) then
       unit =42
       open(unit,file=trim(outfile),action="write")
    else if(comp==2) then
       filename = "Q_"//trim(outfile)
       open(43,file=trim(filename), action="write")
       filename = "U_"//trim(outfile)
       open(44,file=trim(filename), action="write")
    end if
    do i= imin,imax
       if ( nops .AND. i.GE.15 .AND. i.LE.22) cycle
       if ( gal .AND. (i.LE.12 .OR. (i.GE.15 .AND. i.LE.22))) cycle
       if (i<10) then
          call int2string(i,i_char_1)
          prefilename="s"// i_char_1 //"/"
          if(option=="sim") prefilename="Aflat/k"// i_char_1 //"/"
          if(option=="t") prefilename="t"// i_char_1 //"/"
       else
          call int2string(i,i_char_2)
          prefilename="s"// i_char_2 //"/"
          if(option=="sim") prefilename="Aflat/k"// i_char_2 //"/"
          if(option=="t") prefilename="t"// i_char_2 //"/"
       end if
       filename = trim(prefilename) // trim(mapfile)
       call read_map(map, pixels,   nside, order, trim(filename))
       filename = trim(prefilename) // trim(uncertfile)
       call read_map(uncert,pixels2,nside, order, trim(filename))
    
       if (comp==1) then
          write(unit,*) i, map(1,1), uncert(1,1)
       else if(comp==2) then
          write(43,*) i, map(1,2), uncert(1,2)
          write(44,*) real(i+0.2, dp), map(1,3), uncert(1,3)          
       end if

       deallocate(map, uncert, pixels, pixels2)

    end do
    if(comp==1) write(*,*) "Wrote to file ", trim(outfile)
    if(comp==2) write(*,*) "Wrote to file ", "Q_"//trim(outfile)
    close(unit)

  end subroutine command_regiongraph

  subroutine command_gauss
    implicit none
    real(dp)     :: xmin,xmax,range,dx,x
    integer(i4b) :: i
    open(46,file="gauss.txt", action="write")

    xmin=-5.d0
    xmax=5.d0
    range=xmax-xmin
    dx=range/999.
    write(46,*) xmin, exp(-xmin*xmin/2.) / sqrt(2.*pi)
    do i=2,1000
       x=xmin+(i-1)*dx
       write(46,*) x, exp(-x*x/2.) / sqrt(2.*pi)
    end do
    close(46)
  end subroutine command_gauss

  subroutine command_finalregion
    implicit none
    character(len=512) :: prefix, maxlikefile, scatterfile, outfile, prefilename, filename
    character(len=10)  :: option, intext
    character(len=1)   :: i_char_1
    character(len=2)   :: i_char_2
    integer(i4b)       :: unit, nside, order, npix, i, j, imax, imin, comp, reg,ii
    real(dp), dimension(:), allocatable :: maxlike, scatter, uncertmax, uncertscatter
    real(dp),     dimension(:,:), allocatable :: map, uncert
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp)           :: healnan=-1.6375d30, limit, difference, newvalue, newuncert
    logical (lgt)      :: nops, gal

    reg=24
    allocate(maxlike(reg), scatter(reg), uncertmax(reg), uncertscatter(reg))

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, intext)
    read(intext,*) limit
    call getarg(3, maxlikefile)
    call getarg(4, scatterfile)
    call getarg(5, outfile)
    
    open(44,file=trim(maxlikefile), action="read")
    open(45,file=trim(scatterfile), action="read")
    open(46,file=trim(outfile), action="write")
    do i = 1,reg
       read(44,*) ii, maxlike(i),uncertmax(i)
       read(45,*) ii, scatter(i),uncertscatter(i)

       difference = abs(abs(maxlike(i))-abs(scatter(i)))
       if(difference >limit)then 
          write(*,*)i
          cycle
       else
          newvalue = (maxlike(i)+scatter(i)) /2.d0
          newuncert = uncertmax(i) + difference
          write(46,*) i, newvalue, newuncert
       end if
    end do
    close(44)
    close(45)
    close(46)

  end subroutine command_finalregion

  subroutine command_normscatter
    implicit none
    character(len=512) :: prefix, infile, outfile, prefilename, filename
    character(len=10)  :: option, intext
    character(len=1)   :: i_char_1
    character(len=2)   :: i_char_2
    integer(i4b)       :: unit, nside, order, npix, i, j, reg,ii, io, tall, imax
    real(dp), dimension(:), allocatable :: xvar, yvar, newxvar, newyvar, x, y
    real(dp),     dimension(:,:), allocatable :: map, uncert
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp)           :: healnan=-1.6375d30, factor, difference, newvalue, newuncert, xmean, ymean
    logical (lgt)      :: nops, gal

    tall=10000
    allocate(xvar(tall), yvar(tall))
    xvar=healnan
    yvar=healnan
    newxvar=healnan
    newyvar=healnan

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, intext)
    read(intext,*) factor
    call getarg(3, infile)
    call getarg(4, outfile)
    
    open(44,file=trim(infile), action="read")
    i=1
    do while(.true.)
       read(44,*,IOSTAT=io) xvar(i), yvar(i)
       if(io<0) exit
       i=i+1
    end do
    i=i-1
    imax=i
    close(44)
    allocate(x(imax),y(imax), newxvar(imax), newyvar(imax))

    x=xvar(1:imax)
    y=yvar(1:imax)
    deallocate(xvar,yvar)
    write(*,*) 1, x(1), y(1)
    write(*,*) i, x(imax), y(imax)

    xmean = sum(x)/imax
    ymean = sum(y)/imax
    write(*,*) "xmean =",xmean,"ymean =",ymean
    newxvar=x-xmean
    newyvar=y-ymean
!    factor = abs(MAXVAL(x)-MINVAL(x))
    newxvar=newxvar/factor
    newyvar=newyvar/factor

    open(45,file=trim(outfile), action="write")
    do i=1,imax
       write(45,*) newxvar(i), newyvar(i)
    end do
    close(45)
    deallocate(x, y, newxvar, newyvar)

  end subroutine command_normscatter

  subroutine command_minmaxscatter
    implicit none
    character(len=512) :: infile, outfile
    character(len=10)  :: option, intext
    integer(i4b)       :: i, io, tall, imax
    real(dp), dimension(:), allocatable :: xvar, yvar
    real(dp)           :: healnan=-1.6375d30,minxvar,maxxvar

    tall = 20000
    allocate(xvar(tall), yvar(tall))
    xvar=0.
    yvar=0.

    if(iargc() /= 1) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if

    infile="ut_paper_scatter_a00000.txt"
    open(44,file=trim(infile), action="read")
    i=1
    do while(.true.)
       read(44,*,IOSTAT=io) xvar(i), yvar(i)
       if(io<0) exit
       i=i+1
    end do
    i=i-1
    close(44)
    infile="ut_paper_scatter_a00045.txt"
    open(44,file=trim(infile), action="read")
    do while(.true.)
       read(44,*,IOSTAT=io) xvar(i), yvar(i)
       if(io<0) exit
       i=i+1
    end do
    i=i-1
    close(44)

    imax=i

    minxvar=MINVAL(xvar)
    maxxvar=MAXVAL(xvar)

    outfile="ut_paper_scatter_a00_a45_minxvalue.txt"
    open(45,file=trim(outfile), action="write")
    write(45,*) minxvar
    close(45)
    outfile="ut_paper_scatter_a00_a45_maxxvalue.txt"
    open(45,file=trim(outfile), action="write")
    write(45,*) maxxvar
    close(45)
    deallocate(xvar,yvar)

  end subroutine command_minmaxscatter

  subroutine command_makemask
    implicit none

    character(len=512) :: maskfile, outfile
    integer(i4b)       :: nside, order, npix, i, j, n, isnest, teller, ngal, t1
    real(dp),     dimension(:,:), allocatable :: map, map2, uncert, outmap
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2, rpix
    real(dp)           :: theta, phi, lat, lon, vec(3), r, lon2,dlon, val
    real(dp)           :: healnan=-1.6375d30

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, maskfile)
    call getarg(3, outfile)

    call read_map   (map, order, trim(maskfile))
    isnest = 1
    npix=size(map,1)
    nside = npix2nside(npix)
    allocate(outmap(0:npix-1,1), rpix(0:npix-1))
    outmap=healnan

    !point sources

    !mark pixels in Virgo A
    lat=74.49 ; lon=283.78
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 15.d0

    !mark pixels in another region
    lat=90.-25.8 ; lon=290.
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 16.d0

    !mark pixels in the Cen-A
    lat=19.42 ; lon=309.51            !(l,b) = (309.51 deg, 19.42 deg)
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 17.d0

    !mark pixels in another region
    lat=10. ; lon=296.5
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 18.d0

    !mark pixels in 3C273
    lat=4.36 ; lon=289.95
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 19.d0

    !mark pixels in the gc
    lat=0.d0 ; lon=0.d0
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 20.d0

    !mark pixels in the tau-a
    lat=-5.784d0 ; lon=184.557d0
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 21.d0

    !mark pixels in formax A
    lat=-56.69 ; lon=240.16
    r=3.d0
    call ang2vec(pi/2 - lat*DEG2RAD, lon*DEG2RAD, vec)
    call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=isnest)
    outmap(rpix(0:n-1),1) = 22.d0

    ngal=0
    t1=0
    do i=0,npix-1

       if(outmap(i,1).GE.15.d0 .and. outmap(i,1).LE.22.d0) cycle

       call pix2ang_nest(nside, i, theta, phi)
       theta=theta*RAD2DEG; phi = phi * RAD2DEG

       if (map(i,1)>0.5) then
          if (theta<=90. .and. phi<=180. .and. phi>7. ) then  !left/upper quadrant
             if(theta<=30. ) then 
                outmap(i,1)=1.d0
             else if(theta>=30. .and. phi>=90. ) then 
                outmap(i,1)=2.d0
             else if(theta>=30. .and. phi<=90. ) then 
                outmap(i,1)=3.d0
             end if
          else if (theta<=90. .and. ((phi>=180. .and. phi<=360.) .or. phi<=7.)) then  !right/upper quadrant
             if(theta<=30. ) then 
                outmap(i,1)=4.d0
             else if(theta>=30. .and. (phi>=270. .or.phi<=8.) ) then 
                outmap(i,1)=5.d0
             else if(theta>=30. .and. phi<=270. .and. phi>50. ) then 
                outmap(i,1)=6.d0
             end if
          else if (theta>=90. .and. ((phi<=180. .and. phi>=0.) .or. phi>=353.) ) then  !left/lower quadrant
             if(theta>=150. ) then 
                outmap(i,1)=7.d0
             else if(theta<=150. .and. phi>=90. .and.phi<200. ) then 
                outmap(i,1)=8.d0
             else if(theta<=150. .and. (phi<=90. .or. phi>300.)) then 
                outmap(i,1)=9.d0
             end if
          else if (theta>=90. .and. ((phi>=180. .and. phi<=353.) .or. phi<=7. )) then  !right/lower quadrant
             if(theta>=150. ) then 
                outmap(i,1)=10.d0
             else if(theta<=150. .and. phi>=270. ) then 
                outmap(i,1)=11.d0
             else if(theta<=150. .and. phi<=270. ) then 
                outmap(i,1)=12.d0
             end if
          end if
       else
          if (theta<=72. .and. (phi<50. .or. phi>310.) ) then   !north gal spur
             outmap(i,1)=13.d0
             t1=t1+1
          else if (theta>=104. .and. (phi<50. .or. phi>310.) ) then   !south gal spur
             outmap(i,1)=14.d0
          else
             outmap(i,1)=23.d0
             ngal=ngal+1
             rpix(ngal)=i
             
          end if
       end if

    end do

!
    val=23.d0
    teller=0
    dlon=2.
    lon = 0.d0
    lon2= 0.d0
    do while(lon2.LT.360)
       lon= lon2
       lon2 = lon + dlon
       do i=1,ngal
          call pix2ang_nest(nside,rpix(i),theta,phi) !phi=[0,2pi]
          phi = phi*RAD2DEG
          if ((phi.GE.lon).AND.(phi.LE.lon2)) then
             outmap(rpix(i),1) = val
             teller = teller + 1
          end if
       end do
       if(teller.GT.1050) then  !this is the criterium which decides the region size, in pixels. 2280
          val = val+1
          teller = 0
       end if
    end do
    if (teller.LT.500) then  
       !if last region is too small, we combine it with the one before.
       do i = 1, ngal
          if (map(rpix(i),1) == val) then
             map(rpix(i),1) = val - 1
          end if
       end do
    else if (teller.NE.0) then
       val = val + 1 
    end if
    teller = 0
       


    write(*,*),"ngal,t1", ngal, t1
    call system("echo hei")
    call write_map(outmap, order, outfile)
    deallocate(map, outmap)
  end subroutine command_makemask


  subroutine command_medianfilter
    implicit none

    character(len=512) :: mapfile, outfile, intext
    integer(i4b)       :: nside, order, npix, i, j, n, nmaps
    real(dp),     dimension(:,:), allocatable :: map, outmap
    real(dp),     dimension(:), allocatable   :: diskmap
    integer(i4b), dimension(:),   allocatable :: pixels, rpix
    real(dp)           :: vec(3),r, lat, lon, val

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, mapfile)
    call getarg(3, outfile)
    call getarg(4, intext)
    read(intext,*) r

    call read_map   (map, order, trim(mapfile))
    npix=size(map,1)
    nmaps=size(map,2)
    nside = npix2nside(npix)
    allocate(outmap(0:npix-1,nmaps), rpix(0:npix-1))
    outmap=-2

    where (map /= 1.d0)
       map = 0.d0
    end where
    do j=1,nmaps
       if (.not. any(map(:,j)==1.d0)) cycle

       do i=0,npix-1
          if (mod(i,1000)==0) write(*,*) j, i, npix
          call pix2vec_nest(nside, i, vec)
          call query_disc(nside, vec, r*DEG2RAD, rpix, n, nest=1)
          
          if(any(map(rpix(0:n-1),j) .gt. 0.5) .and. any(map(rpix(0:n-1),j) .lt. 0.5)) then
             allocate(diskmap(n))
             diskmap = map(rpix(0:n-1),j)
             call QuickSort_real(diskmap)
             if (mod(n,2) ==0) n=n+1
             val = diskmap(n/2)
             outmap(i,j) = val
             deallocate(diskmap)
          else
             outmap(i,j) = map(i,j)
          end if
       end do
    end do
    call write_map(outmap, order, outfile)
    deallocate(map, outmap,rpix)
  end subroutine command_medianfilter

  subroutine command_makescatterfile
    implicit none

    character(len=512) :: mapfile1, mapfile2, outfile, out1,out2
    integer(i4b)       :: nside, order, npix, i
    real(dp),     dimension(:,:), allocatable :: map1, map2
    integer(i4b), dimension(:),   allocatable :: pixels1, pixels2

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, mapfile1)
    call getarg(3, mapfile2)
    call getarg(4, outfile)

    call read_map   (map1, pixels1,   nside, order, trim(mapfile1))
    call read_map   (map2, pixels2,   nside, order, trim(mapfile2))
    npix=size(map1,1)
    out1="Q_"//outfile
    out2="U_"//outfile
    open(44,file=trim(out1), action="write")
    open(45,file=trim(out2), action="write")

    do i=1,npix
       write(44,*)map1(i,2), map2(i,2)
       write(45,*)map1(i,3), map2(i,3)
    end do
    close(44)
    close(45)
    deallocate(map1,map2)
  end subroutine command_makescatterfile


  subroutine command_makeQUmap
    implicit none

    character(len=512) :: mapfile, uncertfile, outfile
    integer(i4b)       :: nside, order, npix, i, j
    real(dp),     dimension(:,:), allocatable :: map, uncert, totmap
    integer(i4b), dimension(:),   allocatable :: pixels1, pixels2

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, mapfile)
    call getarg(3, uncertfile)
    call getarg(4, outfile)

    call read_map   (map, pixels1,   nside, order, trim(mapfile))
    call read_map   (uncert, pixels2,   nside, order, trim(uncertfile))
    npix=size(map,1)
    allocate(totmap(1:npix,1))
    do i=1,npix
          totmap(i,1)= ( map(i,2)-map(i,3) )/ sqrt( uncert(i,2)**2 + uncert(i,3)**2 )
    end do
    call write_map(totmap, pixels1, nside, order, outfile)
    deallocate(map)
  end subroutine command_makeQUmap

  subroutine command_makebetamap
    implicit none

    character(len=512) :: mapfile, uncertfile, outfile
    integer(i4b)       :: nside, order, npix, i, j
    real(dp),     dimension(:,:), allocatable :: map, uncert
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, mapfile)
    call getarg(3, uncertfile)
    call getarg(4, outfile)

    call read_map   (map, pixels,   nside, order, trim(mapfile))
    call read_map   (uncert, pixels2,   nside, order, trim(uncertfile))
    npix=size(map,1)
    do i=1,npix
       do j=1,1  !change this to 1,3 if you want map 2 and 3 also!
!          map(i,j)= (map(i,j)+3.d0 )/uncert(i,j)    !this is for betamap
           map(i,j)= (map(i,j)- uncert(i,j) )    !this is for diffmap
!          map(i,j) = 10.d0   !for a constant value
       end do
    end do
    call write_map(map, pixels, nside, order, outfile)
    deallocate(map)
  end subroutine command_makebetamap

  subroutine command_removeregions
    implicit none

    character(len=512) :: mapfile, uncertfile, outfile, intext
    integer(i4b)       :: nside, order, npix, i, j
    real(dp),     dimension(:,:), allocatable :: map, uncert
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp)           :: threshold
    real(dp)           :: healnan=-1.6375d30

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, mapfile)
    call getarg(3, uncertfile)
    call getarg(4, outfile)
    call getarg(5, intext)
    read(intext,*) threshold

    call read_map   (map, pixels,   nside, order, trim(mapfile))
    call read_map   (uncert, pixels2,   nside, order, trim(uncertfile))
    npix=size(map,1)
    do i=1,npix
       do j=1,3
          if(uncert(i,j)> threshold ) map(i,j) = healnan
       end do
    end do
    call write_map(map, pixels, nside, order, outfile)
    deallocate(map)
  end subroutine command_removeregions

  subroutine command_removepixels
    implicit none

    character(len=512) :: mapfile, outfile, maskfile
    integer(i4b)       :: nside, order, npix, i, index
    real(dp),     dimension(:,:), allocatable :: map, mask
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp)           :: healnan=-1.6375d30

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help                                                       
       return
    end if
    call getarg(2, mapfile)
    call getarg(3, maskfile)
    call getarg(4, outfile)

    call read_map   (map, pixels,   nside, order, trim(mapfile))
    call read_map  (mask, pixels2,  nside, order, trim(maskfile))

    do i=1,size(pixels2)
       index=pixels2(i)
       !          write(*,*)i,index
       where(pixels == index)
          map(:,1)=healnan
          map(:,2)=healnan
          map(:,3)=healnan
       end where
    end do

    call write_map(map, pixels, nside, order, outfile)
    deallocate(map, mask)
  end subroutine command_removepixels

  subroutine command_removepts
    implicit none

    character(len=512) :: mapfile, outfile
    integer(i4b)       :: nside, order, npix, i, j
    real(dp),     dimension(:,:), allocatable :: map
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp)           :: healnan=-1.6375d30

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help                                                       
       return
    end if
    call getarg(2, mapfile)
    call getarg(3, outfile)

    call read_map   (map, pixels,   nside, order, trim(mapfile))
    npix=size(map,1)
    do i=1,npix
       if(map(i,1).GT.14) then
          if(map(i,1).LT.23) then
             map(i,1)=healnan
          else
             map(i,1)=map(i,1)-8.d0
          endif
       end if
    end do
    call write_map(map, pixels, nside, order, outfile)
    deallocate(map)
  end subroutine command_removepts



  subroutine command_check_m
    implicit none
    character(len=1)         :: i_char
    character(len=512)        :: mapfile, covfile, dir
    real(dp)                  :: offset(3), sigma, p(1), uncert, std, tall
    real(dp)                  :: amp, betain, beta_mean, beta_dev, off_mean, off_dev
    real(dp),     allocatable :: map(:,:),imap(:,:,:,:), incov(:,:,:,:), weights(:,:)
    real(dp),     allocatable :: smap(:,:,:,:), w1(:),w2(:), L(:,:,:,:), A(:,:), m(:,:),loops_beta(:), loops_offset(:,:,:)
    real(dp),     allocatable :: tmap(:,:)
    integer(i4b), allocatable :: pixels(:), pixels2(:)
    integer(i4b)              :: nside, ncomp, order, status, seed, n, i, unit, k, npix, j
    integer(i4b)              :: loops, nreg, r, jj, nband, i2, bi, coff, moff, cilen
    type(planck_rng)     :: rng_handle
    type(pixinfo):: pinfo, pinfo2
    type(banddata), allocatable, dimension(:) :: data

!!!!!! change this
    nside   = 256
    npix    = 5
    nreg    = 1
    loops   = 10
    nband   = 2
    k       = 4
    seed    = -91300
    amp     = 50
    betain  = -3
    offset  = [0d0,100d0,-500d0]
    sigma   = 10000.d0
!!!!!! end change this
    nside = 256
    call rand_init(rng_handle,seed)
    n = npix
    select case(k)
       case(1); coff=0; moff=0; cilen=  n; ncomp=1
       case(2); coff=0; moff=1; cilen=2*n; ncomp=1
       case(3); coff=n; moff=2; cilen=2*n; ncomp=1
       case(4); coff=0; moff=1; cilen=2*n; ncomp=2
    end select

    ! ========= Begin input =========
    allocate(imap(1:n,1:3,nband,1:nreg), incov(cilen,cilen, nband, 1:nreg)) 

    !read in from files
!!$    do i=1,nreg
!!$       call int2string(i,i_char)
!!$       dir='dir'//trim(i_char)
!!$       call read_map   (map, pixels,   nside, order, trim(dir)//"/map_raw0.fits")
!!$       call read_covmat(cov,  pixels2, nside, order, ncomp0, trim(dir)//"/cov0.hdf")
!!$!    npix=size(imap0,1)
!!$!       incov(:,:,i) = cov(1:npix*2,1:npix*2)
!!$       incov(1:npix,1:npix,i) = cov(1:npix,1:npix)
!!$       incov(1+npix:2*npix,1+npix:2*npix,i) = cov(size(map,1):size(map,1)+npix,size(map,1):size(map,1)+npix)
!!$       incov(1+npix:2*npix,1:npix,i)        = cov(size(map,1):size(map,1)+npix, 1:npix)
!!$       incov(1:npix,1+npix:2*npix,i)        = cov(1:npix, size(map,1):size(map,1)+npix)
!!$       imap0(:,:,i) = map(1:npix,:)
!!$       deallocate(cov, pixels,pixels2)
!!$       deallocate(map)
!!$    end do


!    !or create
    incov = 0.d0
    do r=1,nreg
       do bi = 1, nband
          do i=1, ncomp
             do j = 1, n
                imap(j,i+moff,bi,r) = amp*rand_gauss(rng_handle)   !+?
                i2 = coff +(i-1)*n +j           !(i-1)*n+j
                incov(i2,i2,bi,r) = sigma*abs(rand_gauss(rng_handle))
             end do
          end do
write(*,*) "incov", bi, r
call dump_matrix(incov(:,:,bi,r))
write(*,*) "imap", bi, r
call dump_matrix(imap(:,:,bi,r))
       end do
    end do

         
    allocate(A(npix*ncomp,1:nreg), m(ncomp,1:nreg))
    !create data-type
    allocate(data(1:nreg))    
    do r=1,nreg
       allocate(data(r)%freq(nband))
       allocate(data(r)%map(npix,3,nband))
       allocate(data(r)%cov(cilen,cilen,nband))    
       do i = 1, nband
          data(r)%freq(i) = 13+10*i
       end do
       data(r)%pol   = cilen/n
       data(r)%nmaps = 3
       data(r)%npix  = npix
       data(r)%n     = npix
    end do

    do r = 1, nreg
       do bi = 2, nband
          imap(:,:,bi,r) = imap(:,:,1,r) * ant2thermo(data(r)%freq(bi))/ant2thermo(data(r)%freq(1)) * (data(r)%freq(bi)/data(r)%freq(1))**betain
          do i = 1, n
             imap(i,:,bi,r) = imap(i,:,bi,r) + offset
          end do
       end do
    end do

    ! ========== Fake input done ============

    allocate(w1(n*ncomp))
    allocate(L(n*ncomp,n*ncomp,nband,nreg))
    allocate(smap(size(imap,1),size(imap,2),size(imap,3),size(imap,4)))
    allocate(loops_beta(1:loops), loops_offset(1:loops,nreg,ncomp))
    
    do r = 1, nreg
       do bi = 1, nband
          call cholesky_decompose(incov(coff+1:coff+n*ncomp,coff+1:coff+n*ncomp,bi,r),L(:,:,bi,r))
write(*,*) "F", coff+1, coff+n*ncomp, bi, r, shape(incov)
write(*,*) "G", shape(L)
       end do
    end do

!loop here
    do i=1,loops
       !add covnoise
       do bi = 1, nband
          do r=1,nreg
             do j = 1, n*ncomp
                w1(j) = rand_gauss(rng_handle)
             end do
             smap(:,1+moff:ncomp+moff,bi,r) = imap(:,1+moff:ncomp+moff,bi,r) + &
                  reshape(matmul(L(:,:,bi,r),w1),[n,ncomp])
          end do
       end do

       !smooth maps
!!$       nside=256
!!$       order=nest
!!$       if (.not. allocated(pixels)) allocate(pixels(1:npix))
!!$       pixels = irange(npix)-1
!!$       call setup2("cut_beam0.txt", nside, order, pixels, pinfo, weights)
!!$       if (.not. allocated(tmap)) allocate(tmap(1:npix,1:2))
!!$       do r=1,nreg
!!$          call smooth_map  (smap0(:,2:3,r), pinfo, weights, tmap)
!!$          smap0(:,2:3,r)=tmap(:,:)
!!$          call smooth_map  (smap1(:,2:3,r), pinfo, weights, tmap)
!!$          smap1(:,2:3,r)=tmap(:,:)
!!$!          write(*,*) weights, smap0(:,2:3,r)
!!$       end do

       !put in data-type
       do r=1,nreg
          do bi = 1, nband
             data(r)%map(:,:,bi) = smap(:,:,bi,r)
             data(r)%cov(:,:,bi) = incov(:,:,bi,r)
write(*,*) "incov", bi, r
call dump_matrix(data(r)%cov(:,:,bi))
write(*,*) "smap", bi, r
call dump_matrix(data(r)%map(:,:,bi))
          end do
       end do
 
       !overwrite map
!!$       data(1)%map(:,2,1) = [0.,1.]
!!$       data(1)%map(:,3,1) = [1.,0.]
!!$       data(1)%map(:,2,2) = [0.,2.]
!!$       data(1)%map(:,3,2) = [2.,0.]

!!$       data(1)%map(:,2,1) = [1.,1.]
!!$       data(1)%map(:,3,1) = [1.,1.]
!!$       data(1)%map(:,2,2) = [2.,2.]
!!$       data(1)%map(:,3,2) = [2.,2.]

       !call the main functions
       call setup_maxlike_mod(data,k)
       p(1)=-1.
       call gridmin(p, uncert, dir)
!       call get_region_A(p,offsets,Areg)
       do r=1,nreg
          call calculate_Am(p,A(:,r),loops_offset(i,r,:), r)
       end do
       call terminate_maxlike_mod

       loops_beta(i)=p(1)
    end do


!loop end here
!    write(*,*) "offset, m"
!    call dump_matrix(m)
!    write(*,*) "amplitude"
!    call dump_matrix(A)
!    write(*,*)"maps",smap0,smap1
!    write(*,*) "offsets", offsets(1)
!    call dump_matrix(offsets)
!    write(*,*) "amplitude2"
!    call dump_matrix(Areg)

    call dump_matrix(loops_beta,"spectind.txt",fmt="(e15.7)")
    call dump_matrix(reshape(loops_offset,[size(loops_offset,1),size(loops_offset(1,:,:))]),"offset.txt",fmt="(e15.7)")
    
    beta_mean = mean(loops_beta)
    beta_dev  = sqrt(sum((loops_beta-beta_mean)**2))
    write(*,*) "beta = ", beta_mean, " +- ",beta_dev

    do r = 1, nreg
       do i = 1, ncomp
          off_mean = mean(loops_offset(:,r,i))
          off_dev  = sqrt(sum((loops_offset(:,r,i)-off_mean)**2))
          write(*,'(a,i4,a,i1,a,e15.7,a,e15.7)') "offset reg ", r, " comp ", i, " : ", off_mean, " +- ", off_dev
       end do
    end do
    write(*,*) "ferdig check_m"

  end subroutine command_check_m


  subroutine setup2(bfile, nside, order, pixels, pinfo, weights)
    implicit none
    character(len=*) :: bfile
    integer(i4b)     :: nside, order, pixels(:), i
    type(pixinfo)    :: pinfo
    real(dp), dimension(:,:), allocatable :: weights, beam
    real(dp)         :: maxrad
    call read_beam_real(bfile, beam)
    maxrad = radial_to_rmax(beam(:,1), beam(:,2))
    call init_pixinfo(pinfo, nside, order, maxrad, pixels)
    call alloc_weights(pinfo, weights)
    call calc_weights_radial(pinfo, beam(:,1), beam(:,2), weights, normalize=.true.)
    deallocate(beam)
  end subroutine


  subroutine command_nsidemap
    implicit none
    character(len=512) :: nside_char, outfile
    integer(i4b)       :: nside, order,npix,i
    real(dp),     dimension(:,:), allocatable :: map

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, nside_char); read (unit=nside_char,fmt=*) nside
    call getarg(3, outfile)

    npix = 12*nside**2
    allocate(map(0:npix-1,1))
    do i=0,npix-1
       map(i,1)=i+1
    end do
    order=1
    call write_map(map, order, outfile)

  end subroutine command_nsidemap



  subroutine command_fitsmake1map
    implicit none

    character(len=512) :: infile, outfile
    integer(i4b)       :: nside, order, i, tall
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp),     dimension(:,:), allocatable :: map, map2

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)

    call read_map(map, pixels, nside, order, infile)

!    call write_map(map(:,2), pixels, nside, order, outfile)
    map=map*1000d0
    map=map*1000d0
    call write_map(map(:,1:3), pixels, nside, order, outfile)

    deallocate(map, pixels)
  end subroutine


  subroutine command_readwritemap
    implicit none

    character(len=512) :: infile, outfile
    integer(i4b)       :: nside, order, npix
    real(dp),     dimension(:,:), allocatable :: map
    integer(i4b), dimension(:),   allocatable :: pixels
    real (dp)           :: nullval
    logical (lgt)       :: anynull

!    npix = getsize_fits("sim_K23_30.fits", nside=nside, ordering=ordering)
!    allocate(map(0:npix-1,1:nmap))
!  call read_bintab( trim(wmapfile), mapwmap, npix, nmap, nullval, anynull )

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)

    npix = getsize_fits(trim(infile), nside=nside, ordering=order)
    allocate(map(0:npix-1,1:3))
    call read_bintab( trim(infile), map, npix, 3, nullval, anynull )
    pixels = irange(npix)-1
!    call read_map(map, order, infile)
    call write_map(map, order, outfile)
    deallocate(map)


  end subroutine command_readwritemap


  subroutine command_covmat
    implicit none

    character(len=512) :: infile, outfile, intext, covfile
    integer(i4b)       :: nside, order, i, j, seed, unit, k
    real(dp)           :: nu, sigma0
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: map, matrix, L, Vn, eigen_mat
    real(dp),     dimension(:), allocatable   :: w, eigen
    real(dp),     dimension(:,:,:,:), allocatable :: cov
    logical(lgt)                                        :: inv

!    real(dp)           :: n
    !external random
    type(planck_rng)     :: rng_handle

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, covfile)

    call read_map(map, pixels, nside, order, infile)
    allocate(cov(size(map,1),2,size(map,1),2))  ! 2 is pol
    call read_covmat(cov, pixels, nside, order, covfile, inv)
    allocate(L(size(map,1), size(map,1)))
    allocate(matrix(size(map,1),size(map,1)))
    allocate(w(size(map,1)))
    allocate(eigen_mat(size(map,1),size(map,1)))
    allocate(eigen(0:size(map,1)-1))
    allocate(Vn(0:size(map,1)-1, 0:size(map,1)-1))
    do j = 1,2
       matrix(:,:) = cov(:,j,:,j) 
!       call dump_matrix(matrix)
!write(*,*) "cov diag", matrix(1,1), matrix(2,2), matrix(3,3), matrix(4,4), matrix(5,5)

       !check symmetric matrix
       do i = 1,size(matrix,1)
write(*,*) "cov diag",matrix(i,i)
          do k = 1, size(matrix,2) 
             if (abs((matrix(i,k) - matrix(k,i))/sqrt(matrix(i,i)*matrix(k,k)) ) > 1e-12 ) then
                write(*,*)"ikke sym matrise", matrix(i,k), matrix(k,i)
             end if
          end do
       end do
       !add factor on the diagonal
!!$       do i = 1,size(matrix,1)
!!$          matrix(i,i) = matrix(i,i) + 0.1d0
!!$       end do
       !get the eigenvalues
       call get_eigen_decomposition(34, matrix, eigen, Vn)
       unit = getlun()
       open(unit, file="eigen.txt")
       eigen_mat = 0.d0
       do i=0, size(map,1)-1
          write(unit,*) eigen(i)
          if(eigen(i).LT.10**-12) eigen(i) = 0.d0
          eigen_mat(i,i) = eigen(i)
       end do
       close(unit)

!!$       matrix = matmul(Vn,matmul(eigen_mat,transpose(Vn)))
!!$       do i=0, size(map,1)-1
!!$          matrix(i,i) = matrix(i,i) + 0.1d0
!!$          w(i)=  rand_gauss(rng_handle)
!!$       end do
       call cholesky_decompose(matrix, L)
       map(:,j) = map(:,j) + matmul(L,w)
       write(*,*) size(map,1)
!       cov(:,j,:,j) = matrix(:,:) 
    end do

!    call write_covmat(cov,pixels,nside,order,cov_out.)
    deallocate(cov,L,matrix,w, eigen_mat, eigen, Vn)

    !!

    write(*,*)"ferdig m loop", nu, sigma0
    write(*,*)"ferdig m loop", ant2thermo(nu)/ant2thermo(23.d0) * (nu/23.d0)**-3.d0 
!    call write_map(map, pixels, nside, order, outfile)
    deallocate(map, pixels)
  end subroutine command_covmat

  subroutine command_addcovnoise
    implicit none

    character(len=512) :: infile, outfile, intext, covfile
    integer(i4b)       :: nside, order, i, j, seed, n, ncomp
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp),     dimension(:,:), allocatable :: map, L
    real(dp),     dimension(:), allocatable   :: w, W_out
    real(dp),     dimension(:,:), allocatable :: cov

    !external random
    type(planck_rng)     :: rng_handle

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, covfile)
    call getarg(5, intext)
    read(intext,*) seed


    !read in map
    call read_map(map, pixels, nside, order, infile)
    map = map !/ 100000.
    n = size(map,1)
    call rand_init(rng_handle,seed)
write(*,*)"map",map(1:20,2)
    !gaussian noise from covar matrix
    call read_covmat(cov,  pixels2, nside, order, ncomp, trim(covfile))
    allocate(W_out(n))
!    call get_eigenvalues(cov, W_out)
!    call dump_matrix(W_out,"w1.dat")
    !call invert_singular_matrix(cov,1.d-12)
!    call get_eigenvalues(cov, W_out)
!    call dump_matrix(W_out,"w2.dat")
    allocate(L(2*n,2*n))
    allocate(w(2*n))
    do i=1, 2*n
       w(i)=  rand_gauss(rng_handle)
    end do
!    call compute_hermitian_power(0.5d0, 1d-12, cov,L, W_out)
    call cholesky_decompose(cov, L)
    !       if (j==2) call cholesky_decompose(cov(n+1:n*2,n+1:n*2), L)
    w=matmul(L,w)
    map(:,2) = map(:,2) + w(1:n)
    map(:,3) = map(:,3) + w(n+1:2*n)
    write(*,*) size(map,1)!,W_out
!    call dump_matrix(abs(W_out),"w.dat")

!    call write_covmat(cov,pixels,nside,order,cov_out.)
    deallocate(cov,L,w, pixels2)
write(*,*)"map",map(1:20,2)
    !!

    write(*,*)"ferdig"
    call write_map(map, pixels, nside, order, outfile)
    deallocate(map, pixels)
  end subroutine command_addcovnoise




  subroutine command_simulate
    implicit none

    character(len=512) :: infile, outfile, intext, covfile
    integer(i4b)       :: nside, order, i, j, seed, unit, k, n
    real(dp)           :: nu, sigma0, chi
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: map, matrix, imatrix, L, Vn, eigen_mat
    real(dp),     dimension(:), allocatable   :: w, eigen
    real(dp),     dimension(:,:,:,:), allocatable :: cov, icov
    logical(lgt)                                        :: inv

!    real(dp)           :: n
    !external random
    type(planck_rng)     :: rng_handle

    if(iargc() /= 7) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, intext)
    read(intext,*) nu
    call getarg(5, intext)
    read(intext,*) sigma0
    call getarg(6, intext)
    read(intext,*) seed
    call getarg(7, covfile)

    !read in map
    call read_map(map, pixels, nside, order, infile)
    n = size(map,1)

    call rand_init(rng_handle,seed)
    map = map * ant2thermo(nu)/ant2thermo(23.d0) * (nu/23.d0)**-3.d0 

!!$    !gaussian uniform noise
!!$    do i = 1, size(map,1)
!!$       do j = 1, size(map,2)
!!$          map(i,j) = map(i,j) + sigma0 * rand_gauss(rng_handle)
!!$       end do
!!$!    imatrix(i,i) = 1.0 / (sigma0*sigma0)
!!$    end do

    !gaussian noise from covar matrix
    allocate(cov(size(map,1),2,size(map,1),2))  ! 2 is pol
    call read_covmat(cov, pixels, nside, order, covfile, inv)
    allocate(L(size(map,1), size(map,1)))
    allocate(matrix(size(map,1),size(map,1)))
    allocate(w(size(map,1)))
    allocate(eigen_mat(size(map,1),size(map,1)))
    allocate(eigen(0:size(map,1)-1))
    allocate(Vn(0:size(map,1)-1, 0:size(map,1)-1))
    do j = 1,2
       matrix(:,:) = cov(:,j,:,j) 
!       call dump_matrix(matrix)

!!$       call get_eigen_decomposition(34, matrix, eigen, Vn)
!!$       unit = getlun()
!!$       open(unit, file="eigen.txt")
!!$       eigen_mat = 0.d0
!!$       do i=0, size(map,1)-1
!!$          write(unit,*) eigen(i)
!!$          if(eigen(i).LT.10**-12) eigen(i) = 0.d0
!!$          eigen_mat(i,i) = eigen(i)
!!$       end do
!!$       close(unit)
!!$       matrix = matmul(Vn,matmul(eigen_mat,transpose(Vn)))
       do i=1, size(map,1)
!          matrix(i,i) = matrix(i,i) + 0.1d0
          w(i)=  rand_gauss(rng_handle)
       end do
       call cholesky_decompose(matrix, L)
       map(:,j+1) = map(:,j+1) + matmul(L,w) + 30.  !offset
!       map(:,j+1) =  matmul(L,w)
       write(*,*) size(map,1)
!       cov(:,j,:,j) = matrix(:,:) 
    end do

!    call write_covmat(cov,pixels,nside,order,cov_out.)
    deallocate(cov,L,matrix,w, eigen_mat, eigen, Vn)

    !!

    write(*,*)"ferdig m loop", nu, sigma0
    write(*,*)"ferdig m loop", ant2thermo(nu)/ant2thermo(23.d0) * (nu/23.d0)**-3.d0 
    call write_map(map, pixels, nside, order, outfile)
    deallocate(map, pixels)
  end subroutine command_simulate



  subroutine command_simulate_rms
    implicit none

    character(len=512) :: infile, outfile, intext, rmsfile
    integer(i4b)       :: nside, order, i, j, seed, unit, k, n
    real(dp)           :: nu, sigma0, chi
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp),     dimension(:,:), allocatable :: rmsmap, map, matrix, imatrix, L, Vn, eigen_mat
    real(dp),     dimension(:), allocatable   :: w, eigen
    real(dp),     dimension(:,:,:,:), allocatable :: cov, icov
    logical(lgt)                                        :: inv

!    real(dp)           :: n
    !external random
    type(planck_rng)     :: rng_handle

    if(iargc() /= 7) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, intext)
    read(intext,*) nu
    call getarg(5, intext)
    read(intext,*) sigma0
    call getarg(6, intext)
    read(intext,*) seed
    call getarg(7, rmsfile)

    !read in map
    call read_map(map, pixels, nside, order, infile)
    n = size(map,1)

    call rand_init(rng_handle,seed)
    map = map * ant2thermo(nu)/ant2thermo(23.d0) * (nu/23.d0)**-3.d0 

!!$    !gaussian uniform noise
!!$    do i = 1, size(map,1)
!!$       do j = 1, size(map,2)
!!$          map(i,j) = map(i,j) + sigma0 * rand_gauss(rng_handle)
!!$       end do
!!$!    imatrix(i,i) = 1.0 / (sigma0*sigma0)
!!$    end do

    call read_map(rmsmap, pixels2, nside, order, rmsfile)
    allocate(L(size(map,1), size(map,1)))
    allocate(matrix(size(map,1),size(map,1)))
    allocate(imatrix(size(map,1),size(map,1)))
    allocate(w(size(map,1)))
    allocate(eigen_mat(size(map,1),size(map,1)))
    allocate(eigen(0:size(map,1)-1))
    allocate(Vn(0:size(map,1)-1, 0:size(map,1)-1))
    do j = 1,2
       matrix = 0.d0
       imatrix = 0.d0
       do i=1,size(rmsmap,1)
          matrix(i,i) = rmsmap(i,j+1)**2
          imatrix(i,i) = 1.d0 / (rmsmap(i,j+1)**2)
       end do
!       call dump_matrix(matrix)

!!$       call get_eigen_decomposition(34, matrix, eigen, Vn)
!!$       unit = getlun()
!!$       open(unit, file="eigen.txt")
!!$       eigen_mat = 0.d0
!!$       do i=0, size(map,1)-1
!!$          write(unit,*) eigen(i)
!!$          if(eigen(i).LT.10**-12) eigen(i) = 0.d0
!!$          eigen_mat(i,i) = eigen(i)
!!$       end do
!!$       close(unit)
!!$       matrix = matmul(Vn,matmul(eigen_mat,transpose(Vn)))
       do i=1, size(map,1)
!          matrix(i,i) = matrix(i,i) + 0.1d0
          w(i)=  rand_gauss(rng_handle)
       end do
       call cholesky_decompose(matrix, L)
       map(:,j+1) = map(:,j+1) + matmul(L,w)
!       map(:,j+1) =  matmul(L,w)
       write(*,*) size(map,1)
!       cov(:,j,:,j) = matrix(:,:) 

       !calculate and write chisquare
       chi = dot_product( matmul(map(:,j+1), imatrix) ,map(:,j+1) )
       write(*,*) "chi", chi

    end do

!    call write_covmat(cov,pixels,nside,order,cov_out.)
    deallocate(L,matrix,w, eigen_mat, eigen, Vn)

    !!

    write(*,*)"ferdig m loop", nu, sigma0
    write(*,*)"ferdig m loop", ant2thermo(nu)/ant2thermo(23.d0) * (nu/23.d0)**-3.d0 
    call write_map(map, pixels, nside, order, outfile)
    deallocate(map, pixels)
  end subroutine command_simulate_rms


  subroutine command_simulate_map
    implicit none

    character(len=512) :: infile, outfile, intext
    integer(i4b)       :: nside, order, i, j
    real(dp)           :: nu, offset, nu0
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: map

    if(iargc() /= 5) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, intext)
    read(intext,*) nu
    call getarg(5, intext)
    read(intext,*) offset

    nu0=22.45d0

    !read in map
    call read_map(map, pixels, nside, order, infile)

    !scale with the freq factor (for example to go from K to Ka)
    map = map * ant2thermo(nu)/ant2thermo(nu0) * (nu/nu0)**-3.d0 

    !add offset to map
    map = map +offset
    write(*,*) size(map,1)

    write(*,*)"scale factor", ant2thermo(nu)/ant2thermo(nu0) * (nu/nu0)**-3.d0 
    write(*,*)"offset", offset
    call write_map(map, pixels, nside, order, outfile)
    write(*,*)"done!"
    deallocate(map, pixels)
  end subroutine command_simulate_map


 subroutine command_check2
    implicit none
    character(len=512)        :: mapfile, covfile
    integer(i4b)              :: n, nside, order, ncomp
    real(dp)                  :: chisq
    real(dp),     allocatable :: map(:,:), cov(:,:), icov(:,:), qumap(:)
    integer(i4b), allocatable :: pixels(:), pixels2(:)
    call getarg(2, mapfile)
    call getarg(3, covfile)
    call read_map(map, pixels, nside, order, mapfile)
    call read_covmat(cov, pixels2, nside, order, ncomp, covfile)
    n = size(map,1)
    write(*,*) "A", n, shape(cov)
    allocate(icov(2*n,2*n), qumap(2*n))
    icov = cov
    write(*,*) cov(1,1), cov(2,2), cov(3,3), cov(4,4)
    call invert_matrix(icov)
    write(*,*) icov(1,1), icov(2,2), icov(3,3), icov(4,4)
    qumap = reshape(map(:,2:3),(/2*n/))
    chisq = dot_product(matmul(qumap, icov), qumap)
    write(*,*) "chiQ+U", chisq
  end subroutine

  subroutine command_check
    implicit none

    character(len=512) :: infile, outfile, intext, covfile
    integer(i4b)       :: nside, order, i, j, seed, unit, k, n, ncomp
    real(dp)           :: nu, sigma0, chi
    integer(i4b), dimension(:),   allocatable :: pixels, qumap
    real(dp),     dimension(:,:), allocatable :: map, matrix, imatrix, L, Vn, one
    real(dp),     dimension(:),   allocatable :: w, eigen
    real(dp),     dimension(:,:), allocatable :: cov, icov
    logical(lgt)                                        :: inv

!    real(dp)           :: n
    !external random
    type(planck_rng)     :: rng_handle

    if(iargc() /= 7) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)
    call getarg(4, intext)
    read(intext,*) nu
    call getarg(5, intext)
    read(intext,*) sigma0
    call getarg(6, intext)
    read(intext,*) seed
    call getarg(7, covfile)

    !read in map
    call read_map(map, pixels, nside, order, infile)
    n = size(map,1)

    !read in covmat
!    allocate(cov(size(map,1),2,size(map,1),2))  ! 2 is pol

    deallocate(pixels)
    call read_covmat(cov, pixels, nside, order, ncomp, covfile, inv)
    write(*,*) "A", n, shape(cov)

    !Q
    allocate(matrix(n,n))
    allocate(imatrix(n,n))
    matrix = cov(:n,:n)
    imatrix = matrix
    call invert_matrix(imatrix)
    chi = dot_product(matmul(map(:,2), imatrix), map(:,2) )
    write(*,*) "chiQ", chi
    deallocate(matrix,imatrix)

    !U
    allocate(matrix(n,n))
    allocate(imatrix(n,n))
    matrix = cov(n+1:,n+1:)
    imatrix = matrix
    call invert_matrix(imatrix)
    chi = dot_product(matmul(map(:,3), imatrix), map(:,3) )
    write(*,*) "chiU", chi
    deallocate(matrix,imatrix)

    !Q+U
    allocate(matrix(2*n,2*n))
    allocate(imatrix(2*n,2*n))
    matrix = cov
    imatrix = matrix
    call invert_matrix(imatrix)
!    call invert_singular_matrix(imatrix,0.d0)
    !create a stacked map
    allocate(qumap(2*n))
    qumap = reshape(map(:,2:3),(/2*n/))
    chi = dot_product(matmul(qumap, imatrix), qumap)
    write(*,*) "chiQ+U", chi
    deallocate(matrix,imatrix)

    deallocate(map, pixels)
  end subroutine command_check


  subroutine command_plothistogram
    implicit none

    character(len=512) :: infile, uncertfile, outfile, out2file, intext
    integer(i4b)       :: nside, order, i, j, k, unit, unit2, jtot
    real(dp)           :: beta, theta, phi, meanbeta
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2
    real(dp),     dimension(:,:), allocatable :: bmap, umap
    real(dp),     dimension(:),   allocatable :: gmlbeta
    logical (lgt)       :: funn

    if(iargc() /= 6) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, uncertfile)
    call getarg(4, outfile)
    call getarg(5, out2file)
    call getarg(6, intext)
    read(intext,*) beta
    call read_map(bmap, pixels, nside, order, infile)
    call read_map(umap, pixels2, nside, order, uncertfile)
    k = 1
    jtot = 1
    unit = getlun()
    open(unit, file=trim(outfile))
    unit2 = getlun()
    open(unit2, file=trim(out2file))
    allocate(gmlbeta(size(bmap,1)))
    gmlbeta = 0.d0
    do i = 1, size(bmap,1)
       do j=1,jtot
          if(bmap(i,k).EQ.gmlbeta(j)) funn=.true.
       end do
       if (funn == .false.) then
          !make histogram
          write(unit,*) (bmap(i,k)-beta ) / umap(i,k)
!          write(*,*) (bmap(i,k)-beta ) / umap(i,k)
          gmlbeta(jtot) = bmap(i,k)
          jtot=jtot+1
          !make the beta vs lat plot
          call pix2ang_nest(nside,pixels(i),theta,phi) !phi=[0,2pi]
          phi = phi*RAD2DEG
          theta = theta*RAD2DEG
          write(unit2,*) bmap(i,k),  theta - 90.d0
!          write(*,*) bmap(i,k), theta - 90.d0
       end if
       funn=.false.
    end do
    close(unit)
    close(unit2)

    !some calculations on the beta-values
    meanbeta =  sum(gmlbeta(1:jtot-1))/(jtot-1)
    write(*,*) "mean value of spectral index", meanbeta, jtot-1
    write(*,*) "standard deviation", sqrt( (jtot-1)/(jtot-2) * (meanbeta-beta)**2 )

    deallocate(bmap, umap, pixels, pixels2, gmlbeta)
  end subroutine command_plothistogram

  subroutine command_histogram
    implicit none

    character(len=512) :: infile, outfile
    integer(i4b)       :: nside, order, i, unit
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: bmap

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, outfile)

    call read_map(bmap, pixels, nside, order, infile)
    unit = getlun()
    open(unit, file=trim(outfile))
    do i = 1, size(bmap,1)
       !make histogram
       ! Dont print masked pixels (allow for some inaccuracy of the badval)
       if (bmap(i,2).gt.-1d30) then
          write(unit,*) bmap(i,1), bmap(i,2), bmap(i,3)
       end if
    end do
    close(unit)

    deallocate(bmap, pixels)
  end subroutine command_histogram

 subroutine command_makedatfile
    implicit none

    character(len=512) :: infile1, infile2, outfile
    integer(i4b)       :: nside, order, i, j, k, unit
    real(dp)           :: lat,lon,vec0(3),vec(3)
    integer(i4b), dimension(:),   allocatable :: pixels, pixels2,sort_mask
    real(dp),     dimension(:,:), allocatable :: amap, bmap
    real(dp),     dimension(:),   allocatable :: y,dist,tmpdist

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile1)
    call getarg(3, infile2)
    call getarg(4, outfile)

    call read_map(amap, pixels, nside, order, infile1)
    call read_map(bmap, pixels2, nside, order, infile2)

    lon =(-9.09968) *DEG2RAD
    lat =(90.d0 +32.63895) *DEG2RAD
    write(*,*)DEG2RAD
    unit = getlun()
    open(unit, file=trim(outfile))
    allocate(y(size(amap,1)),dist(size(amap,1)), sort_mask(size(amap,1)),tmpdist(size(amap,1)))
    y=amap(:,1)

    call ang2vec(lat,lon,vec0)
    do i = 1, size(amap,1)
       call pix2vec_ring(nside,pixels(i),vec(:))
       call angdist(vec,vec0,dist(i))
    end do

    tmpdist=dist
    do i=1,size(amap,1)
       sort_mask(i)=minloc(tmpdist,1)
       tmpdist(sort_mask(i))=10.d0
    end do
    dist=dist(sort_mask)
    y=y(sort_mask)
    
    do i=1,size(amap,1)
       write(unit,*) dist(i), y(i)       
    end do
    close(unit)

    deallocate(amap, bmap, pixels, pixels2, y, dist)
  end subroutine command_makedatfile

 subroutine command_findlatitude
    implicit none

    character(len=512) :: infile, indatfile, outfile
    integer(i4b)       :: nside, order, i, j, k, unit, unit2, teller, reg
    real(dp)           :: lat,lon,reg_vec(3),vec(3), region_lat(24), region_lon(24),beta,uncert
    integer(i4b), dimension(:),   allocatable :: pixels, sort_mask
    real(dp),     dimension(:,:), allocatable :: amap
    real(dp),     dimension(:),   allocatable :: y,dist,tmpdist

    if(iargc() /= 4) then
       write(*,*) "Wrong number of arguments!"
!       call command_help
       return
    end if
    call getarg(2, infile)
    call getarg(3, indatfile)
    call getarg(4, outfile)

    call read_map(amap, pixels, nside, order, infile)

    do j=1,24
       reg_vec=0.d0
       teller=0
       do i=1,size(amap,1)
          if (amap(i,1).eq.j) then
             call pix2vec_nest(nside,pixels(i),vec(:))
             reg_vec=reg_vec+vec
             teller=teller+1
          end if
       end do
       vec=reg_vec/teller
       call vec2ang(vec,lat,lon)
       region_lat(j)=lat*RAD2DEG
       region_lon(j)=lon*RAD2DEG
       if(region_lon(j).gt.180.) region_lon(j)=region_lon(j)-360.d0 
       write(*,*) j,region_lat(j),region_lon(j)
    end do
    
    unit = getlun()
    open(unit,file=trim(indatfile),action="read",status="old")
    unit2= getlun()
    open(unit2,file=trim(outfile))
    do j=1,24
       read(unit,*)reg,beta,uncert
       write(unit2,*)reg,beta,uncert,region_lat(reg),region_lon(reg)
    end do
    close(unit)
    close(unit2)

!!$
!!$    lon =(-9.09968) *DEG2RAD
!!$    lat =(90.d0 +32.63895) *DEG2RAD
!!$    write(*,*)DEG2RAD
!!$    unit = getlun()
!!$    open(unit, file=trim(outfile))
!!$    allocate(y(size(amap,1)),dist(size(amap,1)), sort_mask(size(amap,1)),tmpdist(size(amap,1)))
!!$    y=amap(:,1)
!!$
!!$    call ang2vec(lat,lon,vec0)
!!$    do i = 1, size(amap,1)
!!$       call pix2vec_ring(nside,pixels(i),vec(:))
!!$       call angdist(vec,vec0,dist(i))
!!$    end do
!!$
!!$    tmpdist=dist
!!$    do i=1,size(amap,1)
!!$       sort_mask(i)=minloc(tmpdist,1)
!!$       tmpdist(sort_mask(i))=10.d0
!!$    end do
!!$    dist=dist(sort_mask)
!!$    y=y(sort_mask)
!!$    
!!$    do i=1,size(amap,1)
!!$       write(unit,*) dist(i), y(i)       
!!$    end do
!!$    close(unit)

    deallocate(amap, pixels)
  end subroutine command_findlatitude

 subroutine command_calculate_rms_value
    implicit none

    character(len=512) :: infile, intext
    integer(i4b)       :: nside, order, i, unit, comp
    real(dp)           :: rms, mu
    integer(i4b), dimension(:),   allocatable :: pixels
    real(dp),     dimension(:,:), allocatable :: amap
    real(dp),     dimension(:),   allocatable :: y2

    if(iargc() /= 3) then
       write(*,*) "Wrong number of arguments!"
!       call command_help                                                                                                        
       return
    end if
    call getarg(2, infile)
    call getarg(3, intext)
    read(intext,*) comp

    call read_map(amap, pixels, nside, order, infile)
    allocate(y2(size(amap,1)))

!    y2=amap(:,comp)*amap(:,comp)
!    rms=sqrt( (SUM(y2))/ real(size(amap,1),dp) )
!    write(*,*) "The RMS value of ",trim(infile),"component ",comp," is: ",rms

    mu = sum(amap(:,comp)) / real(size(amap,1),dp)
    rms = sqrt(sum((amap(:,comp)-mu)**2)/real(size(amap,1),dp))
    write(*,*) "The RMS error value of ",trim(infile),"component ",comp," is: ",rms

    deallocate(amap, pixels,y2)
  end subroutine command_calculate_rms_value

 subroutine command_calculate_sensitivity
    implicit none

    character(len=512) :: covfile
    integer(i4b)       :: nside, ordering, i, unit, polarization
    real(dp)           :: rms, mu, sensitivity, diag
    integer(i8b)                              :: n
    real(dp),     dimension(:,:), allocatable :: cov
    logical(lgt)                              :: inv

    if(iargc() /= 2) then
       write(*,*) "Wrong number of arguments!"
!       call command_help                                                                                                        
       return
    end if
    call getarg(2, covfile)

    unit=42
    call read_covmatrix(unit, trim(covfile), ordering, polarization, cov, inv, n)
    sensitivity=sum(cov(:,:))
    diag=0.d0
    do i=1,n
       diag=diag+cov(i,i)
    end do

    write(*,*) "The overall sum of covmatrix ",trim(covfile)," is: ",sensitivity
    write(*,*) "The sum of the diagonal elements: ",diag, sqrt(1.d0/diag)

!    !optional write out of cov matrix
!    do i=1,n
!       write(*,*) cov(i,i)
!    end do

    deallocate(cov)
  end subroutine command_calculate_sensitivity


end program sindex_utils
