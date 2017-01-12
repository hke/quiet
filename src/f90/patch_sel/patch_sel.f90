program patch_sel
  use rms_mod
  use find_pts_mod
  use plot_patch_mod
  use healpix_types
  use quiet_ephem_mod
  implicit none

  character(len=512)                :: file_in, file_mask, file_rms, file_out, db_path, arg, command
  logical(lgt)                      :: make_map, rms_default

  integer(i4b)                      :: switch, npix_in, ndays, argi
  real(dp)                          :: t1, t2, elevation_limit

  make_map = .true.
  rms_default = .true.

  ! Loop to read optional arguments:
  argi = 0
  do while (argi < iargc())
     argi = argi+1
     call getarg(argi,arg)
     select case(arg)
     case("-db")
        argi = argi+1
        call getarg(argi,db_path)
        call ephem_set_db(trim(db_path))
     case("-h")
        call help
     case("-altrms")
        argi = argi+1
        call getarg(argi, file_rms)
        rms_default = .false.
     case("-nomap")
        make_map = .false.
     case default
        if(arg(1:1) == "-") then
           call help
        end if
        ! If there are no more optional arguments
        exit
     end select
  end do

  ! Loop to read mandatory argument:
  if(argi > iargc()) call help
  call getarg(argi, command)
  select case(command)
  case("rms")
     if(iargc()-argi+1 < 3) call help
     call getarg(argi+1, file_in)
     call getarg(argi+2, file_mask)
      call wall_time(t1)
     call make_rmsmap(file_in, file_mask, file_rms, rms_default)
     call wall_time(t2)
     write(*,*) "Made rms map in ", t2-t1, " seconds"
  
  case("list")
     if(iargc()-argi+1 < 4) call help
     call getarg(argi+1, file_in)
     call getarg(argi+2, arg)
     npix_in = atoi(arg)
     call getarg(argi+3, arg)
     elevation_limit = atof(arg)
     
     call wall_time(t1)
     call find_low_fg_pts(file_in, npix_in, elevation_limit)
     call wall_time(t2)
     write(*,*) "Listed", npix_in, " low-foreground pixels in", t2-t1, " seconds"
  
  case("plot")
     if(iargc()-argi+1 < 4) call help
     call getarg(argi+1, file_in)
     call getarg(argi+2, file_out)
     call getarg(argi+3, arg)
     ndays = atoi(arg)
     
     call wall_time(t1)
     call plot_patch_and_elevation(file_in, file_out, file_rms, ndays, rms_default, make_map)
     call wall_time(t2)
     write(*,*) "Patches plotted and observation site elevation computed in", t2-t1, " seconds"
  
  case default
     call help
  end select

contains

  subroutine help
    implicit none
    write(*,*) "patch_sel: A program for selecting patches"
    write(*,*) " "
    write(*,*) "SYNTAX:"
    write(*,*) "Option 1: ./patch_sel rms fg.fits mask.fits"
    write(*,*) "Provide foreground map and corresponding mask. Program outputs rms of foreground polarization in file rms.fits (unless alternative output filename is given)."
    write(*,*) " "
    write(*,*) "Option 2: ./patch_sel list rms.fits N elevation_limit"
    write(*,*) "Provide the map rms.fits. Program outputs a list of the N pixels with lower foreground that is theoretically visible from the QUIET observation site given a telescope with a lower elevation limit of elevation_limit"
    write(*,*) " "
    write(*,*) "Option 3: ./patch_sel plot list_in.dat list_out.dat N_days"
    write(*,*) "Provide a list of manually chosen pixels in a file providing galactic longitude and latitude in the first two columns."
    write(*,*) "Program outputs files of elevation and azimuth as function of mjd, starting at 1.jan 2012 00:00, for each chosen patch, and (optionally) the chosen patches plotted atop the rms map. It also writes the local RMS at each patch to the patch list file - if input file should not be overwritten, make sure output file has a different name."
    write(*,*) " "
    write(*,*) "Optional arguments (enter before mandatory arguments): "
    write(*,*) "-db filename provides database file for ephem_mod (to compute sun/moon positions in Option 3)"
    write(*,*) "-nomap turns off writing patches to map in option 3"
    write(*,*) "-altrms [filename] specifies a filename for the rms-map which is read or written - default is rms.fits"
    write(*,*) "-h gives program help"
    stop
  end subroutine help

end program patch_sel
