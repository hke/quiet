program descartes_quiet
use ds_types
use ds_fitstools
!use ds_global_pointing
!use ds_focalplane
use ds_multidetector
use ds_solver
use inifile
!use mpi
implicit none

!include 'mpif.h'

type(ds_detpointing), dimension(:), allocatable :: pointing
type(ds_timestream), dimension(:), allocatable :: timestream, tempTimestream
type(ds_correlator) :: correlator
type(ds_offsets), dimension(:), allocatable :: offsetSolution, offsetTarget
type(ds_map) :: qmap, umap

integer na,npix, ierror, status(MPI_STATUS_SIZE)
integer(dp) ::  totPix
integer nside, i, d, nd, p
integer,allocatable,dimension(:) :: nt

!Parameters read fromt he file
integer :: offsetLength, ndet
character(128) :: correlatorFilenameRoot, todName, todroot, outputName, pixelFile, input, destriped_tod_name, pointroot, pointfile, noisefile, rmsnoisemap

!For reading parameters
integer param_lun, narg
logical :: param_status, ignore_correlations
integer ntod
character(128) :: parameter_filename
type(TIniFile) :: param_file
logical, parameter :: FAIL_IF_NO_PARAM = .true.
logical, parameter :: OKAY_IF_NO_PARAM = .false.
real(dp),allocatable,dimension(:) :: whiteNoise
integer :: lun
integer,allocatable,dimension(:) :: originalIndices
integer rank
real(dp) :: datarate
logical :: naive_mode, load_correlation
integer arg

narg = command_argument_count()
call ds_assert(narg > 0, "Program Usage: descarts_quiet parameter_file [more parameter files ... ]")

call MPI_Init(ierror)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierror)
call ds_assert(ierror==0, 'Failed to init mpi in descarts_quiet')
do arg=1,narg





!-------------Read parameters--------------------------
call get_command_argument(arg,parameter_filename)

param_lun = ds_get_lun()
call Ini_Open_File(param_file, parameter_filename, param_lun, param_status, .false.)
call ds_assert(.not. param_status, "Could not open parameter file " // parameter_filename)
todroot = trim( Ini_Read_String_File(param_file, 'tod_filename', FAIL_IF_NO_PARAM) )
outputName = trim( Ini_Read_String_File(param_file, 'output_filename', FAIL_IF_NO_PARAM) )
pixelFile = trim( Ini_Read_String_File(param_file, 'pixel_filename', FAIL_IF_NO_PARAM) )
noiseFile = trim( Ini_Read_String_File(param_file, 'noise_filename', FAIL_IF_NO_PARAM) )
rmsnoisemap= ""
rmsnoisemap = trim( Ini_Read_String_File(param_file, 'rms_noise_map') )

pointroot = trim(ini_read_string_file(param_file, 'pointing_root', FAIL_IF_NO_PARAM) )
destriped_tod_name = trim( Ini_Read_String_File(param_file, 'destriped_tod', OKAY_IF_NO_PARAM) )
correlatorFilenameRoot = Ini_Read_String_File(param_file, 'correlator_filename_root', FAIL_IF_NO_PARAM)
offsetLength = Ini_Read_Int_File(param_file,'offset_length')
naive_mode = Ini_Read_Logical_File(param_file,'naive_mode')
ignore_correlations = Ini_Read_Logical_File(param_file,'ignore_correlations', .false.)
load_correlation = Ini_Read_Logical_File(param_file,'load_correlations', .false.)
ndet = Ini_Read_Int_File(param_file,'number_detectors')
nside = Ini_Read_Int_File(param_file,'nside')
ds_global_feedback = Ini_Read_Int_File(param_file,'feedback',2)
allocate(whiteNoise(0:ndet-1))
!do i=0,ndet-1
!	whiteNoise(i) = Ini_Read_Double_Array_File(param_file,'white_noise', i)
!enddo
if (rank==0) call ds_log("Done reading parameters",ds_feedback_quiet)
!----------------------------------------------------
if (ds_global_feedback >= 2 .and. rank==0) then
!do i=1, param_file%ReadValues%Count
!    write (*,'(a)') trim(param_file%ReadValues%Items(i)%P%Name) // ' = ' //trim(param_file%ReadValues%Items(i)%P%Value)
!end do
endif
call Ini_Close_File(param_file)



 correlator%comm = MPI_COMM_WORLD
call assign_detectors(correlator,ndet)
nd = correlator%my_ndet

!Each proc now has a ready correlator%ndet number and a correlator%my_det array, and they all have the same correlator%owner array

allocate(pointing(0:nd/2-1))
allocate(timestream(0:nd-1))
allocate(tempTimestream(0:nd-1))
allocate(offsetSolution(0:nd-1))
allocate(offsetTarget(0:nd-1))

!read the timestream
do d=0,nd-1
p = correlator%my_det(d)
todname = trim(todroot) // trim(intToStr(p)) // ".fits"
call read_simple_tod_1col(todname, timestream(d))
enddo

if (correlator%proc==0) call ds_log("Read TOD", ds_feedback_quiet)

!read the poinitng.  assume these are numbered so that, e.g. :
!if Q0,U0 are detectors 0 and 1 then pointing is in pointroot.pair_0.fits
!if Q1,U1 are detectors 2 and 3 then pointing is in pointroot.pair_1.fits
!if Q2,U2 are detectors 4 and 5 then pointing is in pointroot.pair_2.fits
do d=0,nd/2-1
	p = correlator%my_det(det_qindex(d)) / 2
	pointfile = trim(pointroot) // "pair_" // trim(intToStr(p)) // ".fits"
	call readDetPointing(pointfile,nside,pointing(d),datarate)
enddo


!We do not load a global pointing as we are only using one detector
!Load detector pointings amd timestream and bits of correlator
!call  READ_QUIET_TOD(todName,pointing,timestream,correlator,noise_only,nside,offsetLength)
if (rank==0) call ds_log("Read Pointing.",ds_feedback_quiet)

!pointing%theta = pointing%theta*DEG2RAD !/ (180.0_8 / pi)

!nt = timestream(0)%nt
allocate(nt(0:ndet-1))
call read_ntod(todname,nt(0))
nt(1:ndet-1)= nt(0)
na = nt(0)/offsetLength

if (load_correlation) then
!	call prepareFFTQ(correlator, na, nt, ndet, correlator%comm, trim(correlatorFilenameRoot), whitenoise, noisefile, ignore_correlations)
else
!	call prepareFFTQ2(correlator, offsetLength, nt, ndet, correlator%comm, trim(correlatorFilenameRoot), whitenoise, noisefile, ignore_correlations)
!  subroutine prepareFFTQ2(correlator,lc,nt,nd,comm,correlatorFilenameRoot,whiteNoise,noisefile,no_correlation_in)
endif
if (rank==0) call ds_log("Multidetector preparation complete.",ds_feedback_quiet)

!this is where tod read should be moved to

totPix = 12*nside*nside
call makeSharedPointingContinuous(correlator,pointing,totPix,pixelFile,originalIndices)
npix = pointing(0)%npix
if (rank==0) call ds_log("Shared pointing complete.",ds_feedback_quiet)


!Output fiducial rms noise map
if(rank==0 .and. rmsnoisemap .ne. '') then
   call writeNoiseRmsFile(rmsnoisemap,pointing,originalIndices,whitenoise,nside)
endif


do d=0,nd - 1
	call prepareOffsets(offsetTarget(d), nt(d)/offsetLength, offsetLength)
	call prepareOffsets(offsetSolution(d), nt(d)/offsetLength, offsetLength)
enddo
call prepareMap(qmap,npix)
call prepareMap(umap,npix)



if(rank==0) then
	call ds_assert(size(qmap%indices)==size(originalIndices),"map size not equal to hit pixels size")
	qmap%indices=originalIndices
	umap%indices=originalIndices
	deallocate(originalIndices)
endif
call MPI_Barrier(correlator%comm,ierror)
if (rank==0) call ds_log("Building RHS.",ds_feedback_quiet)

do d=0,nd-1
	call copyTimestream(timestream(d),tempTimestream(d))
enddo



if (rank==0) call ds_log("Removing Naive Signal",ds_feedback_quiet)
call removeSignalNaive(tempTimestream,pointing,qmap,umap,correlator)

do d=0,nd-1
	call deprojectTimestreamOntoOffset(tempTimestream(d),offsetTarget(d))
	offsetTarget(d)%values= offsetTarget(d)%values / whiteNoise(correlator%my_det(d))
	call destroyTimestream(tempTimestream(d))
enddo

if (naive_mode) then
	if (rank==0) call ds_log("IN NAIVE MODE! NO PCG.",ds_feedback_quiet)
else  !Normal operation mode - PCG
	if (rank==0) call ds_log("Ready for PCG.", ds_feedback_quiet)
!	write(*,*) allocated(pointing), size(pointing), size(whiteNoise), correlator%ndet, rank

	call PCG(correlator, pointing,regularisationParameter * whiteNoise, npix, nt(0), offsetTarget, offsetSolution)

	if (rank==0) call ds_log("Completed PCG",ds_feedback_quiet)

!JZ We only need one timestream for this
	call prepareTimestream(tempTimestream(0),nt(0))
	do d=0,nd-1
		call projectOffsetOntoTimestream(offsetSolution(d),tempTimestream(0))
		do i= 1,timestream(0)%nt
			timestream(d)%timestream(i) = timestream(d)%timestream(i) - tempTimestream(0)%timestream(i)
		enddo
	enddo
	call destroyTimestream(tempTimestream(0))

	if (rank==0) call ds_log("subtracted offsets", ds_feedback_quiet)

endif

!lun = ds_get_lun()
!if (correlator%det == 0) then
!	if (destriped_tod_name .ne. '') then
!		open(unit=lun,file=trim(destriped_tod_name)//'.0')
!		do i=1,timestream%nt
!			write(lun,*) timestream%timestream(i)
!		enddo
!		close(lun)
!	endif
!endif
!if (correlator%det == 1) then
!	if (destriped_tod_name .ne. '') then
!		open(unit=lun,file=trim(destriped_tod_name)//'.1')
!		do i=1,timestream%nt
!			write(lun,*) timestream%timestream(i)
!		enddo
!		close(lun)
!	endif
!endif

!#warning !added code to output 0th detector destriped TOD
!if (rank==0) then
!   call write_simple_tod("/Data/asosx78/point_sources/destriped_tod.fits",0.01_8,512,timestream(0))
!endif


call makeNaiveMap(timestream,qmap,umap,pointing,correlator)

!The master node is a Q node, so it needs to receive the U map from node 1
if (correlator%proc == 0) then
	call savePolMaps(outputName, qmap,umap,nside)
endif


deallocate(whiteNoise)

do d=0,nd-1
	call destroyTimestream(timestream(d))
        call destroyTimestream(tempTimestream(d))
	call destroyOffsets(offsetTarget(d))
	call destroyOffsets(offsetSolution(d))
enddo
do d=0,nd/2-1
	call destroyPointing(pointing(d))
enddo
call destroyMap(qmap)
call destroyMap(umap)
deallocate(timestream)
deallocate(tempTimestream)
deallocate(offsetTarget)
deallocate(offsetSolution)
deallocate(pointing)

deallocate(correlator%owner)
deallocate(correlator%my_det)
	
call destroy_correlator(correlator)

call MPI_Barrier(MPI_COMM_WORLD, ierror)
call ds_assert(ierror == 0, "Error in final mpi_barrier of driver.")




enddo !End loop over parameter files



call MPI_Finalize(ierror)


end program descartes_quiet
