program test_oslo

use ds_types
use ds_fitstools
use ds_multidetector
use ds_solver
use inifile
use mpi
use ds_types
use ds_utils
use quiet_utils
use l2_fileutils
use quiet_db_mod
use quiet_fileutils
use quiet_calib_mod
use quiet_noise_mod
use ds_oslo_fitstools
implicit none

type(ds_detpointing), dimension(:), allocatable :: pointing
type(ds_timestream), dimension(:), allocatable :: timestream
integer rank,nproc
logical CES
character(len=128) :: level2_dir, target_name, target_type
character(len=128) :: parfile_tod2map, parfile_tod_proc, database_file, parfile_common

character(len=256) ::  module_list_file
integer ierror, unit,nside
integer(i8b), dimension(:), allocatable :: originalIndices
real(dp), allocatable, dimension(:) :: ntod


call MPI_Init(ierror)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierror)
call mpi_comm_size(MPI_COMM_WORLD,nproc,ierror)


  call getarg(1,parfile_common)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD2MAP',   par_string=parfile_tod2map)
  call get_parameter(unit, parfile_common, 'PARFILE_TOD_PROC',  par_string=parfile_tod_proc)
  call get_parameter(unit, parfile_common,  'SCAN_DATABASE',        par_string=database_file)
  call get_parameter(unit, parfile_common,  'LEVEL2_DIR',           par_string=level2_dir)
  call get_parameter(unit, parfile_tod2map, 'TARGET_TYPE',          par_string=target_type)
  call get_parameter(unit, parfile_tod2map, 'TARGET_NAME',          par_string=target_name)
  call get_parameter(unit, parfile_common,  'MODULE_LIST',          par_string=module_list_file)
  call get_parameter(unit, parfile_descart, 'NSIDE_OUT',            par_int=nside)

  ! Get parameters

  call initialize_quiet_noise_mod(unit, parfile_tod2map, parfile_tod_proc)
  call initialize_quiet_calib_mod(unit, parfile_common)
  call initialize_quiet_db_mod(unit, database_file, level2_dir)


!JZ this should also 
call readDataAssignWork(pointing,timestream,rank,nproc,CES,level2_dir, target_name, target_type, module_list_file,nside,correlator)

my_nd = size(timestream)
my_np = size(pointing)
call ds_assert(my_nd==my_np*2, "Wrong pair-det relation")
call prepareFFTQ(correlator,offsetLength,ntod,ndet,MPI_COMM_WORLD)


call makeSharedPointingContinuous(correlator,pointing,12*nside*nside,pixelFile,originalIndices)
npix = pointing(0)%npix
if (rank==0) call ds_log("Shared pointing complete.",ds_feedback_quiet)





!Prepare the maps.  We use these both as dummies and as the final output maps
call prepareMap(qmap,npix)
call prepareMap(umap,npix)
if(rank==0) then
	call ds_assert(size(qmap%indices)==size(originalIndices),"map size not equal to hit pixels size")
	qmap%indices=originalIndices
	umap%indices=originalIndices
	deallocate(originalIndices)
endif
call MPI_Barrier(correlator%comm,ierror)


allocate(tempTimestream(0:my_nd-1))
allocate(offsetSolution(0:my_nd-1))
allocate(offsetTarget(0:my_nd-1))




!Build the RHS of the equation.
if (rank==0) call ds_log("Building RHS.",ds_feedback_quiet)
do d=0,my_nd-1
   call copyTimestream(timestream(d),tempTimestream(d))
   call prepareOffsets(offsetSolution(d),offsetLength)
   call prepareOffsets(offsetTarget(d),offsetLength)

enddo
if (rank==0) call ds_log("Removing Naive Signal",ds_feedback_quiet)
call removeSignalNaive(tempTimestream,pointing,qmap,umap,correlator)
do d=0,my_nd-1
   call deprojectTimestreamOntoOffset(tempTimestream(d),offsetTarget(d))
   offsetTarget(d)%values= offsetTarget(d)%values / whiteNoise(correlator%my_det(d))
   call destroyTimestream(tempTimestream(d))
enddo


!Run the PCG
if (rank==0) call ds_log("Ready for PCG.", ds_feedback_quiet)
call PCG(correlator, pointing,regularisationParameter * whiteNoise, npix, nt, offsetTarget, offsetSolution)
if (rank==0) call ds_log("Completed PCG",ds_feedback_quiet)


!Deduct the offsets we have computed from the timestreams
call prepareTimestream(tempTimestream(0),nt)
do d=0,my_nd-1
   call projectOffsetOntoTimestream(offsetSolution(d),tempTimestream(0))
   do i= 1,timestream(0)%nt
      timestream(d)%timestream(i) = timestream(d)%timestream(i) - tempTimestream(0)%timestream(i)
   enddo
enddo
call destroyTimestream(tempTimestream(0))

if (rank==0) call ds_log("subtracted offsets", ds_feedback_quiet)


!Make naive maps from the new white timestreams
call makeNaiveMap(timestream,qmap,umap,pointing,correlator)


!Save the output maps
if (correlator%proc == 0) then
   call savePolMaps(outputName, qmap,umap,nside)
endif


end program test_oslo
