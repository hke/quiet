module ds_oslo_fitstools
use ds_types
use ds_utils
use ds_azimuth
use quiet_utils
use tod2map_utils
use l2_fileutils
use quiet_l2_db_mod
use quiet_fileutils
use quiet_calib_mod
use quiet_noise_mod
use healpix_types
use pix_tools
use ds_multidetector
implicit none
integer, parameter :: ndiodemax = 4
integer, parameter ::npairdiodesmax = 2
integer, parameter :: nmodulemax = 19

type ds_exlist
	integer, allocatable, dimension(:,:,:,:) :: excluded
	integer segmax, runmax
end type

type ds_timestreamInfo
!This concatenates information for PAIRS OF DIODES (Q1,U1) and (Q2,U2)
!The field "pair" says which pair it is.
integer, allocatable, dimension(:) :: run,scan,mod,pair
logical, allocatable, dimension(:) :: needs_azfilter
character(len=128), allocatable, dimension(:) :: filename
integer n
end type

character(len=125),private :: message

contains

subroutine readDataAssignWork(pointing,timestream,correlator,   &
							CES,level2_dir, target_name, target_type, module_list_in,exclusionFile,   &
							nside,alpha,sigma,fknee,   &
							offsetLength,offsets,corrlength,n_bins,   &
							cutfraction)
!JZ This subroutine is getting a bit monolithic.
!Must break it up a little bit.

implicit none
real(dp), optional :: cutfraction
integer offsetLength
integer corrlength
integer n_bins
type(ds_offsets), dimension(:), allocatable :: offsets
type(ds_correlator) :: correlator
type(ds_detpointing), dimension(:), allocatable :: pointing
type(ds_timestream), dimension(:), allocatable :: timestream
character(*) :: exclusionFile
integer rank,nproc
logical(lgt)  CES
character(len=*) :: level2_dir, target_name, target_type
integer nside
real(dp), allocatable, dimension(:) :: alpha,sigma,fknee
real(dp),            pointer, dimension(:,:)            :: az
integer nmodules
integer status
integer ntod
integer ndet
real(dp) :: mu
integer(i4b), pointer,     dimension(:)       :: module_list, module_list_in
integer my_npair, my_ndet, nmin, npair
integer dets_for_this_proc, nextra
integer current_proc
integer modnum
integer f,m,p, pair, q, u, my_pair,i,t, my_q, my_u, my_det,d
character(128) :: current_filename
integer, dimension(:), allocatable ::  ndet_proc
type(module_struct), pointer :: module_data
type(l2_filelist),   pointer, dimension(:) :: target_list
type(module_struct), dimension(:), pointer :: tod_data, tod_data_in
real(dp) :: gain
integer(i4b) :: one, ierr
integer na
real(dp), allocatable, dimension(:) :: azimuth, elevation
real(dp) :: phi,theta,hour_angle
real(dp), parameter :: quiet_latitude = NaN
integer dummy
type(ds_timestreamInfo) timestream_info
integer qdiode, udiode
integer(i4b),        pointer, dimension(:,:,:)          :: pix
real(dp),            pointer, dimension(:,:,:)          :: psi
integer ntarget
integer unit
real(dp) :: delta_el,delta_az,delta_deck
logical goodCalib
logical haveBadScans, thisScanBad
character(*), parameter :: ecliptic='ecliptic'
! Define which modules are temperature and which are polarization; should be moved elsewhere

if (present(cutfraction)) then
	call ds_assert(cutfraction>0 .and. cutfraction<=1,"Invalid cutfraction specified")
endif

unit = ds_get_lun()
call MPI_Comm_size(correlator%comm, nproc,ierr)
call MPI_Comm_rank(correlator%comm, rank, ierr)

dummy=1
call MPI_Allreduce(MPI_IN_PLACE, dummy, 1, MPI_INTEGER,MPI_SUM,correlator%comm,ierr)
call MPI_Allreduce(MPI_IN_PLACE, dummy, 1, MPI_INTEGER,MPI_SUM,correlator%comm,ierr)
call MPI_Allreduce(MPI_IN_PLACE, dummy, 1, MPI_INTEGER,MPI_SUM,correlator%comm,ierr)
!JZ This line is here because otherwise on Titan some nodes don't start until the others reach the first mpi call.

correlator%nproc = nproc
correlator%proc = rank

!Get the data info from the database and then turn it into a format more useful for us
call get_target_list(one,'MJD',target_type, target_name, target_list)
call ds_assert(size(target_list)>0,"No targets found with type:" // trim(target_type) // " name:"//trim(target_name) )
call removeTemperatureFromModuleList(module_list_in, module_list, nmodules)
call buildTimestreamInfo(exclusionFile,target_list,timestream_info,module_list)

npair = timestream_info%n
ndet = npair*2

correlator%ndet = ndet
allocate(correlator%owner(0:ndet-1))
allocate(alpha(0:ndet-1))
allocate(sigma(0:ndet-1))
allocate(fknee(0:ndet-1))
allocate(correlator%ntod(0:ndet-1))
correlator%ntod=0


call create_azimuth_cache(ndet) !creates cache structure with no allocated memory

allocate(ndet_proc(0:nproc-1))


!Work out the number of detectors each proc is responsible for
nmin=floor((1.0*ndet)/nproc)
if (mod(nmin,2)==1) nmin=nmin-1
nextra = ndet - nmin*nproc
ndet_proc = nmin
ndet_proc(0:nextra/2-1) = nmin + 2 

call ds_assert(sum(ndet_proc) == ndet, "Mistake in assignment of dets - not all assigned.")

my_ndet = ndet_proc(rank)
allocate(correlator%my_det(0:my_ndet-1))

my_npair = my_ndet/2
if(rank==0) then
	write(message,*) "Approx number of detectors per proc:", my_ndet
	call ds_log(message,ds_feedback_quiet)
endif

correlator%my_ndet = my_ndet

allocate(timestream(0:my_ndet-1))
allocate(pointing(0:my_npair-1))
allocate(offsets(0:my_ndet-1))

!Assign responsibility for the detector pairs
dets_for_this_proc = 0
current_proc=0
do p=0,npair-1
	q = det_qindex(p)
    u = det_uindex(p)
	correlator%owner(q)=current_proc
	correlator%owner(u)=current_proc	
    dets_for_this_proc = dets_for_this_proc + 2
    if (dets_for_this_proc == ndet_proc(current_proc)) then
       dets_for_this_proc=0
       current_proc = current_proc + 1
    endif
enddo


do d=0,ndet-1
   call ds_assert(correlator%owner(d) >= 0, "Unassigned detector")
enddo

!Recored which detectors this proc owns
my_det=0
do d=0,ndet-1
   if (correlator%owner(d)==correlator%proc) then
      correlator%my_det(my_det) = d
      my_det = my_det + 1
   endif
enddo


fknee=0
alpha=0
sigma=0
correlator%ntod=0
my_det=0
my_pair=0
haveBadScans=.false.
if(rank==0) then
	write(message,*) "Beginning Data Load"
	call ds_log(message,ds_feedback_quiet)
endif

current_filename=""
!Flat is better than nested...
do p=0,npair-1
	q=det_qindex(p)
	u=det_uindex(p)
	if (.not. correlator%owner(q)==rank) cycle   !This is a nice way of avoiding nested if statements

    call ds_assert(correlator%owner(u)==rank,"Something went wrong with pair ownership")

	if (my_pair>0 .and. mod(my_pair,50)==0 .and. rank==0) then
		write(message,*) "Rank 0 completed ",my_det, " of ", my_ndet, "detectors"
		call ds_log(message,ds_feedback_noisy)
	endif

    !Check to see if we already have loaded the data from the file containing this data
    !do this by comparing the loaded filename to the current filename
    if (.not. current_filename==timestream_info%filename(p)  )     then
		thisScanBad=.false.
       if (.not. (current_filename=="")) then
			call deallocate_module_struct(tod_data)
			call deallocate_module_struct(tod_data_in)
			deallocate(psi)
			deallocate(pix)
			deallocate(az)
		endif	
		current_filename = timestream_info%filename(p)
	!	write(*,*) rank, " Reading ", trim(current_filename),my_pair,'of',my_npair-1
		call l2_read(unit,current_filename,data=tod_data_in,status=status)
		if (.not. (status==0)) then 
			write(*,*) "L2 read failed on " // trim(current_filename)
			haveBadScans=.true.
			thisScanBad=.false.
			cycle
		endif
!		call ds_assert(status==0, "L2 read failed on "//trim(current_filename))
		call reorganize_data(.false.,ecliptic,module_list, nside, tod_data_in, tod_data, pix, az,psi)
	endif

	!The parameters of this timestream: module number, pair number, and which the diodes are.
	modnum = timestream_info%mod(p)
	if (timestream_info%pair(p)==0) then
		qdiode = 0
		udiode = 1
	else
		qdiode = 3
		udiode = 2
	endif

	!Find the current data
	module_data=>tod_data(modnum+1)
	ntod = size(module_data%time)
	
	if (present(cutfraction)) ntod = ntod * cutfraction
	
	!JZ Cut the length down to a multiple of the offsetLength
	na = ntod/offsetLength
	ntod = na * offsetLength
	correlator%ntod(q)=ntod
	correlator%ntod(u)=ntod


!    if (my_det==0) write(*,*) "WARNING! Dividing timestream by ntod to compensate for simulation error."
    !Load TOD
    call prepareTimestream(timestream(my_det), ntod)
	my_pair = my_det/2
    my_q = my_det
    my_det = my_det+1
    call prepareTimestream(timestream(my_det), ntod)
    my_u = my_det
    my_det = my_det+1
!#warning dividing by module size
    if (timestream_info%pair(p)==0) then
       timestream(my_q)%timestream = module_data%tod(0,1:ntod)!/size(module_data%time)
       timestream(my_u)%timestream = module_data%tod(1,1:ntod)!/size(module_data%time)
    else
       call ds_assert(timestream_info%pair(p)==1,"The pair is neither 0 nor 1")
       timestream(my_q)%timestream = module_data%tod(3,1:ntod)!/size(module_data%time)
       timestream(my_u)%timestream = module_data%tod(2,1:ntod)!/size(module_data%time)
    endif

	!Prepare offsets.
	call prepareOffsets(offsets(my_q),na,offsetLength)
	call prepareOffsets(offsets(my_u),na,offsetLength)

	!If needed, note that azimuth filter is used.
	if (timestream_info%needs_azfilter(p)) then
		!turn on the azimuth flags in the ds_offset structures
		!offsets(my_q)%azimuth_flag= .true.
		!offsets(my_u)%azimuth_flag= .true.
		!call az_create_pointing(q,ntod,corrlength,n_bins,az(1:ntod,modnum),offsets(my_q)%az)
		!call az_create_pointing(u,ntod,corrlength,n_bins,az(1:ntod,modnum),offsets(my_u)%az)
	endif
				     
	!Load pointing
	call preparePointing(pointing(my_pair),ntod)
	do t=1,ntod
		pointing(my_pair)%pixel(t) = pix(1,t,modnum+1)
		pointing(my_pair)%theta(t) = psi(qdiode,t,modnum+1)
	enddo


	!Load the noise parameters
	if (lbound(module_data%time,1) .gt. ubound(module_data%time,1) ) then
		write(*,*) "Data failed on " // trim(current_filename)

		thisScanBad=.true.
		haveBadScans=.true.
		cycle
	endif

	call get_calibration(run=timestream_info%run(p),segment=timestream_info%scan(p),  &
		mjd=module_data%time(1),module=modnum,diode=qdiode,  &
		sigma0=sigma(q),gain=gain,f_knee=fknee(q),alpha=alpha(q))  
	
		timestream(my_q)%timestream = timestream(my_q)%timestream / gain

	call get_calibration(run=timestream_info%run(p),segment=timestream_info%scan(p),  &
		mjd=module_data%time(1),module=modnum,diode=udiode,  &
		sigma0=sigma(u),gain=gain,f_knee=fknee(u),alpha=alpha(u) )  
		
	timestream(my_u)%timestream = timestream(my_u)%timestream / gain


enddo

call deallocate_module_struct(tod_data)
call deallocate_module_struct(tod_data_in)
deallocate(psi)
deallocate(pix)
deallocate(az)
deallocate(ndet_proc)

write(message,*) "Rank ", rank," loaded all data."
call ds_log(message,ds_feedback_noisy)

do d=0,my_ndet-1
!	write(*,*) rank, "subtracting mean", d, " of ", my_ndet
	call ds_assert(timestream(d)%nt > 0, "Zero length timestream located")
	mu=0.0
	do i=1,timestream(d)%nt
		mu = mu + timestream(d)%timestream(i)
	enddo
	mu=mu/timestream(d)%nt
	do i=1,timestream(d)%nt
		timestream(d)%timestream(i) = timestream(d)%timestream(i) - mu
	enddo
enddo
!write(*,*) "Rank ", rank," done mean removal allreduce."

call MPI_Allreduce(MPI_IN_PLACE, correlator%ntod, ndet, MPI_INTEGER,MPI_SUM,correlator%comm,ierr)
!write(*,*) "Rank ", rank," done first allreduce."
if (rank==0) then
	write(message,*) "Mean scan length = ", sum(correlator%ntod)/size(correlator%ntod)
	call ds_log(message,ds_feedback_noisy)
	write(message,*) "Max scan length = ", maxval(correlator%ntod)
	call ds_log(message,ds_feedback_noisy)
	write(message,*) "Min scan length = ", minval(correlator%ntod)
	call ds_log(message,ds_feedback_noisy)
endif 

call ds_assert(.not. haveBadScans, "Bad scans reported")

call MPI_Allreduce(MPI_IN_PLACE, alpha, ndet, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
!write(*,*) "Rank ", rank," done second allreduce."
call MPI_Allreduce(MPI_IN_PLACE, sigma, ndet, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
!write(*,*) "Rank ", rank," done third allreduce."
call MPI_Allreduce(MPI_IN_PLACE, fknee, ndet, MPI_REAL8,MPI_SUM,correlator%comm,ierr)
!write(*,*) "Rank ", rank," done fourth allreduce."

do d=0,ndet-1
    call ds_assert(correlator%ntod(d) > 0, "Unassigned ntod")
enddo




end subroutine readDataAssignWork

subroutine init_timestream_info(tsinfo, n)
type(ds_timestreamInfo) :: tsinfo
integer n
tsinfo%n = n
allocate(tsinfo%run(0:n-1))
allocate(tsinfo%scan(0:n-1))
allocate(tsinfo%mod(0:n-1))
allocate(tsinfo%pair(0:n-1))
allocate(tsinfo%filename(0:n-1))
allocate(tsinfo%needs_azfilter(0:n-1))
tsinfo%filename = ""
tsinfo%run=-1
tsinfo%scan=-1
tsinfo%mod=-1
tsinfo%pair=-1
tsinfo%needs_azfilter=.false.
end subroutine init_timestream_info

subroutine buildTimestreamInfo(exfile,target_list,timestream_info,module_list)
	type(ds_timestreamInfo) :: timestream_info
	type(l2_filelist),   pointer, dimension(:) :: target_list
	integer(i4b), pointer,     dimension(:)       :: module_list
	character(*) exfile
	character(128) :: filename
	type(ds_exlist) :: exlist
	integer run,scan,m,p,t,ntimestreams,f,n

	call buildExclusions(exfile,exlist,module_list)
	ntimestreams = 0
	
	
	!Loop through the timestreams once to find how many we are actually using
	do t=1,size(target_list)
		do f=1,target_list(t)%numfiles
			run = target_list(t)%id(f,1)   !The ID number, typically in the hundreds/thousands
			scan = target_list(t)%id(f,2)  !The segment number, 1 - 6 or so (or a bit higher)
			do m=0,nmodulemax-1
				do p=0,1
					if (.not. is_excluded(exlist,run,scan,m,p)) then
						ntimestreams = ntimestreams + 1
					endif
				enddo
			enddo
		enddo
	enddo

	call init_timestream_info(timestream_info,ntimestreams)

n=0
!Now loop through again to populate the structure.
do t=1,size(target_list)
	do f=1,target_list(t)%numfiles
		filename = target_list(t)%filenames(f)
		run = target_list(t)%id(f,1)   !The ID number, typically in the hundreds/thousands
		scan = target_list(t)%id(f,2)  !The segment number, 1 - 6 or so
		do m=0,nmodulemax-1
			do p=0,1
				if (.not. is_excluded(exlist,run,scan,m,p)) then
					timestream_info%run(n)=run
					timestream_info%filename(n)=filename
					timestream_info%scan(n)=scan
					timestream_info%mod(n)=m
					timestream_info%pair(n)=p
		     		n=n+1
		 		endif
			enddo
		enddo
	enddo
enddo
	
	
end subroutine buildTimestreamInfo






subroutine init_exclusion(exlist, runmax, segmax)
type(ds_exlist) :: exlist
integer :: runmax, segmax

exlist%runmax = runmax
exlist%segmax=segmax
allocate(exlist%excluded(runmax,segmax,0:nmodulemax-1, 0:npairdiodesmax-1))
exlist%excluded=-1
end subroutine

subroutine set_accepted(exlist,run,seg)
	type(ds_exlist) :: exlist
	integer :: run,seg
	call ds_assert(run>=0, "Negative run number found (acc)")
	call ds_assert(run<=exlist%runmax, "Run number larger than maximum (acc)")
	call ds_assert(seg>0, "Negative segment number found (acc)")
	call ds_assert(seg<=exlist%segmax, "Segment number larger than maximum (acc)")
	
	exlist%excluded(run,seg,:,:)=0
end subroutine set_accepted

subroutine set_excluded(exlist, run,seg,mod,pair)
type(ds_exlist) :: exlist
integer :: run,seg, mod,pair
integer accepted

	!This might seem a bit paranoid, but I am having a bad day.
	call ds_assert(run>=0, "Negative run number found (set)")
	call ds_assert(run<=exlist%runmax, "Run number larger than maximum (set)")
	call ds_assert(seg>0, "Negative segment number found (set)")
	call ds_assert(seg<=exlist%segmax, "Segment number larger than maximum (set)")
	call ds_assert(mod>=0, "Negative segment number found (set)")
	call ds_assert(mod<nmodulemax, "Segment number larger than maximum (set)")
	call ds_assert(pair>=0, "Negative pair number found (set)")
	call ds_assert(pair<npairdiodesmax, "Pair too largs (set)")
	accepted = exlist%excluded(run,seg,mod,pair)
	call ds_assert((accepted .ne. -1), "Tried to excluded diodes from non-accepted scan")

	exlist%excluded(run,seg,mod,pair)=1
end subroutine set_excluded

subroutine destroy_exlist(exlist)
type(ds_exlist) :: exlist
deallocate(exlist%excluded)
exlist%runmax=-1
exlist%segmax=-1
end subroutine destroy_exlist

function is_excluded(exlist,run,seg,mod,pair)
type(ds_exlist) exlist
integer run,seg,mod,pair
logical is_excluded
	call ds_assert(run>=0, "Negative run number found (is_)")
	call ds_assert(run<=exlist%runmax, "Run number larger than maximum (is_)")
	call ds_assert(seg>0, "Negative segment number found (is_)")
	call ds_assert(seg<=exlist%segmax, "Segment number larger than maximum (is_)")
	call ds_assert(mod>=0, "module segment number found (is_)")
	call ds_assert(mod<nmodulemax, "module number larger than maximum (is_)")
	call ds_assert(pair>=0, "Negative pair number found (is_)")
	call ds_assert(pair<npairdiodesmax, "Pair number larger than maximum (is_)")

	call ds_assert(exlist%excluded(run,seg,mod,pair) .ne. -1, "Tried to excluded diodes from non-accepted scan")

	is_excluded=.true.
	if (exlist%excluded(run,seg,mod,pair)==0) is_excluded=.false.

end function is_excluded


subroutine buildExclusions(exfile,exlist,module_list)
	integer(i4b), pointer,     dimension(:)       :: module_list
	character(*) exfile
	type(ds_exlist) :: exlist
	integer segment,segmax, nexclude
	integer run, runmax
	integer m,d,i
	logical module_included
	integer unit
	integer bad(0:1000)
	
	unit = ds_get_lun()
	segmax=100
	runmax=20000
	open(unit=unit,file=exfile)
!	do
!		read(unit,*, end=10) run, segment
!		if (segment .gt. segmax) segment = segmax
!		if (run .lt. runmax) run = runmax
!	enddo
!10 rewind(unit)

	call init_exclusion(exlist, runmax, segmax)

	do
		read(unit,*, end=20) run, segment, nexclude
		call set_accepted(exlist,run,segment)
		!Exclude all the full module exclusions, in case they are not in the list, though they should be

		!Loop through each module ID
		do m=0,nmodulemax-1
			!check to see if it is included
			module_included=.false.
			do i=1,size(module_list)
				if (module_list(i)==m) then
					module_included=.true.
					exit
				endif
			enddo
			!If it is not included exclude both pairs in it.
			if (.not. module_included) then
				call set_excluded(exlist,run,segment,m,0)
				call set_excluded(exlist,run,segment,m,1)				
			endif
		enddo

		!Now re-read the line and remove all the scan-specific exclusions
		backspace(unit)
		read(unit,*) run, segment, nexclude, bad(0:2*nexclude-1)


		do i=0,nexclude-1
			m=bad(2*i)
			d=bad(2*i+1)
			if (d==0 .or. d==3) then
				call set_excluded(exlist,run,segment,m,0)
			else
				call set_excluded(exlist,run,segment,m,1)
			endif
		enddo
	enddo
	
20	 close(unit)
end subroutine buildExclusions

subroutine removeTemperatureFromModuleList(list_in, list_out,nmod)
integer(i4b), pointer,     dimension(:)       :: list_in, list_out
integer n,m,i
integer, intent(out) :: nmod
logical(lgt),              dimension(0:18) :: temperature_module
temperature_module(0:16)  = .false.
temperature_module(17:18) = .true.
nmod = 0
n=size(list_in)

do i=1,n
	if (.not. temperature_module(list_in(i))) nmod = nmod+1
enddo
allocate(list_out(nmod))
m=1
do i=1,n
	if (.not. temperature_module(list_in(i))) then
		list_out(m) = list_in(i)
		m=m+1
	endif
enddo

end subroutine removeTemperatureFromModuleList

end module ds_oslo_fitstools
