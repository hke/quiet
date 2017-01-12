module ds_montecarlo
use random, only : random_normal
use ds_maptools
use ds_mapping
use ds_simple_prior
use ds_types
implicit none

public :: simulate_signal_scans
public :: simulate_noise_scans
public :: simulate_signal_noise_scans


private 

contains

subroutine simulate_signal_scans(scans,naive_q,naive_u,comm,nside,originalIndices,signal_q,signal_u)
	type(ds_modulescan), dimension(:) :: scans
	type(ds_map) :: naive_q,naive_u    !naive maps, to be accumulated
	type(ds_map),optional :: signal_q,signal_u  !Healpix maps
	integer,optional :: nside
	integer comm
	integer(i8b), dimension(:), optional :: originalIndices

	call generate_simulated_scans(scans,naive_q,naive_u,comm,nside,originalIndices,signal_q,signal_u)

end subroutine simulate_signal_scans

subroutine simulate_noise_scans(scans,naive_q,naive_u,comm,noise)
	type(ds_modulescan), dimension(:) :: scans
	type(ds_noiseinfo), optional, dimension(:)  :: noise		
	type(ds_map) :: naive_q,naive_u    !naive maps, to be accumulated
	integer comm
	
	call generate_simulated_scans(scans,naive_q,naive_u,comm,noise=noise)

end subroutine simulate_noise_scans

subroutine simulate_signal_noise_scans(scans,naive_q,naive_u,comm,nside,originalIndices,signal_q,signal_u,noise)
	type(ds_modulescan), dimension(:) :: scans
	type(ds_map) :: naive_q,naive_u    !naive maps, to be accumulated
	type(ds_map),optional :: signal_q,signal_u  !Healpix maps
	type(ds_noiseinfo), optional, dimension(:)  :: noise		
	integer,optional :: nside
	integer comm
	integer(i8b), dimension(:), optional :: originalIndices

	call generate_simulated_scans(scans,naive_q,naive_u,comm,nside,originalIndices,signal_q,signal_u,noise)

end subroutine simulate_signal_noise_scans


subroutine generate_simulated_scans(scans,naive_q,naive_u,comm,nside,originalIndices,signal_q,signal_u,noise)
	use ds_mapping
	type(ds_modulescan), dimension(:) :: scans
	type(ds_map) :: naive_q,naive_u    !naive maps, to be accumulated
	type(ds_map),optional :: signal_q,signal_u  !Healpix maps
	type(ds_noiseinfo), optional, dimension(:)  :: noise		
	integer,optional :: nside
	integer comm
	integer i
	integer npix
	integer ierror
	integer(i8b), dimension(:), optional :: originalIndices
	logical :: do_signal, do_noise
	  include "mpif.h"
	
	npix = p_npix
	
	call prepareMap(naive_q,npix)
	call prepareMap(naive_u,npix)
	
	do_signal = present(signal_q)
	do_noise = present(noise)

	call ds_assert(do_signal .eqv. present(signal_u), "Must specify either both or neither signal Q and U maps (1)" )
	call ds_assert(do_signal .eqv. present(originalIndices), "Must specify original indices if simulating signal" )
	call ds_assert(do_signal .eqv. present(nside), "Must specify original indices if simulating signal" )

	do i=lbound(scans,1),ubound(scans,1)
		if (do_signal .and. do_noise ) then  
			call generate_simulated_scan(scans(i),naive_q,naive_u, nside,originalIndices,signal_q=signal_q,signal_u=signal_u,noise=noise(i))
		else if (do_signal) then
			call generate_simulated_scan(scans(i),naive_q,naive_u, nside,originalIndices,signal_q=signal_q,signal_u=signal_u)
		else if (do_noise) then
			call generate_simulated_scan(scans(i),naive_q,naive_u,noise=noise(i))
		else
			call ds_assert(.false., "Neither signal nor noise specified in simulation (1)")
		endif
	enddo

	call MPI_Allreduce(MPI_IN_PLACE, naive_q%map, naive_q%npix, MPI_REAL8,MPI_SUM,comm,ierror)
	call MPI_Allreduce(MPI_IN_PLACE, naive_u%map, naive_u%npix, MPI_REAL8,MPI_SUM,comm,ierror)
	call ds_assert(ierror==0,"MPI Fail in scan simulation")
	call cov_mult(naive_q,naive_u)

end subroutine generate_simulated_scans




subroutine set_scan_zero(scan)
	type(ds_modulescan) :: scan
	integer i
	logical :: do_o, do_t
	
	do_o = associated(scan%offsets)
	do_t = associated(scan%timestream)
	
	do i=1,ndiodes_max
		if (scan%diode_flag(i)==0) cycle
		if (do_t) scan%timestream(i)%timestream(:)=0
		if (do_o) scan%offsets(i)%values(:)=0
	enddo
	
end subroutine set_scan_zero

subroutine noise_simulation(scan,noise)
	type(ds_modulescan) :: scan
	type(ds_noiseinfo) :: noise
	type(ds_timestream), dimension(ndiodes_max) :: white_noise
	integer i,j,t
	real(dp) :: cij
	real(dp), dimension(ndiodes_max,ndiodes_max) :: correlation_sqrt

	
	
	!Set up modulescan timstreams
	if (.not. associated(scan%timestream)) then
		allocate(scan%timestream(ndiodes_max))
		do i=1,ndiodes_max
			if (scan%diode_flag(i)==1) call prepareTimestream(scan%timestream(i), scan%ntod)
		enddo
	endif
	
	
	!Set up white noise timestreams and fill them with white noise
	do i=1,ndiodes_max
		if (scan%diode_flag(i)==0) cycle
		call prepareTimestream(white_noise(i),scan%ntod)
		do t=1,scan%ntod
			white_noise(i)%timestream(t)=random_normal()*noise%sigma(i)
		enddo
	enddo	
	
	!Set up square root of correlation matrix
	
	correlation_sqrt = noise%corrs
	call matrix_sqrt(ndiodes_max,correlation_sqrt)	
	
	
	!Add white noise according to correlation matrix.
	do i=1,ndiodes_max
		if (scan%diode_flag(i)==0) cycle
		do j=1,ndiodes_max
			cij = correlation_sqrt(i,j)
			if (cij .ne. 0 ) scan%timestream(i)%timestream = scan%timestream(i)%timestream + cij*white_noise(i)%timestream
		enddo
	enddo
	
	!Clean up
	do i=1,ndiodes_max
		if (scan%diode_flag(i)==0) cycle
		call destroyTimestream(white_noise(i))
	enddo		

	
	
end subroutine noise_simulation
	
subroutine one_over_f_simulation(x,sample_frequency,sigma,fknee,alpha)
	real(dp), dimension(0:) :: x
	real(dp) :: sample_frequency
	complex(dpc), dimension(:), allocatable :: noise_ft
	real(dp), dimension(:), allocatable :: noise
	real(dp) :: sigma,fknee,alpha
	real(dp) :: f, power, df, nyquist
	integer nx
	integer nf
	integer i
	
	nx = size(x,1)
	nf = nx/2+1
	allocate(noise_ft(0:nf-1))
	noise_ft(0)=0.0
	nyquist=sample_frequency/2.
	df = 2.0_dp * nyquist/nx
	do i=1,nf-1
		f=i*df
		power = sigma * (1. + (fknee/f)**alpha ) 
		noise_ft(i) = random_normal() * sqrt(power) * random_phase()
	enddo
	allocate(noise(nx))
	noise=0
	call ds_ifft(noise_ft,noise)
	x=x+noise
	deallocate(noise)
	deallocate(noise_ft)
	
end subroutine

function random_phase() result(z)
	real(dp) :: phi
	complex(dpc) :: z
	complex(dpc), parameter :: I = (0.0,1.0)
	call random_number(phi)
	phi=phi*2*PI
	z = exp(I*phi)
end function

subroutine signal_simulation(scan,qmap,umap,nside,originalIndices)
	use pix_tools, only : nside2npix
	type(ds_modulescan) :: scan
	integer nside,npix
	type(ds_map) :: qmap,umap  !Healpix maps
	integer(i8b), dimension(:) :: originalIndices
	integer i,t,p,h
	real(dp) :: q,u,theta
	
	npix=nside2npix(nside)
	call ds_assert(qmap%npix==umap%npix, "maps have difference npixes in signal simulation")	
	call ds_assert(qmap%npix==npix, "Q map has wrong number of indices in signal simulation")

	if (.not. associated(scan%timestream)) then
		allocate(scan%timestream(ndiodes_max))
		do i=1,ndiodes_max
			if (scan%diode_flag(i)==1) call prepareTimestream(scan%timestream(i), scan%ntod)
		enddo
	endif


	do i=1,ndiodes_max
		if (scan%diode_flag(i)==0) cycle
		do t=1,scan%ntod
			theta=scan%theta(t)+scan%dpsi(i)
			p=scan%pointing(t)
			h=originalIndices(p)
			q=qmap%map(h)
			u=umap%map(h)
			scan%timestream(i)%timestream(t) = scan%timestream(i)%timestream(t) +  &
				cos(2*theta)*q + sin(2*theta)*u
		enddo
	enddo
end subroutine signal_simulation




subroutine generate_simulated_scan(scan,naive_q,naive_u,nside,originalIndices,noise,signal_q,signal_u)
	type(ds_modulescan) :: scan
	
	type(ds_map) :: naive_q,naive_u    !naive maps, to be accumulated
	type(ds_map), optional  :: signal_q,signal_u  !Healpix maps
	type(ds_noiseinfo), optional  :: noise	
	integer, optional :: nside
	integer i
	integer(i8b), dimension(:), optional :: originalIndices
	
	
	
	call ds_assert(present(signal_q) .eqv. present(signal_u), "Must specify either both or neither signal Q and U maps" )
	call ds_assert(present(signal_q) .or. present(noise), "Neither signal nor noise specified in simulation")
	call ds_assert(present(signal_q) .eqv. present(originalIndices), "Must specify original indices if simulating signal" )
	call ds_assert(present(signal_q) .eqv. present(nside), "Must specify original indices if simulating signal" )
	

	if (.not. associated(scan%timestream)) allocate(scan%timestream(ndiodes_max))
	do i=1,ndiodes_max
		if (scan%diode_flag(i) .ne. 1) cycle
		if (.not. allocated(scan%timestream(i)%timestream )) call prepareTimestream(scan%timestream(i),scan%ntod)
	enddo

	!Zero the scan and add signal,noise, or both as desired.
	call set_scan_zero(scan)
	if (present(signal_q)) call signal_simulation(scan,signal_q,signal_u,nside,originalIndices)
	if (present(noise)) call noise_simulation(scan,noise)

	!Accumulate naive map
	call add2rhs(scan,naive_q,naive_u)
	
	!Convert to offset space and kill timestreams
	do i=1,ndiodes_max
        if (scan%diode_flag(i)==0) cycle !ignore if bad
        call deprojectTimestreamOntoOffset(scan%timestream(i),scan%offsets(i))
	    call destroyTimestream(scan%timestream(i))
	enddo
	deallocate(scan%timestream)
	
end subroutine generate_simulated_scan
end module 