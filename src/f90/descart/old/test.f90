program testMap
  use healpix_types
  use random, only : random_normal
  use pix_tools
  use fitstools
  use quiet_utils
  use l2_fileutils
  use quiet_fileutils
  use quiet_calib_mod
  use quiet_noise_mod
  use ds_types
  use ds_azimuth
  use ds_multidetector
  use ds_fitstools
	implicit none
	integer err
	integer npix,n,ntotal,size,rank,nt,l,na
	real(dp), allocatable, dimension(:) :: noise_std,noise_var
	type(ds_map) :: TrueQ, TrueU
	type(ds_detpointing), allocatable, dimension(:) ::  p
	type(ds_timestream),  allocatable, dimension(:) ::  t
	type(ds_offsets),  allocatable, dimension(:) ::  a,b
	type(ds_correlator) :: correlator
	integer i,j
	integer, dimension(:), allocatable :: seed
	integer :: seed_size
	character(128) :: filename
	integer :: clock1, clock2, clock3
	character(*), parameter :: filename_format = "('offset.',I1,'.',I1)"

	call MPI_Init(err)	
	call mpi_comm_rank(MPI_COMM_WORLD, rank, err)
	call mpi_comm_size(MPI_COMM_WORLD, size, err)	
	call system_clock(clock1,clock2,clock3)
	call random_seed
	call random_seed(size=seed_size)
	allocate(seed(seed_size))
	seed = clock1
	call random_seed(put=seed)

	npix=500
	n=2
	ntotal=n*size
	nt=10000
	l=100
	na=nt/l

	allocate(noise_std(0:ntotal-1))
	allocate(noise_var(0:ntotal-1))

	noise_std = 100.0
	noise_var = noise_std * noise_std

	call prepareMap(TrueQ,npix)
	call prepareMap(TrueU,npix)
	do i=1,npix
		TrueQ%map(i)=1.0*i
		TrueQ%indices(i)=i
		TrueU%map(i)=1.0*i
		TrueU%indices(i)=i
	enddo



	allocate(a(0:n-1))
	allocate(b(0:n-1))
	allocate(t(0:n-1))
	allocate(p(0:n/2-1))

	do i=0,n/2-1
		call preparePointing(p(i),nt,npix,.true.)
		do j=1,nt
			p(i)%pixel(j)=mod(j,npix)+1
			p(i)%theta(j)=mod(j,2)*pi/4.0
		enddo
	enddo
	
	

	do i=0,n-1
		call prepareTimestream(t(i),nt)
	
		!Add signal
		if (mod(i,2)==0) then
			call projectMapOntoTimestream(TrueQ,t(i),p(i/2))
		else
			call projectMapOntoTimestream(TrueU,t(i),p(i/2))		
		endif
	
		!Add noise
		do j=1,nt
			t(i)%timestream(j) = t(i)%timestream(j) + noise_std(rank*size+i) * random_normal()
		enddo
		!Convert to offsets
		call prepareOffsets(a(i),na,l)
		call prepareOffsets(b(i),na,l)
		call deprojectTimestreamOntoOffset(t(i),a(i))
	enddo

	correlator%my_ndet=n
	correlator%ndet=ntotal
	correlator%comm=MPI_COMM_WORLD
	allocate(correlator%my_det(0:n-1))
	do i=0,n-1
		correlator%my_det(i) = rank*size+i
	enddo
	

	call apply_FtZF(correlator,p,noise_var,npix,a,b)
	!The b offsets should now have no signal and have been divided by the white noise.
	
	
	do i=0,n-1
		write(filename,filename_format) rank,i
		open(file=trim(filename),unit=12)
		do j=1,na
			write(12,*) b(i)%values(j)
		enddo
		close(12)
	enddo

end program testMap
