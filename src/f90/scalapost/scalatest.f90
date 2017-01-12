program scalatest
  use quiet_utils
  use quiet_mpi_mod
  !use scalawrap
  !use scalautils
  implicit none

  integer(i4b)       :: unit, root, ierr, comm, myid, nproc, gridsize(2), context, gridpos(2)
  character(len=30)  :: kommando

  write(*,*) 'q1'
  comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  write(*,*) 'q2'
  
  call mpi_comm_rank(comm, myid, ierr)
  write(*,*) 'q2'
  call mpi_comm_size(comm, nproc, ierr)
  write(*,*) 'q3'
  call proc2grid(nproc, gridsize(1), gridsize(2))
  write(*,*) 'q4', gridsize(1), gridsize(2)
  call sl_init(context, gridsize(1), gridsize(2))
  write(*,*) 'q5', context
  call blacs_gridinfo(context, gridsize(1), gridsize(2), gridpos(1), gridpos(2))
  write(*,*) 'q6'
  
  call mpi_finalize(ierr)
  
contains
  
  subroutine proc2grid(procn, nprow, npcol)
    implicit none
    integer :: procn, nprow, npcol, bestnrow, bestprocn, i, n
    nprow = int(sqrt(real(procn,dp)))
    bestprocn = 0
    do i = nprow, max(1,int(sqrt(real(procn,dp))/2)), -1
       n = procn/i*i
       if(n > bestprocn) then; bestnrow = i; bestprocn = n; end if
    end do
    nprow = bestnrow
    npcol = procn/nprow
  end subroutine proc2grid

end program scalatest

