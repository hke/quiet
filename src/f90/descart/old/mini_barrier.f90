program minibarrier


implicit none
include "mpif.h"

integer ierr,rank,nproc

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

write(*,*) rank, "begun the code"
call MPI_Barrier(MPI_COMM_WORLD,ierr)
write(*,*) rank, "crossed the barrier"

do while (.true.)

enddo

call MPI_Finalize(ierr)
end program minibarrier
