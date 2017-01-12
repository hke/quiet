program gen_covar_mat
  use healpix_types
  implicit none

  integer(i4b) :: i, j, numpar, numsamples, num_burn_in, numchain
  integer(i4b) :: unit
  character(len=2)   :: chain_text
  character(len=128) :: param, filename, prefix

  real(dp), allocatable, dimension(:,:,:) :: samples
  real(dp), allocatable, dimension(:,:)   :: covar, L
  real(dp), allocatable, dimension(:)     :: vals, avg

  unit = 38
  
  call getarg(1,prefix)

  call getarg(2,param)
  read(param,*) numchain

  call getarg(3,param)
  read(param,*) numpar

  call getarg(4,param)
  read(param,*) numsamples

  call getarg(5,param)
  read(param,*) num_burn_in

  ! Read samples
  allocate(samples(numpar, numsamples, numchain))
  allocate(vals(numpar+2))

  do i = 1, numchain

     call int2string(i, chain_text)
     filename = trim(prefix) // '_no' // chain_text // '.dat'
     
     open(unit, file=trim(filename))

     do j = 1, numsamples
        read(unit,*) vals
        samples(:,j,i) = vals(2:1+numpar)
     end do

     close(unit)

  end do

  ! Compute covariance matrix
  allocate(avg(numpar))
  allocate(covar(numpar, numpar))
  allocate(L(numpar, numpar))
  
  do i = 1, numpar
     avg(i) = sum(samples(i,num_burn_in+1:numsamples,:)) / &
          & real((numsamples-num_burn_in)*numchain,dp)
  end do

  do i = 1, numpar
     do j = 1, numpar
        covar(i,j) = sum((samples(i,num_burn_in+1:numsamples,:)-avg(i)) * &
             & (samples(j,num_burn_in+1:numsamples,:)-avg(j)))
        covar(i,j) = covar(i,j) / real((numsamples-num_burn_in)*numchain,dp)
        covar(j,i) = covar(i,j)
     end do
     if (covar(i,i) == 0.d0) then
        covar(:,i) = 0.d0
        covar(i,:) = 0.d0
        covar(i,i) = 1.d0
     end if
  end do

  call cholesky_decompose(covar, L)
  

  ! Output file
  open(unit,file='param_covar_mat.unf', form='unformatted')
  write(unit) L
  close(unit)

  deallocate(avg)
  deallocate(covar)
  deallocate(samples)
  deallocate(L)


contains

  ! Small utility for converting an integer to a string
  subroutine int2string(integer, string)
    implicit none

    integer(i4b),     intent(in)  :: integer
    character(len=*), intent(out) :: string

    integer(i4b)               :: temp_int, i, k

    temp_int = integer
    do i = 1, len(string)
       k = temp_int / 10**(len(string)-i)
       write(string(i:i),'(I1)') k
       temp_int = temp_int - k * 10**(len(string)-i)
    end do

  end subroutine int2string

  subroutine cholesky_decompose(A, L)
    implicit none

    real(dp), dimension(:,:), intent(in)  :: A
    real(dp), dimension(:,:), intent(out) :: L

    integer(i4b) :: N, i, j, k
    real(dp) :: temp
    real(dp), allocatable, dimension(:) :: temp_row

    N = size(A(1,:))

    L = 0.d0

    do j = 1, N
       do i = j, N

          temp = 0.d0

          do k = 1, i-1
             temp = temp + L(i,k) * L(j,k)
          end do

          if (i == j) then
             L(j,j) = sqrt(A(j,j)-temp)
          else
             L(i,j) = (A(i,j) - temp) / L(j,j)
          end if

       end do
    end do

  end subroutine cholesky_decompose


end program gen_covar_mat
