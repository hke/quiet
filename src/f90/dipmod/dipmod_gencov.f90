program gen_covar_mat
  use quiet_utils
  use healpix_types
  implicit none

  integer(i4b) :: i, j, numpar, numsamples, num_burn_in, numchain
  integer(i4b) :: unit
  character(len=2)   :: chain_text
  character(len=128) :: param, filename

  real(dp), allocatable, dimension(:,:,:) :: samples
  real(dp), allocatable, dimension(:,:)   :: covar, L
  real(dp), allocatable, dimension(:)     :: vals, avg

  unit = 38
  
  call getarg(1,param)
  read(param,*) numchain

  call getarg(2,param)
  read(param,*) numpar

  call getarg(3,param)
  read(param,*) numsamples

  call getarg(4,param)
  read(param,*) num_burn_in

  ! Read samples
  allocate(samples(numpar, numsamples, numchain))
  allocate(vals(numpar+2))

  do i = 1, numchain

     call int2string(i, chain_text)
     filename = 'chain_no' // chain_text // '.dat'
     
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


end program gen_covar_mat
