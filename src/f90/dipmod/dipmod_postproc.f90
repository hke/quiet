program postproc
  use healpix_types
  use quiet_utils
  implicit none

  integer(i4b)       :: i, j, num_cl_bin, numbin, numpoint, unit, k, l, ind
  character(len=128) :: binfile, paramfile, resultfile
  
  real(dp), allocatable, dimension(:,:) :: bins, binlim, lnL, grid
  real(dp), allocatable, dimension(:)   :: def_value
  real(dp), allocatable, dimension(:)   :: post_dist, dt

  unit = 24
  
  call getarg(1,paramfile)

  call get_parameter(unit, paramfile, 'NUMBIN',         par_int=num_cl_bin)
  call get_parameter(unit, paramfile, 'NUMBIN_PR_CL',   par_int=numbin)
  call get_parameter(unit, paramfile, 'BINFILE',        par_string=binfile)
  call get_parameter(unit, paramfile, 'RESULTFILE',     par_string=resultfile)
  call get_parameter(unit, paramfile, 'NUMPOINT',       par_int=numpoint)
 
  allocate(bins(num_cl_bin,2))
  allocate(binlim(num_cl_bin,2))
  allocate(def_value(num_cl_bin))
  allocate(dt(num_cl_bin))
  open(unit,file=trim(binfile))
  do i = 1, num_cl_bin
     read(unit,*) bins(i,1), bins(i,2), binlim(i,1), binlim(i,2), def_value(i)
     dt(i) = (binlim(i,2)-binlim(i,1))/real(numbin)
  end do
  close(unit)


  ! Read results
  allocate(lnL(num_cl_bin+1,numpoint))
  open(unit, file=resultfile)
  do i = 1, numpoint
     read(unit,*) lnL(:,i)
  end do
  close(unit)

  ! Normalize results
  ind = num_cl_bin+1
  lnL(ind,:) = lnL(ind,:) - maxval(lnL(ind,:))

  ! Compute unnormalized distribution
  do i = 1, numpoint
     if (lnL(ind,i) < -40.d0) then
        lnL(ind,i) = 0.d0
     else
        lnL(ind,i) = exp(lnL(ind,i))
     end if
  end do

  ! Normalize distribution
  lnL(ind,:) = lnL(ind,:) / sum(lnL(ind,:))
  

  ! Marginalize 
  allocate(post_dist(numbin))
  open(unit,file='marg_dist.dat')

  allocate(grid(numbin,numbin))
  do i = 1, num_cl_bin
     
     post_dist = 0.d0

     do j = 1, numpoint

        if (abs(lnL(i,j)-def_value(i)) > 1.d-6) then
           k = int(lnL(i,j)/dt(i))+1

           if (i == 1) then

              if (abs(lnL(3,j)-def_value(3)) < 1.d-6) then
                 post_dist(k) = post_dist(k) + lnL(num_cl_bin+1,j) * dt(2)

                 l = int(lnL(2,j)/dt(2))+1                 
                 grid(k,l) = lnL(num_cl_bin+1,j)
              end if

           else if (i == num_cl_bin) then

              if (abs(lnL(num_cl_bin-2,j)-def_value(num_cl_bin-2)) < 1.d-6) then
                 post_dist(k) = post_dist(k) + lnL(num_cl_bin+1,j) * dt(num_cl_bin-1)
              end if

           else 

              if ( (abs(lnL(i-1,j)-def_value(i-1)) > 1.d-6) .and. &
                   (abs(lnL(i+1,j)-def_value(i+1)) > 1.d-6)) then
                 post_dist(k) = post_dist(k) + lnL(num_cl_bin+1,j) * dt(i-1) * dt(i+1)
              end if

           end if
           
        end if

     end do

     if (sum(post_dist) > 0.d0) then
        post_dist = post_dist / sum(post_dist) / dt(i)
     end if

     do j = 1, numbin
        write(unit,*) binlim(i,1) + (real(j,sp)-0.5)*dt(i), post_dist(j)
     end do
     write(unit,*) 

  end do
  close(unit)

  open(unit,file='dist_b1_b2.unf', form='unformatted')
  write(unit) grid
  close(unit)
  
  deallocate(grid)
  deallocate(lnL)
  deallocate(bins)
  deallocate(binlim)
  deallocate(def_value)

end program postproc
