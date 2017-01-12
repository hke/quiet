program tod_gibbs
   use quiet_constrained_mod
   use quiet_hdf_mod
   implicit none

   type data_struct
      real(dp),     allocatable :: tod(:), mask(:), templates(:,:), point(:,:)
      integer(i4b), allocatable :: t2p(:)
      real(dp)              :: srate, sangle
      integer(i4b)          :: n, nt, npix
   end type

   type gibbs_params
      real(dp),     allocatable :: noise(:), junk(:), cmb(:), a(:)
      real(dp)                  :: sigma0, fknee, alpha
   end type

   call main

contains

   subroutine main
     implicit none
     character(len=512) :: arg, fname
     type(data_struct)  :: data
     type(gibbs_params) :: gibbs
     integer(i4b)       :: i, j, k, nside
     call getarg(1, fname)
     nside = 512


     call read_data_file(fname, data)
     call setup(data, gibbs)
   end subroutine

   subroutine read_data_file(fname, data)
     implicit none
     character(len=*)   :: fname
     type(data_struct)  :: data
     type(hdf_file)     :: file
     integer(i4b)       :: ext(7)
     call free_data_struct(data)
     call open_hdf_file(fname, file "r")
     call get_size_hdf(file, "templates", ext)
     data%n  = ext(1)
     data%nt = ext(2)
     allocate(data%tod(n), data%templates(n, nt), data%mask(n), data%point(3,n))
     call read_hdf(file, "data",      data%tod)
     call read_hdf(file, "templates", data%templates)
     call read_hdf(file, "mask",      data%mask)
     call read_hdf(file, "samprate",  data%srate)
   end subroutine

   subroutine free_data_struct(data)
     implicit none
     type(data_struct) :: data
     if(allocated(data%tod))       deallocate(data%tod)
     if(allocated(data%templates)) deallocate(data%templates)
     if(allocated(data%mask))      deallocate(data%mask)
     if(allocated(data%point))     deallocate(data%point)
     if(allocated(data%t2p))       deallocate(data%t2p)
   end subroutine

   subroutine free_gibbs_params(p)
     implicit none
     type(gibbs_params) :: p
     if(allocated(p%noise)) deallocate(p%noise)
     if(allocated(p%junk))  deallocate(p%junk)
     if(allocated(p%cmb))   deallocate(p%cmb)
     if(allocated(p%a))     deallocate(p%a)
   end subroutine

end program
